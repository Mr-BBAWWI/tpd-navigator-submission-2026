"""Source-bound medicinal chemistry pipeline shared by M2 jobs and portable CLI."""
import copy
import hashlib
import io
import json
from pathlib import Path
import tempfile

from packages.contracts import encoded
from packages.platform.saved_results import AUTHORITY
from packages.science.dual_e3 import SOURCE, catalog, validate, assemble, attachment_options, verify_sources

VERSION = "design-panel/20260930.4"


def input_binding():
    files = verify_sources()
    return {"format": VERSION, "source_files": files, "digest": hashlib.sha256(encoded(files)).hexdigest()}


def _stereoisomers(records):
    from rdkit import Chem
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions, GetStereoisomerCount
    from packages.science.molecules import canonical
    variants = []
    seen = set()
    for record in records:
        source = Chem.MolFromSmiles(record["mapped_smiles"])
        if source is None:
            continue
        maps = [a.GetAtomMapNum() for a in source.GetAtoms()]
        clean = Chem.Mol(source)
        for atom in clean.GetAtoms():
            atom.SetAtomMapNum(0)
        options = StereoEnumerationOptions(onlyUnassigned=True, unique=True, maxIsomers=16)
        count = GetStereoisomerCount(clean, options=options)
        for ordinal, iso in enumerate(EnumerateStereoisomers(clean, options=options)):
            smi = canonical(iso)
            key = (record.get("rule_id"), tuple(record.get("modified_atom_maps", [])), smi)
            if key in seen:
                continue
            seen.add(key)
            for atom, mapping in zip(iso.GetAtoms(), maps):
                atom.SetAtomMapNum(mapping)
            row = copy.deepcopy(record)
            row.update(mapped_smiles=Chem.MolToSmiles(iso, isomericSmiles=True), canonical_smiles=smi)
            row["id"] = "W-" + hashlib.sha256((str(key)).encode()).hexdigest()[:12]
            unresolved = sum(str(info.specified) == "Unspecified" for info in Chem.FindPotentialStereo(Chem.MolFromSmiles(smi)))
            row["stereochemistry"] = {
                "source": "explicit enumeration of unassigned stereocenters, not measured configuration",
                "isomer_index": ordinal,
                "possible_isomer_count": count,
                "enumeration_truncated": count > 16,
                "unresolved_stereo_count": unresolved,
                "status": "requires_review" if unresolved else "specified_or_not_stereogenic",
            }
            row["risk_flags"] = ["azide_reactivity_requires_review"] if row.get("rule_id") == "LH_AZIDE" else []
            variants.append(row)
    return variants


def _files(mol, identifier, put, highlight=()):
    from rdkit import Chem
    from rdkit.Chem import Draw, rdDepictor
    pic = Chem.Mol(mol)
    pic.RemoveAllConformers()
    rdDepictor.Compute2DCoords(pic)
    sdf = put(identifier + ".sdf", (Chem.MolToMolBlock(pic, forceV3000=True) + "\n$$$$\n").encode(), "chemical/x-mdl-sdfile")
    changed = [a.GetIdx() for a in pic.GetAtoms() if a.GetAtomMapNum() in highlight]
    for atom in pic.GetAtoms():
        atom.SetAtomMapNum(0)
    out = io.BytesIO()
    Draw.MolToImage(pic, size=(900, 440), highlightAtoms=changed).save(out, format="PNG")
    return {"sdf": sdf, "image": put(identifier + ".png", out.getvalue(), "image/png")}


def _ai_review(provider, context, put):
    request = {
        "model": "gpt-5.6-sol",
        "max_output_tokens": 3500,
        "instructions": "당신은 TPD 설계 결과의 검토 보조자다. 제공한 실제 계산과 출처만 사용한다. 생성한 분자를 승인하거나 효능을 주장하지 않는다. 답변은 JSON {\"notes\":[짧은 한국어 문장],\"review_questions\":[검수자 질문]}이다. CRBN과 VHL raw 점수는 교차 순위화하지 않는다.",
        "context": context,
    }
    put("design-agent-input.json", encoded(request))
    response = provider.complete(**request)
    ref = put("design-agent-response.json", encoded({"text": response.text, "metadata": response.metadata}))
    try:
        parsed = json.loads(response.text)
        assert set(parsed) == {"notes", "review_questions"}
        assert all(isinstance(parsed[k], list) and all(type(x) is str and len(x) < 4000 for x in parsed[k]) for k in parsed)
    except (ValueError, AssertionError, TypeError):
        return {"status": "response_requires_format_review", "response_ref": ref, "notes": [], "usage": response.metadata["usage"]}
    return {"status": "live_api_review_received_not_approval", "response_ref": ref, **parsed, "usage": response.metadata["usage"]}


def _scientific_context(parent, protein, source_parent_profile=None):
    from packages.science.chemical_states import (
        before_after_report, enumerate_microstates, interaction_profile,
        optimize_ligand_hydrogens,
    )
    result = {"status": "unknown", "microstates": {}, "ligand_hydrogen_preparation": {}, "parent_interactions": {}, "source_parent_interactions": source_parent_profile}
    try:
        states = enumerate_microstates(parent)
        result["microstates"] = {
            "status": "computed",
            "count": len(states),
            "records": [{"smiles": __import__("rdkit").Chem.MolToSmiles(m, isomericSmiles=True), "origin": m.GetProp("chemical_state_origin")} for m in states],
            "selection": "none; populations and pKa are unknown",
        }
    except Exception as exc:
        result["microstates"] = {"status": "unknown", "error": str(exc)}
    try:
        structures, receipt = optimize_ligand_hydrogens(parent)
        result["ligand_hydrogen_preparation"] = {
            "status": "computed_fixed_heavy_parent_reference",
            "receipt": receipt,
            "docking_input_used": False,
        }
        prepared_profile = interaction_profile(structures["hydrogens_optimized"], protein, None)
        prepared_profile["protein_hydrogen_status"] = "missing_requires_directional_review"
        prepared_profile["chemical_geometry_pass"] = None
        result["parent_interactions"] = prepared_profile
        if isinstance(source_parent_profile, dict) and isinstance(source_parent_profile.get("interactions"), list):
            result["parent_before_after"] = before_after_report(
                source_parent_profile, prepared_profile
            )
            result["parent_before_after"]["hbond_loss_claimed"] = False
            result["parent_before_after"]["review_reason"] = (
                "Protein hydrogens are unavailable; differences are descriptive "
                "contact-profile changes, not negative hydrogen-bond findings."
            )
    except Exception as exc:
        result["ligand_hydrogen_preparation"] = {"status": "unknown", "error": str(exc), "docking_input_used": False}
        result["parent_interactions"] = {"status": "unknown", "error": str(exc)}
    result["status"] = "computed_with_visible_unknowns"
    return result


def _source_parent_profile(parent, protein):
    from packages.science.chemical_states import interaction_profile
    try:
        profile = interaction_profile(parent, protein, None)
        profile["profile_role"] = "source_parent_before_chemical_preparation"
        profile["protein_hydrogen_status"] = "missing_requires_directional_review"
        profile["chemical_geometry_pass"] = None
        return profile
    except Exception as exc:
        return {
            "status": "unknown",
            "profile_role": "source_parent_before_chemical_preparation",
            "protein_hydrogen_status": "missing_requires_directional_review",
            "chemical_geometry_pass": None,
            "error": str(exc),
        }


def _pose_directional_profiles(folder, diagnostics, protein, parent_profile):
    """Attach an actual coordinate profile to every pose in poses.sdf."""
    from rdkit import Chem
    from packages.science.chemical_states import (
        before_after_report, interaction_profile, optimize_ligand_hydrogens,
    )

    path = Path(folder) / "poses.sdf"
    if not path.is_file():
        for diagnostic in diagnostics:
            diagnostic["directional_contacts"] = {
                "status": "unknown",
                "error": "DOCKING_POSES_SDF_MISSING",
                "chemical_geometry_pass": None,
            }
        return

    try:
        with path.open("rb") as stream:
            poses = list(Chem.ForwardSDMolSupplier(stream, removeHs=False))
    except Exception as exc:
        for diagnostic in diagnostics:
            diagnostic["directional_contacts"] = {
                "status": "unknown",
                "error": str(exc),
                "chemical_geometry_pass": None,
            }
        return

    while len(diagnostics) < len(poses):
        diagnostics.append({
            "mode": len(diagnostics) + 1,
            "status": "review",
            "reason": "POSE_DIAGNOSTIC_MISSING_FROM_DOCKING_RESULT",
        })

    for index, pose in enumerate(poses):
        diagnostic = diagnostics[index]
        if pose is None:
            diagnostic["directional_contacts"] = {
                "status": "unknown",
                "error": "POSE_SDF_RECORD_INVALID",
                "chemical_geometry_pass": None,
            }
            continue
        try:
            structures, receipt = optimize_ligand_hydrogens(pose)
            profile = interaction_profile(
                structures["hydrogens_optimized"], protein, None
            )
            delta = (
                before_after_report(parent_profile, profile)
                if isinstance(parent_profile, dict)
                and isinstance(parent_profile.get("interactions"), list)
                else {
                    "status": "unknown",
                    "error": "SOURCE_PARENT_PROFILE_UNAVAILABLE",
                }
            )
            if isinstance(delta, dict):
                delta["hbond_loss_claimed"] = False
                delta["review_reason"] = (
                    "Protein hydrogens are unavailable; lost/gained identifiers "
                    "are coordinate annotations and not chemical-pass findings."
                )
            diagnostic["directional_contacts"] = {
                "status": "computed_from_actual_pose_sdf",
                "pose_record_index": index,
                "ligand_hydrogen_preparation": receipt,
                "heavy_coordinates_exactly_preserved": receipt[
                    "heavy_coordinates_exactly_preserved"
                ],
                "protein_hydrogen_status": "missing_requires_directional_review",
                "profile": profile,
                "parent_before_after": delta,
                "chemical_geometry_pass": None,
                "interpretation": (
                    "Directional geometry is diagnostic only and does not turn a "
                    "pose-preservation pass into a chemical or efficacy pass."
                ),
            }
        except Exception as exc:
            diagnostic["directional_contacts"] = {
                "status": "unknown",
                "error": str(exc),
                "protein_hydrogen_status": "missing_requires_directional_review",
                "chemical_geometry_pass": None,
            }


def _catalog_card_count(source):
    if source.get("status") != "configured_hash_verified":
        return None
    data = source.get("data")
    if isinstance(data, list):
        return len(data)
    if isinstance(data, dict):
        parent_catalog = data.get("parent_catalog")
        if isinstance(parent_catalog, dict) and isinstance(parent_catalog.get("records"), list):
            return len(parent_catalog["records"])
        if isinstance(data.get("records"), list):
            return len(data["records"])
        if isinstance(data.get("cards"), list):
            return len(data["cards"])
    return None


def _broad_family(value):
    return "ring_modification" if value in {"ring_expansion", "ring_contraction"} else value


def _calibration(evidence, reference_packet):
    source = evidence["sources"]["crbn_calibration.json"]
    data = source["data"]
    if source["status"] != "configured_hash_verified":
        crbn = {"status": "not_configured", "seed_receipts": [], "metrics": None, "reported_reproduction_status": None}
    else:
        crbn = {
            "status": "measured_source_exposed_without_cross_e3_ranking",
            "seed_receipts": data.get("seed_receipts", data.get("receipts", [])) if isinstance(data, dict) else [],
            "metrics": data.get("metrics") if isinstance(data, dict) else None,
            "reported_reproduction_status": data.get("reproduction_status") if isinstance(data, dict) else None,
            "source_record": data,
        }
    return {
        "VHL": {"status": "existing_C01_C02_saved_benchmark_separate_from_new_panel"},
        "CRBN": crbn,
        "new_MSA_baseline": data.get("new_MSA_baseline", data.get("new_msa_baseline")) if isinstance(data, dict) else None,
        "cross_e3_raw_score_ranking": False,
        "reference": reference_packet(),
    }


def run_panel(parameters, put, check_active=lambda: None, provider=None, progress=lambda stage: None,
              scientific_policy=None):
    from rdkit import Chem
    from packages.science.structures import atom_sites
    from packages.science.warhead_sites import analyze_sites, funnel
    from packages.science.analog_generation import generate, select_diverse
    from packages.science.analog_filters import cheap_filter, count_by_site_and_broad_family
    from packages.science.design_docking import dock
    from packages.science.molecules import identity
    from packages.science.e3_benchmark import reference_packet
    from packages.science.linker_assessment import assess_linker, route_evidence

    parameters = validate(parameters)
    binding = input_binding()
    c = catalog()
    meta = c["warhead"]
    evidence = c["acceptance_evidence"]
    check_active()
    if parameters["use_api"] and provider is None:
        raise ValueError("DESIGN_API_NOT_CONFIGURED")

    stage_order = []
    stage_counts = {}
    def stage(name, **counts):
        stage_order.append(name)
        stage_counts[name] = counts
        progress(name)

    selected_parent_id = parameters.get("parent_id", "SMARCA2-FX5")
    if selected_parent_id == "SMARCA2-FX5":
        # Deliberately retain the original neutral SDF and current curated SAR.
        parent = Chem.MolFromMolBlock((SOURCE / "SMARCA2-neutral-design.sdf").read_text())
        meta = c["warhead"]
        parent_input_role = "original_neutral_default"
    else:
        from packages.science.dual_e3 import REFERENCE_SOURCE
        from packages.science.reference_parents import load_reference_parent
        parent, meta = load_reference_parent(selected_parent_id, REFERENCE_SOURCE)
        parent_input_role = "hash_verified_literature_pose_aligned_to_6HAZ"
    if parent is None or parent.GetNumConformers() != 1:
        raise ValueError("SELECTED_PARENT_STRUCTURE_INVALID")
    protein = atom_sites(SOURCE / "6HAZ.cif", ["A"])
    source_parent_profile = _source_parent_profile(parent, protein)
    stage("sites")
    sites = analyze_sites(parent, protein, meta["sar"], meta["protected_maps"])
    scientific_policy_binding=None
    policy_parent_funnel=None
    if scientific_policy is not None:
        from packages.science.dual_e3 import validate_scientific_policy,apply_scientific_policy
        scientific_policy=validate_scientific_policy(scientific_policy,parent_id=selected_parent_id)
        sites["atoms"]=apply_scientific_policy(sites["atoms"],scientific_policy,meta["protected_maps"])
        scientific_policy_binding={"id":scientific_policy["scope"]["assessment_id"],
            "revision":scientific_policy["scope"]["policy_revision"],
            "digest":scientific_policy["scope"]["policy_digest"]}
        policy_parent_funnel=copy.deepcopy(scientific_policy["parent_funnel"])
    protected = [r["atom_map"] for r in sites["atoms"] if r["state"] == "PROTECTED"]

    stage("generate")
    generated = generate(parent, sites["atoms"], parameters["exploratory"])
    stage_counts["generate"] = {"raw": len(generated["analogs"]), "generation_rejections": len(generated["rejections"])}

    stage("cheap_filter")
    cheap_pass = []
    cheap_rejections = []
    for source_row in generated["analogs"]:
        row = copy.deepcopy(source_row)
        row["cheap_filter"] = cheap_filter(parent, row)
        if row["cheap_filter"]["valid"]:
            cheap_pass.append(row)
        else:
            cheap_rejections.append({"stage": "cheap_filter", "record": row, "reasons": row["cheap_filter"]["hard_reasons"]})
    stage_counts["cheap_filter"] = {"input": len(generated["analogs"]), "passed": len(cheap_pass), "rejected": len(cheap_rejections)}

    stage("stereoisomer_expansion")
    variants = _stereoisomers(cheap_pass)
    stage_counts["stereoisomer_expansion"] = {"input": len(cheap_pass), "output": len(variants)}

    stage("eligibility")
    eligible = []
    eligibility_rejections = []
    for row in variants:
        if row["stereochemistry"]["unresolved_stereo_count"]:
            row["pipeline_status"] = "rejected_unresolved_stereochemistry"
            eligibility_rejections.append({"id": row["id"], "reason": "UNRESOLVED_STEREOCHEMISTRY", "record": copy.deepcopy(row)})
            continue
        handles = attachment_options(row)
        row["attachment_options"] = copy.deepcopy(handles)
        row["attachment_feasibility"] = "direct_graph_handle" if handles else "derivatization_required"
        if not handles:
            row["pipeline_status"] = "rejected_no_supported_attachment"
            eligibility_rejections.append({"id": row["id"], "reason": "NO_SUPPORTED_NH_OH_HANDLE", "record": copy.deepcopy(row)})
            continue
        row["pipeline_status"] = "eligible_for_docking" if parameters["dock"] else "preview_only_not_docked"
        eligible.append(row)
    stage_counts["eligibility"] = {"eligible": len(eligible), "rejected": len(eligibility_rejections)}

    def dock_one(molecule, name, folder):
        check_active()
        molecule = Chem.Mol(molecule)
        molecule.RemoveAllConformers()
        docking_options = {"check_active": check_active}
        if selected_parent_id != "SMARCA2-FX5":
            docking_options["selected_parent_id"] = selected_parent_id
        result = dock(
            molecule,
            parent,
            protein,
            protected,
            folder,
            **docking_options,
        )
        diagnostics = result.get("results", {}).get(
            "pose_preservation", {}
        ).get("all_poses_diagnostics", [])
        _pose_directional_profiles(folder, diagnostics, protein, source_parent_profile)
        result["files"] = {
            key: put(name + "-" + filename, (folder / filename).read_bytes(), "chemical/x-mdl-sdfile" if filename.endswith(".sdf") else "text/plain")
            for key, filename in result.get("files", {}).items() if (folder / filename).is_file()
        }
        put(name + "-docking.json", encoded(result))
        return result

    stage("dock_all_eligible")
    with tempfile.TemporaryDirectory(prefix="tpd-dock-") as work_text:
        work = Path(work_text)
        if parameters["dock"]:
            for ordinal, row in enumerate(eligible, 1):
                check_active()
                progress("dock_all_eligible " + str(ordinal) + "/" + str(len(eligible)))
                mol = Chem.MolFromSmiles(row["mapped_smiles"])
                full = dock_one(mol, row["id"], work / row["id"])
                preservation = full.get("results", {}).get("pose_preservation", {})
                diagnostics = preservation.get("all_poses_diagnostics", [])
                row["docking"] = {
                    "status": full["status"],
                    "pose_preserved": preservation.get("docking_pose_preserved", "review"),
                    "best_core_rmsd_A": min((d.get("core_RMSD_A_in_receptor_frame") for d in diagnostics if d.get("core_RMSD_A_in_receptor_frame") is not None), default=None),
                    "pose_count": len(diagnostics),
                    "passing_pose_count": sum(d.get("status") == "pass" for d in diagnostics),
                    "error": full.get("results", {}).get("error"),
                    "source_ref": put(row["id"] + "-pose-evaluation.json", encoded(full)),
                    "files": full["files"],
                }
            stage_counts["dock_all_eligible"] = {"attempted": len(eligible), "completed": sum(r["docking"]["status"] == "completed_with_limits" for r in eligible)}
        else:
            for row in eligible:
                row["docking"] = {"status": "not_run", "pose_preserved": "review", "pose_count": 0, "passing_pose_count": 0}
            stage_counts["dock_all_eligible"] = {"attempted": 0, "qualified": 0, "preview_selection": True}

        stage("parent_redock_and_pose_filter")
        if parameters["dock"]:
            parent_2d = Chem.Mol(parent)
            parent_2d.RemoveAllConformers()
            reference_docking = dock_one(parent_2d, "parent-redocking", work / "parent")
            reference_docking["parent_redocking"] = True
        else:
            reference_docking = {"status": "not_run", "real_docking": False, "parent_redocking": True}

    parent_pass = reference_docking.get("results", {}).get("pose_preservation", {}).get("docking_pose_preserved") is True
    qualified = []
    pose_failures = []
    for row in eligible:
        row["parent_redocking_supported"] = parent_pass
        row["docking_pose_preserved"] = row["docking"]["pose_preserved"]
        row["assembly_eligible"] = bool(parameters["dock"] and parent_pass and row["docking_pose_preserved"] is True and row["docking"]["passing_pose_count"] > 0)
        if row["assembly_eligible"]:
            row["pipeline_status"] = "qualified"
            qualified.append(row)
        else:
            row["pipeline_status"] = "preview_not_qualified" if not parameters["dock"] else "failed_pose_filter"
            pose_failures.append({"id": row["id"], "reason": "DOCKING_DISABLED" if not parameters["dock"] else "POSE_OR_PARENT_REDOCK_FAILED", "docking": copy.deepcopy(row["docking"])})
    stage_counts["parent_redock_and_pose_filter"] = {"parent_pass": parent_pass, "qualified": len(qualified), "failed_or_preview": len(pose_failures)}

    stage("diverse_final_selection")
    selection_pool = qualified if parameters["dock"] else eligible
    selected = select_diverse(selection_pool, parameters["panel_size"])
    selected_ids = {row["id"] for row in selected}
    selection_fields = {
        row["id"]: {
            key: copy.deepcopy(row[key])
            for key in ("cluster_id", "cluster_size")
            if key in row
        }
        for row in selected
    }
    for row in variants:
        row["selected"] = row["id"] in selected_ids
        row.update(selection_fields.get(row["id"], {}))
        if "docking" not in row:
            row["docking"] = {"status": "not_eligible", "pose_preserved": "review", "pose_count": 0, "passing_pose_count": 0}
    selection_mode = "qualified_final" if parameters["dock"] else "explicit_preview_not_qualified"
    qualified_count = len(qualified) if parameters["dock"] else 0
    shortfall = max(0, parameters["panel_size"] - len(selected))
    stage_counts["diverse_final_selection"] = {"mode": selection_mode, "selected": len(selected), "qualified": qualified_count, "shortfall": shortfall}

    stage("assemble_qualified")
    linker_cache = {}
    for linker in c["linkers"]["templates"]:
        if linker["id"] in parameters["linker_ids"]:
            linker_cache[linker["id"]] = assess_linker(linker)
    route_source = evidence["sources"]["route_evidence.json"]
    route_records = route_source["data"] if route_source["status"] == "configured_hash_verified" else []
    if isinstance(route_records, dict):
        route_records = route_records.get("records", route_records.get("evidence_records", []))
    if not isinstance(route_records, list):
        route_records = []

    candidates = []
    assembly_rejections = []
    seen = set()
    if parameters["dock"]:
        for row in selected:
            for recruiter in c["recruiters"]["recruiters"]:
                for linker in c["linkers"]["templates"]:
                    if linker["id"] not in parameters["linker_ids"]:
                        continue
                    check_active()
                    try:
                        candidate = assemble(copy.deepcopy(row), copy.deepcopy(recruiter), copy.deepcopy(linker))
                    except ValueError as exc:
                        assembly_rejections.append({"analog_id": row["id"], "e3_type": recruiter["e3_type"], "linker_id": linker["id"], "reason": str(exc)})
                        continue
                    if candidate["candidate_id"] in seen:
                        continue
                    seen.add(candidate["candidate_id"])
                    candidate["assembly_mode"] = "pose_supported_hypothesis"
                    candidate["linker_assessment"] = copy.deepcopy(linker_cache[linker["id"]])
                    candidate["route_evidence"] = route_evidence(route_records, candidate["candidate_id"], {"novel_transformations": [row["rule_id"]]})
                    candidate["route_evidence"]["C01_C02_policy"] = "reusable precedent only; never exact validation for this novel graph"
                    candidate["files"] = _files(Chem.MolFromSmiles(candidate["mapped_smiles"]), candidate["candidate_id"], put, candidate["atom_roles"]["linker_maps"])
                    candidates.append(candidate)
    stage_counts["assemble_qualified"] = {
        "assembled": len(candidates),
        "rejected": len(assembly_rejections),
        "linker_assessments_cached": len(linker_cache),
        "status": (
            "not_run_docking_disabled_explicit_no_assembly"
            if not parameters["dock"] else "completed"
        ),
    }

    for row in variants:
        mol = Chem.MolFromSmiles(row["mapped_smiles"])
        row["properties"] = identity(mol)["properties"]
        row["files"] = _files(mol, row["id"], put, row.get("modified_atom_maps", []) + row.get("added_atom_maps", []))

    coverage = count_by_site_and_broad_family(cheap_pass, sites["atoms"])
    catalog_families = sorted({_broad_family(r["transformation_class"]) for r in generated["rule_catalog"]})
    actual_families = sorted({_broad_family(r["transformation_class"]) for r in cheap_pass})
    counts = {
        "raw_graph_analogs": len(generated["analogs"]),
        "cheap_filter_passed": len(cheap_pass),
        "valid_analogs": len(variants),
        "docking_eligible": len(eligible),
        "docking_attempted": len(eligible) if parameters["dock"] else 0,
        "qualified_analogs": qualified_count,
        "selected_analogs": len(selected),
        "panel_size_requested": parameters["panel_size"],
        "shortfall": shortfall,
        "selection_mode": selection_mode,
        "protac_count": len(candidates),
        "branch_counts": {e3: sum(r["e3_type"] == e3 for r in candidates) for e3 in ("CRBN", "VHL")},
        "broad_family_coverage": {"catalog": catalog_families, "actual": actual_families},
        "per_site_distinct_graph_counts": coverage,
        "stage_order": stage_order,
        "stage_counts": stage_counts,
        "modifiable_sites": sum(s["state"] == "MODIFIABLE" for s in sites["atoms"]),
        "unknown_sites": sum(s["state"] == "UNKNOWN" for s in sites["atoms"]),
    }

    scientific = _scientific_context(parent, protein, source_parent_profile)
    preparation_source = evidence["sources"]["preparation.json"]
    scientific["optional_reference_preparation"] = {
        "status": preparation_source["status"],
        "source_record": preparation_source["data"],
        "applied_to_design_ligand": False,
        "policy": "reference-specific preparation is not transferred to a different ligand",
    }
    warhead_source = evidence["sources"]["warhead_catalog.json"]
    medchem_source = evidence["sources"]["medchem_evidence.json"]
    from packages.science.medchem_sar import link_selected_parent
    selected_parent_measured = link_selected_parent(parent, selected_parent_id, medchem_source)
    collected_parent_count = _catalog_card_count(warhead_source)
    medchem_parent_count = _catalog_card_count(medchem_source)

    ledger = {
        "functional_implementation": {"status": "implemented", "pipeline_order": stage_order, "version": VERSION},
        "measured_data": {"status": evidence["status"], "source_version": evidence["manifest"]["version"], "sources": {k: v["status"] for k, v in evidence["sources"].items()}},
        "scientific_approval": {"status": "not_approved", "reviewer": None, "approved": False, "cancel_behavior": "workbench cancellation and source-version checks remain authoritative"},
    }
    questions = [
        "현재 원자별 SAR와 UNKNOWN 위치 중 후속 합성 검토할 구조를 골라 주세요.",
        "도킹 통과가 활성 유지 증거가 아니라는 범위와 부모 재도킹 결과를 검토해 주세요.",
        "단백질 수소가 없는 방향성 상호작용은 REVIEW로 유지했습니다. 실제 protonation/tautomer 상태를 검토해 주세요.",
        "신규 그래프의 실제 시약·반응 경로·보호기·SI 선례를 검토해 주세요.",
    ]
    limitations = [
        f"이번 실행의 실제 설계 범위는 선택한 {selected_parent_id} 단일 parent입니다.",
        "문헌 parent의 ready 표시는 hash 검증과 6HAZ frame 정렬이 끝난 기술 입력이라는 뜻이며 과학적 승인이 아닙니다.",
        "고정 6HAZ receptor frame은 parent 간 pose 비교와 설계 가설에만 사용하며 결합 모드 동등성을 주장하지 않습니다.",
        "선택 parent에 명시적으로 연결된 SAR만 사용하며 다른 문헌 parent의 SAR를 자동 전이하지 않습니다.",
        "보호 코어 검사와 도킹은 결합·분해 효능·선택성 또는 합성 성공을 검증하지 않습니다.",
        "CRBN과 VHL raw score를 합치거나 교차 순위화하지 않습니다.",
    ]
    result = {
        "format": VERSION,
        "candidate_id": "PANEL-" + hashlib.sha256(encoded({"parameters": parameters, "binding": binding})).hexdigest()[:12],
        "status": "hypotheses_pending_review",
        "parameters": parameters,
        "input_binding": binding,
        "source_catalog": evidence,
        "parent_scope": {
            "collected_parent_catalog_count": collected_parent_count,
            "medchem_parent_catalog_count": medchem_parent_count,
            "actual_design_parent_count": 1,
            "actual_design_parent_id": selected_parent_id,
            "selected_parent_title": selected_parent_id,
            "selected_parent_input_role": parent_input_role,
            "receptor_frame": "6HAZ chain A",
            "fixed_receptor_alignment_hypothesis": "The selected aligned ligand is compared and docked in the same fixed 6HAZ chain-A receptor frame; this is a pose-comparison hypothesis, not scientific approval.",
            "parent_redock_required_for_qualification": True,
            "technical_ready_is_scientific_approval": False,
            "medchem_measurement_join": selected_parent_measured["status"],
        },
        "selected_parent_measured_evidence": selected_parent_measured,
        "warhead": meta,
        "sites": sites,
        "funnel": policy_parent_funnel if policy_parent_funnel is not None else funnel([dict(meta, state="registered_reference", SAR_status="one_attachment_precedent", exit_vector_confidence="medium")]),
        "scientific_policy_binding": scientific_policy_binding,
        "analogs": variants,
        "rejections": generated["rejections"] + cheap_rejections + eligibility_rejections + pose_failures,
        "rule_catalog": generated["rule_catalog"],
        "protac_candidates": candidates,
        "assembly_rejections": assembly_rejections,
        "summary": counts,
        "parent_redocking": reference_docking,
        "chemical_states": scientific,
        "calibration": _calibration(evidence, reference_packet),
        "acceptance_ledger": ledger,
        "review_questions": questions,
        "limitations": limitations,
        "authority": AUTHORITY.copy(),
        "ai": {"status": "not_requested", "notes": []},
        "files": {},
    }
    if parameters["use_api"]:
        put("design-panel-before-api.json", encoded(result))
        check_active()
        result["ai"] = _ai_review(provider, {"summary": counts, "calibration": result["calibration"], "ledger": ledger, "limitations": limitations}, put)

    for key, records, name in [("sdf", candidates, "PROTAC-panel.sdf"), ("warheads", variants, "warhead-panel.sdf")]:
        text = []
        for row in records:
            mol = Chem.MolFromSmiles(row["mapped_smiles"])
            mol.SetProp("_Name", row.get("candidate_id", row.get("id")))
            from rdkit.Chem import rdDepictor
            rdDepictor.Compute2DCoords(mol)
            text.append(Chem.MolToMolBlock(mol, forceV3000=True) + "\n$$$$\n")
        result["files"][key] = put(name, "".join(text).encode(), "chemical/x-mdl-sdfile")
    import jsonschema
    schema = json.loads((SOURCE.parents[1] / "contracts/drafts/design_panel.schema.json").read_text())
    jsonschema.Draft202012Validator(schema).validate(result)
    result["files"]["json"] = put("design-panel.json", encoded(result))
    result["files"]["report"] = put("설계결과_검수요청.md", report_markdown(result).encode(), "text/markdown")
    return result


def report_markdown(r):
    s = r["summary"]
    scope = r["parent_scope"]
    lines = [
        "# Warhead SAR·CRBN/VHL 통합 설계 결과", "",
        "상태: 실제 계산으로 생성된 설계 가설이며 사람 승인·효능·합성 검증은 없습니다.", "",
        f"- 수집 parent catalog: {scope['collected_parent_catalog_count']} / 실제 설계 범위: {scope['actual_design_parent_count']}개 parent",
        f"- 선택 parent: {scope.get('selected_parent_title', scope.get('actual_design_parent_id'))} / 입력 역할: {scope.get('selected_parent_input_role')}",
        f"- receptor frame: {scope.get('receptor_frame')} / parent redock 필수: {scope.get('parent_redock_required_for_qualification')}",
        f"- medchem.parent_catalog: {scope.get('medchem_parent_catalog_count')}",
        f"- 선택 parent 측정 join: {r.get('selected_parent_measured_evidence', {}).get('status', '미보고')}; 신규 analog 측정 증명: {r.get('selected_parent_measured_evidence', {}).get('new_analog_measurement_status', '없음')}",
        f"- 생성 graph {s['raw_graph_analogs']} / cheap-filter 통과 {s['cheap_filter_passed']} / 입체 상태 {s['valid_analogs']}",
        f"- 도킹 적격 {s['docking_eligible']} / 실제 qualified {s['qualified_analogs']} / 선택 {s['selected_analogs']} / shortfall {s['shortfall']}",
        f"- 선택 모드: {s['selection_mode']}; CRBN 조립 {s['branch_counts']['CRBN']}, VHL 조립 {s['branch_counts']['VHL']}",
        f"- 실제 broad family: {', '.join(s['broad_family_coverage']['actual']) or '없음'}", "",
        "## Linker·route·수소·CRBN 상태", "",
        f"- linker 평가 cache: {s['stage_counts']['assemble_qualified']['linker_assessments_cached']}개; descriptor와 sampled 3D는 후보별 JSON에 기록했습니다.",
        "- ternary geometry 입력이 없으면 미지이며 exclusion으로 사용하지 않았습니다.",
        "- C01/C02 route는 재사용 선례로만 취급하고 신규 graph의 정확한 합성 검증으로 사용하지 않았습니다.",
        f"- ligand 수소 준비: {r['chemical_states']['ligand_hydrogen_preparation']['status']}; 이 구조는 docking metric 입력으로 사용하지 않았습니다.",
        f"- CRBN calibration: {r['calibration']['CRBN']['status']}; seed receipt와 metric은 source에 존재하는 값만 노출합니다.",
        "- new MSA baseline은 기존 E3 raw score와 별도로 기록하며 교차 순위화하지 않습니다.", "",
        "## 대표 Warhead", "", "| ID | 변형 | pipeline | 도킹 | 통과 pose |", "|---|---|---|---|---|",
    ]
    for analog in r["analogs"]:
        if analog.get("selected"):
            docking = analog["docking"]
            lines.append(f"| {analog['id']} | {analog['transformation_class']} | {analog.get('pipeline_status')} | {docking['status']} | {docking.get('passing_pose_count', 0)} |")
    lines += ["", "## 확인이 필요한 부분", ""] + ["- " + q for q in r["review_questions"]]
    lines += ["", "## 해석 범위", ""] + ["- " + value for value in r["limitations"]]
    lines += ["", "## API 검토", "", r["ai"]["status"]] + ["- " + value for value in r["ai"].get("notes", [])]
    return "\n".join(lines) + "\n"
