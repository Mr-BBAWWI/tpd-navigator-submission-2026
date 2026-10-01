"""Native 9D12 parent expansion v7.

This module runs a predeclared 20-graph exploratory N3 panel and a separate
installed-rule applicability audit. It does not promote UNKNOWN or PROTECTED
sites, alter strict policy, or make potency, binding, approval, or route claims.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
import os
import shutil
from pathlib import Path
from typing import Any, Callable

from rdkit import Chem

from packages.science.analog_filters import cheap_filter
from packages.science.analog_generation import RULE_CATALOG, generate
from packages.science.native_attachment_experiment import (
    EXHAUSTIVENESS,
    FRAGMENTS as ORIGINAL_FRAGMENTS,
    PARENT_ID,
    SEED,
    _assemble_branches,
    _graph_hash,
    _manifest_files,
    _original_maps,
    _pose_metrics,
    _remote_exclusions,
    _source_hash,
    _write_sdf,
    build_n3_candidate,
    find_source_linkage,
    validate_parent_baseline,
    verify_native_probe,
)
from packages.science.native_parent_docking import (
    dock_native,
    load_native_context,
    verify_source_pair,
)

NEW_FRAGMENTS = (
    ("N3-V7-13", "COCCN"),
    ("N3-V7-14", "CCNCCO"),
    ("N3-V7-15", "CCNCCN"),
    ("N3-V7-16", "CCCOCCO"),
    ("N3-V7-17", "CCCOCCN"),
    ("N3-V7-18", "CCOCCCO"),
    ("N3-V7-19", "CCOCCCN"),
    ("N3-V7-20", "CC(=O)NCCO"),
)
FRAGMENTS = tuple(ORIGINAL_FRAGMENTS) + NEW_FRAGMENTS
PANEL_VERSION = "native-9D12-N3-parent-expansion-v7"


def _finite(value: Any) -> Any:
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("Non-finite JSON number is forbidden")
        return value
    if isinstance(value, dict):
        return {str(key): _finite(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_finite(item) for item in value]
    if hasattr(value, "item"):
        return _finite(value.item())
    return value


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            _finite(value), indent=2, sort_keys=True, ensure_ascii=False,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("xb") as stream:
            stream.write(_json_bytes(value))
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(value)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _append_case_log(path: Path, event: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(_finite(event), sort_keys=True, ensure_ascii=False, allow_nan=False)
    with path.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(line + "\n")
        stream.flush()
        os.fsync(stream.fileno())


def _map_free_graph(mol: Chem.Mol, *, isomeric: bool) -> str:
    copied = Chem.Mol(mol)
    for atom in copied.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(copied, canonical=True, isomericSmiles=isomeric)


def enumerate_expansion_candidates(parent: Chem.Mol) -> list[dict[str, Any]]:
    """Build and constitutionally deduplicate all 20 graphs before any docking."""
    records = [
        build_n3_candidate(parent, fragment_id, fragment_smiles)
        for fragment_id, fragment_smiles in FRAGMENTS
    ]
    constitutional: set[str] = set()
    candidate_ids: set[str] = set()
    for record in records:
        mol = Chem.MolFromSmiles(record["mapped_smiles"])
        if mol is None:
            raise ValueError(f"Candidate graph cannot be parsed: {record['id']}")
        graph = _map_free_graph(mol, isomeric=False)
        if graph in constitutional:
            raise ValueError(f"Duplicate constitutional graph in predeclared panel: {record['id']}")
        constitutional.add(graph)
        if record["candidate_id"] in candidate_ids:
            raise ValueError("Duplicate candidate identity in predeclared panel")
        candidate_ids.add(record["candidate_id"])
        record["panel_version"] = "v7"
        record["broad_family"] = "linker_handle_introduction"
        record["fragment_motif_counts_as_new_family"] = False
    if len(records) != 20 or len(constitutional) != 20:
        raise ValueError("The v7 panel must contain exactly 20 unique constitutional graphs")
    return records


def _merged_family(value: Any) -> str:
    family = str(value or "unknown")
    if family in {"ring_expansion", "ring_contraction"}:
        return "ring_modification"
    return family


def _reason_category(reason: Any) -> str:
    text = str(reason or "")
    if text == "NO_APPLICABLE_SITE":
        return "chemical_pattern_missing"
    if text in {"UNKNOWN_REQUIRES_EXPLORATORY"} or "EXACT" in text.upper():
        return "source_missing"
    if text == "PROTECTED_TOUCHED":
        return "protected"
    if "QUOTA" in text or "LIMIT_REACHED" in text:
        return "quota"
    return "generation_or_validation_rejection"


def _canonical_record_graph(record: dict[str, Any]) -> str | None:
    smiles = record.get("mapped_smiles") or record.get("canonical_smiles")
    mol = Chem.MolFromSmiles(str(smiles)) if smiles else None
    return _map_free_graph(mol, isomeric=False) if mol is not None else None


def audit_rule_applicability(
    parent: Chem.Mol,
    sites: list[dict[str, Any]],
    *,
    generate_fn: Callable[..., dict[str, Any]] | None = None,
    cheap_filter_fn: Callable[[Chem.Mol, dict[str, Any]], dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Audit the full installed catalog under exact strict and exploratory states.

    Catalog availability is reported separately from actual emitted, graph-valid,
    cheap-filter-passing, and strict-qualified counts.
    """
    generator = generate if generate_fn is None else generate_fn
    filterer = cheap_filter if cheap_filter_fn is None else cheap_filter_fn
    site_snapshot = copy.deepcopy(sites)
    strict_result = generator(parent, copy.deepcopy(site_snapshot), exploratory=False)
    exploratory_result = generator(parent, copy.deepcopy(site_snapshot), exploratory=True)
    if sites != site_snapshot:
        raise ValueError("Rule audit mutated the original site declarations")

    catalog = {str(row["rule_id"]): copy.deepcopy(row) for row in RULE_CATALOG}
    modes: dict[str, dict[str, Any]] = {}
    for mode, generated in (("strict", strict_result), ("exploratory", exploratory_result)):
        emitted_by_rule: dict[str, list[dict[str, Any]]] = {}
        for record in generated.get("analogs", []):
            if not isinstance(record, dict):
                continue
            row = copy.deepcopy(record)
            row["cheap_filter"] = filterer(parent, row)
            emitted_by_rule.setdefault(str(row.get("rule_id")), []).append(row)
        rejection_by_rule: dict[str, list[dict[str, Any]]] = {}
        for rejection in generated.get("rejections", []):
            if not isinstance(rejection, dict):
                continue
            row = copy.deepcopy(rejection)
            row["category"] = _reason_category(row.get("reason"))
            rejection_by_rule.setdefault(str(row.get("rule_id")), []).append(row)

        rules = []
        for rule_id, catalog_row in catalog.items():
            emitted = emitted_by_rule.get(rule_id, [])
            graphs = {
                graph for graph in (_canonical_record_graph(row) for row in emitted)
                if graph is not None
            }
            cheap_graphs = {
                graph for row in emitted
                if row.get("cheap_filter", {}).get("valid") is True
                for graph in [_canonical_record_graph(row)] if graph is not None
            }
            strict_graphs = {
                graph for row in emitted
                if row.get("qualified_for_counts") is True
                and row.get("cheap_filter", {}).get("valid") is True
                for graph in [_canonical_record_graph(row)] if graph is not None
            }
            rejections = rejection_by_rule.get(rule_id, [])
            rules.append({
                "rule_id": rule_id,
                "scope": catalog_row.get("scope"),
                "transformation_class": catalog_row.get("transformation_class"),
                "broad_family": _merged_family(catalog_row.get("transformation_class")),
                "catalog_available": True,
                "actual_emitted_count": len(emitted),
                "actual_graph_valid_count": len(graphs),
                "actual_cheap_pass_count": len(cheap_graphs),
                "actual_strict_qualified_count": len(strict_graphs),
                "exact_touched_maps": sorted({
                    int(atom_map)
                    for row in emitted
                    for key in (
                        "attachment_site_atom_maps", "modified_atom_maps",
                        "removed_atom_maps", "touched_atom_maps",
                    )
                    for atom_map in (row.get(key, []) or [])
                }),
                "cheap_filter_failures": [{
                    "candidate_id": row.get("candidate_id") or row.get("id"),
                    "reasons": copy.deepcopy(row.get("cheap_filter", {}).get("hard_reasons", [])),
                    "modified_atom_maps": sorted(int(atom_map) for atom_map in (row.get("modified_atom_maps", []) or [])),
                    "diagnostics": copy.deepcopy(row.get("cheap_filter", {})),
                } for row in emitted if row.get("cheap_filter", {}).get("valid") is not True],
                "rejections": [{
                    "reason": row.get("reason"),
                    "category": row["category"],
                    "touched_atom_maps": row.get("touched_atom_maps", []),
                    "protected_atom_maps": row.get("protected_atom_maps", []),
                    "unknown_atom_maps": row.get("unknown_atom_maps", []),
                } for row in rejections],
                "rejection_reason_counts": {
                    category: sum(row["category"] == category for row in rejections)
                    for category in sorted({row["category"] for row in rejections})
                },
            })
        actual_families = sorted({
            row["broad_family"] for row in rules if row["actual_cheap_pass_count"] > 0
        })
        modes[mode] = {
            "rules": rules,
            "actual_graph_valid_count": sum(row["actual_graph_valid_count"] for row in rules),
            "actual_cheap_pass_count": sum(row["actual_cheap_pass_count"] for row in rules),
            "actual_strict_qualified_count": sum(row["actual_strict_qualified_count"] for row in rules),
            "actual_broad_families": actual_families,
            "actual_broad_family_count": len(actual_families),
        }

    catalog_families = sorted({_merged_family(row.get("transformation_class")) for row in catalog.values()})
    return {
        "format_version": "installed-rule-applicability-audit-v1",
        "site_states_preserved_exactly": site_snapshot,
        "catalog_rule_count": len(catalog),
        "catalog_broad_families_available": catalog_families,
        "catalog_availability_does_not_count_as_actual": True,
        "strict": modes["strict"],
        "exploratory": modes["exploratory"],
        "ring_expansion_and_contraction_merged_as": "ring_modification",
        "unknown_or_protected_conversion_performed": False,
    }


def _copy_all_pose_sdf(case_dir: Path) -> str | None:
    source = case_dir / "dock" / "poses.sdf"
    if not source.is_file():
        return None
    destination = case_dir / "all-poses.sdf"
    shutil.copyfile(source, destination)
    return destination.relative_to(case_dir).as_posix()


def _execute_case(
    parent: Chem.Mol,
    protein: list[dict[str, Any]],
    protected_maps: list[int],
    receptor_pdbqt: Path,
    record: dict[str, Any],
    case_dir: Path,
    *,
    audit_only: bool,
    docker: Callable[..., dict[str, Any]] = dock_native,
    assembler: Callable[[dict[str, Any], Path], list[dict[str, Any]]] = _assemble_branches,
) -> dict[str, Any]:
    """Execute one retained case; every post-preflight failure becomes an outcome."""
    case_dir.mkdir(parents=True, exist_ok=False)
    log_path = case_dir / "case.log.jsonl"
    outcome: dict[str, Any] = {
        "id": record["id"],
        "candidate_id": record["candidate_id"],
        "mapped_smiles": record["mapped_smiles"],
        "raw_candidate_sdf": "raw-candidate.sdf",
        "all_pose_sdf": None,
        "raw_dock_receipt": None,
        "assembly_sdfs": [],
        "cheap_filter": record["cheap_filter"],
        "broad_family": "linker_handle_introduction",
        "strict_qualified_for_counts": False,
        "scientific_pending": True,
        "docking_attempted": False,
        "assemblies": [],
        "claims": {"binding": False, "potency": False, "approval": False, "synthetic_route": False},
    }
    try:
        mol = Chem.MolFromSmiles(record["mapped_smiles"])
        if mol is None:
            raise ValueError("Generated candidate cannot be parsed")
        _write_sdf(case_dir / "raw-candidate.sdf", mol)
        _atomic_json(case_dir / "candidate.json", record)
        _append_case_log(log_path, {"event": "candidate_written", "id": record["id"]})
        if record["cheap_filter"].get("valid") is not True:
            outcome["status"] = "filtered_not_replaced"
            outcome["docking"] = None
        elif audit_only:
            outcome["status"] = "audit_only_no_docking"
            outcome["docking"] = None
        else:
            ligand = Chem.Mol(mol)
            ligand.RemoveAllConformers()
            _append_case_log(log_path, {"event": "docking_started", "real_docking_requested": True})
            outcome["docking_attempted"] = True
            docking = docker(
                ligand, parent, protein, protected_maps, receptor_pdbqt,
                case_dir / "dock", seed=23, exhaustiveness=16, timeout=180,
            )
            _atomic_json(case_dir / "raw-dock-receipt.json", docking)
            outcome["raw_dock_receipt"] = "raw-dock-receipt.json"
            outcome["docking"] = docking
            outcome["all_pose_sdf"] = _copy_all_pose_sdf(case_dir)
            if docking.get("status") != "completed_with_limits" or docking.get("real_docking") is not True:
                outcome["status"] = "docking_failed_not_replaced"
            else:
                if outcome["all_pose_sdf"] is None:
                    raise FileNotFoundError("Completed real docking did not produce dock/poses.sdf")
                outcome["pose_metrics"] = _pose_metrics(docking)
                if outcome["pose_metrics"]["pose_preservation_pass_count"] < 1:
                    outcome["status"] = "docked_pose_failed_not_replaced"
                else:
                    try:
                        outcome["assemblies"] = assembler(record, case_dir)
                        outcome["assembly_sdfs"] = [
                            str(row["sdf"]) for row in outcome["assemblies"] if row.get("sdf")
                        ]
                        outcome["status"] = "docked_and_assembled_hypotheses"
                    except Exception as exc:
                        outcome["status"] = "assembly_failed_not_replaced"
                        outcome["assembly_error"] = str(exc)
        _append_case_log(log_path, {"event": "case_finished", "status": outcome["status"]})
    except Exception as exc:
        outcome["status"] = "case_error_not_replaced"
        outcome["error"] = str(exc)
        outcome.setdefault("docking", None)
        _append_case_log(log_path, {"event": "case_error", "error": str(exc)})
    _atomic_json(case_dir / "outcome.json", outcome)
    return outcome


def _actual_counts(cases: list[dict[str, Any]], candidates: list[dict[str, Any]]) -> dict[str, int]:
    attempted = [row for row in cases if row.get("docking_attempted") is True]
    completed = [
        row for row in attempted
        if (row.get("docking") or {}).get("status") == "completed_with_limits"
        and (row.get("docking") or {}).get("real_docking") is True
        and row.get("all_pose_sdf") is not None
    ]
    return {
        "graph_valid": len(candidates),
        "cheap_pass": sum(row["cheap_filter"].get("valid") is True for row in candidates),
        "dock_attempt": len(attempted),
        "dock_completed": len(completed),
        "dock_failed": len(attempted) - len(completed),
        "pose_pass": sum(row.get("pose_metrics", {}).get("pose_preservation_pass_count", 0) > 0 for row in cases),
        "assembled_cases": sum(bool(row.get("assemblies")) for row in cases),
        "assembled_branches": sum(len(row.get("assemblies", [])) for row in cases),
    }


def _report(cases: list[dict[str, Any]], counts: dict[str, int], audit_only: bool) -> str:
    mode = "감사 전용 — 도킹 실행 안 함" if audit_only else "실제 native Vina 도킹 실행"
    return "\n".join([
        "# Native 9D12 N3 parent expansion v7",
        "",
        f"실행 모드: **{mode}**",
        "",
        "## 추가 개발 필요",
        "",
        "- strict broad-family 요구 6개는 충족되지 않았다. catalog availability는 actual family count가 아니다.",
        "- 20개 panel은 UNKNOWN map 3의 exploratory source-extrapolation이며 strict qualified panel 10개로 계산하지 않는다.",
        "- microstate, exact synthesis, calibration, ternary geometry 및 human review 요구사항은 계속 PENDING이다.",
        "- UNKNOWN 또는 PROTECTED site를 자동 승격하지 않았고 기존 threshold를 변경하지 않았다.",
        "",
        "## 실제 exploratory 결과",
        "",
        f"- graph valid: {counts['graph_valid']}",
        f"- cheap-filter pass: {counts['cheap_pass']}",
        f"- docking attempted: {counts['dock_attempt']}",
        f"- docking completed: {counts['dock_completed']}",
        f"- docking failed: {counts['dock_failed']}",
        f"- pose-pass cases: {counts['pose_pass']}",
        f"- assembled cases/branches: {counts['assembled_cases']}/{counts['assembled_branches']}",
        "- strict counts: 0",
        "",
        "All 20 predeclared outcomes, including failures, were retained without replacement.",
        "Docking scores are not affinity measurements; assemblies are unvalidated graph hypotheses.",
        "",
        "## Strict 실행 준비 요건",
        "",
        "M2 evidence가 검증된 뒤 실행할 명령:",
        "`python scripts/run_verified_parent_panel.py --data-dir <M2_STORE> --project <PROJECT> --parent-id SMARCA2-9D12-A1A1P --scientific-policy-id <ASSESSMENT_ID> --panel-size 20 --receipt <NEW_RECEIPT>`",
        "",
        "주의: matching native-science context M2 worker는 아직 구현되지 않았으며, 기존 strictDesign worker는 fixed 6HAZ context를 사용한다.",
        "",
    ])


def run_expansion(
    native_probe: Path,
    native_cif: Path,
    source_export: Path,
    output: Path,
    *,
    audit_only: bool = False,
) -> dict[str, Any]:
    """Run the v7 audit and, unless explicitly disabled, real native docking."""
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")

    native_receipt = verify_native_probe(native_probe)
    baseline = validate_parent_baseline(native_receipt)
    verified, selected, source_object, medchem_ref = verify_source_pair(source_export, PARENT_ID)
    source_before = _source_hash(source_object)
    linkage = find_source_linkage(selected, source_object.get("data"))
    if _source_hash(source_object) != source_before:
        raise ValueError("Verified source object was mutated")

    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir(exist_ok=False)
    exclusions = _remote_exclusions(native_receipt["protocol"])
    parent, protein, context = load_native_context(
        PARENT_ID, native_cif, output / "receptorfiles",
        exclude_remote_incomplete_residues=exclusions,
    )
    baseline_cif_hash = native_receipt["protocol"].get("native_context", {}).get("native_cif_sha256")
    if not baseline_cif_hash or context.get("native_cif_sha256") != baseline_cif_hash:
        raise ValueError("The native CIF SHA256 does not exactly match the verified 9D12 baseline")

    from packages.science.medchem_sar import link_selected_parent
    recomputed = link_selected_parent(parent, PARENT_ID, source_object)
    recomputed_linkage = find_source_linkage(recomputed, source_object.get("data"))
    for key in (
        "pair_id", "requested_change_label", "selected_parent_pair_side",
        "other_compound_id", "source_affected_atom_maps",
        "selected_parent_affected_atom_maps", "source_site_consensus", "source_row_locators",
    ):
        if recomputed_linkage["source_record"].get(key) != linkage["source_record"].get(key):
            raise ValueError("Fresh source-pair join differs from the verified linkage")

    from packages.science.warhead_sites import analyze_sites
    native_sites = analyze_sites(parent, protein, {}, [])
    site_rows = native_sites.get("atoms", [])
    n3_rows = [row for row in site_rows if isinstance(row, dict) and row.get("atom_map") == 3]
    if len(n3_rows) != 1 or str(n3_rows[0].get("state", "UNKNOWN")).upper() != "UNKNOWN":
        raise ValueError("Native parent map 3 must remain exactly UNKNOWN")
    protected_maps = sorted(
        int(row["atom_map"]) for row in site_rows
        if isinstance(row, dict) and row.get("state") == "PROTECTED"
    )
    if len(protected_maps) != 14 or 3 in protected_maps:
        raise ValueError("Native docking requires exactly 14 PROTECTED maps excluding UNKNOWN map 3")

    candidates = enumerate_expansion_candidates(parent)
    rule_audit = audit_rule_applicability(parent, site_rows)
    _atomic_json(output / "rule-audit.json", rule_audit)

    requirements_to_run_strict = {
        "status": "not_satisfied_by_this_exploratory_run",
        "requirements": [
            "M2 source-bound trusted site policy",
            "exact allowed rule IDs tied to source DOI, source locator, and selected-parent atom map",
        ],
        "human_approval_created": False,
        "human_review_status": "pending",
        "automatic_strict_policy_import_performed": False,
        "runnable_after_verified_M2_evidence": "python scripts/run_verified_parent_panel.py --data-dir <M2_STORE> --project <PROJECT> --parent-id SMARCA2-9D12-A1A1P --scientific-policy-id <ASSESSMENT_ID> --panel-size 20 --receipt <NEW_RECEIPT>",
        "cli_caution": "A matching native-science context M2 worker is not yet implemented; the existing strictDesign worker uses a fixed 6HAZ context.",
    }
    protocol = {
        "format_version": PANEL_VERSION,
        "status": "predeclared_before_first_dock",
        "mode": "audit_only_no_docking" if audit_only else "real_native_docking",
        "parent_id": PARENT_ID,
        "requested_count": 20,
        "predeclared_fragments": [{"id": row[0], "smiles": row[1]} for row in FRAGMENTS],
        "constitutional_graph_dedup_completed_before_any_vina": True,
        "failure_replacement_policy": "none",
        "source_export_receipt": verified,
        "source_medchem_evidence": medchem_ref,
        "source_linkage": linkage,
        "source_object_sha256": source_before,
        "native_probe": {key: native_receipt[key] for key in ("manifest_sha256", "result_sha256", "protocol_sha256")},
        "parent_baseline": baseline,
        "native_context": context,
        "native_cif_sha256_fixed_to_verified_baseline": baseline_cif_hash,
        "explicit_remote_exclusions_reused_exactly": exclusions,
        "site_analysis": native_sites,
        "masks": {"native_PROTECTED_14_excluding_UNKNOWN_map_3": protected_maps},
        "docking": {"seed": SEED, "exhaustiveness": EXHAUSTIVENESS, "num_modes": 5, "timeout_seconds": 180},
        "actual_matched_core_thresholds": "unchanged existing pose_preservation implementation",
        "panel_broad_family": "linker_handle_introduction",
        "fragment_motifs_are_not_new_families": True,
        "requirements_to_run_strict": requirements_to_run_strict,
        "claims": {"binding": False, "potency": False, "approval": False, "synthetic_route": False},
        "cases": [{
            "id": row["id"], "candidate_id": row["candidate_id"],
            "fragment_smiles": row["fragment_smiles"],
            "family": "linker_handle_introduction", "seed": 23,
        } for row in candidates],
    }
    _atomic_json(output / "protocol.json", protocol)

    receptor_pdbqt = Path(context["receptor_preparation"]["pdbqt"])
    outcomes: list[dict[str, Any]] = []
    for record in candidates:
        outcome = _execute_case(
            parent, protein, protected_maps, receptor_pdbqt, record,
            output / "cases" / record["id"], audit_only=audit_only,
        )
        outcomes.append(outcome)
        print(f"[{len(outcomes)}/20] {record['id']}: {outcome['status']}", flush=True)
        _atomic_json(output / "panel.incremental.json", {
            "format_version": PANEL_VERSION,
            "requested_count": 20,
            "completed_outcome_count": len(outcomes),
            "replacement_count": 0,
            "cases": outcomes,
        })

    counts = _actual_counts(outcomes, candidates)
    achieved_families = sorted({
        row["broad_family"] for row in candidates
        if row.get("cheap_filter", {}).get("valid") is True
    })
    final = {
        "format_version": PANEL_VERSION,
        "status": "audit_only_completed_no_docking" if audit_only else "completed_with_limits",
        "parent_id": PARENT_ID,
        "requested_count": 20,
        "retained_count": len(outcomes),
        "replacement_count": 0,
        "source_state_map_3": "UNKNOWN",
        "strict_qualified_for_counts": False,
        "strict_counts": {"families": 0, "panel": 0, "total": 0},
        "actual_counts": counts,
        "achieved_exploratory_broad_families": achieved_families,
        "achieved_exploratory_broad_family_count": len(achieved_families),
        "requirements_to_run_strict": requirements_to_run_strict,
        "cases": outcomes,
        "claims": {"binding": False, "potency": False, "approval": False, "synthetic_route": False},
    }
    compact = {
        "format_version": "native-parent-expansion-v7-compact-v1",
        "mode": protocol["mode"],
        "requested": 20,
        "retained": len(outcomes),
        "actual_counts": counts,
        "strict_counts": {"families": 0, "panel": 0, "total": 0},
        "all_claims_false": True,
        "candidate_rows": [{
            "id": row["id"],
            "candidate_id": row["candidate_id"],
            "status": row["status"],
            "global_best_raw_score": row.get("pose_metrics", {}).get("global_best_raw_score"),
            "best_geometry_core_RMSD_A": row.get("pose_metrics", {}).get("best_geometry_core_RMSD_A"),
            "best_geometry_contact_retention": row.get("pose_metrics", {}).get("best_geometry_contact_retention"),
            "assembly_count": len(row.get("assemblies", [])),
        } for row in outcomes],
    }
    _atomic_json(output / "panel.json", final)
    _atomic_json(output / "compact-summary.json", compact)
    _atomic_text(output / "report.md", _report(outcomes, counts, audit_only))
    _atomic_json(output / "manifest.json", {
        "format_version": PANEL_VERSION,
        "hash_algorithm": "sha256",
        "files": _manifest_files(output),
        "manifest_self_hash_policy": "manifest.json is excluded",
    })
    return final


__all__ = [
    "FRAGMENTS",
    "NEW_FRAGMENTS",
    "audit_rule_applicability",
    "enumerate_expansion_candidates",
    "run_expansion",
]
