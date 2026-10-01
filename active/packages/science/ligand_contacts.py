"""Offline execution/archive reader for ligand contacts, with explicit chemical limits."""
import copy
from datetime import datetime, timezone
import math
import os
from pathlib import Path
import subprocess
import time

from .evidence_common import file_index, safe_relative, schema_check, seal, verify_seal, write_new_directory
from .handoff import ACTIVE, check, encoded, sha
from .ligand_preparation import read_preparation
from .quality_evidence import tree_files
from .structure_quality import probe_pairs, read_quality, strict_json

INPUT = "tpd-ligand-contact-input/0.1.0-draft"
REPORT = "tpd-ligand-contact-report/0.1.0-draft"
ARCHIVE = "tpd-ligand-contact-archive/0.1.0-draft"
LIMITS = ["INPUT_PROTONATION_NOT_PREDICTED", "LIGAND_H_ORIENTATION_UNOPTIMIZED",
          "ATOM_TYPING_REQUIRES_REVIEW", "PROTEIN_H_REUSED_FROM_PRIOR_REDUCE2",
          "PROBE_INTERNAL_COORDINATE_ROUNDING", "NO_WHOLE_STRUCTURE_OR_EFFICACY_APPROVAL"]


def archive_files(directory, expected_format):
    directory = Path(directory)
    manifest = strict_json((directory / "manifest.json").read_bytes()); verify_seal(manifest)
    check(manifest["format"] == expected_format, "CONTACT_ARCHIVE_FORMAT")
    files = {}
    check(not any(p.is_symlink() for p in directory.rglob("*")), "CONTACT_ARCHIVE_SYMLINK")
    for name, ref in manifest["files"].items():
        check(name != "manifest.json", "CONTACT_RESERVED_PATH")
        data = safe_relative(directory, name).read_bytes()
        check(sha(data) == ref["sha256"] and len(data) == ref["bytes"], "CONTACT_ARCHIVE_HASH")
        files[name] = data
    check({p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file()} == set(files) | {"manifest.json"}, "CONTACT_UNINDEXED_FILE")
    return manifest, files


def read_contact_input(directory, *, allow_synthetic=False):
    directory = Path(directory)
    manifest, files = archive_files(directory, INPUT)
    context = strict_json(files["context.json"]); verify_seal(context)
    check(context["format"] == INPUT and context["digest"] == manifest["context_digest"], "CONTACT_CONTEXT_DIGEST")
    ligand = read_preparation(directory / "source/ligand", allow_synthetic=allow_synthetic)
    quality = read_quality(directory / "source/quality", allow_synthetic=allow_synthetic)
    check(context["binding"] == ligand["binding"] == quality["binding"], "CONTACT_SAMPLE_BINDING")
    check(context["data_mode"] == ligand["data_mode"] == quality["data_mode"], "CONTACT_DATA_MODE")
    check(context["ligand_preparation_digest"] == ligand["digest"] and context["prior_quality_digest"] == quality["digest"], "CONTACT_SOURCE_DIGEST")
    check(sha(files["raw/restraints.cif"]) == context["dictionary_sha256"], "CONTACT_DICTIONARY_HASH")
    check(sha(files["raw/protein-normalization.json"]) == context["protein_normalization_sha256"], "CONTACT_NORMALIZATION_HASH")
    check(set(context["models"]) == {"xray", "nuclear"}, "CONTACT_CONVENTIONS")
    for model in context["models"].values():
        check(sha(files[model["path"]]) == model["sha256"], "CONTACT_MODEL_HASH")
    return context


def analyze_atom_audit(context, mode, audit):
    model = context["models"][mode]
    check(audit["format"] == "tpd-probe-atom-audit/0.1.0-draft" and audit["cctbx_version"] == "2025.11", "CONTACT_AUDIT_VERSION")
    check(audit["input_sha256"] == model["sha256"] and audit["hydrogen_convention"] == mode, "CONTACT_AUDIT_BINDING")
    expected = {tuple(a["atom"]): a for a in model["atoms"]}
    actual = {tuple(a["atom"]): a for a in audit["atoms"]}
    check(len(expected) == len(model["atoms"]) and len(actual) == len(audit["atoms"]) and set(actual) == set(expected), "CONTACT_AUDIT_ATOM_SET")
    coordinate_deltas = []
    for key, atom in actual.items():
        check(atom["element"] == expected[key]["element"] and len(atom["xyz_A"]) == 3
              and all(math.isfinite(v) for v in atom["xyz_A"]), "CONTACT_PROBE_CHANGED_MODEL")
        check(all(abs(v-w) < 1e-8 or (abs(v-round(v, 3)) < 1e-8 and abs(v-w) <= 0.00050001)
                  for v, w in zip(atom["xyz_A"], expected[key]["xyz_A"])), "CONTACT_PROBE_CHANGED_MODEL")
        coordinate_deltas.append(math.dist(atom["xyz_A"], expected[key]["xyz_A"]))
        check(type(atom["probe_acceptor"]) is bool and type(atom["probe_donor"]) is bool, "CONTACT_AUDIT_ROLE_TYPE")
    prefix = context["ligand_prefix"]
    ligand = {k[4]: a for k, a in actual.items() if list(k[:4]) == prefix}
    check(set(ligand) == set(context["dictionary_atoms"]), "CONTACT_AUDIT_LIGAND_SET")
    neighbors = {n: set() for n in ligand}
    for bond in context["dictionary_bonds"]:
        a, b = bond["atoms"]; neighbors[a].add(b); neighbors[b].add(a)
    for name, atom in ligand.items():
        check(all(n[:4] == prefix for n in atom["neighbors"]), "CONTACT_UNEXPECTED_COVALENT_LINK")
        check({n[4] for n in atom["neighbors"]} == neighbors[name], "CONTACT_PROBE_BOND_GRAPH")
    missing = [a["atom"] for a in actual.values() if a["typing_error"] is not None
               or a["restraint_h_bond_type"] not in ("A", "B", "D", "N", "H")
               or not isinstance(a["probe_radius_A"], (int, float)) or not math.isfinite(a["probe_radius_A"])
               or a["probe_radius_A"] <= 0 or not a["energy_type"]]
    disagreements = []
    for role in context["heavy_atom_roles"]:
        atom = ligand[role["name"]]
        check(atom["probe_charge"] == role["formal_charge"], "CONTACT_PROBE_CHARGE_MISMATCH")
        if role["element"] not in ("N", "O", "S"):
            continue
        donor = any(ligand[n]["element"] == "H" and ligand[n]["probe_donor"] for n in neighbors[role["name"]])
        if bool(atom["probe_acceptor"]) != role["rdkit_acceptor"] or donor != role["rdkit_donor"]:
            disagreements.append({"name": role["name"], "atom_map": role["atom_map"],
                                  "rdkit_acceptor": role["rdkit_acceptor"], "probe_acceptor": atom["probe_acceptor"],
                                  "rdkit_donor": role["rdkit_donor"], "probe_attached_h_donor": donor,
                                  "restraint_h_bond_type": atom["restraint_h_bond_type"]})
    return {"status": "incomplete" if missing else "review_required", "all_atom_count": len(actual),
            "ligand_atom_count": len(ligand), "missing_atom_types": missing, "role_disagreements": disagreements,
            "probe_coordinate_rounding_max_A": max(coordinate_deltas),
            "meaning": "RDKit and Probe use different rules; discrepancies require review, not automatic correction."}


def project_contacts(context, files):
    models = []; all_completed = True
    for mode in ("xray", "nuclear"):
        model = context["models"][mode]; directions = {}; merged = {}; audits = []
        for direction in ("forward", "reverse"):
            folder = f"runs/{mode}/{direction}"
            record = strict_json(files[folder + "/execution.json"])
            check(type(record["exit_code"]) is int, "CONTACT_EXECUTION_EXIT_TYPE")
            source = [context["ligand_prefix"][0]] if direction == "forward" else model["target_chains"]
            target = model["target_chains"] if direction == "forward" else [context["ligand_prefix"][0]]
            check(record["input_sha256"] == model["sha256"] and record["mode"] == mode
                  and record["source_chains"] == source and record["target_chains"] == target, "CONTACT_EXECUTION_BINDING")
            check(sha(files[folder + "/probe.log"]) == record["log_sha256"], "CONTACT_EXECUTION_LOG_HASH")
            error = record["error"]; typing = None; pairs = {}
            if record["exit_code"] == 0 and error is None:
                try:
                    output = files[folder + "/probe.json"]; raw_audit = files[folder + "/atom-audit.json"]
                    check(sha(output) == record["output_sha256"] and sha(raw_audit) == record["audit_sha256"], "CONTACT_EXECUTION_OUTPUT_HASH")
                    audit = strict_json(raw_audit)
                    check(audit["source_chains"] == source and audit["target_chains"] == target, "CONTACT_AUDIT_SCOPE")
                    typing = analyze_atom_audit(context, mode, audit)
                    pairs = probe_pairs(output, source, target)
                    known = {tuple(a["atom"]) for a in model["atoms"]}
                    check(all(set(pair) <= known for pair in pairs), "CONTACT_PROBE_UNKNOWN_ATOM")
                    audits.append(typing)
                except (ValueError, KeyError, TypeError) as exc:
                    error = "CONTACT_RESULT_REJECTED: " + str(exc)
            else:
                error = error or "CONTACT_WORKER_FAILED"
            complete = error is None; all_completed &= complete
            directions[direction] = {"execution_status": "completed" if complete else "failed", "error": error,
                                     "typing": typing, "elapsed_seconds": record["elapsed_seconds"]}
            if complete:
                for key, pair in pairs.items():
                    item = merged.setdefault(key, {"atoms": pair["atoms"], "types": set(), "gap_A": pair["gap_A"], "directions": []})
                    check(abs(item["gap_A"] - pair["gap_A"]) < 0.002, "CONTACT_PAIR_GAP_DISAGREEMENT")
                    item["types"].update(pair["types"]); item["directions"].append(direction)
        if len(audits) == 2:
            check(audits[0] == audits[1], "CONTACT_DIRECTION_TYPING_DISAGREEMENT")
        bad = [{**v, "types": sorted(v["types"])} for _, v in sorted(merged.items()) if "bo" in v["types"]]
        models.append({"hydrogen_convention": mode, "model_sha256": model["sha256"], "directions": directions,
                       "bad_overlap_pair_count": len(bad) if audits else None, "bad_overlap_pairs": bad,
                       "both_directions_completed": len(audits) == 2, "coverage": "incomplete",
                       "ligand_hydrogen_count": len(model["h_bond_lengths"])})
    report = seal({"format": REPORT, "data_mode": context["data_mode"], "binding": copy.deepcopy(context["binding"]),
                   "input_digest": context["digest"], "prior_quality_digest": context["prior_quality_digest"],
                   "status": "completed_with_limits" if all_completed else "failed_or_partial", "models": models,
                   "method": context["method"], "protein_name_alias_count": len(context["protein_name_aliases"]),
                   "limitations": LIMITS, "hydrogen_orientation": "unoptimized", "probe_atom_typing": "review_required",
                   "human_review": "pending", "automatic_acceptance": False, "whole_structure_validated": False})
    schema_check(report, "b_ligand_contacts.schema.json")
    return report


def validate_runtime(library_manifest_path, ccd_manifest_path):
    library = strict_json(Path(library_manifest_path).read_bytes()); ccd = strict_json(Path(ccd_manifest_path).read_bytes())
    for ref in library["files"]:
        check(sha(safe_relative(library["library"], ref["relative"]).read_bytes()) == ref["sha256"], "CONTACT_LIBRARY_HASH")
    for ref in ccd["files"]:
        path = Path(ref["path"]).resolve()
        check(path.is_relative_to(Path(ccd["repository"]).resolve()) and sha(path.read_bytes()) == ref["sha256"], "CONTACT_CCD_HASH")
    return library, ccd


def execute_contacts(input_directory, destination, validation_python, library_manifest_path, ccd_manifest_path, *, allow_synthetic=False, timeout_seconds=120):
    check(1 <= timeout_seconds <= 120, "CONTACT_TIMEOUT_LIMIT")
    input_directory = Path(input_directory); destination = Path(destination)
    context = read_contact_input(input_directory, allow_synthetic=allow_synthetic)
    library, ccd = validate_runtime(library_manifest_path, ccd_manifest_path)
    files = {"input/" + n: d for n, d in tree_files(input_directory).items()}
    files["runtime/library.json"] = encoded(library); files["runtime/ccd.json"] = encoded(ccd)
    write_new_directory(destination, files)
    env = dict(os.environ, MMTBX_CCP4_MONOMER_LIB=library["library"])
    for mode in ("xray", "nuclear"):
        model = context["models"][mode]
        for direction in ("forward", "reverse"):
            folder = destination / "runs" / mode / direction; folder.mkdir(parents=True)
            source = [context["ligand_prefix"][0]] if direction == "forward" else model["target_chains"]
            target = model["target_chains"] if direction == "forward" else [context["ligand_prefix"][0]]
            command = [str(validation_python), "-m", "packages.science.probe_worker", "--model", str((destination / "input" / model["path"]).resolve()),
                       "--input-sha256", model["sha256"], "--restraints", str((destination / "input/raw/restraints.cif").resolve()),
                       "--repository", ccd["repository"], "--out", str((folder / "probe.json").resolve()), "--audit", str((folder / "atom-audit.json").resolve()),
                       "--mode", mode, "--source-chains", *source, "--target-chains", *target]
            started = datetime.now(timezone.utc).isoformat(); clock = time.monotonic(); error = None
            with (folder / "probe.log").open("xb") as log:
                try:
                    result = subprocess.run(command, cwd=ACTIVE, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=timeout_seconds)
                    code = result.returncode
                except subprocess.TimeoutExpired:
                    code = 124; error = "CONTACT_TOOL_TIMEOUT"
                except OSError as exc:
                    code = 127; error = "CONTACT_TOOL_START_FAILED: " + str(exc)
            record = {"started_at": started, "elapsed_seconds": time.monotonic()-clock, "command": command, "exit_code": code, "error": error,
                      "mode": mode, "source_chains": source, "target_chains": target, "input_sha256": model["sha256"],
                      "output_sha256": sha((folder / "probe.json").read_bytes()) if (folder / "probe.json").exists() else None,
                      "audit_sha256": sha((folder / "atom-audit.json").read_bytes()) if (folder / "atom-audit.json").exists() else None,
                      "log_sha256": sha((folder / "probe.log").read_bytes())}
            (folder / "execution.json").write_bytes(encoded(record))
    files = tree_files(destination)
    report = project_contacts(context, files)
    (destination / "report.json").write_bytes(encoded(report))
    files["report.json"] = encoded(report)
    (destination / "manifest.json").write_bytes(encoded(seal({"format": ARCHIVE, "files": file_index(files), "report_digest": report["digest"]})))
    return report


def read_contacts(directory, *, allow_synthetic=False):
    directory = Path(directory)
    manifest, files = archive_files(directory, ARCHIVE)
    context = read_contact_input(directory / "input", allow_synthetic=allow_synthetic)
    report = project_contacts(context, files)
    check(encoded(report) == files["report.json"] and manifest["report_digest"] == report["digest"], "CONTACT_REPORT_PROJECTION")
    return report
