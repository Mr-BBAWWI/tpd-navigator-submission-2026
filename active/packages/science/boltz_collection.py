"""Collect all saved samples and retain missing/invalid/failed attempts. Never run Boltz."""
from __future__ import annotations

import copy
from pathlib import Path
import re
import tempfile

from .evidence_common import file_index, safe_relative, schema_check, seal, verify_seal, write_new_directory
from .handoff import check, encoded, parse, sha, build_material


def candidate_from_receipt(receipt, cid):
    matches = [c for c in receipt["candidates"] if c["compound_id"] == cid]
    check(len(matches) == 1, "RUN_CANDIDATE_NOT_IN_RECEIPT")
    return matches[0]


def validate_run_info(info, receipt):
    schema_check(info, "b_boltz_run.schema.json")
    check(info["environment"]["rdkit"] == receipt["runtime"]["rdkit"], "RUN_RDKIT_DIFFERS_FROM_ATOM_NAMING_RECEIPT")
    return candidate_from_receipt(receipt, info["compound_id"])


def collect_run(predictions: Path, input_path: Path, receipt, info, destination: Path):
    """Archive a single record's output folder, including malformed and auxiliary files.

    Run metadata is a caller report, not an authenticated execution or approval record.
    """
    verify_seal(receipt)
    schema_check(receipt, "b_boltz_parser.schema.json")
    c = validate_run_info(info, receipt)
    input_bytes = input_path.read_bytes()
    check(sha(input_bytes) == c["input_sha256"], "RUN_INPUT_HASH_MISMATCH")
    check(predictions.is_dir(), "PREDICTION_DIRECTORY_REQUIRED")
    files = {"run-info.json": encoded(info), "input.yaml": input_bytes}
    for path in sorted(predictions.iterdir()):
        check(path.is_file() and not path.is_symlink(), "SCOPED_PREDICTION_FILES_ONLY")
        files["raw/" + path.name] = path.read_bytes()
    manifest = seal({"format": "tpd-boltz-run-archive/0.1.0-draft", "run_id": info["run_id"],
                     "compound_id": c["compound_id"], "molecule_id": c["molecule_id"],
                     "receipt_digest": receipt["digest"], "input_sha256": c["input_sha256"],
                     "files": file_index(files), "execution_provenance": "caller_reported_not_authenticated"})
    files["archive.json"] = encoded(manifest)
    write_new_directory(destination, files)
    return manifest


def load_archive(directory, receipt):
    manifest = parse((directory / "archive.json").read_bytes())
    schema_check(manifest, "b_boltz_run_archive.schema.json")
    verify_seal(manifest)
    check(manifest["receipt_digest"] == receipt["digest"], "RUN_STALE_RECEIPT")
    files = {}
    for name, item in manifest["files"].items():
        data = safe_relative(directory, name).read_bytes()
        check(sha(data) == item["sha256"] and len(data) == item["bytes"], "RUN_ARCHIVE_HASH_MISMATCH")
        files[name] = data
    actual = {p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file()}
    check(actual == set(files) | {"archive.json"}, "RUN_ARCHIVE_UNINDEXED_FILES")
    info = parse(files["run-info.json"])
    c = validate_run_info(info, receipt)
    check(info["run_id"] == manifest["run_id"] and c["compound_id"] == manifest["compound_id"]
          and c["molecule_id"] == manifest["molecule_id"] and sha(files["input.yaml"]) == c["input_sha256"] == manifest["input_sha256"],
          "RUN_ARCHIVE_IDENTITY")
    return manifest, info, c, files


def conditions(info):
    # Seed and sample count intentionally remain part of each separately grouped attempt.
    return {k: copy.deepcopy(info[k]) for k in ("environment", "settings", "weights_sha256", "msa_files_sha256", "requested_samples")}


def summarize_result(result):
    comparison = result["comparison"]
    mappings = comparison["ligand"]["all_mappings"]
    warnings = []
    if not all(m["coordinates_match_input_stereochemistry"] for m in mappings):
        warnings.append("STEREOCHEMISTRY_REQUIRES_REVIEW")
    if result["confidence"]["status"] == "missing":
        warnings.append("CONFIDENCE_MISSING")
    return {"confidence_score": None if result["confidence"]["raw"] is None else result["confidence"]["raw"]["confidence_score"],
            "target_rmsd_A": comparison["target_alignment"]["rmsd_A"],
            "vhl_rmsd_after_target_alignment_A": comparison["proteins"][comparison["vhl_chain"]]["rmsd_after_target_alignment_A"],
            "ligand_rmsd_range_A": comparison["ligand"]["whole_ligand_rmsd_range_A"],
            "warhead_rmsd_range_A": [min(m["parts_rmsd_A"]["warhead"] for m in mappings), max(m["parts_rmsd_A"]["warhead"] for m in mappings)],
            "symmetry_mapping_count": len(mappings), "warnings": warnings}


def build_collection(snapshot, receipt, directories, *, allow_synthetic=False):
    from .boltz_parser import validate_receipt
    from .boltz_results import read_result

    validate_receipt(receipt, snapshot)
    material = build_material(snapshot)
    runs, exported, seen = [], {}, set()
    for directory in directories:
        directory = Path(directory)
        manifest, info, c, files = load_archive(directory, receipt)
        check(info["run_id"] not in seen, "DUPLICATE_RUN_ID")
        seen.add(info["run_id"])
        check(info["origin"] != "synthetic_test" or allow_synthetic, "SYNTHETIC_RUN_NOT_ALLOWED")
        prefix = "runs/" + info["run_id"] + "/"
        exported.update({prefix + n: d for n, d in files.items()})
        exported[prefix + "archive.json"] = encoded(manifest)
        # Include every rank encountered, including orphan confidence/npz files.
        pattern = re.compile(r"^(?:confidence_|pae_|pde_|plddt_)?" + re.escape(c["record_id"]) + r"_model_(\d+)\.(?:cif|json|npz|pdb)$")
        found = {int(match.group(1)) for n in files if n.startswith("raw/") and (match := pattern.fullmatch(n[4:]))}
        expected = set(range(info["requested_samples"]))
        unrecognized = [n for n in files if n.startswith("raw/") and not pattern.fullmatch(n[4:])]
        warnings = []
        if info["exit_code"] != 0:
            warnings.append("PROCESS_FAILED_OR_EXIT_UNKNOWN")
        if not info["msa_files_sha256"]:
            warnings.append("MSA_PROVENANCE_MISSING")
        if info["log_sha256"] is None:
            warnings.append("EXECUTION_LOG_PROVENANCE_MISSING")
        if info["weights_sha256"] is None:
            warnings.append("WEIGHTS_PROVENANCE_MISSING")
        if found - expected:
            warnings.append("UNEXPECTED_MODEL_RANKS")
        if unrecognized:
            warnings.append("UNRECOGNIZED_OUTPUTS_PRESERVED")
        models = []
        for rank in sorted(expected | found):
            row = {"rank": rank, "expected": rank in expected, "status": "missing", "reason": "STRUCTURE_NOT_PROVIDED",
                   "summary": None, "comparison_file": None}
            cif_name = f'{c["record_id"]}_model_{rank}.cif'
            conf_name = f'confidence_{c["record_id"]}_model_{rank}.json'
            if "raw/" + cif_name in files:
                data = {c["input_name"]: files["input.yaml"], cif_name: files["raw/" + cif_name]}
                if "raw/" + conf_name in files:
                    data[conf_name] = files["raw/" + conf_name]
                single = {"format": "tpd-boltz-output-manifest/0.1.0-draft", "origin": info["origin"],
                          "compound_id": c["compound_id"], "molecule_id": c["molecule_id"], "receipt_digest": receipt["digest"],
                          "input_sha256": c["input_sha256"], "record_id": c["record_id"], "model_rank": rank,
                          "runtime_reported": {"boltz_version": info["environment"]["boltz"], "weights_sha256": info["weights_sha256"], "run_id": info["run_id"]},
                          "files": {"input": {"name": c["input_name"], "sha256": sha(data[c["input_name"]])},
                                    "structure": {"name": cif_name, "sha256": sha(data[cif_name])},
                                    "confidence": {"name": conf_name, "sha256": sha(data[conf_name])} if conf_name in data else None}}
                data["manifest.json"] = encoded(single)
                try:
                    with tempfile.TemporaryDirectory(prefix="tpd-boltz-read-") as tmp:
                        model_dir = Path(tmp) / "model"
                        write_new_directory(model_dir, data)
                        result = read_result(model_dir, receipt, snapshot, allow_synthetic=allow_synthetic)
                    result_name = f"comparisons/{info['run_id']}-rank-{rank}.json"
                    result_bytes = encoded(result)
                    exported[result_name] = result_bytes
                    row.update(status="parsed", reason=None, summary=summarize_result(result),
                               comparison_file={"path": result_name, "sha256": sha(result_bytes)})
                except (ValueError, RuntimeError, KeyError, TypeError, IndexError) as error:
                    # Retain the raw bytes and failure; continue parsing the other samples.
                    row.update(status="invalid", reason=f"{type(error).__name__}: {error}")
            models.append(row)
        all_parsed = all(r["status"] == "parsed" for r in models if r["expected"])
        all_confidence = all(not r["summary"]["warnings"] for r in models if r["status"] == "parsed")
        runs.append({"run_id": info["run_id"], "compound_id": c["compound_id"], "molecule_id": c["molecule_id"],
                     "origin": info["origin"], "reported_execution": info,
                     "condition_digest": sha(encoded(conditions(info))), "input_sha256": c["input_sha256"],
                     "archive_file": {"path": prefix + "archive.json", "sha256": sha(encoded(manifest))},
                     "readiness": "available_unreviewed" if all_parsed and all_confidence and not warnings else "needs_review",
                     "warnings": warnings, "unrecognized_files": unrecognized, "models": models})
    collection = seal({"format": "tpd-boltz-collection/0.1.0-draft", "case_id": material["case_id"],
                       "material_digest": material["material_digest"], "receipt_digest": receipt["digest"],
                       "data_mode": "synthetic_test" if any(r["origin"] == "synthetic_test" for r in runs) else "real",
                       "runs": runs, "aggregation": "all_samples_kept_no_pooled_score_or_automatic_selection",
                       "human_review": "pending", "dispatch_authorized": False})
    validate_collection(collection, material, allow_synthetic=allow_synthetic)
    return collection, exported


def empty_collection(material):
    return seal({"format": "tpd-boltz-collection/0.1.0-draft", "case_id": material["case_id"],
                 "material_digest": material["material_digest"], "receipt_digest": None, "data_mode": "real", "runs": [],
                 "aggregation": "all_samples_kept_no_pooled_score_or_automatic_selection", "human_review": "pending", "dispatch_authorized": False})


def validate_collection(collection, material, *, allow_synthetic=False):
    schema_check(collection, "b_boltz_collection.schema.json")
    verify_seal(collection)
    check(collection["material_digest"] == material["material_digest"] and collection["case_id"] == material["case_id"], "COLLECTION_STALE_MATERIAL")
    check(allow_synthetic or collection["data_mode"] != "synthetic_test", "SYNTHETIC_COLLECTION_NOT_ALLOWED")
    check(len({r["run_id"] for r in collection["runs"]}) == len(collection["runs"]), "DUPLICATE_RUN_ID")
    molecules = {c["compound_id"]: c["molecule_id"] for c in material["h1"]["compounds"] if c["role"] == "known_reference_protac"}
    for run in collection["runs"]:
        info = run["reported_execution"]
        schema_check(info, "b_boltz_run.schema.json")
        check(run["molecule_id"] == molecules.get(run["compound_id"]), "COLLECTION_MOLECULE_MISMATCH")
        check(all(run[k] == info[k] for k in ("run_id", "compound_id", "origin")), "COLLECTION_RUN_IDENTITY")
        check(run["condition_digest"] == sha(encoded(conditions(info))), "COLLECTION_CONDITION_MISMATCH")
        expected = set(range(info["requested_samples"]))
        models = run["models"]
        check(len({m["rank"] for m in models}) == len(models) and {m["rank"] for m in models if m["expected"]} == expected
              and all(m["expected"] == (m["rank"] in expected) for m in models), "COLLECTION_SAMPLE_COVERAGE")
        check(collection["data_mode"] == "synthetic_test" or run["origin"] != "synthetic_test", "COLLECTION_ORIGIN_MISMATCH")
        ready = (info["exit_code"] == 0 and bool(info["msa_files_sha256"]) and info["log_sha256"] is not None
                 and info["weights_sha256"] is not None
                 and not run["warnings"] and not run["unrecognized_files"]
                 and {m["rank"] for m in models} == expected
                 and all(m["status"] == "parsed" and not m["summary"]["warnings"] for m in models))
        check(run["readiness"] == ("available_unreviewed" if ready else "needs_review"), "COLLECTION_READINESS_MISMATCH")
    return collection
