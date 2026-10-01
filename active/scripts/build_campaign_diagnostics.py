#!/usr/bin/env python3
"""Build a portable, read-only campaign diagnostics evidence pack.

The builder validates the original reports and raw model artifacts, copies only
an explicit allowlist, writes a compact derived index, and never changes source
scientific state or approval fields.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import html
import importlib
import json
import math
import os
import re
import shutil
import statistics
import sys
import uuid
from pathlib import Path, PurePosixPath
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

ALLOWED_SUFFIXES = {".json", ".cif", ".sdf", ".yaml", ".html", ".py", ".md"}
SEEDS = [23, 41, 61, 79, 97]
COMPACT_LIMIT = 20 * 1024


class BuildError(ValueError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise BuildError("JSON_READ_FAILED:" + str(path)) from exc


def is_relative_to(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def reject_symlink_components(path: Path, stop: Path | None = None) -> None:
    current = path
    while True:
        if current.is_symlink():
            raise BuildError("SYMLINK_REJECTED:" + str(path))
        if stop is not None and current == stop:
            return
        parent = current.parent
        if parent == current:
            return
        current = parent


def safe_source(path: Path, active_root: Path, suffixes: set[str] = ALLOWED_SUFFIXES) -> Path:
    """Return a confined regular source path without following symlinks."""
    active_root = active_root.resolve()
    raw = str(path).replace("\\", "/")
    raw_parts = [part for part in raw.split("/") if part]
    if ".." in raw_parts:
        raise BuildError("SOURCE_TRAVERSAL_REJECTED:" + str(path))
    lowered = [part.lower() for part in raw_parts]
    protected = {".git", "boltz-cache", "secrets", ".secrets", "auth"}
    if any(part in protected or part.startswith(".venv") for part in lowered):
        raise BuildError("PROTECTED_SOURCE_PATH_REJECTED:" + str(path))
    filename = lowered[-1] if lowered else ""
    if (filename in {"id_rsa", "id_ed25519", ".env", "credentials", "credentials.json"}
            or re.search(r"(?:credential|password|passwd|api[_-]?key|access[_-]?token|secret[_-]?key)", filename)):
        raise BuildError("CREDENTIAL_SOURCE_REJECTED:" + str(path))
    candidate = path if path.is_absolute() else active_root / path
    reject_symlink_components(candidate.absolute(), active_root)
    resolved = candidate.resolve(strict=True)
    if not is_relative_to(resolved, active_root):
        raise BuildError("SOURCE_OUTSIDE_ACTIVE_ROOT:" + str(path))
    if resolved.is_symlink() or not resolved.is_file():
        raise BuildError("REGULAR_SOURCE_REQUIRED:" + str(path))
    if resolved.suffix.lower() not in suffixes:
        raise BuildError("SOURCE_TYPE_NOT_ALLOWED:" + str(path))
    return resolved


def safe_pack_relative(value: str | Path) -> Path:
    original = str(value)
    if re.match(r"^[A-Za-z]:", original) or original.startswith(("\\\\", "//")):
        raise BuildError("INVALID_PACK_PATH:" + original)
    raw = original.replace("\\", "/")
    segments = raw.split("/")
    if not raw or any(part in {"", ".", ".."} for part in segments):
        raise BuildError("INVALID_PACK_PATH:" + raw)
    pure = PurePosixPath(raw)
    if pure.is_absolute() or not pure.parts:
        raise BuildError("INVALID_PACK_PATH:" + raw)
    path = Path(*pure.parts)
    if path.suffix.lower() not in ALLOWED_SUFFIXES:
        raise BuildError("PACK_TYPE_NOT_ALLOWED:" + raw)
    return path


def source_record(source: Path, destination: Path, active_root: Path) -> dict[str, str]:
    return {
        "source_sha256": sha256_file(source),
        "active_root_relative_path": source.relative_to(active_root).as_posix(),
        "portable_pack_relative_path": destination.as_posix(),
    }


def copy_verified(source: Path, stage: Path, destination: str | Path,
                  active_root: Path, expected_hash: str | None = None) -> dict[str, str]:
    source = safe_source(source, active_root)
    relative = safe_pack_relative(destination)
    actual = sha256_file(source)
    if expected_hash is not None and actual != str(expected_hash).lower():
        raise BuildError("SOURCE_HASH_MISMATCH:" + str(source))
    target = stage / relative
    if target.exists():
        raise BuildError("DUPLICATE_PACK_PATH:" + relative.as_posix())
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as incoming, target.open("xb") as outgoing:
        shutil.copyfileobj(incoming, outgoing, 1024 * 1024)
    if sha256_file(target) != actual:
        raise BuildError("COPY_HASH_MISMATCH:" + relative.as_posix())
    return source_record(source, relative, active_root)


def write_json_exclusive(path: Path, value: Any, compact: bool = False) -> None:
    with path.open("xb") as stream:
        if compact:
            data = json.dumps(value, ensure_ascii=False, sort_keys=True,
                              separators=(",", ":")).encode("utf-8") + b"\n"
        else:
            data = json.dumps(value, ensure_ascii=False, sort_keys=True,
                              indent=2).encode("utf-8") + b"\n"
        stream.write(data)


def write_compact(path: Path, value: Any, limit: int = COMPACT_LIMIT) -> int:
    data = json.dumps(value, ensure_ascii=False, sort_keys=True,
                      separators=(",", ":")).encode("utf-8") + b"\n"
    if len(data) >= limit:
        raise BuildError(f"COMPACT_JSON_TOO_LARGE:{len(data)}")
    with path.open("xb") as stream:
        stream.write(data)
    return len(data)


def verify_artifact_manifest(report_root: Path, manifest: dict, active_root: Path) -> None:
    report_root = report_root.resolve(strict=True)
    if not is_relative_to(report_root, active_root.resolve()):
        raise BuildError("REPORT_ROOT_OUTSIDE_ACTIVE_ROOT:" + str(report_root))
    entries = manifest.get("artifacts")
    if isinstance(entries, list):
        pairs = [(item.get("manifest_relativepath", item.get("path")), item.get("sha256"))
                 for item in entries if isinstance(item, dict)]
    else:
        files = manifest.get("hashed_files")
        if not isinstance(files, dict):
            files = manifest.get("files", manifest.get("output_files_sha256"))
        pairs = list(files.items()) if isinstance(files, dict) else []
    if not pairs:
        raise BuildError("MANIFEST_ARTIFACTS_MISSING:" + str(report_root))
    for relative, expected in pairs:
        if not isinstance(relative, str) or not isinstance(expected, str) or len(expected) != 64:
            raise BuildError("INVALID_MANIFEST_ENTRY:" + str(report_root))
        confined = safe_pack_relative(relative)
        source = safe_source(report_root / confined, report_root)
        if sha256_file(source) != expected.lower():
            raise BuildError("MANIFEST_ARTIFACT_HASH_MISMATCH:" + relative)
    exact = manifest.get("exact_output_file_set")
    if exact is not None:
        if not isinstance(exact, list) or not all(isinstance(item, str) for item in exact):
            raise BuildError("INVALID_EXACT_OUTPUT_FILE_SET:" + str(report_root))
        expected_set = {safe_pack_relative(item).as_posix() for item in exact}
        actual_set = {
            path.relative_to(report_root).as_posix()
            for path in report_root.rglob("*") if path.is_file()
        }
        if actual_set != expected_set:
            raise BuildError("MANIFEST_EXACT_OUTPUT_FILE_SET_MISMATCH:" + str(report_root))


def _walk_values(value: Any) -> Iterable[dict]:
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from _walk_values(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk_values(child)


def _rows(summary: dict) -> list[dict]:
    for key in ("rows", "state_rows", "results"):
        value = summary.get(key)
        if isinstance(value, list) and all(isinstance(row, dict) for row in value):
            return value
    raise BuildError("CORE_ROWS_MISSING")


def build_core_compact(summary: dict) -> dict:
    rows = _rows(summary)
    if len(rows) != 40:
        raise BuildError("CORE_ROW_COUNT_NOT_40")
    true_count = sum(row.get("required_contact_evidence") is True for row in rows)
    false_count = sum(row.get("required_contact_evidence") is False for row in rows)
    if (true_count, false_count) != (16, 24):
        raise BuildError("CORE_CONTACT_DENOMINATOR_INVALID")
    failure_count = sum(row.get("failure") is not None for row in rows)
    null_failure_count = sum(row.get("failure") is None for row in rows)

    pka_rows = []
    for item in _walk_values(summary):
        keys = {"chain", "residue", "sequence_number", "pKa"}
        if keys.issubset(item):
            pka_rows.append(item)
    if len(pka_rows) != 88:
        raise BuildError("PROTEIN_PKA_RAW_ROW_COUNT_NOT_88")
    unique = {
        (str(row["chain"]), str(row["residue"]), str(row["sequence_number"]),
         json.dumps(row["pKa"], sort_keys=True))
        for row in pka_rows
    }
    if len(unique) != 44:
        raise BuildError("PROTEIN_PKA_UNIQUE_COUNT_NOT_44")

    analogs = summary.get("analogs")
    if not isinstance(analogs, list):
        analogs = []
        for analog_id in sorted({str(row.get("analog_id")) for row in rows}):
            group = [row for row in rows if str(row.get("analog_id")) == analog_id]
            analogs.append({
                "analog_id": analog_id,
                "retained_row_count": len(group),
                "supported_pose_indices": sorted({row.get("pose_index") for row in group
                                                   if row.get("required_contact_evidence") is True}),
            })
    analog_summary = [{key: item.get(key) for key in (
        "analog_id", "pose_count", "requested_state_count_per_pose",
        "retained_row_count", "pose0_policy", "supported_pose_indices",
        "all_states_robust_at_same_pose") if key in item} for item in analogs]
    return {
        "row_count": 40,
        "required_contact_true_count": true_count,
        "required_contact_false_count": false_count,
        "failure_count": failure_count,
        "failure_null_count": null_failure_count,
        "analogs": analog_summary,
        "ligand_pKa_prediction_performed": False,
        "population_prediction_performed": False,
        "protein_pKa_raw_row_count": 88,
        "protein_pKa_unique_chain_residue_sequence_pKa_count": 44,
        "protein_pKa_duplicate_explanation": "88 raw rows repeat 44 unique entries across two log sources; they are not 88 unique pKa values.",
        "microstate_gate": "pending",
        "protein_hydrogen_orientation": "pending_review",
        "state_selection_performed": False,
        "tautomer_selection_performed": False,
        "best_pose_selected": False,
        "diagnostic_only": True,
        "scientific_approved": False,
    }


def _finite(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        raise BuildError("FINITE_NUMBER_REQUIRED:" + label)
    return float(value)


def _stats(values: list[Any], label: str) -> dict[str, float | None]:
    numbers = [_finite(value, label) for value in values if value is not None]
    if not numbers:
        return {"min": None, "median": None, "max": None}
    return {"min": min(numbers), "median": statistics.median(numbers), "max": max(numbers)}


def _candidate_rows(topology: dict) -> list[dict]:
    value = topology.get("candidates", topology.get("candidate_diagnostics"))
    if not isinstance(value, list):
        raise BuildError("TOPOLOGY_CANDIDATES_MISSING")
    return value


def build_topology_compact(summary: dict, comparison: dict | None = None) -> dict:
    candidates = _candidate_rows(summary)
    if len(candidates) != 6 or summary.get("raw_seed_result_count") != 30:
        raise BuildError("TOPOLOGY_EXACT_6_BY_5_REQUIRED")
    compact_rows = []
    for candidate in candidates:
        cid = candidate.get("candidate_id")
        seeds = candidate.get("seeds")
        if not isinstance(seeds, list) or len(seeds) != 5:
            raise BuildError("TOPOLOGY_ALL_FIVE_SEEDS_REQUIRED:" + str(cid))
        under2 = []
        vdw = []
        for seed in seeds:
            clashes = seed.get("interface_clashes") if isinstance(seed, dict) else None
            if not isinstance(clashes, dict):
                raise BuildError("TOPOLOGY_INTERFACE_CLASHES_MISSING:" + str(cid))
            under2.append(clashes.get("raw_under_2A_pair_count"))
            vdw.append(clashes.get("vdw_pair_count"))
        if any(isinstance(value, bool) or not isinstance(value, int) for value in under2 + vdw):
            raise BuildError("TOPOLOGY_CLASH_COUNT_INVALID:" + str(cid))
        pairs = candidate.get("pairwise_contact_jaccard")
        if not isinstance(pairs, list) or len(pairs) != 10:
            raise BuildError("TOPOLOGY_TEN_PAIRWISE_ROWS_REQUIRED:" + str(cid))
        def values(key: str) -> list[Any]:
            return [item.get(key) for item in pairs if isinstance(item, dict)]
        sensitivity = candidate.get("cluster_sensitivity",
                                    candidate.get("complete_link_contact_jaccard_sensitivity"))
        if (not isinstance(sensitivity, list) or len(sensitivity) != 3
                or [item.get("threshold") for item in sensitivity if isinstance(item, dict)] != [0.3, 0.5, 0.7]):
            raise BuildError("TOPOLOGY_CLUSTER_THRESHOLDS_CHANGED:" + str(cid))
        compact_rows.append({
            "candidate_id": cid,
            "e3_type": candidate.get("e3_type"),
            "priority_diagnostic": candidate.get("priority_diagnostic"),
            "protein_protein_contact_jaccard": _stats(values("protein_protein_contact_jaccard"), key_name("pp")),
            "target_warhead_site_jaccard": _stats(values("target_warhead_site_jaccard"), key_name("target")),
            "e3_recruiter_site_jaccard": _stats(values("e3_recruiter_site_jaccard"), key_name("e3")),
            "clashes_under2_A_all_five": under2,
            "vdw_clashes_all_five": vdw,
            "cluster_sensitivity": sensitivity,
        })
    return {
        "candidate_count": 6,
        "raw_seed_result_count": 30,
        "seeds": summary.get("seeds", SEEDS),
        "candidates": compact_rows,
        "cluster_thresholds": [0.3, 0.5, 0.7],
        "cluster_labels_are_descriptive_not_acceptance": True,
        "priority_interpretation": "Priority labels identify poor-convergence diagnostics from the descriptive narrative only; no scientifically accepted cutoff is introduced.",
        "reference_free": True,
        "scientific_approved": False,
    }


def key_name(value: str) -> str:
    return "topology_" + value


def _selected_index(receipt: dict) -> int:
    selection = receipt.get("selection")
    candidates = []
    if isinstance(selection, dict):
        candidates += [selection.get("model_index"), selection.get("selected_model_index")]
    candidates += [receipt.get("selected_model_index")]
    for value in candidates:
        if isinstance(value, int) and not isinstance(value, bool):
            return value
    models = receipt.get("models", [])
    selected = [model.get("model_index") for model in models
                if isinstance(model, dict) and model.get("selected") is True]
    if len(selected) == 1:
        return selected[0]
    raise BuildError("CALIBRATION_SELECTED_MODEL_INDEX_MISSING")


def build_calibration_compact(plan: dict, receipts: list[dict], metrics_pass: Any) -> dict:
    if len(receipts) != 5 or {item.get("seed") for item in receipts} != set(SEEDS):
        raise BuildError("CALIBRATION_EXACT_FIVE_SEEDS_REQUIRED")
    all_models = []
    selected = []
    per_seed = []
    for receipt in sorted(receipts, key=lambda item: item["seed"]):
        models = receipt.get("models")
        if not isinstance(models, list) or len(models) != 5:
            raise BuildError("CALIBRATION_FIVE_MODELS_PER_SEED_REQUIRED")
        all_models.extend(models)
        index = _selected_index(receipt)
        model = next((item for item in models if item.get("model_index") == index), None)
        if model is None:
            raise BuildError("CALIBRATION_SELECTION_NOT_IN_MODELS")
        passed = bool(metrics_pass(model.get("metrics", {})))
        selected.append(passed)
        per_seed.append({
            "seed": receipt["seed"],
            "selected_model_index": index,
            "selected_pass": passed,
            "selected_e3_CA_RMSD_after_target_alignment_A": model.get("metrics", {}).get(
                "e3_CA_RMSD_after_target_alignment_A"),
        })
    raw_pass = sum(bool(metrics_pass(model.get("metrics", {}))) for model in all_models)
    if len(all_models) != 25 or raw_pass != 16 or sum(selected) != 2:
        raise BuildError("CALIBRATION_RAW_OR_SELECTED_DENOMINATOR_CHANGED")
    baseline = plan.get("preserved_original_baseline")
    if baseline != {"geometric_successes": 2, "seed_count": 5}:
        raise BuildError("CALIBRATION_ORIGINAL_BASELINE_CHANGED")
    return {
        "case_id": "6BOY_CRBN",
        "original_baseline": "2_of_5",
        "seed_count": 5,
        "raw_model_count": 25,
        "raw_all_four_pass_count": 16,
        "selected_pass_count": 2,
        "per_seed_selected": per_seed,
        "criteria_unchanged": True,
        "cutoffs": plan.get("cutoffs"),
        "effective_protocol_settings": plan.get("effective_protocol_settings"),
        "failure_cause": "E3_relative_placement",
        "scientific_approved": False,
    }


def validate_provenance(value: Any, active_root: Path) -> None:
    for item in _walk_values(value):
        if "path" in item and "sha256" in item:
            path, expected = item["path"], item["sha256"]
            if not isinstance(path, str) or not isinstance(expected, str) or len(expected) != 64:
                raise BuildError("RANKING_PROVENANCE_INVALID")
            source = safe_source(Path(path), active_root)
            if sha256_file(source) != expected.lower():
                raise BuildError("RANKING_PROVENANCE_HASH_MISMATCH:" + path)


def build_ranking_compact(diagnosis: dict) -> dict:
    evaluations = diagnosis.get("evaluations")
    if not isinstance(evaluations, dict) or len(evaluations) != 7:
        raise BuildError("RANKING_EXACT_SEVEN_RULES_REQUIRED")
    rules = {}
    for name, result in evaluations.items():
        if not isinstance(result, dict):
            raise BuildError("RANKING_RESULT_INVALID")
        count = result.get("selected_pass_count")
        valid = result.get("valid", count is not None)
        error = result.get("error")
        if not isinstance(valid, bool) or (error is not None and not isinstance(error, str)):
            raise BuildError("RANKING_VALIDITY_INVALID:" + name)
        if count is not None and (isinstance(count, bool) or not isinstance(count, int) or not 0 <= count <= 5):
            raise BuildError("RANKING_PASS_COUNT_INVALID:" + name)
        if valid and count is None:
            raise BuildError("RANKING_VALID_RESULT_WITHOUT_COUNT:" + name)
        rules[name] = {"selected_pass_count": count, "valid": valid, "error": error}
    return {
        "rules": rules,
        "development_only": True,
        "does_not_replace_original_2_of_5": True,
        "scientific_approved": False,
    }


def _slim_novel_candidate(candidate: dict) -> dict:
    seeds = candidate.get("seeds")
    if not isinstance(seeds, list) or len(seeds) != 5:
        raise BuildError("NOVEL_ALL_FIVE_SEEDS_REQUIRED")
    selections = []
    for row in seeds:
        selection = row.get("selection") if isinstance(row, dict) else None
        selected = selection.get("selected_model_index") if isinstance(selection, dict) else None
        if isinstance(selected, bool) or not isinstance(selected, int):
            raise BuildError("NOVEL_SELECTED_MODEL_INDEX_MISSING")
        selections.append({"seed": row.get("seed"), "selected_model_index": selected})
    quantitative = candidate.get("old_expert_quantitative_metrics")
    if not isinstance(quantitative, dict):
        raise BuildError("NOVEL_QUANTITATIVE_METRICS_MISSING")
    for key in ("finite_selected_ipTM_values", "finite_selected_endpoint_values_A"):
        values = quantitative.get(key)
        if not isinstance(values, list) or len(values) != 5:
            raise BuildError("NOVEL_QUANTITATIVE_ARRAY_INVALID:" + key)
        for index, value in enumerate(values):
            _finite(value, "novel_" + key + "_" + str(index))
    for key in ("ipTM_IQR", "ipTM_IQR_maximum", "endpoint_IQR_A", "endpoint_IQR_maximum_A"):
        _finite(quantitative.get(key), "novel_" + key)
    contacts = quantitative.get("positive_target_warhead_and_e3_recruiter_contacts")
    if not isinstance(contacts, list) or len(contacts) != 5 or not all(isinstance(value, bool) for value in contacts):
        raise BuildError("NOVEL_QUANTITATIVE_CONTACT_FLAGS_INVALID")
    for key in ("quantitative_criteria_met", "automatic_approval",
                "qualitative_topology_or_clash_acceptance_used"):
        if not isinstance(quantitative.get(key), bool):
            raise BuildError("NOVEL_QUANTITATIVE_BOOLEAN_INVALID:" + key)
    diagnostics = candidate.get("geometric_diagnostics")
    if not isinstance(diagnostics, dict):
        raise BuildError("NOVEL_GEOMETRIC_DIAGNOSTICS_MISSING")
    clashes = diagnostics.get("all_five_seeds_including_failures")
    pairs = diagnostics.get("pairwise_contact_jaccard")
    clusters = diagnostics.get("cluster_sensitivity")
    if not isinstance(clashes, list) or len(clashes) != 5:
        raise BuildError("NOVEL_FIVE_CLASH_ROWS_REQUIRED")
    if not isinstance(pairs, list) or len(pairs) != 10:
        raise BuildError("NOVEL_TEN_PAIRWISE_ROWS_REQUIRED")
    if (not isinstance(clusters, list) or len(clusters) != 3
            or [row.get("threshold") for row in clusters if isinstance(row, dict)] != [0.3, 0.5, 0.7]):
        raise BuildError("NOVEL_CLUSTER_THRESHOLDS_CHANGED")
    def pair_stats(key: str) -> dict[str, float | None]:
        return _stats([row.get(key) for row in pairs if isinstance(row, dict)], "novel_" + key)
    return {
        "candidate_id": candidate.get("candidate_id"),
        "e3_type": candidate.get("e3_type"),
        "selected_top1_per_seed": selections,
        "all_five_quantitative_metrics": quantitative,
        "geometric_diagnostics": {
            "protein_protein_jaccard_stats": pair_stats("protein_protein_jaccard"),
            "target_site_jaccard_stats": pair_stats("target_site_jaccard"),
            "e3_site_jaccard_stats": pair_stats("e3_site_jaccard"),
            "raw_under_2A_pair_count_all_five": [row.get("interface_clashes", {}).get("raw_under_2A_pair_count") for row in clashes],
            "vdw_pair_count_all_five": [row.get("interface_clashes", {}).get("vdw_pair_count") for row in clashes],
            "cluster_sensitivity": clusters,
        },
    }


def build_novel_compact(summary: dict) -> dict:
    if summary.get("protocol_complete") is not True:
        raise BuildError("NOVEL_PROTOCOL_NOT_COMPLETE")
    required = {
        "candidate_count": 2, "seed_count_per_candidate": 5,
        "models_per_seed": 5, "raw_model_count": 50,
        "raw_model_artifact_file_count": 100, "actual_completion_count": 10,
    }
    if any(summary.get(key) != value for key, value in required.items()):
        raise BuildError("NOVEL_PROTOCOL_COUNTS_INVALID")
    candidates = summary.get("candidates")
    if not isinstance(candidates, list) or len(candidates) != 2:
        raise BuildError("NOVEL_TWO_CANDIDATES_REQUIRED")
    return {
        "format": "campaign-diagnostics-novel-msa/2",
        "protocol_complete": True,
        "completed_seed_count": 10,
        "raw_model_count": 50,
        "candidates": [_slim_novel_candidate(item) for item in candidates],
        "development_set": True,
        "scientific_approved": False,
    }


def _find_input(receipt: dict, root: Path, seed: int) -> Path:
    for key in ("input_path", "input_yaml", "input"):
        value = receipt.get(key)
        if isinstance(value, str):
            return Path(value)
        if isinstance(value, dict):
            for subkey in ("path", "input_path", "yaml_path"):
                if isinstance(value.get(subkey), str):
                    return Path(value[subkey])
    for candidate in (root / f"seed-{seed}" / "input.yaml",
                      root / f"seed-{seed}" / "packet" / "input.yaml"):
        if candidate.is_file():
            return candidate
    raise BuildError("INPUT_YAML_MISSING:" + str(seed))


def _copy_calibration(protocol_root: Path, receipts: list[dict], stage: Path,
                      active_root: Path, provenance: dict) -> None:
    for receipt in receipts:
        seed = receipt["seed"]
        base = Path("evidence/calibration/raw") / f"seed-{seed}"
        receipt_path = protocol_root / f"seed-{seed}" / "receipt.json"
        provenance[f"calibration_receipt_{seed}"] = copy_verified(
            receipt_path, stage, base / "receipt.json", active_root)
        provenance[f"calibration_input_{seed}"] = copy_verified(
            _find_input(receipt, protocol_root, seed), stage, base / "input.yaml", active_root)
        for model in receipt["models"]:
            index = model["model_index"]
            provenance[f"calibration_{seed}_{index}_cif"] = copy_verified(
                Path(model["prediction_path"]), stage, base / f"model-{index}.cif",
                active_root, model["prediction_sha256"])
            provenance[f"calibration_{seed}_{index}_confidence"] = copy_verified(
                Path(model["confidence_path"]), stage, base / f"confidence-{index}.json",
                active_root, model["confidence_sha256"])


def _resolve_output_file(entry: dict, relative: str, active_root: Path) -> Path:
    rel = safe_pack_relative(relative)
    receipt_path = entry.get("receipt_path")
    if not isinstance(receipt_path, (str, Path)):
        raise BuildError("PREFLIGHT_RECEIPT_PATH_MISSING")
    output_root = Path(receipt_path).parent / "boltz_output"
    candidate = output_root / rel
    source = safe_source(candidate, active_root)
    if not is_relative_to(source, output_root.resolve(strict=True)):
        raise BuildError("PREFLIGHT_OUTPUT_PATH_INVALID:" + relative)
    return source


def _copy_baseline(preflight: list[dict], stage: Path, active_root: Path,
                   provenance: dict) -> None:
    if len(preflight) != 6:
        raise BuildError("BASELINE_EXACT_SIX_PLANS_REQUIRED")
    total = 0
    for candidate in preflight:
        identity = candidate.get("identity")
        cid = identity[0] if isinstance(identity, (list, tuple)) else identity.get("id")
        base = Path("evidence/baseline30") / str(cid)
        provenance[f"baseline_plan_{cid}"] = copy_verified(
            Path(candidate["plan_path"]), stage, base / "plan.json", active_root)
        seeds = candidate.get("seeds")
        if not isinstance(seeds, list) or len(seeds) != 5:
            raise BuildError("BASELINE_EXACT_FIVE_SEEDS_REQUIRED:" + str(cid))
        for seed_entry in seeds:
            seed = seed_entry["seed"]
            seed_base = base / f"seed-{seed}"
            provenance[f"baseline_receipt_{cid}_{seed}"] = copy_verified(
                Path(seed_entry["receipt_path"]), stage, seed_base / "receipt.json", active_root)
            provenance[f"baseline_input_{cid}_{seed}"] = copy_verified(
                Path(seed_entry["input_path"]), stage, seed_base / "input.yaml", active_root)
            output_hashes = seed_entry.get("output_hashes")
            if not isinstance(output_hashes, dict):
                raise BuildError("BASELINE_OUTPUT_HASHES_MISSING")
            prediction = safe_source(Path(seed_entry["prediction"]), active_root)
            pred_matches = [(path, digest) for path, digest in output_hashes.items()
                            if str(path).lower().endswith(".cif") and digest == sha256_file(prediction)]
            confidence = [(path, digest) for path, digest in output_hashes.items()
                          if re.fullmatch(r"confidence_[^/\\]+\.json", Path(str(path)).name)]
            if len(pred_matches) != 1 or len(confidence) != 1:
                raise BuildError("BASELINE_EXACT_CIF_AND_CONFIDENCE_REQUIRED")
            provenance[f"baseline_cif_{cid}_{seed}"] = copy_verified(
                prediction, stage, seed_base / "prediction.cif", active_root, pred_matches[0][1])
            conf_source = _resolve_output_file(seed_entry, confidence[0][0], active_root)
            provenance[f"baseline_confidence_{cid}_{seed}"] = copy_verified(
                conf_source, stage, seed_base / "confidence.json", active_root, confidence[0][1])
            total += 1
    if total != 30:
        raise BuildError("BASELINE_RAW_COUNT_NOT_30")


def _copy_novel(root: Path, summary: dict, stage: Path, active_root: Path,
                 provenance: dict) -> None:
    protocol_module = importlib.import_module("scripts.run_novel_msa_protocol")
    final_path = safe_source(root / "final-protocol.json", active_root)
    summary_path = safe_source(root / "protocol-summary.json", active_root)
    final = read_json(final_path)
    final_without_digest = copy.deepcopy(final)
    recorded_digest = final_without_digest.pop("protocol_digest", None)
    if (not isinstance(recorded_digest, str)
            or protocol_module.digest(final_without_digest) != recorded_digest):
        raise BuildError("NOVEL_FINAL_PROTOCOL_DIGEST_MISMATCH")
    if summary.get("protocol_digest") != recorded_digest:
        raise BuildError("NOVEL_SUMMARY_PROTOCOL_DIGEST_MISMATCH")

    final_candidates = final.get("candidates")
    if not isinstance(final_candidates, list) or len(final_candidates) != 2:
        raise BuildError("NOVEL_FINAL_PROTOCOL_CANDIDATES_INVALID")
    final_by_identity = {}
    for candidate in final_candidates:
        if not isinstance(candidate, dict):
            raise BuildError("NOVEL_FINAL_PROTOCOL_CANDIDATE_INVALID")
        identity = (candidate.get("candidate_id"), candidate.get("e3_type"))
        if identity in final_by_identity:
            raise BuildError("NOVEL_FINAL_PROTOCOL_CANDIDATE_DUPLICATE")
        final_by_identity[identity] = candidate

    provenance["novel_final-protocol.json"] = copy_verified(
        final_path, stage, "evidence/novel-msa/final-protocol.json", active_root)
    provenance["novel_protocol-summary.json"] = copy_verified(
        summary_path, stage, "evidence/novel-msa/protocol-summary.json", active_root)

    count = 0
    summary_identities = set()
    for candidate in summary["candidates"]:
        cid = candidate.get("candidate_id")
        e3_type = candidate.get("e3_type")
        if not isinstance(cid, str) or not cid or not isinstance(e3_type, str) or not e3_type:
            raise BuildError("NOVEL_CANDIDATE_IDENTITY_INVALID")
        identity = (cid, e3_type)
        if identity in summary_identities or identity not in final_by_identity:
            raise BuildError("NOVEL_FINAL_SUMMARY_CANDIDATE_MISMATCH")
        summary_identities.add(identity)
        final_candidate = final_by_identity[identity]

        plan_path = final_candidate.get("new_plan_path")
        plan_hash = final_candidate.get("new_plan_sha256")
        if (not isinstance(plan_path, str)
                or not isinstance(plan_hash, str)
                or re.fullmatch(r"[0-9a-fA-F]{64}", plan_hash) is None):
            raise BuildError("NOVEL_PLAN_PROVENANCE_INVALID:" + cid)
        candidate_dir_name = cid + "--" + e3_type
        provenance[f"novel_plan_{cid}_{e3_type}"] = copy_verified(
            safe_source(Path(plan_path), active_root), stage,
            Path("evidence/novel-msa/plans") / candidate_dir_name / "plan.json",
            active_root, plan_hash)

        receipts = candidate.get("seeds")
        if (not isinstance(receipts, list) or len(receipts) != 5
                or {receipt.get("seed") for receipt in receipts if isinstance(receipt, dict)} != set(SEEDS)):
            raise BuildError("NOVEL_EXACT_SEED_SET_REQUIRED:" + cid)
        candidate_root = root / candidate_dir_name
        for receipt in receipts:
            if not isinstance(receipt, dict):
                raise BuildError("NOVEL_RECEIPT_INVALID")
            seed = receipt["seed"]
            seed_root = candidate_root / f"seed-{seed}"
            receipt_path = safe_source(seed_root / "receipt.json", active_root)
            if read_json(receipt_path) != receipt:
                raise BuildError("NOVEL_SUMMARY_RECEIPT_MISMATCH:" + cid + ":" + str(seed))
            base = Path("evidence/novel-msa/raw") / candidate_dir_name / f"seed-{seed}"
            provenance[f"novel_receipt_{cid}_{seed}"] = copy_verified(
                receipt_path, stage, base / "receipt.json", active_root)
            provenance[f"novel_input_{cid}_{seed}"] = copy_verified(
                safe_source(seed_root / "novel_ternary.yaml", active_root), stage,
                base / "novel_ternary.yaml", active_root)
            for name in ("settings.json", "selection-receipt.json"):
                provenance[f"novel_{name}_{cid}_{seed}"] = copy_verified(
                    safe_source(seed_root / name, active_root), stage, base / name, active_root)
            for name in ("failure-receipt.json", "cancellation-receipt.json"):
                source = seed_root / name
                if source.exists():
                    provenance[f"novel_{name}_{cid}_{seed}"] = copy_verified(
                        safe_source(source, active_root), stage, base / name, active_root)

            models = receipt.get("models")
            if (not isinstance(models, list) or len(models) != 5
                    or {model.get("model_index") for model in models if isinstance(model, dict)} != set(range(5))):
                raise BuildError("NOVEL_EXACT_MODEL_INDEX_SET_REQUIRED:" + cid + ":" + str(seed))
            for model in models:
                pred_path = model.get("raw_prediction_path", model.get("prediction_path"))
                conf_path = model.get("raw_confidence_path", model.get("confidence_path"))
                pred_hash = model.get("raw_prediction_sha256", model.get("prediction_sha256"))
                conf_hash = model.get("raw_confidence_sha256", model.get("confidence_sha256"))
                if (not isinstance(pred_path, str) or not isinstance(conf_path, str)
                        or not isinstance(pred_hash, str)
                        or re.fullmatch(r"[0-9a-fA-F]{64}", pred_hash) is None
                        or not isinstance(conf_hash, str)
                        or re.fullmatch(r"[0-9a-fA-F]{64}", conf_hash) is None):
                    raise BuildError("NOVEL_MODEL_PROVENANCE_INVALID:" + cid + ":" + str(seed))
                index = model["model_index"]
                provenance[f"novel_cif_{cid}_{seed}_{index}"] = copy_verified(
                    safe_source(Path(pred_path), active_root), stage,
                    base / f"model-{index}.cif", active_root, pred_hash)
                provenance[f"novel_conf_{cid}_{seed}_{index}"] = copy_verified(
                    safe_source(Path(conf_path), active_root), stage,
                    base / f"confidence-{index}.json", active_root, conf_hash)
                count += 1
    if summary_identities != set(final_by_identity):
        raise BuildError("NOVEL_FINAL_SUMMARY_CANDIDATE_MISMATCH")
    if count != 50:
        raise BuildError("NOVEL_RAW_COUNT_NOT_50")


def _copy_exact_report(root: Path, files: list[str], destination: str, stage: Path,
                       active_root: Path, provenance: dict, label: str) -> None:
    for name in files:
        source = root / name
        if not source.is_file():
            raise BuildError("REQUIRED_REPORT_FILE_MISSING:" + str(source))
        provenance[label + "_" + name] = copy_verified(
            source, stage, Path(destination) / name, active_root)


def _copy_owned_sources(stage: Path, active_root: Path, provenance: dict) -> None:
    owned = {
        "run_calibration_protocol": "test_calibration_protocol.py",
        "inspect_ternary_ensemble": "test_ternary_ensemble.py",
        "inspect_core_state_robustness": "test_core_state_robustness.py",
        "run_novel_msa_protocol": "test_novel_msa_protocol.py",
        "diagnose_calibration_ranking": "test_calibration_ranking.py",
        "build_campaign_diagnostics": "test_campaign_diagnostics.py",
    }
    for stem, test_name in owned.items():
        script = active_root / "scripts" / f"{stem}.py"
        test = active_root / "tests" / test_name
        if not script.is_file():
            raise BuildError("OWNED_SCRIPT_MISSING:" + str(script))
        if not test.is_file():
            raise BuildError("OWNED_TEST_MISSING:" + str(test))
        provenance["owned_script_" + stem] = copy_verified(
            script, stage, Path("sources/scripts") / script.name, active_root)
        provenance["owned_test_" + stem] = copy_verified(
            test, stage, Path("sources/tests") / test.name, active_root)
    docs = ["science_diagnostics_replay_20261001.md"]
    for name in docs:
        source = active_root / "docs" / name
        if not source.is_file():
            raise BuildError("OWNED_DOC_MISSING:" + str(source))
        provenance["owned_doc_" + name] = copy_verified(
            source, stage, Path("sources/docs") / name, active_root)


def _html_value(value: Any) -> str:
    if value is None:
        return "미계산"
    if value is True:
        return "참"
    if value is False:
        return "거짓"
    if isinstance(value, float):
        return f"{value:.6f}".rstrip("0").rstrip(".")
    return str(value)


def _h(value: Any) -> str:
    return html.escape(_html_value(value), quote=True)


def _local_link_path(path: str) -> str | None:
    raw = path.replace("\\", "/")
    if re.fullmatch(r"[A-Za-z0-9._/-]+", raw) is None:
        return None
    try:
        relative = safe_pack_relative(raw)
    except BuildError:
        return None
    return relative.as_posix()


def _html_link(path: str, label: str | None = None) -> str:
    escaped_label = html.escape(label if label is not None else path, quote=True)
    safe_path = _local_link_path(path)
    if safe_path is None:
        return escaped_label
    return f'<a href="{html.escape(safe_path, quote=True)}">{escaped_label}</a>'


def _markdown_value(value: Any) -> str:
    text = _html_value(value)
    return text.replace("\\", "\\\\").replace("`", "\\`").replace("|", "\\|")


def _markdown_link(path: str, label: str | None = None) -> str:
    safe_label = _markdown_value(label if label is not None else path).replace("[", "\\[").replace("]", "\\]")
    safe_path = _local_link_path(path)
    return f"[{safe_label}]({safe_path})" if safe_path is not None else safe_label


def _finite_metric_values(metrics: Any, key: str) -> list[float]:
    if not isinstance(metrics, dict) or not isinstance(metrics.get(key), list):
        return []
    return [float(value) for value in metrics[key]
            if isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(float(value))]


def _iqr(values: list[float]) -> tuple[float, float] | None:
    if len(values) < 2:
        return None
    quartiles = statistics.quantiles(values, n=4, method="inclusive")
    return quartiles[0], quartiles[2]


def _old_pp_value(candidate: dict, topology: dict) -> float | None:
    candidate_id = candidate.get("candidate_id")
    e3_type = candidate.get("e3_type")
    for baseline in topology.get("candidates", []):
        if (isinstance(baseline, dict)
                and baseline.get("candidate_id") == candidate_id
                and baseline.get("e3_type") == e3_type):
            value = baseline.get("protein_protein_contact_jaccard", {}).get("median")
            if (isinstance(value, (int, float)) and not isinstance(value, bool)
                    and math.isfinite(float(value))):
                return float(value)
            return None
    return None


def _raw_links(compact: dict) -> list[str]:
    links: list[str] = []
    calibration = compact.get("known_calibration", {})
    for row in calibration.get("per_seed_selected", []):
        seed = row.get("seed")
        base = f"evidence/calibration/raw/seed-{seed}"
        links.extend([f"{base}/receipt.json", f"{base}/input.yaml"])
        for index in range(5):
            links.extend([f"{base}/model-{index}.cif", f"{base}/confidence-{index}.json"])
    topology = compact.get("topology_single_sequence", {})
    for candidate in topology.get("candidates", []):
        cid = candidate.get("candidate_id")
        links.append(f"evidence/baseline30/{cid}/plan.json")
        for seed in topology.get("seeds", SEEDS):
            base = f"evidence/baseline30/{cid}/seed-{seed}"
            links.extend([f"{base}/receipt.json", f"{base}/input.yaml",
                          f"{base}/prediction.cif", f"{base}/confidence.json"])
    novel = compact.get("novel_msa")
    if isinstance(novel, dict):
        links.extend(["evidence/novel-msa/final-protocol.json",
                      "evidence/novel-msa/protocol-summary.json"])
        for candidate in novel.get("candidates", []):
            cid = candidate.get("candidate_id")
            e3 = candidate.get("e3_type")
            directory = f"{cid}--{e3}"
            links.append(f"evidence/novel-msa/plans/{directory}/plan.json")
            for selected in candidate.get("selected_top1_per_seed", []):
                seed = selected.get("seed")
                base = f"evidence/novel-msa/raw/{directory}/seed-{seed}"
                links.extend([f"{base}/receipt.json", f"{base}/novel_ternary.yaml",
                              f"{base}/settings.json", f"{base}/selection-receipt.json"])
                for index in range(5):
                    links.extend([f"{base}/model-{index}.cif",
                                  f"{base}/confidence-{index}.json"])
    return links


def render_index(compact: dict) -> str:
    calibration = compact.get("known_calibration", {})
    ranking = compact.get("ranking_development", {})
    core = compact.get("core_state", {})
    topology = compact.get("topology_single_sequence", {})
    novel = compact.get("novel_msa")

    calibration_rows = "".join(
        "<tr>"
        f"<td>{_h(row.get('seed'))}</td>"
        f"<td>{_h(row.get('selected_model_index'))}</td>"
        f"<td>{_h(row.get('selected_e3_CA_RMSD_after_target_alignment_A'))}</td>"
        f"<td>{_h(row.get('selected_pass'))}</td>"
        "</tr>"
        for row in calibration.get("per_seed_selected", [])
    )
    ranking_rows = "".join(
        "<tr>"
        f"<td>{_h(name)}</td><td>{_h(result.get('selected_pass_count'))}</td>"
        f"<td>{_h(result.get('valid'))}</td><td>{_h(result.get('error'))}</td>"
        "</tr>"
        for name, result in ranking.get("rules", {}).items()
    )
    topology_rows = "".join(
        "<tr>"
        f"<td>{_h(candidate.get('priority_diagnostic'))}</td>"
        f"<td><code>{_h(candidate.get('candidate_id'))}</code></td>"
        f"<td>{_h(candidate.get('e3_type'))}</td>"
        f"<td>{_h(candidate.get('protein_protein_contact_jaccard', {}).get('median'))}</td>"
        f"<td>{_h(candidate.get('target_warhead_site_jaccard', {}).get('median'))}</td>"
        f"<td>{_h(candidate.get('e3_recruiter_site_jaccard', {}).get('median'))}</td>"
        f"<td><code>{_h(candidate.get('clashes_under2_A_all_five'))}</code></td>"
        f"<td><code>{_h(candidate.get('vdw_clashes_all_five'))}</code></td>"
        "</tr>"
        for candidate in topology.get("candidates", [])
    )

    if isinstance(novel, dict):
        novel_parts = [
            '<section id="novel"><h2>선택적 신규 공동 프로토콜</h2>',
            '<div class="notice"><strong>정량 진단이며 승인이 아닙니다.</strong> '
            '기존 30개 결과는 단일 서열/r3/d1이고 신규 50개 원시 모델은 고정 MSA/r10/d5입니다. '
            '여러 설정이 함께 바뀐 공동 프로토콜 비교이므로 관찰된 차이를 MSA 단독 효과로 해석할 수 없습니다.</div>',
            f"<p>완료 seed 결과: <strong>{_h(novel.get('completed_seed_count'))}</strong>; "
            f"원시 모델: <strong>{_h(novel.get('raw_model_count'))}</strong>. "
            '후보별 5개 seed에서 각각 5개 모델을 계산하고, 기하 평가 전에 고정된 최고 confidence 규칙으로 top-1을 선택했습니다. '
            '따라서 두 후보 합계는 50개 원시 모델에서 선택한 10개이며, 50개 독립 seed도 아니고 총 5개 선택도 아닙니다.</p>'
        ]
        for candidate in novel.get("candidates", []):
            cid = candidate.get("candidate_id")
            e3 = candidate.get("e3_type")
            selections = candidate.get("selected_top1_per_seed", [])
            novel_parts.append(f"<h3><code>{_h(cid)}</code> · {_h(e3)}</h3>")
            novel_parts.append("<p>완료 5 / 원시 25 / top-1 선택 5.</p>")
            novel_parts.append('<div class="table-wrap"><table><caption>선택 인덱스 5쌍</caption>'
                               '<thead><tr><th>seed</th><th>선택 model index</th></tr></thead><tbody>')
            for row in selections:
                novel_parts.append(f"<tr><td>{_h(row.get('seed'))}</td>"
                                   f"<td>{_h(row.get('selected_model_index'))}</td></tr>")
            novel_parts.append("</tbody></table></div>")

            metrics = candidate.get("all_five_quantitative_metrics", {})
            novel_parts.append('<div class="table-wrap"><table><caption>제공된 선택 정량 배열 '
                               '(포함 분위수법)</caption><thead><tr><th>지표</th><th>n</th>'
                               '<th>중앙값</th><th>Q1</th><th>Q3</th><th>계산 IQR</th>'
                               '<th>저장 IQR</th><th>원래 maximum</th></tr></thead><tbody>')
            metric_rows = (
                ("finite_selected_ipTM_values", "ipTM_IQR", "ipTM_IQR_maximum"),
                ("finite_selected_endpoint_values_A", "endpoint_IQR_A", "endpoint_IQR_maximum_A"),
            )
            for values_key, stored_key, maximum_key in metric_rows:
                values = _finite_metric_values(metrics, values_key)
                interval = _iqr(values)
                q1, q3 = interval if interval is not None else (None, None)
                calculated_iqr = q3 - q1 if q1 is not None and q3 is not None else None
                median = statistics.median(values) if values else None
                novel_parts.append(
                    f"<tr><td><code>{_h(values_key)}</code></td><td>{_h(len(values))}</td>"
                    f"<td>{_h(median)}</td><td>{_h(q1)}</td><td>{_h(q3)}</td>"
                    f"<td>{_h(calculated_iqr)}</td><td>{_h(metrics.get(stored_key))}</td>"
                    f"<td>{_h(metrics.get(maximum_key))}</td></tr>")
            novel_parts.append("</tbody></table></div>")
            novel_parts.append(
                '<p><code>quantitative_criteria_met</code>: '
                f"<strong>{_h(metrics.get('quantitative_criteria_met'))}</strong> — "
                '<strong>quantitative only, not scientific approval</strong>. '
                '<code>automatic_approval</code>: '
                f"<strong>{_h(metrics.get('automatic_approval'))}</strong>; "
                '<code>qualitative_topology_or_clash_acceptance_used</code>: '
                f"<strong>{_h(metrics.get('qualitative_topology_or_clash_acceptance_used'))}</strong>; "
                '<code>positive_target_warhead_and_e3_recruiter_contacts</code>: '
                f"<code>{_h(metrics.get('positive_target_warhead_and_e3_recruiter_contacts'))}</code>.</p>")

            diagnostics = candidate.get("geometric_diagnostics", {})
            new_pp = diagnostics.get("protein_protein_jaccard_stats", {}).get("median")
            old_pp = _old_pp_value(candidate, topology)
            delta = (float(new_pp) - old_pp
                     if old_pp is not None and isinstance(new_pp, (int, float))
                     and not isinstance(new_pp, bool) and math.isfinite(float(new_pp)) else None)
            novel_parts.append('<div class="table-wrap"><table><caption>관찰된 PP Jaccard 비교</caption>'
                               '<thead><tr><th>기존 중앙값</th><th>신규 중앙값</th><th>신규−기존</th>'
                               '</tr></thead><tbody><tr>'
                               f"<td>{_h(old_pp)}</td><td>{_h(new_pp)}</td><td>{_h(delta)}</td>"
                               '</tr></tbody></table></div>')
            novel_parts.append('<p class="muted">위 delta는 제공된 값의 관찰 차이이며 인과 효과가 아닙니다. '
                               f"2 Å 미만 충돌 배열: <code>{_h(diagnostics.get('raw_under_2A_pair_count_all_five'))}</code>; "
                               f"vdW 충돌 배열: <code>{_h(diagnostics.get('vdw_pair_count_all_five'))}</code>.</p>")
        novel_parts.append("</section>")
        novel_html = "".join(novel_parts)
    else:
        novel_html = (
            '<section id="novel"><h2>선택적 신규 공동 프로토콜</h2>'
            '<div class="notice neutral"><strong>계산되지 않음.</strong> 이 팩에는 완료가 검증된 '
            '신규 MSA/r10/d5 프로토콜이 없습니다. 이를 성공, 실패 또는 승인으로 해석하지 않습니다.</div></section>'
        )

    source_rows = "".join(
        "<tr>"
        f"<td>{_h(name)}</td>"
        f"<td>{_html_link(str(record.get('portable_pack_relative_path')), str(record.get('portable_pack_relative_path')))}</td>"
        f"<td><code>{_h(record.get('source_sha256'))}</code></td>"
        f"<td><code>{_h(record.get('active_root_relative_path'))}</code></td>"
        "</tr>"
        for name, record in compact.get("sources", {}).items()
    )
    provenance = compact.get("source_provenance", {})
    raw_links = "".join(f"<li>{_html_link(path)}</li>" for path in _raw_links(compact))

    return f'''<!doctype html>
<html lang="ko">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>TPD Navigator · 과학 계산 검토</title>
<style>
:root{{--bg:#f4f7fb;--panel:#fff;--ink:#172033;--muted:#58657a;--line:#cbd5e1;--accent:#174ea6;--warn:#8a4b08;--ok:#17633a}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--bg);color:var(--ink);font:16px/1.6 system-ui,-apple-system,"Segoe UI",sans-serif}}
a{{color:var(--accent);text-underline-offset:.18em}} a:focus{{outline:3px solid #f5b700;outline-offset:2px}}
.skip{{position:absolute;left:-9999px}} .skip:focus{{left:1rem;top:1rem;background:#fff;padding:.6rem;z-index:2}}
header,main,footer{{width:min(1180px,calc(100% - 2rem));margin:auto}} header{{padding:2.5rem 0 1rem}}
h1{{font-size:clamp(1.8rem,5vw,3rem);line-height:1.15;margin:.2rem 0}} h2{{margin-top:0}} h3{{margin-top:2rem}}
section{{background:var(--panel);border:1px solid var(--line);border-radius:14px;padding:clamp(1rem,3vw,1.8rem);margin:1rem 0;box-shadow:0 4px 16px #1831530d}}
.badge{{display:inline-block;border:1px solid #b45309;background:#fff7ed;color:#7c2d12;border-radius:999px;padding:.18rem .65rem;font-weight:700}}
.notice{{border-left:5px solid var(--warn);background:#fff8ed;padding:.8rem 1rem;margin:1rem 0}} .neutral{{border-color:#64748b;background:#f8fafc}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:.8rem}} .card{{border:1px solid var(--line);border-radius:10px;padding:1rem}}
.card strong{{display:block;font-size:1.5rem}} .table-wrap{{overflow-x:auto;margin:1rem 0}} table{{border-collapse:collapse;width:100%;min-width:650px}}
caption{{font-weight:700;text-align:left;padding:.5rem 0}} th,td{{border:1px solid var(--line);padding:.55rem .65rem;text-align:left;vertical-align:top}} th{{background:#edf2f7}}
code{{overflow-wrap:anywhere}} details{{margin:.8rem 0}} .muted{{color:var(--muted)}} footer{{padding:1rem 0 3rem}}
@media(max-width:600px){{header,main,footer{{width:min(100% - 1rem,1180px)}} section{{border-radius:9px}}}}
</style>
</head>
<body>
<a class="skip" href="#main">본문으로 건너뛰기</a>
<header><span class="badge">과학 승인: {_h(compact.get('scientific_approved'))}</span>
<h1>TPD Navigator · 과학 계산 검토</h1>
<p>오프라인 검토용 읽기 전용 증거 팩입니다. 원본 과학 상태나 승인 필드를 변경하지 않습니다.</p></header>
<main id="main">
<section id="scope"><h2>해석 범위</h2>
<div class="notice"><strong>6BOY/CRBN은 BRD4–dBET6 알려진 보정 사례입니다.</strong>
RCSB 6BOY 구조를 이용한 벤치마크이며,
신규 SMARCA2 후보의 참조 정확도를 입증하지 않습니다. 알려진 입력에는 DDB1 제외가 기록되어 있으나,
그 제외가 오차를 일으켰다는 인과 증거는 아닙니다.</div>
<p>기존 30개 결과(단일 서열/r3/d1)와 선택적 신규 50개 결과(고정 MSA/r10/d5)는
<strong>공동 프로토콜 비교</strong>입니다. MSA만이 변화를 일으켰다고 주장할 수 없습니다.
모든 과학 승인 값은 거짓이며, native evidence는 상위 캠페인 소유로 이 팩에 포함되지 않았습니다.</p></section>
<section id="calibration"><h2>알려진 6BOY/CRBN 보정</h2>
<div class="cards"><div class="card"><span>사전 선언 top-1 선택 통과</span><strong>{_h(calibration.get('selected_pass_count'))} / {_h(calibration.get('seed_count'))}</strong></div>
<div class="card"><span>원시 모델 네 기준 통과</span><strong>{_h(calibration.get('raw_all_four_pass_count'))} / {_h(calibration.get('raw_model_count'))}</strong></div></div>
<p>25개 원시 모델 중 16개 통과와 사전 선언 선택 5개 중 2개 통과는 서로 다른 분모입니다.
선택은 기하 평가 전에 고정된 최고 confidence 규칙을 사용했습니다. 원래 네 cut-off는 보존되었습니다.</p>
<div class="table-wrap"><table><caption>seed별 사전 선언 선택 모델과 E3 Cα RMSD</caption><thead><tr><th>seed</th><th>model index</th><th>E3 Cα RMSD (Å)</th><th>네 기준 통과</th></tr></thead><tbody>{calibration_rows}</tbody></table></div>
<h3>랭킹 규칙 진단</h3><p>아래는 사후 개발 분석이며 원래 2/5 baseline을 바꾸지 않습니다.</p>
<div class="table-wrap"><table><thead><tr><th>규칙</th><th>선택 통과 수 / 5</th><th>유효</th><th>오류</th></tr></thead><tbody>{ranking_rows}</tbody></table></div></section>
<section id="core"><h2>Core-state 진단</h2>
<div class="cards"><div class="card"><span>전체 행</span><strong>{_h(core.get('row_count'))}</strong></div>
<div class="card"><span>contact 양성</span><strong>{_h(core.get('required_contact_true_count'))}</strong></div>
<div class="card"><span>contact 음성</span><strong>{_h(core.get('required_contact_false_count'))}</strong></div>
<div class="card"><span>failure / null</span><strong>{_h(core.get('failure_count'))} / {_h(core.get('failure_null_count'))}</strong></div></div>
<p>단백질 pKa 로그는 원시 {_h(core.get('protein_pKa_raw_row_count'))}행, 고유 chain/residue/sequence/pKa 조합 {_h(core.get('protein_pKa_unique_chain_residue_sequence_pKa_count'))}개입니다.
리간드 pKa 및 population은 추정되지 않았고 microstate gate는 <strong>{_h(core.get('microstate_gate'))}</strong>입니다.
Native evidence는 상위 캠페인에 속하며 여기에는 포함되지 않습니다.</p></section>
<section id="baseline"><h2>기존 6후보 · 30개 단일 서열 결과</h2>
<p>seed 순서는 <code>{_h(topology.get('seeds'))}</code>입니다. Jaccard, 충돌 수, priority는 설명적 진단입니다.
과학적 보정의 원래 네 cut-off를 대체하지 않습니다.</p>
<div class="table-wrap"><table><caption>후보별 5-seed 진단</caption><thead><tr><th>priority 진단</th><th>ID</th><th>E3</th><th>PP Jaccard 중앙값</th><th>target site 중앙값</th><th>E3 site 중앙값</th><th>&lt;2 Å 충돌 5개</th><th>vdW 충돌 5개</th></tr></thead><tbody>{topology_rows}</tbody></table></div></section>
{novel_html}
<section id="sources"><h2>원본과 출처</h2>
<p>원본 보고서는 byte-for-byte 복사되었습니다. JSON 안의 원래 로컬 경로 문자열은 출처 기록일 뿐이며 자동으로 읽거나 가져오지 않습니다.</p>
<div class="table-wrap"><table><caption>compact.sources의 정확한 SHA-256 및 경로</caption><thead><tr><th>항목</th><th>팩 경로</th><th>원본 SHA-256</th><th>원래 workspace 경로</th></tr></thead><tbody>{source_rows}</tbody></table></div>
<ul><li>{_html_link(str(provenance.get('portable_pack_relative_path', 'source-provenance.json')))} — SHA-256 <code>{_h(provenance.get('sha256'))}</code></li>
<li>{_html_link('output-manifest.json')} — 이 페이지를 포함한 출력 파일 SHA-256 목록</li>
<li>{_html_link('sources/docs/science_diagnostics_replay_20261001.md')}</li></ul>
<details><summary>포함된 모든 원시 artifact 경로 링크</summary><ul>{raw_links}</ul></details></section>
</main><footer><p>로컬 파일로 열 수 있으며 GPU, 네트워크, 원격 스크립트, 웹폰트가 필요하지 않습니다.</p></footer>
</body></html>
'''


def render_readme(compact: dict) -> str:
    calibration = compact.get("known_calibration", {})
    core = compact.get("core_state", {})
    topology = compact.get("topology_single_sequence", {})
    novel = compact.get("novel_msa")
    optional_flag = (
        " --novel-protocol-root .localdata/final-sprint-20261001/science/novel-msa-r10-d5"
        if isinstance(novel, dict) else ""
    )
    sources = "\n".join(
        f"- {_markdown_link(str(record.get('portable_pack_relative_path')), str(name))} — "
        f"SHA-256 `{_markdown_value(record.get('source_sha256'))}`"
        for name, record in compact.get("sources", {}).items()
    )
    if isinstance(novel, dict):
        novel_text = (
            f"신규 공동 프로토콜은 완료 seed {_markdown_value(novel.get('completed_seed_count'))}개, "
            f"원시 모델 {_markdown_value(novel.get('raw_model_count'))}개를 포함합니다. "
            "후보별 25개 원시 모델에서 5개 top-1을 선택하여 합계 50개 원시 모델에서 10개를 선택했습니다. "
            "이는 50개 독립 seed가 아니며 MSA 단독 인과 효과를 보여주지 않습니다."
        )
    else:
        novel_text = (
            "완료가 검증된 신규 MSA/r10/d5 결과는 이 팩에서 계산되지 않았습니다. "
            "이는 성공, 실패 또는 승인을 뜻하지 않습니다."
        )
    return f"""# TPD Navigator · 과학 계산 검토

읽기 시작: [오프라인 HTML 인덱스](index.html)

이 팩은 GPU와 네트워크 없이 읽을 수 있는 읽기 전용 증거 사본입니다. 모든 과학 승인 값은 `false`이며 빌더는 승인 상태를 쓰지 않습니다.

## 핵심 분모와 해석

- 알려진 6BOY/CRBN 사례는 [RCSB 6BOY](https://www.rcsb.org/structure/6BOY)의 BRD4–dBET6 보정 벤치마크이며, 신규 SMARCA2 참조 정확도 입증이 아닙니다.
- 사전 선언 선택 통과는 `{_markdown_value(calibration.get('selected_pass_count'))}/{_markdown_value(calibration.get('seed_count'))}`, 원시 네 기준 통과는 `{_markdown_value(calibration.get('raw_all_four_pass_count'))}/{_markdown_value(calibration.get('raw_model_count'))}`입니다. 두 분모는 다릅니다.
- Core-state는 양성 `{_markdown_value(core.get('required_contact_true_count'))}`, 음성 `{_markdown_value(core.get('required_contact_false_count'))}`, 전체 `{_markdown_value(core.get('row_count'))}`행입니다.
- 기존 topology는 후보 `{_markdown_value(topology.get('candidate_count'))}`개, 원시 seed 결과 `{_markdown_value(topology.get('raw_seed_result_count'))}`개입니다.
- {novel_text}
- 기존 단일 서열/r3/d1과 신규 고정 MSA/r10/d5는 공동 프로토콜 비교입니다. MSA만이 차이를 일으켰다고 주장할 수 없습니다.
- 선택은 기하 평가 전에 고정된 최고 confidence 규칙으로 수행됩니다. Jaccard와 충돌 정의는 설명적이며 알려진 보정의 원래 네 cut-off는 보존됩니다.
- Native evidence는 상위 캠페인 소유이며 포함되지 않았습니다. 알려진 입력의 DDB1 제외 기록은 오차 원인 증명이 아닙니다.

## 출처 JSON

- [compact-campaign-diagnostics.json](compact-campaign-diagnostics.json)
- [source-provenance.json](source-provenance.json)
- [output-manifest.json](output-manifest.json)
- [재현 및 범위 문서](sources/docs/science_diagnostics_replay_20261001.md)
{sources}

원본 JSON에 남은 로컬 경로 문자열은 provenance 전용이며 자동 fetch 위치가 아닙니다.

## 출력 manifest 검증

팩 디렉터리로 이동한 뒤 번들된 표준 라이브러리 helper를 사용합니다.

```powershell
python -B -X utf8 -c "from pathlib import Path; from sources.scripts.build_campaign_diagnostics import verify_manifest; verify_manifest(Path('.')); print('manifest verified')"
```

`-B`는 캐시 파일 생성을 막아 동결된 manifest의 정확한 파일 집합을 유지하며, 이전 실행으로 캐시가 생겼다면 깨끗한 원본 팩을 사용하세요.

## 원래 workspace에서 fresh v2 재빌드

```powershell
.venv/Scripts/python.exe -X utf8 scripts/build_campaign_diagnostics.py --science-root .localdata/final-sprint-20261001/science{optional_flag} --output .localdata/final-sprint-20261001/science-evidence-v2-REPLAY-NEW
```

출력은 정상 UUID staging 디렉터리에서 검증된 뒤 `os.replace`로 게시되며 상위 디렉터리 ACL 상속을 유지합니다. 기존 출력과 원본 source는 변경하지 않습니다.

## 선택적 GPU 재실행

휴대용 팩만으로는 실행할 수 없습니다. 원래 검증 workspace, frozen weights/checkpoint, 처리된 MSA/cache 및 관리되는 Boltz 환경이 선행 조건이며 팩에는 포함되지 않습니다. 공식 설정 의미는 [Boltz prediction 문서](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md)를 함께 확인하십시오.

```powershell
.venv/Scripts/python.exe -X utf8 scripts/run_novel_msa_protocol.py --batch-root .localdata/expert-closure-20260930/priority-batch --output .localdata/novel-msa-r10-d5-rerun-NEW --boltz-executable .venv-boltz/Scripts/boltz.exe --checkpoint .localdata/boltz-cache/boltz2_conf.ckpt --cache .localdata/boltz-cache --timeout-seconds 1800 --msa-timeout-seconds 600
```

프로토콜은 recycling 10, sampling 200, diffusion samples 5, seeds `23 41 61 79 97`, max MSA 256, potentials OFF를 고정합니다. 이들은 CLI step knob가 아닙니다. 재실행 결과 역시 과학 승인을 자동으로 쓰지 않습니다.
"""


def write_text_exclusive(path: Path, value: str) -> None:
    with path.open("xb") as stream:
        stream.write(value.encode("utf-8"))


def make_manifest(stage: Path) -> dict:
    files = []
    for path in sorted(stage.rglob("*"), key=lambda item: item.as_posix()):
        if path.is_symlink():
            raise BuildError("OUTPUT_SYMLINK_REJECTED:" + str(path))
        if not path.is_file() or path.relative_to(stage).as_posix() == "output-manifest.json":
            continue
        if path.suffix.lower() not in ALLOWED_SUFFIXES:
            raise BuildError("OUTPUT_TYPE_NOT_ALLOWED:" + str(path))
        relative = path.relative_to(stage).as_posix()
        files.append({"path": relative, "bytes": path.stat().st_size,
                      "sha256": sha256_file(path)})
    return {"format": "campaign-diagnostics-output-manifest/1",
            "manifest_excludes_itself": True, "files": files}


def verify_manifest(root: Path) -> None:
    manifest = read_json(root / "output-manifest.json")
    expected = {item["path"]: item for item in manifest.get("files", [])}
    actual = {}
    for path in root.rglob("*"):
        if path.is_symlink():
            raise BuildError("OUTPUT_SYMLINK_REJECTED:" + str(path))
        if path.is_file() and path.relative_to(root).as_posix() != "output-manifest.json":
            relative = path.relative_to(root).as_posix()
            actual[relative] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    if set(actual) != set(expected):
        raise BuildError("OUTPUT_MANIFEST_FILE_SET_MISMATCH")
    for name, observed in actual.items():
        if observed["bytes"] != expected[name]["bytes"] or observed["sha256"] != expected[name]["sha256"]:
            raise BuildError("OUTPUT_MANIFEST_HASH_MISMATCH:" + name)


def _import_functions() -> tuple[Any, Any, Any]:
    ranking = importlib.import_module("scripts.diagnose_calibration_ranking")
    calibration = importlib.import_module("scripts.run_calibration_protocol")
    ensemble = importlib.import_module("scripts.inspect_ternary_ensemble")
    return ranking.load_verified_protocol, calibration.metrics_pass, ensemble._preflight


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--science-root", "-S", required=True)
    value.add_argument("--baseline-batch-root",
                       default=".localdata/expert-closure-20260930/priority-batch")
    value.add_argument("--output", required=True)
    value.add_argument("--novel-protocol-root")
    return value


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    active_root = Path.cwd().resolve()
    science = Path(args.science_root).expanduser().resolve()
    baseline = Path(args.baseline_batch_root).expanduser().resolve()
    output_requested = Path(args.output).expanduser()
    output_absolute = output_requested if output_requested.is_absolute() else Path.cwd() / output_requested
    try:
        reject_symlink_components(output_absolute.absolute())
    except BuildError as exc:
        print("campaign diagnostics build error: BuildError: " + str(exc), file=sys.stderr)
        return 2
    output = output_absolute.resolve()
    novel = Path(args.novel_protocol_root).expanduser().resolve() if args.novel_protocol_root else None
    stage: Path | None = None
    try:
        if output.exists():
            raise BuildError("OUTPUT_ALREADY_EXISTS")
        if not is_relative_to(output, active_root):
            raise BuildError("OUTPUT_OUTSIDE_ACTIVE_ROOT")
        roots = [science, baseline] + ([novel] if novel else [])
        for root in roots:
            if root is None or not root.is_dir() or root.is_symlink():
                raise BuildError("INPUT_ROOT_DIRECTORY_REQUIRED:" + str(root))
            if not is_relative_to(root, active_root):
                raise BuildError("INPUT_ROOT_OUTSIDE_ACTIVE_ROOT:" + str(root))
            if is_relative_to(output, root) or is_relative_to(root, output):
                raise BuildError("OUTPUT_INPUT_SUBTREE_OVERLAP:" + str(root))
        if not output.parent.exists() or output.parent.is_symlink():
            raise BuildError("OUTPUT_PARENT_DIRECTORY_REQUIRED")

        load_protocol, metrics_pass, baseline_preflight = _import_functions()
        calibration_root = science / "calibration-r10-d5"
        core_root = science / "core-state-robustness"
        topology_root = science / "ternary-ensemble"
        ranking_root = science / "calibration-ranking-diagnosis"

        core_manifest = read_json(core_root / "manifest.json")
        verify_artifact_manifest(core_root, core_manifest, active_root)
        topology_manifest = read_json(topology_root / "output-manifest.json")
        verify_artifact_manifest(topology_root, topology_manifest, active_root)
        plan, receipts = load_protocol(calibration_root)
        preflight = baseline_preflight(baseline)
        if isinstance(preflight, dict):
            preflight = list(preflight.values())
        if not isinstance(preflight, list):
            raise BuildError("BASELINE_PREFLIGHT_SCHEMA_INVALID")

        ranking_provenance = read_json(ranking_root / "source-provenance.json")
        validate_provenance(ranking_provenance, active_root)
        diagnosis = read_json(ranking_root / "diagnosis.json")
        comparison_path = science / "candidate-comparison-derived.json"
        comparison = read_json(comparison_path) if comparison_path.is_file() else None
        topology_summary_hash = sha256_file(topology_root / "summary.json")
        if comparison is not None and comparison.get("source_sha256") != topology_summary_hash:
            raise BuildError("CANDIDATE_COMPARISON_SOURCE_HASH_MISMATCH")
        novel_summary = None
        if novel is not None:
            summary_path = novel / "protocol-summary.json"
            if not summary_path.is_file():
                raise BuildError("NOVEL_COMPLETED_SUMMARY_REQUIRED")
            novel_summary = read_json(summary_path)
            build_novel_compact(novel_summary)

        stage = output.parent / ("." + output.name + ".staging-" + uuid.uuid4().hex)
        stage.mkdir(parents=False, exist_ok=False)
        provenance: dict[str, Any] = {}
        _copy_exact_report(core_root, ["summary.json", "manifest.json", "all-states.sdf", "index.html"],
                           "evidence/core", stage, active_root, provenance, "core")
        _copy_exact_report(topology_root, ["summary.json", "output-manifest.json", "report.html"],
                           "evidence/topology", stage, active_root, provenance, "topology")
        _copy_exact_report(calibration_root, ["protocol-summary.json", "plan.json"],
                           "evidence/calibration", stage, active_root, provenance, "calibration")
        _copy_exact_report(ranking_root,
                           ["diagnosis.json", "ranking-proposals.json", "source-provenance.json", "diagnosis.html"],
                           "evidence/ranking", stage, active_root, provenance, "ranking")
        for name in ("candidate-comparison-derived.json",
                     "baseline-e3-fold-diagnostic.json",
                     "calibration-protein-interface-derived.json"):
            source = science / name
            if not source.is_file():
                raise BuildError("REQUIRED_DERIVED_REPORT_MISSING:" + str(source))
            provenance["derived_" + name] = copy_verified(
                source, stage, Path("evidence/derived") / name, active_root)
        _copy_calibration(calibration_root, receipts, stage, active_root, provenance)
        _copy_baseline(preflight, stage, active_root, provenance)
        if novel is not None and novel_summary is not None:
            _copy_novel(novel, novel_summary, stage, active_root, provenance)
        _copy_owned_sources(stage, active_root, provenance)

        core_summary = read_json(core_root / "summary.json")
        topology_summary = read_json(topology_root / "summary.json")
        provenance_path = stage / "source-provenance.json"
        write_json_exclusive(provenance_path, provenance)
        concise_source_keys = [
            "core_summary.json", "topology_summary.json",
            "calibration_protocol-summary.json", "calibration_plan.json",
            "ranking_diagnosis.json",
        ]
        if novel_summary is not None:
            concise_source_keys += ["novel_protocol-summary.json", "novel_final-protocol.json"]
        concise_sources = {key: provenance[key] for key in concise_source_keys if key in provenance}
        compact = {
            "format": "compact-campaign-diagnostics/20261001.1",
            "portable_evidence": True,
            "scientific_approved": False,
            "human_review_performed": False,
            "native_docking_evidence_included": False,
            "native_evidence_owner": "parent campaign",
            "known_calibration": build_calibration_compact(plan, receipts, metrics_pass),
            "ranking_development": build_ranking_compact(diagnosis),
            "core_state": build_core_compact(core_summary),
            "topology_single_sequence": build_topology_compact(topology_summary),
            "sources": concise_sources,
            "source_provenance": {
                "portable_pack_relative_path": "source-provenance.json",
                "sha256": sha256_file(provenance_path),
            },
            "offline_reader_formats": ["json", "sdf", "cif", "html", "md"],
            "gpu_rerun_requires_original_validated_workspace": True,
            "frozen_msa_cache_and_checkpoint_intentionally_not_embedded": True,
        }
        if novel_summary is not None:
            compact["format"] = "compact-campaign-diagnostics/20261001.2"
            compact["novel_msa"] = build_novel_compact(novel_summary)
        write_compact(stage / "compact-campaign-diagnostics.json", compact)
        write_text_exclusive(stage / "index.html", render_index(compact))
        write_text_exclusive(stage / "README.md", render_readme(compact))
        write_json_exclusive(stage / "output-manifest.json", make_manifest(stage))
        verify_manifest(stage)
        os.replace(stage, output)
        stage = None
        return 0
    except (OSError, BuildError, ImportError, KeyError, TypeError) as exc:
        print("campaign diagnostics build error: " + type(exc).__name__ + ": " + str(exc),
              file=sys.stderr)
        return 2
    finally:
        if stage is not None and stage.exists():
            expected_prefix = "." + output.name + ".staging-"
            try:
                parent = output.parent.resolve(strict=True)
                candidate = stage.resolve(strict=True)
                if candidate.parent == parent and candidate.name.startswith(expected_prefix):
                    shutil.rmtree(candidate, ignore_errors=True)
            except (OSError, BuildError):
                pass


if __name__ == "__main__":
    raise SystemExit(main())
