"""Strict aggregation of source-bound known-case Boltz receipts.

The resulting document is a protocol proposal, never a scientific approval.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path, PurePath
from typing import Any


KNOWN_BASE_SUFFIX = PurePath("cases/acceptance_sources/crbn_calibration.json")
KNOWN_BASE_SHA256 = "225a4264dc499b87816b11b3c9293ae5f35b8a4392e4dceb40f3652181af6b08"


class BenchmarkDistributionError(ValueError):
    pass


def _check(condition: bool, code: str) -> None:
    if not condition:
        raise BenchmarkDistributionError(code)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load(value: Any) -> Any:
    if isinstance(value, (str, Path)):
        return json.loads(Path(value).read_text(encoding="utf-8"))
    return value


def _path_has_symlink(path: Path) -> bool:
    absolute = path.absolute()
    current = Path(absolute.anchor)
    for part in absolute.parts[1:]:
        current = current / part
        if current.is_symlink():
            return True
    return False


def _source_file(root: Path, reference: Any, expected: PurePath) -> Path:
    _check(isinstance(reference, str) and reference != "", "SOURCE_REFERENCE_REQUIRED")
    normalized = reference.replace("\\", "/")
    relative = PurePath(normalized)
    _check(not relative.is_absolute() and ".." not in relative.parts,
           "SOURCE_REFERENCE_TRAVERSAL")
    _check(relative == expected, "SOURCE_REFERENCE_UNREGISTERED")
    candidate = root.joinpath(*relative.parts)
    _check(candidate.is_file() and not _path_has_symlink(candidate),
           "SOURCE_REGULAR_FILE_REQUIRED")
    resolved_root = root.resolve()
    resolved = candidate.resolve()
    _check(resolved == resolved_root or resolved_root in resolved.parents,
           "SOURCE_REFERENCE_TRAVERSAL")
    return candidate


def _contains_hash(
    value: Any, expected: str, context: tuple[str, ...] = ()
) -> bool:
    if isinstance(value, dict):
        for key, child in value.items():
            child_context = context + (str(key).lower(),)
            binding_name = ".".join(child_context)
            if (
                child == expected
                and "sha256" in binding_name
                and any(
                    token in binding_name
                    for token in ("source", "original", "receipt")
                )
            ):
                return True
            if _contains_hash(child, expected, child_context):
                return True
    elif isinstance(value, list):
        return any(_contains_hash(child, expected, context) for child in value)
    return False


def _verified_base(path_value: str | Path) -> dict:
    path = Path(path_value)
    _check(".." not in path.parts, "BASE_PATH_TRAVERSAL")
    _check(path.is_file() and not _path_has_symlink(path), "BASE_REGULAR_FILE_REQUIRED")
    _check(PurePath(*path.parts[-len(KNOWN_BASE_SUFFIX.parts):]) == KNOWN_BASE_SUFFIX,
           "BASE_SOURCE_UNREGISTERED")
    manifest_hash = _sha256(path)
    _check(manifest_hash == KNOWN_BASE_SHA256, "FROZEN_BASE_MANIFEST_HASH_MISMATCH")
    value = json.loads(path.read_text(encoding="utf-8"))
    _check(isinstance(value, dict), "FROZEN_BASE_OBJECT_REQUIRED")
    source_root = path.parent
    rows = value.get("seed_receipts")
    if not isinstance(rows, list):
        rows = value.get("receipts", value.get("known_receipts"))
    _check(isinstance(rows, list) and all(isinstance(row, dict) for row in rows),
           "RECEIPTS_REQUIRED")

    verified_rows = []
    for source_row in rows:
        row = dict(source_row)
        seed = row.get("seed")
        _check(type(seed) is int and seed >= 0, "SEED_INVALID")
        original_ref = row.get("original_receipt_ref", row.get("original_ref"))
        expected_original_ref = PurePath(
            f"supporting/seeds/seed-{seed}/receipt.json"
        )
        original_path = _source_file(source_root, original_ref, expected_original_ref)
        recorded_hash = row.get("original_receipt_sha256")
        _check(isinstance(recorded_hash, str) and len(recorded_hash) == 64,
               "ORIGINAL_RECEIPT_HASH_REQUIRED")
        _check(_sha256(original_path) == recorded_hash,
               "ORIGINAL_RECEIPT_HASH_MISMATCH")
        original = json.loads(original_path.read_text(encoding="utf-8"))
        _check(isinstance(original, dict), "ORIGINAL_RECEIPT_OBJECT_REQUIRED")

        reinspection = None
        reinspection_ref = row.get("reinspection_ref")
        if reinspection_ref is not None:
            expected_reinspection_ref = PurePath(
                f"supporting/seeds/seed-{seed}/reinspection.json"
            )
            reinspection_path = _source_file(
                source_root, reinspection_ref, expected_reinspection_ref
            )
            reinspection = json.loads(reinspection_path.read_text(encoding="utf-8"))
            _check(isinstance(reinspection, dict), "REINSPECTION_OBJECT_REQUIRED")
            recorded_reinspection_hash = row.get("reinspection_sha256")
            if recorded_reinspection_hash is not None:
                _check(recorded_reinspection_hash == _sha256(reinspection_path),
                       "REINSPECTION_HASH_MISMATCH")
            _check(_contains_hash(reinspection, recorded_hash),
                   "REINSPECTION_SOURCE_BINDING_INVALID")

        original_completed = (
            original.get("exit_code") == 0
            and original.get("process_state") == "completed"
        )
        original_inspected = (
            isinstance(_assessment(original), dict)
            and _assessment(original).get("status") == "success"
        )
        derived_inspected = (
            seed == 23
            and isinstance(reinspection, dict)
            and isinstance(_assessment(reinspection), dict)
            and _assessment(reinspection).get("status") == "success"
        )
        proof_valid = original_completed and (original_inspected or derived_inspected)
        _check(row.get("execution_success") is not True or proof_valid,
               "FROZEN_BASE_EXECUTION_PROVENANCE_INVALID")

        row["original_receipt_ref"] = original_ref
        row["original_receipt"] = original
        row["reinspection_receipt"] = reinspection
        row["original_status"] = original.get("status")
        row["original_failure_reason"] = original.get("failure_reason")
        row["_source_bound_original_verified"] = proof_valid
        row["_derived_reinspection_verified"] = derived_inspected
        row["_frozen_manifest_sha256"] = manifest_hash
        verified_rows.append(row)

    copied = dict(value)
    copied["seed_receipts"] = verified_rows
    return copied


def _receipts(value: Any, *, base: bool = False) -> list[dict]:
    if base and isinstance(value, (str, Path)):
        value = _verified_base(value)
    else:
        value = _load(value)
    container = value if isinstance(value, dict) else None
    if isinstance(value, list):
        result = value
    elif isinstance(value, dict) and isinstance(value.get("seed_receipts"), list):
        result = value["seed_receipts"]
    elif isinstance(value, dict) and isinstance(value.get("receipts"), list):
        result = value["receipts"]
    elif isinstance(value, dict) and isinstance(value.get("known_receipts"), list):
        result = value["known_receipts"]
    elif isinstance(value, dict) and "seed" in value:
        result = [value]
    else:
        raise BenchmarkDistributionError("RECEIPTS_REQUIRED")
    _check(all(isinstance(item, dict) for item in result), "RECEIPT_OBJECT_REQUIRED")
    copied = [dict(item) for item in result]
    if base and container is not None:
        baseline = container.get("new_MSA_baseline", {}).get(
            "processed_msa_hashes_by_seed", {}
        )
        for item in copied:
            seed = item.get("seed")
            if str(seed) in baseline:
                item["_base_processed_msa"] = baseline[str(seed)]
            elif seed in baseline:
                item["_base_processed_msa"] = baseline[seed]
    return copied


def _assessment(receipt: dict) -> dict | None:
    for key in ("assessment", "reinspection", "inspection"):
        value = receipt.get(key)
        if isinstance(value, dict):
            nested = value.get("assessment")
            if isinstance(nested, dict):
                return nested
            if value.get("status") == "success":
                return value
    return None


def _execution(receipt: dict) -> tuple[bool, str]:
    actual = receipt.get("actual_execution_success")
    _check(actual is None or type(actual) is bool, "ACTUAL_EXECUTION_SUCCESS_NOT_BOOLEAN")
    if actual is True:
        return True, "completed_exit0"
    if receipt.get("exit_code") == 0 and receipt.get("process_state") == "completed":
        return True, "completed_exit0"
    if receipt.get("execution_success") is True:
        _check(receipt.get("_source_bound_original_verified") is True,
               "FROZEN_BASE_EXECUTION_PROVENANCE_INVALID")
        if receipt.get("_derived_reinspection_verified") is True:
            return True, "frozen_base_hash_verified_reinspection"
        return True, "frozen_base_hash_verified_original_inspection"
    return False, "process_failure"


def _binding(receipt: dict, key: str) -> Any:
    hashes = receipt.get("hashes", {})
    if key in hashes:
        return hashes[key]
    verified = receipt.get("verified_original_recorded_hashes", {})
    if key in verified:
        return verified[key]
    original = receipt.get("original_receipt", {})
    return original.get("hashes", {}).get(key)


def _environment(receipt: dict) -> Any:
    value = receipt.get("tool_environment")
    if not isinstance(value, dict):
        value = receipt.get("original_receipt", {}).get("tool_environment")
    if not isinstance(value, dict):
        return None
    canonical = {
        "boltz_version": value.get("boltz_version"),
        "torch_version": value.get("torch_version"),
        "cuda_available": value.get("cuda_available"),
        "cuda_runtime": value.get("cuda_runtime"),
        "cuda_devices": value.get("cuda_devices"),
    }
    rdkit = value.get("rdkit_version", value.get("RDKit_version"))
    if rdkit is not None:
        canonical["rdkit_version"] = rdkit
    _check(all(item is not None for item in canonical.values()), "ENVIRONMENT_BINDING_MISSING")
    return canonical


def _settings(receipt: dict) -> Any:
    return (receipt.get("settings") or receipt.get("confirmed_execution_metadata", {}).get("applied_settings")
            or receipt.get("original_receipt", {}).get("settings"))


def _msa(receipt: dict) -> Any:
    baseline = receipt.get("_base_processed_msa")
    if isinstance(baseline, dict):
        normalized = {PurePath(str(key).replace("\\", "/")).name: value
                      for key, value in baseline.items()}
        _check(len(normalized) == len(baseline), "MSA_FILENAME_COLLISION")
        return normalized
    hashes = receipt.get("hashes", {}).get("processed_msa_files_sha256")
    if hashes is None:
        hashes = receipt.get("hashes", {}).get("processed_output_files_sha256")
    if hashes is None:
        hashes = receipt.get("verified_processed_output_files_sha256")
    if isinstance(hashes, dict):
        entries = [
            (PurePath(str(key).replace("\\", "/")).name, value)
            for key, value in hashes.items()
            if "/processed/msa/" in "/" + str(key).replace("\\", "/")
        ]
        filtered = dict(entries)
        _check(len(filtered) == len(entries), "MSA_FILENAME_COLLISION")
        if filtered:
            return filtered
    return None


def _numeric_leaves(value: Any, prefix: str = "") -> dict[str, float]:
    result: dict[str, float] = {}
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return result
    if isinstance(value, (int, float)):
        _check(math.isfinite(float(value)), "METRIC_NONFINITE")
        result[prefix] = float(value)
    elif isinstance(value, dict):
        for key, child in value.items():
            child_prefix = f"{prefix}.{key}" if prefix else str(key)
            result.update(_numeric_leaves(child, child_prefix))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            child_prefix = f"{prefix}[{index}]"
            result.update(_numeric_leaves(child, child_prefix))
    return result


def _linear(values: list[float], q: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * q
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] + (ordered[high] - ordered[low]) * (position - low)


def _summary(values: list[float]) -> dict:
    q1, median, q3 = _linear(values, .25), _linear(values, .5), _linear(values, .75)
    return {"count": len(values), "min": min(values) if values else None,
            "max": max(values) if values else None, "median": median,
            "q1": q1, "q3": q3, "iqr": (q3 - q1) if values else None,
            "quantile_method": "linear interpolation at (n-1)q"}


def _number(value: Any, code: str) -> float:
    _check(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value), code)
    return float(value)


def build_benchmark_distribution(base_receipts: Any, new_receipts: Any) -> dict:
    receipts = _receipts(base_receipts, base=True) + _receipts(new_receipts)
    seeds = [item.get("seed", item.get("confirmed_execution_metadata", {}).get("applied_seed"))
             for item in receipts]
    _check(all(type(seed) is int and seed >= 0 for seed in seeds), "SEED_INVALID")
    _check(len(seeds) == len(set(seeds)), "DUPLICATE_SEED")

    executions = [_execution(item) for item in receipts]
    fields = {
        "input": [_binding(item, "input_yaml_sha256") for item in receipts],
        "checkpoint": [_binding(item, "checkpoint_sha256") for item in receipts],
        "settings": [_settings(item) for item in receipts],
        "environment": [_environment(item) for item in receipts],
        "msa": [_msa(item) for item in receipts],
    }
    unknown_fields: dict[str, list[int]] = {}
    for name, values in fields.items():
        known_computed = [value for value, (executed, _status) in zip(values, executions)
                          if executed and value is not None]
        if known_computed:
            first = known_computed[0]
            _check(all(value == first for value in known_computed), "MIXED_" + name.upper())
        missing = [seed for seed, value in zip(seeds, values) if value is None]
        if missing:
            unknown_fields[name] = missing

    rows: list[dict] = []
    metric_values: dict[str, list[float]] = {}
    technical_successes = 0
    for seed, receipt, (executed, execution_status) in zip(seeds, receipts, executions):
        technical_successes += int(executed)
        assessment = _assessment(receipt)
        inspected = isinstance(assessment, dict) and assessment.get("status") == "success"
        row = {"seed": seed, "technical_execution_success": executed,
               "execution_status": execution_status,
               "inspection_status": "success" if inspected else "failed_or_missing",
               "original_status": receipt.get("original_status", receipt.get("status")),
               "original_failure_reason": receipt.get("original_failure_reason", receipt.get("failure_reason")),
               "provenance": receipt.get("provenance"), "metrics": None}
        if inspected:
            supplied = {}
            supplied.update({"comparison.metrics." + key: value for key, value in
                             _numeric_leaves(assessment.get("comparison", {}).get("metrics", {})).items()})
            supplied.update({"model_confidence." + key: value for key, value in
                             _numeric_leaves(assessment.get("model_confidence", {})).items()})
            row["metrics"] = supplied
            for key, value in supplied.items():
                metric_values.setdefault(key, []).append(value)
        rows.append(row)

    comparability_known = not unknown_fields
    return {
        "format": "tpd-known-case-benchmark-distribution/1",
        "status": "proposal_protocol_ready_not_calibration_approved",
        "seeds": rows, "seed_count": len(rows),
        "technical_success_count": technical_successes,
        "technical_failure_count": len(rows) - technical_successes,
        "scientific_model_quality_success_count": None,
        "expert_geometric_success_cutoff": None,
        "distribution_strength": "weak_single_known_case",
        "poor_replicates_retained": True, "outlier_deletion_performed": False,
        "statistics": {key: _summary(values) for key, values in sorted(metric_values.items())},
        "comparability": {"uniform": True if comparability_known else None,
                          "comparability_unknown": not comparability_known,
                          "unknown_binding_seeds": unknown_fields,
                          "bindings": fields},
        "interpretation": "Known-case technical distribution only; no efficacy, E3 superiority, calibration approval, confidence cutoff, or invented geometric-success cutoff.",
    }


merge_benchmark_distribution = build_benchmark_distribution
