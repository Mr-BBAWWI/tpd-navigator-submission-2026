"""Fail-closed verification of portable campaign diagnostic evidence packs.

Only manifest-declared, validated pack-relative files are opened. Historical
absolute provenance strings are treated as inert data. Bundled Python source is
hashed as evidence but is never imported or executed.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
import os
import re
import shutil
import stat
import tempfile
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Any

import yaml

from . import boltz_worker
from . import novel_ternary

MANIFEST_NAME = "output-manifest.json"
COMPACT_NAME = "compact-campaign-diagnostics.json"
PROVENANCE_NAME = "source-provenance.json"
MANIFEST_FORMAT = "campaign-diagnostics-output-manifest/1"
COMPACT_FORMAT = "compact-campaign-diagnostics/20261001.2"
SEEDS = [23, 41, 61, 79, 97]
MODELS = [0, 1, 2, 3, 4]
MAX_FILES = 1024
MAX_BYTES = 256 * 1024 * 1024
ALLOWED_SUFFIXES = {".json", ".cif", ".sdf", ".yaml", ".html", ".py", ".md"}
SHA_RE = re.compile(r"[0-9a-f]{64}\Z")
CALIBRATION_CUTOFFS = {
    "target_CA_RMSD_A": 1.5,
    "ligand_heavy_atom_RMSD_after_target_alignment_A": 3.0,
    "e3_CA_RMSD_after_target_alignment_A": 10.0,
    "contact_jaccard": 0.4,
}
CALIBRATION_SETTINGS = {
    "recycling_steps": 10,
    "sampling_steps": 200,
    "diffusion_samples": 5,
    "max_parallel_samples": 1,
    "max_msa_seqs": 256,
    "use_potentials": False,
}
NOVEL_SETTINGS = {
    "recycling_steps": 10,
    "sampling_steps": 200,
    "diffusion_samples": 5,
    "max_parallel_samples": 1,
    "max_msa_seqs": 256,
    "potentials": False,
}


class ProtocolEvidenceError(ValueError):
    """The evidence pack is malformed, incomplete, unbound, or inconsistent."""


def _fail(code: str) -> None:
    raise ProtocolEvidenceError(code)


def _check(condition: bool, code: str) -> None:
    if not condition:
        _fail(code)


def _sha_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _pairs_no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            _fail("JSON_DUPLICATE_KEY:" + key)
        result[key] = value
    return result


def _reject_constant(value: str) -> None:
    _fail("JSON_NONFINITE_NUMBER:" + value)


def _json_bytes(data: bytes, label: str) -> Any:
    try:
        text = data.decode("utf-8")
        value = json.loads(
            text,
            object_pairs_hook=_pairs_no_duplicates,
            parse_constant=_reject_constant,
        )
        def finite(item: Any) -> bool:
            if isinstance(item, float):
                return math.isfinite(item)
            if isinstance(item, dict):
                return all(finite(child) for child in item.values())
            if isinstance(item, list):
                return all(finite(child) for child in item)
            return True
        _check(finite(value), "JSON_NONFINITE_NUMBER:" + label)
        return value
    except ProtocolEvidenceError:
        raise
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ProtocolEvidenceError("JSON_INVALID:" + label) from error


def _json_file(path: Path, label: str) -> Any:
    try:
        return _json_bytes(path.read_bytes(), label)
    except OSError as error:
        raise ProtocolEvidenceError("FILE_READ_FAILED:" + label) from error


def _canonical_digest(value: Any) -> str:
    try:
        encoded = json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
            allow_nan=False,
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise ProtocolEvidenceError("CANONICAL_JSON_INVALID") from error
    return _sha_bytes(encoded)


def _safe_relative(value: Any) -> str:
    _check(isinstance(value, str) and value != "", "PATH_INVALID")
    _check("\\" not in value and "\x00" not in value and ":" not in value,
           "PATH_BACKSLASH_NUL_OR_COLON")
    posix = PurePosixPath(value)
    windows = PureWindowsPath(value)
    _check(not posix.is_absolute() and not windows.is_absolute() and not windows.drive,
           "PATH_ABSOLUTE_OR_DRIVE")
    _check(posix.as_posix() == value, "PATH_NOT_NORMALIZED_POSIX")
    _check(all(part not in {"", ".", ".."} for part in posix.parts), "PATH_TRAVERSAL")
    _check(PurePosixPath(value).suffix.lower() in ALLOWED_SUFFIXES, "PATH_SUFFIX_FORBIDDEN")
    return value


def _regular_no_links(path: Path, root: Path) -> None:
    try:
        root_real = root.resolve(strict=True)
        resolved = path.resolve(strict=True)
        resolved.relative_to(root_real)
    except (OSError, ValueError) as error:
        raise ProtocolEvidenceError("PATH_ESCAPES_PACK") from error
    current = path
    while True:
        try:
            info = current.lstat()
        except OSError as error:
            raise ProtocolEvidenceError("PATH_STAT_FAILED") from error
        _check(not stat.S_ISLNK(info.st_mode), "SYMLINK_OR_REPARSE_FORBIDDEN")
        if hasattr(info, "st_file_attributes"):
            reparse = getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400)
            tag = getattr(info, "st_reparse_tag", 0)
            name_surrogate = 0x20000000
            _check(not (info.st_file_attributes & reparse and
                        (tag in {0xA000000C, 0xA0000003} or tag & name_surrogate)),
                   "SYMLINK_OR_REPARSE_FORBIDDEN")
        if current == root:
            break
        current = current.parent
    _check(path.is_file(), "REGULAR_FILE_REQUIRED")


def _walk_files(root: Path) -> dict[str, Path]:
    _check(root.is_dir() and not root.is_symlink(), "PACK_ROOT_DIRECTORY_REQUIRED")
    result: dict[str, Path] = {}
    folded: set[str] = set()
    total = 0
    for current, directories, files in os.walk(root, followlinks=False):
        base = Path(current)
        for name in directories:
            child = base / name
            _check(not child.is_symlink(), "SYMLINK_OR_REPARSE_FORBIDDEN")
            try:
                mode = child.lstat().st_mode
            except OSError as error:
                raise ProtocolEvidenceError("PATH_STAT_FAILED") from error
            _check(stat.S_ISDIR(mode), "DIRECTORY_REQUIRED")
            info = child.lstat()
            if hasattr(info, "st_file_attributes"):
                reparse = getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400)
                tag = getattr(info, "st_reparse_tag", 0)
                _check(not (info.st_file_attributes & reparse and
                            (tag in {0xA000000C, 0xA0000003} or tag & 0x20000000)),
                       "SYMLINK_OR_REPARSE_FORBIDDEN")
        for name in files:
            path = base / name
            relative = _safe_relative(path.relative_to(root).as_posix())
            folded_name = relative.casefold()
            _check(folded_name not in folded, "CASEFOLD_PATH_DUPLICATE")
            folded.add(folded_name)
            _regular_no_links(path, root)
            size = path.stat().st_size
            total += size
            _check(len(result) + 1 <= MAX_FILES, "PACK_FILE_COUNT_LIMIT")
            _check(total <= MAX_BYTES, "PACK_BYTE_LIMIT")
            result[relative] = path
    return result


def _manifest_entries(manifest: dict) -> dict[str, tuple[str, int]]:
    _check(manifest.get("format") == MANIFEST_FORMAT, "MANIFEST_FORMAT_INVALID")
    _check(set(manifest) == {"format", "files", "manifest_excludes_itself"} and
           manifest.get("manifest_excludes_itself") is True,
           "MANIFEST_FIELDS_INVALID")
    raw = manifest.get("files")
    _check(isinstance(raw, list), "MANIFEST_FILES_REQUIRED")
    entries: dict[str, tuple[str, int]] = {}
    folded: set[str] = set()
    for item in raw:
        _check(isinstance(item, dict) and set(item) >= {"path", "sha256", "bytes"},
               "MANIFEST_FILE_RECORD_INVALID")
        name = _safe_relative(item["path"])
        _check(name != MANIFEST_NAME, "MANIFEST_SELF_ENTRY_FORBIDDEN")
        _check(name not in entries and name.casefold() not in folded, "MANIFEST_PATH_DUPLICATE")
        digest, size = item["sha256"], item["bytes"]
        _check(isinstance(digest, str) and SHA_RE.fullmatch(digest) is not None,
               "MANIFEST_SHA256_INVALID")
        _check(type(size) is int and 0 <= size <= MAX_BYTES, "MANIFEST_SIZE_INVALID")
        folded.add(name.casefold())
        entries[name] = (digest, size)
    _check(len(entries) + 1 <= MAX_FILES, "PACK_FILE_COUNT_LIMIT")
    return entries


def _verify_manifest(root: Path, expected_sha: str) -> tuple[dict, dict[str, Path], list[dict]]:
    _check(isinstance(expected_sha, str) and SHA_RE.fullmatch(expected_sha) is not None,
           "EXPECTED_MANIFEST_SHA256_INVALID")
    files = _walk_files(root)
    _check(MANIFEST_NAME in files, "OUTPUT_MANIFEST_MISSING")
    manifest_bytes = files[MANIFEST_NAME].read_bytes()
    _check(_sha_bytes(manifest_bytes) == expected_sha, "MANIFEST_HASH_PIN_MISMATCH")
    manifest = _json_bytes(manifest_bytes, MANIFEST_NAME)
    _check(isinstance(manifest, dict), "MANIFEST_OBJECT_REQUIRED")
    entries = _manifest_entries(manifest)
    _check(set(files) == set(entries) | {MANIFEST_NAME}, "PACK_FILE_SET_MISMATCH")
    normalized = []
    for name in sorted(entries):
        digest, size = entries[name]
        path = files[name]
        _check(path.stat().st_size == size, "FILE_SIZE_MISMATCH:" + name)
        _check(_sha_file(path) == digest, "FILE_HASH_MISMATCH:" + name)
        normalized.append({"path": name, "sha256": digest, "bytes": size})
    normalized.append({"path": MANIFEST_NAME, "sha256": expected_sha,
                       "bytes": len(manifest_bytes)})
    normalized.sort(key=lambda item: item["path"])
    return manifest, files, normalized


def _approval_false(value: Any, location: str = "$.") -> None:
    if isinstance(value, dict):
        for key, child in value.items():
            low = key.lower()
            if ("approval" in low or low in {"approved", "scientific_approved",
                                             "scientifically_approved"}) and child is True:
                _fail("APPROVAL_TRUE_FORBIDDEN:" + location + key)
            _approval_false(child, location + key + ".")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _approval_false(child, f"{location}{index}.")


def _provenance_records(value: Any) -> list[dict]:
    if isinstance(value, list):
        records = value
    elif isinstance(value, dict):
        records = value.get("records") if isinstance(value.get("records"), list) else list(value.values())
    else:
        _fail("SOURCE_PROVENANCE_INVALID")
    _check(all(isinstance(item, dict) for item in records), "SOURCE_PROVENANCE_RECORD_INVALID")
    return records


def _verify_provenance(compact: dict, files: dict[str, Path]) -> None:
    pointer = compact.get("source_provenance")
    _check(isinstance(pointer, dict), "COMPACT_SOURCE_PROVENANCE_REQUIRED")
    name = _safe_relative(pointer.get("portable_pack_relative_path"))
    _check(name == PROVENANCE_NAME and name in files, "SOURCE_PROVENANCE_PATH_INVALID")
    digest = pointer.get("sha256")
    _check(isinstance(digest, str) and _sha_file(files[name]) == digest,
           "SOURCE_PROVENANCE_POINTER_HASH_MISMATCH")
    provenance = _json_file(files[name], name)
    records = _provenance_records(provenance)
    by_pair: dict[tuple[str, str], dict] = {}
    for record in records:
        path = _safe_relative(record.get("portable_pack_relative_path"))
        digest = record.get("source_sha256")
        _check(path in files and isinstance(digest, str) and SHA_RE.fullmatch(digest) is not None,
               "SOURCE_PROVENANCE_REFERENCE_INVALID")
        _check(_sha_file(files[path]) == digest, "SOURCE_PROVENANCE_HASH_MISMATCH:" + path)
        pair = (path, digest)
        _check(pair not in by_pair, "SOURCE_PROVENANCE_DUPLICATE")
        by_pair[pair] = record
    sources = compact.get("sources")
    _check(isinstance(sources, dict), "COMPACT_SOURCES_REQUIRED")
    for source in sources.values():
        _check(isinstance(source, dict), "COMPACT_SOURCE_RECORD_INVALID")
        pair = (_safe_relative(source.get("portable_pack_relative_path")),
                source.get("source_sha256"))
        _check(pair in by_pair, "COMPACT_SOURCE_PROVENANCE_MISMATCH")


def _load_required(files: dict[str, Path], name: str) -> dict:
    _check(name in files, "REQUIRED_FILE_MISSING:" + name)
    value = _json_file(files[name], name)
    _check(isinstance(value, dict), "JSON_OBJECT_REQUIRED:" + name)
    return value


def _same_number(first: Any, second: Any, tolerance: float = 1e-6) -> bool:
    if isinstance(first, bool) or isinstance(second, bool):
        return first is second
    if isinstance(first, (int, float)) and isinstance(second, (int, float)):
        return math.isfinite(float(first)) and math.isfinite(float(second)) and math.isclose(
            float(first), float(second), rel_tol=tolerance, abs_tol=tolerance)
    return first == second


def _strict_scalar(actual: Any, expected: Any, code: str) -> None:
    _check(type(actual) is type(expected) and actual == expected, code)


def _exact_int_list(actual: Any, expected: list[int], code: str) -> None:
    _check(isinstance(actual, list) and len(actual) == len(expected) and
           all(type(item) is int for item in actual) and actual == expected, code)


def _verify_effective_command(receipt: dict, seed: int, code: str) -> None:
    command = receipt.get("effective_command")
    _check(isinstance(command, list) and all(isinstance(item, str) for item in command), code)
    options: dict[str, str | None] = {}
    index = 0
    while index < len(command):
        token = command[index]
        if not token.startswith("--"):
            index += 1
            continue
        _check("=" not in token and token not in options, code)
        value = None
        if index + 1 < len(command) and not command[index + 1].startswith("--"):
            value = command[index + 1]
            index += 1
        options[token] = value
        index += 1
    _check("--use_potentials" not in options and "--use_msa_server" not in options, code)
    expected = {
        "--seed": str(seed),
        "--recycling_steps": "10",
        "--sampling_steps": "200",
        "--diffusion_samples": "5",
        "--max_parallel_samples": "1",
        "--max_msa_seqs": "256",
        "--model": "boltz2",
    }
    _check(all(options.get(key) == value for key, value in expected.items()), code)


def _compare_tree(actual: Any, expected: Any, ignored_keys: set[str] | None = None,
                  location: str = "$", tolerance: float = 1e-6) -> None:
    ignored = ignored_keys or set()
    if isinstance(expected, dict):
        _check(isinstance(actual, dict), "VALUE_TYPE_MISMATCH:" + location)
        for key, value in expected.items():
            if key in ignored:
                continue
            _check(key in actual, "VALUE_MISSING:" + location + "." + key)
            _compare_tree(actual[key], value, ignored, location + "." + key, tolerance)
    elif isinstance(expected, list):
        _check(isinstance(actual, list) and len(actual) == len(expected),
               "LIST_MISMATCH:" + location)
        for index, value in enumerate(expected):
            _compare_tree(actual[index], value, ignored, f"{location}[{index}]", tolerance)
    else:
        _check(_same_number(actual, expected, tolerance), "VALUE_MISMATCH:" + location)


def _finite_score(model: dict) -> float:
    value = model.get("confidence", {}).get("confidence_score")
    _check(not isinstance(value, bool) and isinstance(value, (int, float)) and
           math.isfinite(float(value)), "FINITE_CONFIDENCE_SCORE_REQUIRED")
    return float(value)


def _selected_index(models: list[dict]) -> int:
    _check(all(isinstance(item, dict) and type(item.get("model_index")) is int
               for item in models), "MODEL_INDEX_INVALID")
    _check(len(models) == 5 and {item.get("model_index") for item in models} == set(MODELS),
           "ALL_FIVE_MODELS_REQUIRED")
    return min(((-_finite_score(item), item["model_index"]) for item in models))[1]


def _metrics_pass(metrics: dict) -> bool:
    values: dict[str, float] = {}
    for key in CALIBRATION_CUTOFFS:
        value = metrics.get(key)
        _check(not isinstance(value, bool) and isinstance(value, (int, float)) and
               math.isfinite(float(value)), "CALIBRATION_METRIC_INVALID")
        values[key] = float(value)
    _check(values["target_CA_RMSD_A"] >= 0 and
           values["ligand_heavy_atom_RMSD_after_target_alignment_A"] >= 0 and
           values["e3_CA_RMSD_after_target_alignment_A"] >= 0 and
           0 <= values["contact_jaccard"] <= 1,
           "CALIBRATION_METRIC_RANGE_INVALID")
    return (
        values["target_CA_RMSD_A"] <= 1.5 and
        values["ligand_heavy_atom_RMSD_after_target_alignment_A"] <= 3.0 and
        values["e3_CA_RMSD_after_target_alignment_A"] <= 10.0 and
        values["contact_jaccard"] >= 0.4
    )


def _calibration_raw_paths(seed: int, index: int) -> tuple[str, str]:
    base = f"evidence/calibration/raw/seed-{seed}"
    return f"{base}/model-{index}.cif", f"{base}/confidence-{index}.json"


def _verify_calibration_receipt_msa_settings(settings: dict) -> None:
    _strict_scalar(settings.get("msa_mode"), "cached_processed_offline",
                   "CALIBRATION_RECEIPT_SETTINGS_INVALID:msa_mode")
    if "use_msa_server" in settings:
        _strict_scalar(settings["use_msa_server"], False,
                       "CALIBRATION_RECEIPT_SETTINGS_INVALID:use_msa_server")


def _verify_calibration(files: dict[str, Path], compact: dict, trusted_root: Path) -> dict:
    plan_name = "evidence/calibration/plan.json"
    summary_name = "evidence/calibration/protocol-summary.json"
    plan = _load_required(files, plan_name)
    summary = _load_required(files, summary_name)
    _check(plan.get("format") == "tpd-crbn-preregistered-development-protocol/1.0",
           "CALIBRATION_PLAN_FORMAT_INVALID")
    _check(summary.get("format") == "tpd-crbn-preregistered-development-summary/1.0",
           "CALIBRATION_SUMMARY_FORMAT_INVALID")
    _exact_int_list(plan.get("seeds"), SEEDS, "CALIBRATION_SEED_MODEL_SET_INVALID")
    _exact_int_list(plan.get("expected_models_per_seed"), MODELS,
                    "CALIBRATION_SEED_MODEL_SET_INVALID")
    _check(plan.get("cutoffs") == CALIBRATION_CUTOFFS and
           summary.get("cutoffs") == CALIBRATION_CUTOFFS, "CALIBRATION_CUTOFFS_CHANGED")
    settings = plan.get("effective_protocol_settings", {})
    _check(isinstance(settings, dict), "CALIBRATION_SETTINGS_INVALID")
    for key, expected in CALIBRATION_SETTINGS.items():
        _strict_scalar(settings.get(key), expected, "CALIBRATION_SETTINGS_INVALID:" + key)
    _check(plan.get("execution") == "serial_all_five_no_early_stop",
           "CALIBRATION_NOT_SERIAL_COMPLETE")
    _check(plan.get("mutates_store_or_acceptance") is False,
           "CALIBRATION_ACCEPTANCE_MUTATION_FORBIDDEN")

    compact_cal = compact.get("known_calibration")
    _check(isinstance(compact_cal, dict), "COMPACT_CALIBRATION_REQUIRED")
    _check(compact_cal.get("cutoffs") == CALIBRATION_CUTOFFS,
           "COMPACT_CALIBRATION_CUTOFF_MISMATCH")
    receipts = summary.get("seeds")
    _check(isinstance(receipts, list) and len(receipts) == 5 and
           all(isinstance(item, dict) and type(item.get("seed")) is int for item in receipts),
           "CALIBRATION_RECEIPTS_INCOMPLETE")
    by_seed = {item["seed"]: item for item in receipts}
    _check(set(by_seed) == set(SEEDS), "CALIBRATION_RECEIPTS_INCOMPLETE")

    reference = trusted_root / "cases" / "design_sources" / "6BOY.cif"
    metadata = trusted_root / "cases" / "design_sources" / "crbn_benchmark.json"
    reference_hash = plan.get("sources", {}).get("reference", {}).get("sha256")
    _check(reference.is_file() and not reference.is_symlink() and
           metadata.is_file() and not metadata.is_symlink(),
           "TRUSTED_CALIBRATION_SOURCE_MISSING")
    _check(_sha_file(reference) == reference_hash, "TRUSTED_REFERENCE_HASH_MISMATCH")
    try:
        with tempfile.TemporaryDirectory(prefix="protocol-calibration-") as temporary:
            trusted_yaml = Path(temporary) / "6BOY_CRBN.yaml"
            trusted_packet = copy.deepcopy(
                boltz_worker.prepare_input(reference, metadata, trusted_yaml, "server")
            )
            trusted_input_bytes = trusted_yaml.read_bytes()
    except Exception as error:
        raise ProtocolEvidenceError("TRUSTED_CALIBRATION_INPUT_FAILED") from error
    geometry_recomputed = True

    raw_pass = 0
    selected_pass = 0
    selected_rows = []
    for seed in SEEDS:
        receipt = by_seed[seed]
        portable_receipt_name = f"evidence/calibration/raw/seed-{seed}/receipt.json"
        portable_receipt = _load_required(files, portable_receipt_name)
        _compare_tree(portable_receipt, receipt)
        _check(receipt.get("status") == "success" and
               type(receipt.get("exit_code")) is int and receipt["exit_code"] == 0 and
               receipt.get("process_state") == "completed" and
               receipt.get("scientifically_approved") is False and
               type(receipt.get("seed")) is int and receipt["seed"] == seed,
               "CALIBRATION_SEED_NOT_COMPLETE")
        receipt_settings = receipt.get("effective_protocol_settings")
        _check(isinstance(receipt_settings, dict), "CALIBRATION_RECEIPT_SETTINGS_REQUIRED")
        for key, expected in CALIBRATION_SETTINGS.items():
            _strict_scalar(receipt_settings.get(key), expected,
                           "CALIBRATION_RECEIPT_SETTINGS_INVALID:" + key)
        _verify_calibration_receipt_msa_settings(receipt_settings)
        _verify_effective_command(receipt, seed, "CALIBRATION_EFFECTIVE_COMMAND_INVALID")
        models = receipt.get("models")
        _check(isinstance(models, list), "CALIBRATION_MODELS_REQUIRED")
        selected = _selected_index(models)
        selection = receipt.get("selection", {})
        _check(selection.get("selected_model_index") == selected,
               "CALIBRATION_SELECTION_ORACLE_MISMATCH")
        candidates = selection.get("candidates")
        _check(isinstance(candidates, list) and
               {(x.get("model_index"), float(x.get("confidence_score"))) for x in candidates}
               == {(m["model_index"], _finite_score(m)) for m in models},
               "CALIBRATION_SELECTION_CANDIDATES_MISMATCH")
        packet = receipt.get("input")
        _check(isinstance(packet, dict), "CALIBRATION_INPUT_PACKET_REQUIRED")
        for key in ("protein_sequences", "reference_chain_mapping", "ligand",
                    "input_chain_roles", "explicit_exclusions", "omissions", "msa_mode"):
            _check(packet.get(key) == trusted_packet.get(key),
                   "CALIBRATION_INPUT_PACKET_UNTRUSTED:" + key)
        input_name = f"evidence/calibration/raw/seed-{seed}/input.yaml"
        _check(input_name in files and
               _sha_file(files[input_name]) == packet.get("input_sha256") ==
               trusted_packet.get("input_sha256") == _sha_bytes(trusted_input_bytes),
               "CALIBRATION_INPUT_YAML_HASH_MISMATCH")
        for model in models:
            index = model["model_index"]
            cif_name, confidence_name = _calibration_raw_paths(seed, index)
            _check(cif_name in files and confidence_name in files,
                   "CALIBRATION_RAW_MODEL_MISSING")
            _check(_sha_file(files[cif_name]) == model.get("prediction_sha256"),
                   "CALIBRATION_RAW_CIF_HASH_MISMATCH")
            _check(_sha_file(files[confidence_name]) == model.get("confidence_sha256"),
                   "CALIBRATION_RAW_CONFIDENCE_HASH_MISMATCH")
            raw_confidence = _json_file(files[confidence_name], confidence_name)
            _check(raw_confidence == model.get("confidence") == model.get("confidence_raw"),
                   "CALIBRATION_CONFIDENCE_CONTENT_MISMATCH")
            metrics = model.get("metrics")
            _check(isinstance(metrics, dict), "CALIBRATION_METRICS_REQUIRED")
            claimed = _metrics_pass(metrics)
            _check(model.get("passes_cutoffs") is claimed,
                   "CALIBRATION_MODEL_PASS_MISMATCH")
            raw_pass += int(claimed)
            if geometry_recomputed:
                comparison = boltz_worker.compare_structures(reference, files[cif_name], packet)
                _compare_tree(comparison["metrics"], metrics, tolerance=1e-6)
        chosen = next(item for item in models if item["model_index"] == selected)
        claimed_selected = _metrics_pass(chosen["metrics"])
        _check(receipt.get("selected_model_pass") is claimed_selected,
               "CALIBRATION_SELECTED_PASS_MISMATCH")
        selected_pass += int(claimed_selected)
        selected_rows.append({
            "seed": seed,
            "selected_model_index": selected,
            "selected_e3_CA_RMSD_after_target_alignment_A":
                chosen["metrics"]["e3_CA_RMSD_after_target_alignment_A"],
            "selected_pass": claimed_selected,
        })

    _check(raw_pass == summary.get("raw_all_four_pass_count", raw_pass),
           "CALIBRATION_SUMMARY_RAW_PASS_MISMATCH")
    _check(selected_pass == summary.get("selected_pass_count"),
           "CALIBRATION_SUMMARY_SELECTED_PASS_MISMATCH")
    _check(raw_pass == compact_cal.get("raw_all_four_pass_count") and
           selected_pass == compact_cal.get("selected_pass_count") and
           compact_cal.get("raw_model_count") == 25 and
           compact_cal.get("seed_count") == 5,
           "COMPACT_CALIBRATION_COUNT_MISMATCH")
    _compare_tree(compact_cal.get("per_seed_selected"), selected_rows, tolerance=1e-6)
    claimed_protocol_pass = selected_pass >= 3
    _check(summary.get("protocol_pass") is claimed_protocol_pass,
           "CALIBRATION_PROTOCOL_PASS_MISMATCH")
    protocol_id = "calibration:" + _sha_file(files[plan_name])
    return {
        "protocol_id": protocol_id,
        "kind": "known_structure_calibration",
        "summary": {
            "benchmark_scope": "BRD4-CRBN 6BOY related diagnostic; not the SMARCA2 novel target",
            "ligand_scope": "RN6 benchmark ligand; no claim that it matches a SMARCA2 candidate ligand",
            "raw_pass_count": raw_pass,
            "raw_model_denominator": 25,
            "selected_pass_count": selected_pass,
            "selected_seed_denominator": 5,
            "claimed_protocol_pass": claimed_protocol_pass,
            "failure_count": 0,
            "selected": selected_rows,
            "preserved_original_baseline": "2_of_5",
        },
        "measurement_verification": {
            "verification_level": "raw hashes + confidence top1 + all 25 geometry recomputed",
            "geometry_recomputed": True,
            "stored_metrics_only": False,
            "claimed_pass_reported_separately": True,
        },
        "scientific_approved": False,
        "gates_affected": False,
    }


def _validate_expected_inputs(expected_binding: dict, expected_graphs: dict) -> None:
    required = {"project", "job_id", "input_sha256", "result_sha256", "parent_id",
                "result_binding_kind"}
    _check(isinstance(expected_binding, dict) and set(expected_binding) == required,
           "EXPECTED_BINDING_FIELDS_INVALID")
    _check(all(isinstance(expected_binding[key], str) and expected_binding[key]
               for key in required), "EXPECTED_BINDING_VALUE_INVALID")
    _check(SHA_RE.fullmatch(expected_binding["input_sha256"]) is not None and
           SHA_RE.fullmatch(expected_binding["result_sha256"]) is not None and
           expected_binding["result_binding_kind"] == "exact_archived_design_json_bytes",
           "EXPECTED_BINDING_VALUE_INVALID")
    _check(isinstance(expected_graphs, dict) and len(expected_graphs) >= 2,
           "EXPECTED_CANDIDATE_GRAPHS_INVALID")
    _check(all(isinstance(key, str) and isinstance(value, dict)
               for key, value in expected_graphs.items()), "EXPECTED_CANDIDATE_GRAPHS_INVALID")


def _trusted_regular_source(path: Path, trusted_root: Path) -> Path:
    """Require an installed source to be a regular, non-link file under the trust root."""
    try:
        root_real = trusted_root.resolve(strict=True)
        resolved = path.resolve(strict=True)
        resolved.relative_to(root_real)
    except (OSError, ValueError) as error:
        raise ProtocolEvidenceError("TRUSTED_SOURCE_PATH_INVALID") from error

    current = path
    while True:
        try:
            info = current.lstat()
        except OSError as error:
            raise ProtocolEvidenceError("TRUSTED_SOURCE_STAT_FAILED") from error
        _check(not stat.S_ISLNK(info.st_mode), "TRUSTED_SOURCE_SYMLINK_OR_JUNCTION_FORBIDDEN")
        if hasattr(info, "st_file_attributes"):
            reparse = getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400)
            tag = getattr(info, "st_reparse_tag", 0)
            # Name-surrogate reparse points include symlinks and junctions. Cloud-file
            # reparse points used by ordinary OneDrive installations remain allowed.
            _check(not (info.st_file_attributes & reparse and tag & 0x20000000),
                   "TRUSTED_SOURCE_SYMLINK_OR_JUNCTION_FORBIDDEN")
        if current == root_real:
            break
        current = current.parent
    try:
        mode = path.lstat().st_mode
    except OSError as error:
        raise ProtocolEvidenceError("TRUSTED_SOURCE_STAT_FAILED") from error
    _check(stat.S_ISREG(mode), "TRUSTED_SOURCE_REGULAR_FILE_REQUIRED")
    return resolved


def _validated_trusted_plan_sources(plan: dict, graph: dict) -> dict[str, dict]:
    trusted_root = Path(__file__).resolve().parents[2]
    e3_type = graph.get("e3_type")
    fixed_paths = {
        "target": trusted_root / "cases" / "design_sources" / "6HAZ.cif",
        "CRBN": trusted_root / "cases" / "design_sources" / "6BOY.cif",
        "VHL": trusted_root / "cases" / "acceptance_sources" / "supporting" /
               "catalog" / "source_snapshots" /
               "f9a7ab0aca167f758cad8d77b2173cdfbd7ab5edd2cb91bc6dacd6e9901b5d30" /
               "6HAY.cif",
    }
    _check(e3_type in {"CRBN", "VHL"}, "NOVEL_E3_TYPE_INVALID")
    sources = plan.get("sources")
    _check(isinstance(sources, dict), "NOVEL_SOURCES_REQUIRED")

    required_fields = {
        "expected_sha256", "label_asym_id", "auth_asym_id", "entity_id", "role",
        "description", "description_terms", "canonical_sequence", "sequence_length",
        "sequence_sha256", "sequence_source", "expected_deposit", "observed_entry_id",
        "excluded_proteins",
    }
    validated_sources: dict[str, dict] = {}
    for role, expected, fixed_path in (
        ("target", "target", fixed_paths["target"]),
        ("e3", e3_type, fixed_paths[e3_type]),
    ):
        declared = sources.get(role)
        _check(isinstance(declared, dict), "NOVEL_SOURCE_RECORD_INVALID:" + role)
        trusted_path = _trusted_regular_source(fixed_path, trusted_root)
        rebased = copy.deepcopy(declared)
        # Historical source.path is deliberately inert and is never resolved or opened.
        rebased["path"] = str(trusted_path)
        try:
            validated = novel_ternary._validate_source(rebased, expected)
        except Exception as error:
            raise ProtocolEvidenceError(
                "TRUSTED_NOVEL_SOURCE_VALIDATION_FAILED:" + role
            ) from error
        _check(isinstance(validated, dict) and required_fields <= set(validated),
               "TRUSTED_NOVEL_SOURCE_RESULT_INVALID:" + role)
        for key, value in validated.items():
            if key == "path":
                continue
            _check(key in declared and declared[key] == value,
                   "NOVEL_SOURCE_BINDING_MISMATCH:" + role + ":" + key)
        validated_sources[role] = validated
    return validated_sources


def _verify_plan_portable(plan: dict, expected_binding: dict, expected_graph: dict) -> None:
    _check(plan.get("format") == novel_ternary.PLAN_FORMAT, "NOVEL_PLAN_FORMAT_INVALID")
    supplied = plan.get("plan_digest")
    unsigned = copy.deepcopy(plan)
    unsigned.pop("plan_digest", None)
    _check(supplied == novel_ternary._digest(unsigned), "NOVEL_PLAN_DIGEST_INVALID")
    graph = plan.get("candidate_graph")
    _check(graph == expected_graph, "NOVEL_EXPECTED_GRAPH_MISMATCH")
    _check(novel_ternary._mapped_graph(graph) == graph, "NOVEL_MAPPED_GRAPH_INVALID")
    scope = plan.get("scientific_scope")
    _check(isinstance(scope, dict) and scope.get("reference_free") is True and
           scope.get("experimental_ternary_reference") is None,
           "NOVEL_SCIENTIFIC_SCOPE_INVALID")
    trusted_sources = _validated_trusted_plan_sources(plan, graph)
    job = plan.get("bindings", {}).get("job", {})
    for key in ("project", "job_id", "input_sha256", "result_sha256"):
        _check(job.get(key) == expected_binding[key], "NOVEL_JOB_BINDING_MISMATCH:" + key)
    _check(job.get("candidate_graph_sha256") == graph.get("graph_sha256"),
           "NOVEL_JOB_GRAPH_DIGEST_MISMATCH")
    _exact_int_list(plan.get("seeds"), SEEDS, "NOVEL_PLAN_SEEDS_OR_MSA_INVALID")
    _check(plan.get("msa_mode") == "server", "NOVEL_PLAN_SEEDS_OR_MSA_INVALID")
    boltz = plan.get("boltz_input", {})
    text = boltz.get("yaml")
    _check(isinstance(text, str) and _sha_bytes(text.encode("utf-8")) == boltz.get("sha256"),
           "NOVEL_YAML_HASH_MISMATCH")
    try:
        document = yaml.safe_load(text)
    except yaml.YAMLError as error:
        raise ProtocolEvidenceError("NOVEL_YAML_INVALID") from error
    expected_document = novel_ternary._document(
        trusted_sources["target"]["canonical_sequence"],
        trusted_sources["e3"]["canonical_sequence"],
        graph["actualmapped_smiles"], "server",
    )
    _check(document == expected_document, "NOVEL_YAML_CONTENT_MISMATCH")
    _check(boltz.get("templates") == [] and boltz.get("constraints") == [] and
           boltz.get("potentials") is False, "NOVEL_UNSAFE_PLAN_FEATURE")


def _novel_raw_paths(candidate_id: str, e3: str, seed: int, index: int) -> tuple[str, str]:
    base = f"evidence/novel-msa/raw/{candidate_id}--{e3}/seed-{seed}"
    return f"{base}/model-{index}.cif", f"{base}/confidence-{index}.json"


def _strip_paths(value: Any) -> Any:
    ignored = {"prediction", "confidence", "confidence_file", "selected_prediction",
               "raw_prediction_path", "raw_confidence_path", "view_directory"}
    if isinstance(value, dict):
        return {key: _strip_paths(child) for key, child in value.items() if key not in ignored}
    if isinstance(value, list):
        return [_strip_paths(child) for child in value]
    return value


def _positive_contact(value: Any) -> bool:
    if isinstance(value, bool):
        return False
    if isinstance(value, (int, float)):
        return math.isfinite(float(value)) and float(value) > 0
    if isinstance(value, list):
        return bool(value)
    if isinstance(value, dict):
        return _positive_contact(value.get("count"))
    return False


def _iqr5(values: list[float]) -> float:
    _check(len(values) == 5 and all(math.isfinite(value) for value in values),
           "NOVEL_IQR_INPUT_INVALID")
    ordered = sorted(values)
    return ordered[3] - ordered[1]


def _verify_novel(files: dict[str, Path], compact: dict, expected_binding: dict,
                  expected_graphs: dict) -> dict:
    protocol_name = "evidence/novel-msa/final-protocol.json"
    summary_name = "evidence/novel-msa/protocol-summary.json"
    protocol = _load_required(files, protocol_name)
    summary = _load_required(files, summary_name)
    _check(protocol.get("format") == "tpd-novel-msa-preregistered-protocol/1.0",
           "NOVEL_PROTOCOL_FORMAT_INVALID")
    _check(summary.get("format") == "tpd-novel-msa-protocol-summary/1.0",
           "NOVEL_SUMMARY_FORMAT_INVALID")
    unsigned = copy.deepcopy(protocol)
    supplied_digest = unsigned.pop("protocol_digest", None)
    _check(supplied_digest == _canonical_digest(unsigned), "FINAL_PROTOCOL_DIGEST_MISMATCH")
    _check(summary.get("protocol_digest") == supplied_digest,
           "NOVEL_SUMMARY_PROTOCOL_DIGEST_MISMATCH")
    _exact_int_list(protocol.get("seeds"), SEEDS, "NOVEL_PROTOCOL_SEED_MODEL_SET_INVALID")
    _exact_int_list(protocol.get("model_indices"), MODELS,
                    "NOVEL_PROTOCOL_SEED_MODEL_SET_INVALID")
    settings = protocol.get("settings", {})
    _check(isinstance(settings, dict), "NOVEL_SETTINGS_INVALID")
    for key, expected in NOVEL_SETTINGS.items():
        _strict_scalar(settings.get(key), expected, "NOVEL_SETTINGS_INVALID:" + key)
    _check(protocol.get("execution") == "serial", "NOVEL_EXECUTION_NOT_SERIAL")
    candidates = protocol.get("candidates")
    _check(isinstance(candidates, list) and len(candidates) == 2,
           "NOVEL_CANDIDATE_COUNT_INVALID")
    identities = {(item.get("candidate_id"), item.get("e3_type")) for item in candidates}
    _check(len(identities) == 2 and {item[1] for item in identities} == {"CRBN", "VHL"},
           "NOVEL_CANDIDATE_IDENTITIES_INVALID")
    declared_ids = {item[0] for item in identities}
    _check(declared_ids <= set(expected_graphs),
           "NOVEL_EXPECTED_CANDIDATE_SET_MISMATCH")

    summary_candidates = summary.get("candidates")
    _check(isinstance(summary_candidates, list) and len(summary_candidates) == 2,
           "NOVEL_SUMMARY_CANDIDATES_INVALID")
    _check(all(isinstance(item, dict) for item in summary_candidates),
           "NOVEL_SUMMARY_CANDIDATES_INVALID")
    summary_identities = {(item.get("candidate_id"), item.get("e3_type"))
                          for item in summary_candidates}
    _check(summary_identities == identities and len(summary_identities) == 2,
           "NOVEL_SUMMARY_CANDIDATE_IDENTITY_MISMATCH")
    summary_by_id = {item.get("candidate_id"): item for item in summary_candidates}
    compact_candidates = compact.get("novel_msa", {}).get("candidates")
    _check(isinstance(compact_candidates, list) and len(compact_candidates) == 2,
           "COMPACT_NOVEL_CANDIDATES_INVALID")
    _check(all(isinstance(item, dict) for item in compact_candidates),
           "COMPACT_NOVEL_CANDIDATES_INVALID")
    compact_identities = {(item.get("candidate_id"), item.get("e3_type"))
                          for item in compact_candidates}
    _check(compact_identities == identities and len(compact_identities) == 2,
           "COMPACT_NOVEL_CANDIDATE_IDENTITY_MISMATCH")
    compact_by_id = {item.get("candidate_id"): item for item in compact_candidates}

    result_candidates = []
    with tempfile.TemporaryDirectory(prefix="protocol-evidence-") as temporary:
        temporary_root = Path(temporary)
        for declared in candidates:
            candidate_id, e3 = declared["candidate_id"], declared["e3_type"]
            slug = f"{candidate_id}--{e3}"
            plan_name = f"evidence/novel-msa/plans/{slug}/plan.json"
            plan = _load_required(files, plan_name)
            _check(_sha_file(files[plan_name]) == declared.get("new_plan_sha256"),
                   "NOVEL_PLAN_FILE_HASH_MISMATCH")
            _check(plan.get("plan_digest") == declared.get("new_plan_digest"),
                   "NOVEL_PLAN_DECLARED_DIGEST_MISMATCH")
            _verify_plan_portable(plan, expected_binding, expected_graphs[candidate_id])
            _check(plan["candidate_graph"].get("candidate_id") == candidate_id and
                   plan["candidate_graph"].get("e3_type") == e3,
                   "NOVEL_PLAN_IDENTITY_MISMATCH")
            claimed = summary_by_id.get(candidate_id)
            compact_item = compact_by_id.get(candidate_id)
            _check(isinstance(claimed, dict) and isinstance(compact_item, dict),
                   "NOVEL_CANDIDATE_SUMMARY_MISSING")
            receipts = claimed.get("seeds")
            _check(isinstance(receipts, list) and len(receipts) == 5,
                   "NOVEL_RECEIPTS_INCOMPLETE")
            _check(all(isinstance(item, dict) and type(item.get("seed")) is int
                       for item in receipts), "NOVEL_RECEIPTS_INCOMPLETE")
            by_seed = {item["seed"]: item for item in receipts}
            _check(set(by_seed) == set(SEEDS), "NOVEL_RECEIPTS_INCOMPLETE")
            iptm: list[float] = []
            endpoint: list[float] = []
            positives: list[bool] = []
            selected_rows = []
            diagnostics_receipts = []
            for seed in SEEDS:
                receipt_name = f"evidence/novel-msa/raw/{slug}/seed-{seed}/receipt.json"
                receipt = _load_required(files, receipt_name)
                _compare_tree(receipt, by_seed[seed])
                _check(receipt.get("status") == "completed" and
                       type(receipt.get("exit_code")) is int and receipt["exit_code"] == 0 and
                       receipt.get("process_state") == "completed" and
                       receipt.get("actual_computation") is True and
                       receipt.get("scientific_approved") is False and
                       receipt.get("candidate_id") == candidate_id and
                       receipt.get("e3_type") == e3 and
                       type(receipt.get("seed")) is int and receipt["seed"] == seed,
                       "NOVEL_RECEIPT_BINDING_INVALID")
                receipt_settings = receipt.get("settings")
                _check(isinstance(receipt_settings, dict), "NOVEL_RECEIPT_SETTINGS_INVALID")
                for key, expected in NOVEL_SETTINGS.items():
                    _strict_scalar(receipt_settings.get(key), expected,
                                   "NOVEL_RECEIPT_SETTINGS_INVALID:" + key)
                _strict_scalar(receipt_settings.get("use_msa_server"), False,
                               "NOVEL_RECEIPT_SETTINGS_INVALID:use_msa_server")
                _verify_effective_command(receipt, seed, "NOVEL_EFFECTIVE_COMMAND_INVALID")
                models = receipt.get("models")
                _check(isinstance(models, list), "NOVEL_MODELS_REQUIRED")
                selected = _selected_index(models)
                selection = receipt.get("selection", {})
                _check(selection.get("selected_model_index") == selected,
                       "NOVEL_SELECTION_ORACLE_MISMATCH")
                selection_candidates = selection.get("candidates")
                _check(isinstance(selection_candidates, list) and
                       {(item.get("model_index"), item.get("confidence_score"))
                        for item in selection_candidates if isinstance(item, dict)} ==
                       {(model["model_index"], _finite_score(model)) for model in models},
                       "NOVEL_SELECTION_CANDIDATES_MISMATCH")
                selection_name = f"evidence/novel-msa/raw/{slug}/seed-{seed}/selection-receipt.json"
                selection_file = _load_required(files, selection_name)
                _check(selection_file == receipt.get("selection"),
                       "NOVEL_SELECTION_RECEIPT_MISMATCH")
                yaml_name = f"evidence/novel-msa/raw/{slug}/seed-{seed}/novel_ternary.yaml"
                _check(yaml_name in files and
                       _sha_file(files[yaml_name]) == plan["boltz_input"]["sha256"],
                       "NOVEL_INPUT_YAML_HASH_MISMATCH")
                settings_name = f"evidence/novel-msa/raw/{slug}/seed-{seed}/settings.json"
                settings_file = _load_required(files, settings_name)
                _check(settings_file.get("settings") == receipt.get("settings") and
                       settings_file.get("seed") == seed and
                       settings_file.get("plan_digest") == plan.get("plan_digest"),
                       "NOVEL_SETTINGS_FILE_MISMATCH")
                chosen = None
                for model in models:
                    index = model["model_index"]
                    cif_name, confidence_name = _novel_raw_paths(candidate_id, e3, seed, index)
                    _check(cif_name in files and confidence_name in files,
                           "NOVEL_RAW_MODEL_MISSING")
                    _check(_sha_file(files[cif_name]) == model.get("raw_prediction_sha256"),
                           "NOVEL_RAW_CIF_HASH_MISMATCH")
                    _check(_sha_file(files[confidence_name]) == model.get("raw_confidence_sha256"),
                           "NOVEL_RAW_CONFIDENCE_HASH_MISMATCH")
                    raw_confidence = _json_file(files[confidence_name], confidence_name)
                    _check(raw_confidence == model.get("confidence"),
                           "NOVEL_CONFIDENCE_CONTENT_MISMATCH")
                    if index == selected:
                        chosen = (model, cif_name, confidence_name)
                _check(chosen is not None, "NOVEL_SELECTED_MODEL_MISSING")
                model, cif_name, confidence_name = chosen
                view = temporary_root / slug / f"seed-{seed}" / "predictions" / "novel_ternary"
                view.mkdir(parents=True)
                copied_cif = view / Path(cif_name).name
                copied_confidence = view / Path(confidence_name).name
                shutil.copyfile(files[cif_name], copied_cif)
                shutil.copyfile(files[confidence_name], copied_confidence)
                _check(_sha_file(copied_cif) == _sha_file(files[cif_name]) and
                       _sha_file(copied_confidence) == _sha_file(files[confidence_name]),
                       "NOVEL_INSPECTION_COPY_HASH_MISMATCH")
                inspection = novel_ternary.assess_novel_outputs(view.parents[1], plan)
                _check(inspection.get("status") == "computed_hypothesis" and
                       inspection.get("actual_computation") is True,
                       "NOVEL_SELECTED_GEOMETRY_RECOMPUTE_FAILED")
                stored = receipt.get("selected_model_inspection")
                _check(isinstance(stored, dict), "NOVEL_SELECTED_INSPECTION_MISSING")
                _compare_tree(_strip_paths(inspection), _strip_paths(stored), tolerance=1e-6)
                score = inspection.get("model_confidence", {}).get("iptm")
                distance = inspection.get("descriptive_metrics", {}).get(
                    "linker_endpoint_distance_A")
                _check(not isinstance(score, bool) and isinstance(score, (int, float)) and
                       math.isfinite(float(score)) and
                       not isinstance(distance, bool) and isinstance(distance, (int, float)) and
                       math.isfinite(float(distance)) and float(distance) > 0,
                       "NOVEL_SELECTED_QUANTITATIVE_METRIC_INVALID")
                iptm.append(float(score))
                endpoint.append(float(distance))
                contacts = inspection["descriptive_metrics"][
                    "contacts_by_protein_role_and_ligand_group"]
                positives.append(_positive_contact(contacts.get("target", {}).get("warhead")) and
                                 _positive_contact(contacts.get("e3", {}).get("recruiter")))
                selected_rows.append({"seed": seed, "selected_model_index": selected})
                rebased = copy.deepcopy(receipt)
                rebased["selected_prediction"] = str(copied_cif)
                diagnostics_receipts.append(rebased)

            iptm_iqr, endpoint_iqr = _iqr5(iptm), _iqr5(endpoint)
            quantitative = iptm_iqr <= 0.20 and endpoint_iqr <= 1.5 and all(positives)
            quant = {
                "finite_selected_ipTM_values": iptm,
                "ipTM_IQR": iptm_iqr,
                "ipTM_IQR_maximum": 0.20,
                "finite_selected_endpoint_values_A": endpoint,
                "endpoint_IQR_A": endpoint_iqr,
                "endpoint_IQR_maximum_A": 1.5,
                "positive_target_warhead_and_e3_recruiter_contacts": positives,
                "quantitative_criteria_met": quantitative,
                "automatic_approval": False,
                "qualitative_topology_or_clash_acceptance_used": False,
            }
            stored_quant = claimed.get("old_expert_quantitative_metrics")
            _check(isinstance(stored_quant, dict), "NOVEL_QUANTITATIVE_SUMMARY_MISSING")
            _compare_tree(quant, stored_quant, tolerance=1e-6)
            _compare_tree(quant, compact_item.get("all_five_quantitative_metrics"), tolerance=1e-6)
            _check(compact_item.get("selected_top1_per_seed") == selected_rows,
                   "COMPACT_NOVEL_SELECTION_MISMATCH")

            try:
                from scripts.run_novel_msa_protocol import selected_diagnostics
            except Exception as error:
                raise ProtocolEvidenceError("TRUSTED_DIAGNOSTIC_IMPORT_FAILED") from error
            recomputed_diagnostics = selected_diagnostics(
                {"id": candidate_id, "e3": e3, "slug": slug, "plan": plan},
                diagnostics_receipts,
            )
            stored_diagnostics = claimed.get("geometric_diagnostics")
            _check(isinstance(stored_diagnostics, dict), "NOVEL_DIAGNOSTICS_MISSING")
            _compare_tree(_strip_paths(recomputed_diagnostics), _strip_paths(stored_diagnostics),
                          tolerance=1e-6)
            rows = recomputed_diagnostics.get("all_five_seeds_including_failures")
            pairs = recomputed_diagnostics.get("pairwise_contact_jaccard")
            clusters = recomputed_diagnostics.get("cluster_sensitivity")
            _check(isinstance(rows, list) and len(rows) == 5 and
                   isinstance(pairs, list) and len(pairs) == 10 and
                   isinstance(clusters, list) and len(clusters) == 3,
                   "NOVEL_RECOMPUTED_DIAGNOSTIC_DENOMINATOR_INVALID")
            def stats(key: str) -> dict:
                values = [float(item[key]) for item in pairs
                          if isinstance(item, dict) and
                          not isinstance(item.get(key), bool) and
                          isinstance(item.get(key), (int, float)) and
                          math.isfinite(float(item[key])) and 0 <= float(item[key]) <= 1]
                _check(len(values) == 10, "NOVEL_JACCARD_RANGE_INVALID:" + key)
                ordered = sorted(values)
                return {"min": ordered[0], "median": (ordered[4] + ordered[5]) / 2,
                        "max": ordered[-1]}
            projected_diagnostics = {
                "protein_protein_jaccard_stats": stats("protein_protein_jaccard"),
                "target_site_jaccard_stats": stats("target_site_jaccard"),
                "e3_site_jaccard_stats": stats("e3_site_jaccard"),
                "raw_under_2A_pair_count_all_five": [
                    item["interface_clashes"]["raw_under_2A_pair_count"] for item in rows],
                "vdw_pair_count_all_five": [
                    item["interface_clashes"]["vdw_pair_count"] for item in rows],
                "cluster_sensitivity": clusters,
            }
            _compare_tree(compact_item.get("geometric_diagnostics"), projected_diagnostics,
                          tolerance=1e-6)
            result_candidates.append({
                "candidate_id": candidate_id,
                "e3_type": e3,
                "selected_models": selected_rows,
                "selected_model_denominator": 5,
                "raw_model_denominator": 25,
                "failure_count": 0,
                "ipTM_IQR": iptm_iqr,
                "endpoint_IQR_A": endpoint_iqr,
                "positive_fragment_contacts_all_five": all(positives),
                "quantitative_criteria_met": quantitative,
                "geometric_diagnostics": projected_diagnostics,
            })

    _check(summary.get("raw_model_count") == 50 and
           summary.get("expected_raw_model_count") == 50 and
           summary.get("protocol_complete") is True,
           "NOVEL_SUMMARY_RAW_DENOMINATOR_INVALID")
    _check(compact.get("novel_msa", {}).get("raw_model_count") == 50 and
           compact.get("novel_msa", {}).get("completed_seed_count") == 10,
           "COMPACT_NOVEL_DENOMINATOR_INVALID")
    return {
        "protocol_id": "novel-msa:" + supplied_digest,
        "kind": "novel_reference_free_msa",
        "summary": {
            "candidate_count": 2,
            "candidate_denominator": 2,
            "selected_geometry_recomputed_count": 10,
            "selected_geometry_denominator": 10,
            "raw_model_count": 50,
            "raw_model_denominator": 50,
            "failure_count": 0,
            "candidates": result_candidates,
            "baseline30": "integrity hash verified; reported metrics not recomputed",
            "core_raw40": "integrity hash verified; reported metrics not recomputed",
            "known_calibration_comparison": "BRD4 diagnostic is related to the job but is not the SMARCA2 target",
        },
        "measurement_verification": {
            "verification_level": "raw hashes + confidence top1 + selected geometry recomputed",
            "raw_models_geometry_recomputed": 10,
            "raw_models_integrity_verified": 50,
            "limitations": [
                "MSA cache and checkpoint are not embedded.",
                "Original execution was not independently reproduced.",
                "The unsigned manifest is not expert proof.",
                "The joint MSA plus r10/d5 comparison does not establish MSA causation.",
                "No docking, efficacy, binding, degradation, or mechanism claim is made.",
            ],
        },
        "scientific_approved": False,
        "gates_affected": False,
    }


def verify_pack(pack_root: Path, *, expected_manifest_sha256: str,
                expected_binding: dict, expected_candidate_graphs: dict) -> dict:
    """Verify a portable evidence pack and return a normalized safe report.

    The caller-provided manifest digest is mandatory. The function never follows
    source provenance paths and never imports or executes source code from the
    evidence pack.
    """
    root = Path(pack_root)
    _validate_expected_inputs(expected_binding, expected_candidate_graphs)
    _manifest, files, normalized_files = _verify_manifest(root, expected_manifest_sha256)
    _check(COMPACT_NAME in files, "COMPACT_FILE_MISSING")
    compact = _load_required(files, COMPACT_NAME)
    _check(compact.get("format") == COMPACT_FORMAT, "COMPACT_FORMAT_INVALID")
    _approval_false(compact)
    _verify_provenance(compact, files)
    for name, path in files.items():
        if name.endswith(".json"):
            _approval_false(_json_file(path, name), name + ".")
    trusted_root = Path(__file__).resolve().parents[2]
    calibration = _verify_calibration(files, compact, trusted_root)
    novel = _verify_novel(files, compact, expected_binding, expected_candidate_graphs)
    return {
        "pack_manifest_sha256": expected_manifest_sha256,
        "compact_sha256": _sha_file(files[COMPACT_NAME]),
        "files": normalized_files,
        "protocols": [calibration, novel],
    }
