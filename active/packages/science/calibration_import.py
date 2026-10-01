"""Verification of expanded, source-bound CRBN 6BOY calibration receipts."""
from __future__ import annotations

import hashlib
import json
import tempfile
from pathlib import Path
from typing import Any

from packages.science import benchmark_distribution, boltz_worker

VERSION = "crbn-calibration-import/1"
ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "cases" / "acceptance_sources" / "crbn_calibration.json"
DESIGN_MANIFEST = ROOT / "cases" / "design_sources" / "manifest.json"
REFERENCE = ROOT / "cases" / "design_sources" / "6BOY.cif"
METADATA = ROOT / "cases" / "design_sources" / "crbn_benchmark.json"
EXPECTED_REFERENCE_SHA256 = "18eff4e76d65482f74beeff70d69e375cda0acb966e834a6a8a7f0d26f44b0d7"
EXPECTED_METADATA_SHA256 = "fbd191fb7cd5d716c63f9a1993b4238072f873e46807fb7b25b5dccba326bf70"
EXPECTED_EXTRA_SEEDS = {79, 97}


class CalibrationImportError(ValueError):
    pass


def _check(condition: bool, code: str) -> None:
    if not condition:
        raise CalibrationImportError(code)


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _regular(path: Path, code: str) -> Path:
    path = Path(path).absolute()
    current = path
    while True:
        if current.exists():
            _check(not current.is_symlink(), code + "_SYMLINK")
        if current.parent == current:
            break
        current = current.parent
    _check(path.is_file() and not path.is_symlink(), code)
    return path.resolve(strict=True)


def _inside(root: Path, value: Path, code: str, *, directory: bool = False) -> Path:
    root = root.resolve(strict=True)
    value = Path(value).absolute()
    current = value
    while True:
        if current.exists():
            _check(not current.is_symlink(), code + "_SYMLINK")
        if current == root or current.parent == current:
            break
        current = current.parent
    value = value.resolve(strict=True)
    _check(value == root or root in value.parents, code + "_ESCAPE")
    _check(value.is_dir() if directory else value.is_file(), code)
    _check(not value.is_symlink(), code + "_SYMLINK")
    return value


def _canonical_packet(value: Any) -> Any:
    if not isinstance(value, dict):
        return value
    # input_path is the sole location-dependent field emitted by prepare_input.
    # Preserve its presence while normalizing only its runtime-specific value;
    # every other source-binding field must compare exactly.
    return {
        key: "<runtime-input-path>" if key == "input_path" else child
        for key, child in value.items()
    }


def _manifest_binds_file(manifest: Any, relative_path: str, sha256: str) -> bool:
    if not isinstance(manifest, dict):
        return False
    files = manifest.get("files")
    return isinstance(files, dict) and files.get(relative_path) == sha256


def _verify_design_sources() -> tuple[str, str]:
    manifest_path = _regular(DESIGN_MANIFEST, "CALIBRATION_DESIGN_MANIFEST")
    reference_path = _regular(REFERENCE, "CALIBRATION_REFERENCE_FILE")
    metadata_path = _regular(METADATA, "CALIBRATION_METADATA_FILE")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CalibrationImportError("CALIBRATION_DESIGN_MANIFEST_JSON") from error
    reference_hash = _sha(reference_path)
    metadata_hash = _sha(metadata_path)
    _check(reference_hash == EXPECTED_REFERENCE_SHA256,
           "CALIBRATION_REFERENCE_KNOWN_HASH")
    _check(metadata_hash == EXPECTED_METADATA_SHA256,
           "CALIBRATION_METADATA_KNOWN_HASH")
    _check(_manifest_binds_file(
               manifest, "6BOY.cif", EXPECTED_REFERENCE_SHA256),
           "CALIBRATION_REFERENCE_MANIFEST_BINDING")
    _check(_manifest_binds_file(
               manifest, "crbn_benchmark.json", EXPECTED_METADATA_SHA256),
           "CALIBRATION_METADATA_MANIFEST_BINDING")
    return reference_hash, metadata_hash


def _canonical_inspection(value: Any) -> dict:
    _check(isinstance(value, dict), "CALIBRATION_INSPECTION_REQUIRED")
    comparison = value.get("comparison")
    _check(isinstance(comparison, dict), "CALIBRATION_COMPARISON_REQUIRED")
    mapping = comparison.get("mapping")
    metrics = comparison.get("metrics")
    confidence = value.get("model_confidence")
    _check(isinstance(mapping, dict), "CALIBRATION_CHAIN_MAPPING_REQUIRED")
    _check(isinstance(metrics, dict), "CALIBRATION_METRICS_REQUIRED")
    _check(isinstance(confidence, dict), "CALIBRATION_CONFIDENCE_REQUIRED")
    return {
        "status": value.get("status"),
        "mapping": mapping,
        "metrics": metrics,
        "model_confidence": confidence,
    }


def measurement_fingerprint() -> str:
    digest = hashlib.sha256()
    for relative in (
        "packages/science/calibration_import.py",
        "packages/science/benchmark_distribution.py",
        "packages/science/boltz_worker.py",
    ):
        path = ROOT / relative
        _check(path.is_file() and not path.is_symlink(),
               "CALIBRATION_IMPLEMENTATION_FILE_MISSING")
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def verify_distribution(additional_receipt_paths: list[str | Path]) -> dict:
    _regular(BASE, "CALIBRATION_BASE_FILE")
    reference_hash, metadata_hash = _verify_design_sources()
    base = benchmark_distribution._verified_base(BASE)
    base_rows = base["seed_receipts"]
    base_seeds = {row.get("seed") for row in base_rows}
    _check(base_seeds == {23, 41, 61}, "CALIBRATION_FROZEN_BASE_SEEDS")

    _check(isinstance(additional_receipt_paths, list) and
           len(additional_receipt_paths) == 2,
           "CALIBRATION_EXTRA_RECEIPT_COUNT")
    receipt_paths = [_regular(Path(value), "CALIBRATION_EXTRA_RECEIPT_FILE")
                     for value in additional_receipt_paths]
    _check(len(set(receipt_paths)) == 2, "CALIBRATION_EXTRA_RECEIPT_DUPLICATE")

    extras: list[dict] = []
    artifacts: list[tuple[str, Path]] = []
    seen: set[int] = set()
    for receipt_path in receipt_paths:
        try:
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as error:
            raise CalibrationImportError("CALIBRATION_RECEIPT_JSON") from error
        _check(isinstance(receipt, dict), "CALIBRATION_RECEIPT_OBJECT")
        seed = receipt.get("seed")
        _check(type(seed) is int and seed in EXPECTED_EXTRA_SEEDS and seed not in seen,
               "CALIBRATION_EXTRA_SEED")
        seen.add(seed)
        _check(receipt.get("format") == "tpd-boltz-crbn-calibration/1.0",
               "CALIBRATION_RECEIPT_FORMAT")
        _check(receipt.get("exit_code") == 0 and
               receipt.get("process_state") == "completed" and
               receipt.get("status") == "success",
               "CALIBRATION_PROCESS_NOT_COMPLETED")
        _check(receipt.get("seed_application", {}).get("actually_applied") is True and
               receipt.get("seed_application", {}).get("seed") == seed,
               "CALIBRATION_SEED_NOT_APPLIED")

        run_root = receipt_path.parent.resolve(strict=True)
        yaml_path = _inside(run_root, run_root / "6BOY_CRBN.yaml",
                            "CALIBRATION_INPUT_YAML")
        log_path = _inside(run_root, run_root / "inference.log",
                           "CALIBRATION_INFERENCE_LOG")
        output_root = _inside(run_root, run_root / "boltz_output",
                              "CALIBRATION_OUTPUT_DIRECTORY", directory=True)
        hashes = receipt.get("hashes")
        _check(isinstance(hashes, dict), "CALIBRATION_HASHES_REQUIRED")
        _check(_sha(yaml_path) == hashes.get("input_yaml_sha256") ==
               boltz_worker.FROZEN_CRBN_INPUT_SHA256,
               "CALIBRATION_INPUT_HASH")
        _check(_sha(log_path) == hashes.get("inference_log_sha256"),
               "CALIBRATION_LOG_HASH")
        _check(hashes.get("checkpoint_sha256") ==
               boltz_worker.FROZEN_CRBN_CHECKPOINT_SHA256,
               "CALIBRATION_CHECKPOINT_HASH")
        _check(hashes.get("reference_sha256") == reference_hash,
               "CALIBRATION_REFERENCE_HASH")
        recorded_tree = hashes.get("processed_output_files_sha256")
        _check(isinstance(recorded_tree, dict) and bool(recorded_tree),
               "CALIBRATION_OUTPUT_HASH_TREE_REQUIRED")
        actual_tree = boltz_worker._hash_tree(output_root)
        _check(actual_tree == recorded_tree,
               "CALIBRATION_OUTPUT_HASH_TREE")
        for field in ("msa_files_sha256", "processed_msa_files_sha256"):
            msa_hashes = hashes.get(field)
            _check(isinstance(msa_hashes, dict) and bool(msa_hashes),
                   "CALIBRATION_" + field.upper() + "_REQUIRED")
            _check(all(actual_tree.get(relative) == digest
                       for relative, digest in msa_hashes.items()),
                   "CALIBRATION_" + field.upper() + "_MISMATCH")

        input_record = receipt.get("input")
        _check(isinstance(input_record, dict), "CALIBRATION_INPUT_RECORD")
        msa_mode = input_record.get("msa_mode")
        with tempfile.TemporaryDirectory(prefix="crbn-calibration-") as temporary:
            prepared_path = Path(temporary) / "6BOY_CRBN.yaml"
            prepared = boltz_worker.prepare_input(
                REFERENCE, METADATA, prepared_path, msa_mode)
            _check(prepared_path.read_bytes() == yaml_path.read_bytes(),
                   "CALIBRATION_PREPARED_INPUT_MISMATCH")
            _check(_canonical_packet(prepared) == _canonical_packet(input_record),
                   "CALIBRATION_INPUT_SOURCE_BINDING")

        settings = receipt.get("settings")
        _check(isinstance(settings, dict) and
               settings.get("recycling_steps") == 3 and
               settings.get("sampling_steps") == 200 and
               settings.get("diffusion_samples") == 1 and
               settings.get("max_msa_seqs") == 256 and
               settings.get("max_parallel_samples") == 1 and
               settings.get("num_workers") == 0 and
               settings.get("preprocessing_threads") == 1 and
               settings.get("no_kernels") is True and
               settings.get("use_potentials") is False and
               settings.get("msa_mode") == msa_mode,
               "CALIBRATION_MODEL_SETTINGS")

        rerun = boltz_worker.assess_outputs(output_root, REFERENCE, input_record)
        _check(rerun.get("status") == "success",
               "CALIBRATION_REINSPECTION_FAILED")
        _check(_canonical_inspection(rerun) ==
               _canonical_inspection(receipt.get("inspection")),
               "CALIBRATION_REINSPECTION_MISMATCH")

        discovered = boltz_worker.discover_outputs(output_root)
        _check(discovered.get("status") == "located",
               "CALIBRATION_OUTPUT_DISCOVERY")
        prediction = _inside(output_root, Path(discovered["prediction"]),
                             "CALIBRATION_PREDICTION")
        confidence = _inside(output_root, Path(discovered["confidence"]),
                             "CALIBRATION_CONFIDENCE")
        verified = dict(receipt)
        verified["inspection"] = rerun
        verified["actual_execution_success"] = True
        verified["provenance"] = {"sha256": _sha(receipt_path)}
        extras.append(verified)
        artifacts.extend((
            (f"seed-{seed}-receipt", receipt_path),
            (f"seed-{seed}-yaml", yaml_path),
            (f"seed-{seed}-log", log_path),
            (f"seed-{seed}-prediction", prediction),
            (f"seed-{seed}-confidence", confidence),
        ))

    _check(seen == EXPECTED_EXTRA_SEEDS, "CALIBRATION_EXTRA_SEED_COVERAGE")
    distribution = benchmark_distribution.build_benchmark_distribution(BASE, extras)
    _check(distribution.get("seed_count") == 5 and
           {row.get("seed") for row in distribution.get("seeds", [])} ==
           {23, 41, 61, 79, 97},
           "CALIBRATION_DISTRIBUTION_COVERAGE")
    _check(distribution.get("technical_success_count") == 5 and
           distribution.get("scientific_model_quality_success_count") is None and
           distribution.get("expert_geometric_success_cutoff") is None,
           "CALIBRATION_DISTRIBUTION_AUTHORITY")

    artifacts.extend((
        ("frozen-base-manifest", BASE),
        ("design-source-manifest", DESIGN_MANIFEST),
        ("6BOY-reference", REFERENCE),
        ("6BOY-metadata", METADATA),
    ))
    for row in base_rows:
        artifacts.append((f"seed-{row['seed']}-original-receipt",
                          BASE.parent / row["original_receipt_ref"]))
        if row.get("reinspection_ref"):
            artifacts.append((f"seed-{row['seed']}-reinspection",
                              BASE.parent / row["reinspection_ref"]))
    return {
        "distribution": distribution,
        "artifacts": artifacts,
        "base_manifest_sha256": benchmark_distribution.KNOWN_BASE_SHA256,
        "reference_sha256": reference_hash,
        "metadata_sha256": metadata_hash,
        "measurement_version": VERSION,
        "measurement_fingerprint": measurement_fingerprint(),
    }
