"""Reinspect an immutable completed CRBN calibration without another GPU call."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path, PurePosixPath
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from packages.science.boltz_worker import assess_outputs, sha256_file


def _regular(path: Path, label: str) -> Path:
    path = Path(path)
    if not path.is_file() or path.is_symlink():
        raise ValueError(label + "_REGULAR_FILE_REQUIRED")
    return path


def _source_hash() -> str:
    import packages.science.boltz_worker as worker
    return sha256_file(Path(worker.__file__))


def _verify_hash(path: Path, expected: str, label: str) -> str:
    actual = sha256_file(_regular(path, label))
    if actual != expected:
        raise ValueError(label + "_HASH_MISMATCH")
    return actual


def _input_packet(receipt: dict[str, Any], input_yaml: Path) -> dict[str, Any]:
    document = yaml.safe_load(input_yaml.read_text(encoding="utf-8"))
    sequences = document.get("sequences") if isinstance(document, dict) else None
    if not isinstance(sequences, list):
        raise ValueError("INPUT_YAML_INVALID")
    by_id: dict[str, dict] = {}
    for item in sequences:
        if not isinstance(item, dict) or len(item) != 1:
            raise ValueError("INPUT_YAML_INVALID")
        material = next(iter(item.values()))
        if not isinstance(material, dict) or not isinstance(material.get("id"), str):
            raise ValueError("INPUT_YAML_INVALID")
        by_id[material["id"]] = material
    original = receipt.get("input")
    if not isinstance(original, dict):
        raise ValueError("ORIGINAL_INPUT_PACKET_MISSING")
    roles = original.get("input_chain_roles")
    if not isinstance(roles, dict):
        raise ValueError("ORIGINAL_INPUT_ROLES_MISSING")
    proteins = {}
    for role in ("target", "e3"):
        chain = roles.get(role)
        material = by_id.get(chain)
        if not isinstance(material, dict) or not isinstance(material.get("sequence"), str):
            raise ValueError("INPUT_YAML_PROTEIN_MISSING")
        proteins[role] = material["sequence"]
        if proteins[role] != original.get("protein_sequences", {}).get(role):
            raise ValueError("INPUT_YAML_RECEIPT_SEQUENCE_MISMATCH")
    ligand_chain = roles.get("ligand")
    if by_id.get(ligand_chain, {}).get("ccd") != original.get("ligand", {}).get("ccd"):
        raise ValueError("INPUT_YAML_RECEIPT_LIGAND_MISMATCH")
    packet = dict(original)
    packet["protein_sequences"] = proteins
    return packet


def _verify_processed_outputs(receipt: dict[str, Any], output_dir: Path) -> dict[str, str]:
    hashes = receipt.get("hashes")
    if not isinstance(hashes, dict):
        raise ValueError("ORIGINAL_HASHES_MISSING")
    expected = hashes.get("processed_output_files_sha256")
    if not isinstance(expected, dict) or not expected:
        raise ValueError("ORIGINAL_OUTPUT_HASHES_MISSING")

    unresolved_root = Path(output_dir)
    if unresolved_root.is_symlink():
        raise ValueError("OUTPUT_DIR_SYMLINK_FORBIDDEN")
    if not unresolved_root.is_dir():
        raise ValueError("OUTPUT_DIR_INVALID")
    root = unresolved_root.resolve(strict=True)

    current: dict[str, Path] = {}
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ValueError("OUTPUT_SYMLINK_FORBIDDEN")
        if path.is_file():
            current[path.relative_to(root).as_posix()] = path
        elif not path.is_dir():
            raise ValueError("OUTPUT_TREE_ENTRY_INVALID")

    normalized: dict[str, str] = {}
    for stored_path, digest in expected.items():
        if not isinstance(stored_path, str) or not stored_path:
            raise ValueError("OUTPUT_PATH_INVALID")
        if "\\" in stored_path:
            raise ValueError("OUTPUT_PATH_INVALID")
        relative_path = PurePosixPath(stored_path)
        if (relative_path.is_absolute() or relative_path.as_posix() in {"", "."}
                or ".." in relative_path.parts
                or relative_path.as_posix() != stored_path):
            raise ValueError("OUTPUT_PATH_TRAVERSAL")
        candidate = root.joinpath(*relative_path.parts).resolve(strict=False)
        try:
            candidate.relative_to(root)
        except ValueError as exc:
            raise ValueError("OUTPUT_PATH_TRAVERSAL") from exc
        if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
            raise ValueError("ORIGINAL_OUTPUT_HASH_INVALID")
        normalized[relative_path.as_posix()] = digest

    if set(current) != set(normalized):
        raise ValueError("OUTPUT_FILE_SET_MISMATCH")
    for relative, expected_hash in normalized.items():
        if sha256_file(current[relative]) != expected_hash:
            raise ValueError("OUTPUT_HASH_MISMATCH:" + relative)
    return normalized


def _execution_metadata(receipt: dict[str, Any]) -> dict[str, Any]:
    seed = receipt.get("seed")
    seed_application = receipt.get("seed_application")
    if type(seed) is not int or seed < 0:
        raise ValueError("CONFIRMED_EXECUTION_METADATA_MISSING")
    if (not isinstance(seed_application, dict)
            or seed_application.get("actually_applied") is not True):
        raise ValueError("CONFIRMED_SEED_APPLICATION_MISSING")
    applied_seed = seed_application.get("seed")
    if type(applied_seed) is not int or applied_seed != seed:
        raise ValueError("CONFIRMED_SEED_APPLICATION_MISMATCH")

    settings = receipt.get("settings")
    environment = receipt.get("tool_environment")
    if not isinstance(settings, dict) or not isinstance(environment, dict):
        raise ValueError("CONFIRMED_EXECUTION_METADATA_MISSING")
    model_version = environment.get("boltz_version")
    cuda_available = environment.get("cuda_available")
    cuda_device_count = environment.get("cuda_device_count")
    cuda_devices = environment.get("cuda_devices")
    if (not isinstance(model_version, str) or not model_version
            or type(cuda_available) is not bool
            or type(cuda_device_count) is not int or cuda_device_count < 0
            or not isinstance(cuda_devices, list)
            or not all(isinstance(device, str) for device in cuda_devices)
            or "cuda_runtime" not in environment):
        raise ValueError("CONFIRMED_EXECUTION_METADATA_MISSING")
    return {
        "applied_seed": applied_seed,
        "seed_application_method": seed_application.get("method"),
        "applied_settings": settings,
        "model_version": model_version,
        "torch_version": environment.get("torch_version"),
        "cuda": {
            "available": cuda_available,
            "runtime": environment.get("cuda_runtime"),
            "device_count": cuda_device_count,
            "devices": cuda_devices,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", required=True, type=Path)
    parser.add_argument("--input-yaml", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--ligand-groups", type=Path)
    args = parser.parse_args()

    receipt_path = _regular(args.receipt, "ORIGINAL_RECEIPT")
    evidence_path = Path(args.evidence)
    if evidence_path.exists():
        raise ValueError("EVIDENCE_ALREADY_EXISTS")
    if evidence_path.resolve() == receipt_path.resolve():
        raise ValueError("ORIGINAL_RECEIPT_OVERWRITE_FORBIDDEN")
    receipt_bytes = receipt_path.read_bytes()
    receipt = json.loads(receipt_bytes)
    if receipt.get("exit_code") != 0 or receipt.get("process_state") != "completed":
        raise ValueError("ORIGINAL_EXECUTION_NOT_COMPLETED_EXIT0")
    hashes = receipt.get("hashes")
    if not isinstance(hashes, dict):
        raise ValueError("ORIGINAL_HASHES_MISSING")

    verified = {
        "input_yaml_sha256": _verify_hash(args.input_yaml, hashes["input_yaml_sha256"], "INPUT_YAML"),
        "checkpoint_sha256": _verify_hash(args.checkpoint, hashes["checkpoint_sha256"], "CHECKPOINT"),
        "reference_sha256": _verify_hash(args.reference, hashes["reference_sha256"], "REFERENCE"),
        "provenance_sha256": sha256_file(_regular(args.provenance, "PROVENANCE")),
    }
    provenance = json.loads(args.provenance.read_text(encoding="utf-8"))
    if not isinstance(provenance, dict):
        raise ValueError("PROVENANCE_INVALID")
    packet = _input_packet(receipt, args.input_yaml)
    verified_outputs = _verify_processed_outputs(receipt, args.output_dir)
    output_root = Path(args.output_dir).resolve(strict=True)
    try:
        evidence_path.resolve(strict=False).relative_to(output_root)
    except ValueError:
        pass
    else:
        raise ValueError("EVIDENCE_INSIDE_ORIGINAL_OUTPUT_FORBIDDEN")
    execution_metadata = _execution_metadata(receipt)
    assessment = assess_outputs(args.output_dir, args.reference, packet, args.ligand_groups)
    source_hash = _source_hash()
    evidence = {
        "format": "tpd-boltz-crbn-reinspection/1.0",
        "actual_execution_success": True,
        "structural_assessment_completed": assessment.get("status") == "success",
        "assessment": assessment,
        "original_receipt_path": str(receipt_path),
        "original_receipt_sha256": hashlib.sha256(receipt_bytes).hexdigest(),
        "original_status": receipt.get("status"),
        "original_inspection": receipt.get("inspection"),
        "original_failure_reason": receipt.get("failure_reason"),
        "original_exit_code": receipt.get("exit_code"),
        "original_process_state": receipt.get("process_state"),
        "original_elapsed_seconds": receipt.get("elapsed_seconds"),
        "verified_original_recorded_hashes": verified,
        "verified_processed_output_files_sha256": verified_outputs,
        "confirmed_execution_metadata": execution_metadata,
        "provenance": provenance,
        "provenance_authority": "external_context_only_not_authoritative_for_execution_or_output_integrity",
        "source_code_sha256": source_hash,
        "no_gpu_call_performed": True,
    }
    evidence_path.parent.mkdir(parents=True, exist_ok=True)
    with evidence_path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(evidence, stream, ensure_ascii=False, indent=2, sort_keys=True)
        stream.write("\n")
    return 0 if evidence["structural_assessment_completed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
