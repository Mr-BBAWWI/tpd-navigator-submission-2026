"""Deterministic, verified evidence assembly for research campaign agents.

This module validates exported job receipts and their complete artifact closure. It
never executes artifact content and never grants scientific or human approval.
"""
from __future__ import annotations

import copy
import datetime as _datetime
import hashlib
import json
import math
import os
import re
import stat
from pathlib import Path
from typing import Any, Iterable

from rdkit import Chem

FORMAT = "research-campaign-evidence/1"
_MAX_JSON_BYTES = 64 * 1024 * 1024
_MAX_MANIFEST_BYTES = 16 * 1024 * 1024
_MAX_BLOB_BYTES = 128 * 1024 * 1024
_MAX_EXPORT_BYTES = 1024 * 1024 * 1024
_MAX_JSON_DEPTH = 128
_MAX_JSON_NODES = 2_000_000
_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,191}$")
_ARTIFACT_ID_RE = re.compile(r"^a-[a-f0-9]{32}$")
_JOB_ID_RE = re.compile(r"^job-[A-Za-z0-9][A-Za-z0-9._-]{0,186}$")
_PARENT_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9-]{0,191}$")
_HASH_RE = re.compile(r"^[0-9a-f]{64}$")
_REF_FIELDS = ("artifact_id", "version", "sha256", "media_type", "schema_id", "provenance")
_EXPECTED_CRITERION_IDS = (
    "parent_funnel",
    "expert_parent_selection",
    "modifiable_sites",
    "distinct_constitutional_graphs",
    "actual_broad_families",
    "qualified_panel",
    "both_e3_assembly",
    "core_interaction_preservation",
    "microstates_h_direction",
    "known_crbn_calibration",
    "novel_ternary_repeats",
    "novel_ternary_geometry",
    "exact_synthesis_review",
    "formal_expert_decision",
)
_CRITERION_STATUSES = {"pass", "failed", "pending", "blocked"}


class EvidenceValidationError(ValueError):
    """Raised when an evidence input cannot be safely and completely verified."""


def _reject_constant(value: str) -> None:
    raise EvidenceValidationError(f"nonfinite JSON number is not allowed: {value}")


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise EvidenceValidationError(f"duplicate JSON object key: {key!r}")
        result[key] = value
    return result


def _validate_json_shape(value: Any, label: str) -> None:
    nodes = 0
    stack: list[tuple[Any, int]] = [(value, 0)]
    while stack:
        current, depth = stack.pop()
        nodes += 1
        if nodes > _MAX_JSON_NODES:
            raise EvidenceValidationError(f"{label} exceeds JSON node limit {_MAX_JSON_NODES}")
        if depth > _MAX_JSON_DEPTH:
            raise EvidenceValidationError(f"{label} exceeds JSON depth limit {_MAX_JSON_DEPTH}")
        if isinstance(current, float) and not math.isfinite(current):
            raise EvidenceValidationError(f"{label} contains a nonfinite number")
        if isinstance(current, dict):
            stack.extend((item, depth + 1) for item in current.values())
        elif isinstance(current, list):
            stack.extend((item, depth + 1) for item in current)
        elif current is not None and type(current) not in (str, int, float, bool):
            raise EvidenceValidationError(f"{label} contains unsupported JSON value {type(current).__name__}")


def _loads_json(raw: bytes, label: str, *, limit: int = _MAX_JSON_BYTES) -> Any:
    if len(raw) > limit:
        raise EvidenceValidationError(f"{label} is too large: {len(raw)} bytes (limit {limit})")
    try:
        value = json.loads(
            raw.decode("utf-8"),
            parse_constant=_reject_constant,
            object_pairs_hook=_unique_object,
        )
    except UnicodeDecodeError as exc:
        raise EvidenceValidationError(f"{label} is not valid UTF-8") from exc
    except json.JSONDecodeError as exc:
        raise EvidenceValidationError(f"invalid JSON in {label}: {exc.msg} at line {exc.lineno}") from exc
    _validate_json_shape(value, label)
    return value


def _read_regular(path: Path, label: str, limit: int) -> bytes:
    try:
        info = path.lstat()
    except FileNotFoundError as exc:
        raise EvidenceValidationError(f"missing {label}: {path.name}") from exc
    if stat.S_ISLNK(info.st_mode):
        raise EvidenceValidationError(f"symlink is not allowed for {label}: {path.name}")
    if not stat.S_ISREG(info.st_mode):
        raise EvidenceValidationError(f"{label} is not a regular file: {path.name}")
    if info.st_size > limit:
        raise EvidenceValidationError(f"{label} is too large: {info.st_size} bytes (limit {limit})")
    return path.read_bytes()


def _scan_export_tree(root: Path) -> None:
    try:
        root_info = root.lstat()
    except FileNotFoundError as exc:
        raise EvidenceValidationError(f"parent export directory does not exist: {root}") from exc
    if stat.S_ISLNK(root_info.st_mode) or not stat.S_ISDIR(root_info.st_mode):
        raise EvidenceValidationError(f"parent export must be a real directory, not a symlink: {root}")
    total = 0
    for current, directories, files in os.walk(root, followlinks=False):
        current_path = Path(current)
        for name in directories:
            path = current_path / name
            info = path.lstat()
            if stat.S_ISLNK(info.st_mode):
                raise EvidenceValidationError(f"symlink directory is not allowed in export: {name}")
            if not stat.S_ISDIR(info.st_mode):
                raise EvidenceValidationError(f"non-directory entry encountered as directory: {name}")
        for name in files:
            path = current_path / name
            info = path.lstat()
            if stat.S_ISLNK(info.st_mode):
                raise EvidenceValidationError(f"symlink file is not allowed in export: {name}")
            if not stat.S_ISREG(info.st_mode):
                raise EvidenceValidationError(f"non-regular file is not allowed in export: {name}")
            total += info.st_size
            if total > _MAX_EXPORT_BYTES:
                raise EvidenceValidationError(
                    f"export exceeds aggregate size limit {_MAX_EXPORT_BYTES} bytes"
                )


def _sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def _require_identifier(value: Any, label: str) -> str:
    if not isinstance(value, str) or not _ID_RE.fullmatch(value):
        raise EvidenceValidationError(f"invalid {label}: {value!r}")
    if value in {".", ".."} or "/" in value or "\\" in value or ":" in value:
        raise EvidenceValidationError(f"non-portable or unsafe {label}: {value!r}")
    return value


def _require_typed_identifier(value: Any, label: str, pattern: re.Pattern[str]) -> str:
    if not isinstance(value, str) or not pattern.fullmatch(value):
        raise EvidenceValidationError(f"invalid {label}: {value!r}")
    return value


def _require_hash(value: Any, label: str) -> str:
    if not isinstance(value, str) or not _HASH_RE.fullmatch(value):
        raise EvidenceValidationError(f"invalid lowercase SHA-256 for {label}: {value!r}")
    return value


def _artifact_metadata(value: Any, label: str, *, complete: bool) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise EvidenceValidationError(f"{label} must be an object")
    artifact_id = _require_typed_identifier(
        value.get("artifact_id"), f"{label}.artifact_id", _ARTIFACT_ID_RE
    )
    metadata: dict[str, Any] = {"artifact_id": artifact_id}
    if complete or "version" in value:
        version = value.get("version")
        if type(version) is not int or version <= 0:
            raise EvidenceValidationError(f"{label}.version must be a positive integer")
        metadata["version"] = version
    if complete or "sha256" in value:
        metadata["sha256"] = _require_hash(value.get("sha256"), f"{label}.sha256")
    if complete or "media_type" in value:
        media_type = value.get("media_type")
        if not isinstance(media_type, str) or not media_type or len(media_type) > 255:
            raise EvidenceValidationError(f"{label}.media_type must be a nonempty string")
        metadata["media_type"] = media_type
    for field in ("schema_id", "provenance"):
        if complete and field not in value:
            raise EvidenceValidationError(f"{label}.{field} is required")
        if field in value:
            item = value[field]
            if item is not None and not isinstance(item, (str, dict)):
                raise EvidenceValidationError(f"{label}.{field} has invalid type")
            metadata[field] = item
    return metadata


def _manifest_artifacts(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    value = manifest.get("artifacts")
    if isinstance(value, dict):
        records = []
        for key, item in value.items():
            if not isinstance(item, dict):
                raise EvidenceValidationError(f"manifest artifact {key!r} must be an object")
            item = dict(item)
            item.setdefault("artifact_id", key)
            records.append(item)
        return records
    if isinstance(value, list):
        return value
    raise EvidenceValidationError("manifest.artifacts must be a list or object")


def _extract_ids(value: dict[str, Any], label: str) -> tuple[str, str]:
    job_id = value.get("job_id")
    if job_id is None and label == "receipt":
        job_id = value.get("id")
    parent_id = value.get("parent_id")
    parameters = value.get("parameters")
    if parent_id is None and isinstance(parameters, dict):
        parent_id = parameters.get("parent_id")
    binding = value.get("binding")
    if parent_id is None and isinstance(binding, dict):
        parent_id = binding.get("parent_id")
    return (
        _require_typed_identifier(job_id, f"{label}.job_id", _JOB_ID_RE),
        _require_typed_identifier(parent_id, f"{label}.parent_id", _PARENT_ID_RE),
    )


def _walk(value: Any) -> Iterable[tuple[str | None, Any]]:
    stack: list[tuple[str | None, Any]] = [(None, value)]
    while stack:
        key, current = stack.pop()
        yield key, current
        if isinstance(current, dict):
            for child_key, child in reversed(list(current.items())):
                stack.append((child_key, child))
        elif isinstance(current, list):
            for child in reversed(current):
                stack.append((key, child))


def _artifact_refs(value: Any, label: str) -> list[dict[str, Any]]:
    refs: list[dict[str, Any]] = []
    for _, current in _walk(value):
        if isinstance(current, dict) and "artifact_id" in current:
            refs.append(_artifact_metadata(current, f"{label} artifact reference", complete=True))
    return refs


def _metadata_matches(reference: dict[str, Any], manifest: dict[str, Any]) -> bool:
    for field in _REF_FIELDS:
        if field in reference and reference[field] != manifest.get(field):
            return False
    return True


def _safe_blob_path(root: Path, artifact_id: str) -> Path:
    blobs = root / "blobs"
    path = blobs / artifact_id
    try:
        blobs_resolved = blobs.resolve(strict=True)
        path_parent = path.parent.resolve(strict=True)
    except FileNotFoundError as exc:
        raise EvidenceValidationError("export is missing blobs directory") from exc
    if path_parent != blobs_resolved:
        raise EvidenceValidationError(f"artifact path escapes blobs directory: {artifact_id!r}")
    return path


def _verify_export(root: Path) -> dict[str, Any]:
    _scan_export_tree(root)
    manifest_raw = _read_regular(root / "manifest.json", "manifest.json", _MAX_MANIFEST_BYTES)
    manifest = _loads_json(manifest_raw, "manifest.json", limit=_MAX_MANIFEST_BYTES)
    if not isinstance(manifest, dict):
        raise EvidenceValidationError("manifest.json must contain an object")
    receipt_raw = _read_regular(root / "job-receipt.json", "job-receipt.json", _MAX_JSON_BYTES)
    expected_receipt_hash = _require_hash(
        manifest.get("raw_receipt_sha256"), "manifest.raw_receipt_sha256"
    )
    if _sha256(receipt_raw) != expected_receipt_hash:
        raise EvidenceValidationError("job-receipt.json does not match manifest.raw_receipt_sha256")
    receipt = _loads_json(receipt_raw, "job-receipt.json")
    if not isinstance(receipt, dict):
        raise EvidenceValidationError("job-receipt.json must contain an object")
    manifest_job, manifest_parent = _extract_ids(manifest, "manifest")
    receipt_job, receipt_parent = _extract_ids(receipt, "receipt")
    if (manifest_job, manifest_parent) != (receipt_job, receipt_parent):
        raise EvidenceValidationError(
            "manifest and receipt job_id/parent_id must match exactly: "
            f"{manifest_job!r}/{manifest_parent!r} != {receipt_job!r}/{receipt_parent!r}"
        )

    artifacts: dict[str, dict[str, Any]] = {}
    loaded_json: dict[str, Any] = {}
    for index, record in enumerate(_manifest_artifacts(manifest)):
        metadata = _artifact_metadata(record, f"manifest.artifacts[{index}]", complete=True)
        artifact_id = metadata["artifact_id"]
        if artifact_id in artifacts:
            raise EvidenceValidationError(
                f"duplicate manifest artifact_id/version: {artifact_id}/{metadata['version']}"
            )
        artifacts[artifact_id] = metadata
        path = _safe_blob_path(root, artifact_id)
        raw = _read_regular(path, f"blob {artifact_id}", _MAX_BLOB_BYTES)
        if _sha256(raw) != metadata["sha256"]:
            raise EvidenceValidationError(f"blob hash mismatch for artifact {artifact_id}")
        if metadata["media_type"].split(";", 1)[0].strip().lower() == "application/json":
            loaded_json[artifact_id] = _loads_json(raw, f"blob {artifact_id}")

    seen_refs: dict[str, dict[str, Any]] = {}
    sources = [("job-receipt.json", receipt)] + [
        (f"blob {artifact_id}", loaded_json[artifact_id]) for artifact_id in sorted(loaded_json)
    ]
    for label, document in sources:
        for reference in _artifact_refs(document, label):
            artifact_id = reference["artifact_id"]
            prior = seen_refs.get(artifact_id)
            if prior is not None and any(
                field in prior and field in reference and prior[field] != reference[field]
                for field in _REF_FIELDS
            ):
                raise EvidenceValidationError(
                    f"contradictory references for artifact {artifact_id}"
                )
            seen_refs.setdefault(artifact_id, reference)
            metadata = artifacts.get(artifact_id)
            if metadata is None:
                raise EvidenceValidationError(
                    f"incomplete export: referenced artifact {artifact_id} is absent from manifest"
                )
            if not _metadata_matches(reference, metadata):
                raise EvidenceValidationError(
                    f"artifact reference metadata does not match manifest for {artifact_id}"
                )

    return {
        "root": root,
        "manifest": manifest,
        "manifest_sha256": _sha256(manifest_raw),
        "receipt": receipt,
        "receipt_sha256": expected_receipt_hash,
        "job_id": receipt_job,
        "parent_id": receipt_parent,
        "artifacts": artifacts,
        "json_blobs": loaded_json,
    }


def _canonical_smiles(value: Any) -> str | None:
    if not isinstance(value, str) or not value:
        return None
    molecule = Chem.MolFromSmiles(value)
    if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
        return None
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def _constitutional_smiles(value: Any) -> str | None:
    if not isinstance(value, str) or not value:
        return None
    molecule = Chem.MolFromSmiles(value)
    if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
        return None
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for bond in molecule.GetBonds():
        bond.SetStereo(Chem.BondStereo.STEREONONE)
        bond.SetBondDir(Chem.BondDir.NONE)
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)


def _record_smiles(record: dict[str, Any]) -> str | None:
    identity = record.get("identity")
    values = []
    if isinstance(identity, dict):
        values.append(identity.get("canonical_isomeric_smiles"))
    values.extend((record.get("canonical_smiles"), record.get("constitutional_smiles")))
    for value in values:
        canonical = _canonical_smiles(value)
        if canonical is not None:
            return canonical
    return None


def _record_list(value: Any, label: str) -> list[dict[str, Any]]:
    if not isinstance(value, list):
        raise EvidenceValidationError(f"{label} must be a list")
    if any(not isinstance(item, dict) for item in value):
        raise EvidenceValidationError(f"{label} entries must be objects")
    return value


def _sites(value: Any, label: str) -> set[int]:
    if isinstance(value, dict):
        rows = _record_list(value.get("atoms"), f"{label}.atoms")
    else:
        rows = _record_list(value, label)
    result: set[int] = set()
    for item in rows:
        if item.get("state") != "MODIFIABLE":
            continue
        atom_map = item.get("atom_map")
        if type(atom_map) is int and atom_map > 0:
            result.add(atom_map)
    return result


def _site_maps(record: dict[str, Any]) -> list[int]:
    values = record.get("attachment_site_atom_maps", record.get("modified_atom_maps"))
    if not isinstance(values, list):
        return []
    return sorted({item for item in values if type(item) is int and item > 0})


def _pose_details(record: dict[str, Any]) -> tuple[bool, dict[str, Any], Any]:
    docking = record.get("docking")
    if not isinstance(docking, dict):
        return False, {}, None
    status = docking.get("status")
    passing = docking.get("passing_pose_count")
    pose_qualified = (
        isinstance(status, str)
        and status.startswith("completed")
        and docking.get("pose_preserved") is True
        and type(passing) is int
        and passing > 0
    )
    metrics: dict[str, Any] = {}
    for key in (
        "best_core_rmsd_A",
        "pose_count",
        "passing_pose_count",
        "score",
        "best_score",
        "affinity",
    ):
        value = docking.get(key)
        if type(value) in (int, float) and not isinstance(value, bool):
            if isinstance(value, float) and not math.isfinite(value):
                raise EvidenceValidationError(f"nonfinite docking metric {key}")
            metrics[key] = value
    return pose_qualified, metrics, docking.get("error") or docking.get("failure")


def _broad_family(value: Any) -> str | None:
    if value in {"ring_expansion", "ring_contraction"}:
        return "ring_modification"
    return value if isinstance(value, str) and value else None


def _compact_diagnostic(value: Any, depth: int = 0) -> Any:
    if value is None or type(value) in (str, int, bool):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise EvidenceValidationError("nonfinite supplied summary diagnostic")
        return value
    if depth >= 2:
        if isinstance(value, list):
            return {"omitted_nested_items": len(value)}
        if isinstance(value, dict):
            return {"omitted_nested_keys": len(value)}
    if isinstance(value, list):
        return [_compact_diagnostic(item, depth + 1) for item in value[:100]]
    if isinstance(value, dict):
        return {
            str(key): _compact_diagnostic(item, depth + 1)
            for key, item in sorted(value.items())[:100]
        }
    return None


def _criteria_snapshot(assessment: dict[str, Any]) -> list[dict[str, Any]]:
    criteria = assessment.get("criteria")
    if not isinstance(criteria, list):
        raise EvidenceValidationError("assessment.criteria must be a list")
    if len(criteria) != len(_EXPECTED_CRITERION_IDS):
        raise EvidenceValidationError("assessment.criteria must contain exactly 14 criteria")
    ids: list[str] = []
    for item in criteria:
        if not isinstance(item, dict):
            raise EvidenceValidationError("assessment criterion must be an object")
        criterion_id = item.get("id")
        if not isinstance(criterion_id, str):
            raise EvidenceValidationError("assessment criterion id must be a string")
        ids.append(criterion_id)
        status = item.get("status")
        if not isinstance(status, str) or status.lower() not in _CRITERION_STATUSES:
            raise EvidenceValidationError(
                f"invalid assessment criterion status for {criterion_id!r}: {status!r}"
            )
    if len(set(ids)) != len(ids) or set(ids) != set(_EXPECTED_CRITERION_IDS):
        raise EvidenceValidationError("assessment criteria must have the exact 14 unique known ids")

    by_id = {item["id"]: item for item in criteria}
    result = []
    for criterion_id in _EXPECTED_CRITERION_IDS:
        item = by_id[criterion_id]
        row = {
            "id": item.get("id"),
            "title": item.get("title"),
            "status": item.get("status"),
            "required": _compact_diagnostic(item.get("required")),
            "reason": item.get("reason"),
        }
        observed = item.get("observed")
        if observed is None or type(observed) in (str, int, float, bool):
            row["observed"] = observed
        elif isinstance(observed, dict):
            compact = {
                key: value
                for key, value in observed.items()
                if value is None
                or type(value) in (str, int, bool)
                or isinstance(value, float) and math.isfinite(value)
                or isinstance(value, list)
                and len(value) <= 30
                and all(x is None or type(x) in (str, int, float, bool) for x in value)
            }
            row["observed"] = compact
            if len(compact) != len(observed):
                row["nested_evidence_omitted"] = True
        else:
            row["nested_evidence_omitted"] = True
        result.append(row)
    return result


def _extract_runtime_freshness(receipt: dict[str, Any]) -> Any:
    freshness = receipt.get("freshness")
    if freshness is not None:
        return _compact_diagnostic(freshness)
    runtime = receipt.get("runtime")
    if isinstance(runtime, dict):
        return _compact_diagnostic(runtime.get("freshness", runtime.get("current")))
    return None


def _build_parent(export: dict[str, Any], source_id: str) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    receipt = export["receipt"]
    result = receipt.get("result")
    if not isinstance(result, dict):
        raise EvidenceValidationError(f"receipt.result must be an object for {export['job_id']}")
    documents = [result]
    analogs = _record_list(result.get("analogs", []), "receipt.result.analogs")
    candidates = _record_list(
        result.get("protac_candidates", result.get("candidates", [])),
        "receipt.result.protac_candidates/candidates",
    )
    supported_sites = _sites(result.get("sites", []), "receipt.result.sites")
    parameters = receipt.get("parameters")
    if not isinstance(parameters, dict) or type(parameters.get("exploratory")) is not bool:
        raise EvidenceValidationError(
            f"receipt.parameters.exploratory must be explicitly boolean for {export['job_id']}"
        )
    mode = "exploratory" if parameters["exploratory"] else "strict"

    analog_rows: list[dict[str, Any]] = []
    analog_by_id: dict[str, dict[str, Any]] = {}
    graph_sets: dict[int, set[str]] = {}
    actual_generated_families: set[str] = set()
    supported_site_families: set[str] = set()
    attempted = completed = failed_docking = failed_posefilter = selected_count = skipped = 0
    unknown_statuses: dict[str, int] = {}
    analog_seen: dict[str, bytes] = {}

    for record in analogs:
        analog_id = record.get("id") or record.get("analog_id")
        _require_identifier(analog_id, "analog_id")
        if analog_id in analog_seen:
            raise EvidenceValidationError(
                f"duplicate analog record for {export['job_id']}/{analog_id}"
            )
        analog_seen[analog_id] = _canonical_bytes(record)
        exact_smiles = _record_smiles(record)
        site_maps = _site_maps(record)
        cheap = record.get("cheap_filter")
        cheap_valid = isinstance(cheap, dict) and cheap.get("valid") is True
        explicitly_qualified = record.get("qualified_for_counts") is True
        supported = bool(site_maps) and bool(set(site_maps) & supported_sites)
        graph = _constitutional_smiles(exact_smiles) if exact_smiles else None
        actual_generated = bool(cheap_valid and graph)
        supported_site_qualified = bool(
            actual_generated and explicitly_qualified and supported
        )
        family = _broad_family(record.get("transformation_class"))
        if actual_generated and family:
            actual_generated_families.add(family)
        if supported_site_qualified:
            if family:
                supported_site_families.add(family)
            for atom_map in set(site_maps) & supported_sites:
                graph_sets.setdefault(atom_map, set()).add(graph)
        pose_qualified, docking_metrics, docking_failure = _pose_details(record)
        docking = record.get("docking")
        status = docking.get("status") if isinstance(docking, dict) else None
        was_completed = isinstance(status, str) and status.startswith("completed")
        was_failed_docking = status == "failed_docking"
        was_attempted = was_completed or was_failed_docking
        if status in (None, "not_eligible", "not_run", "preview", "skipped"):
            skipped += 1
        elif not was_attempted:
            key = str(status)
            unknown_statuses[key] = unknown_statuses.get(key, 0) + 1
        selected = record.get("selected") is True
        attempted += int(was_attempted)
        completed += int(was_completed)
        failed_docking += int(was_failed_docking)
        failed_posefilter += int(was_completed and not pose_qualified)
        selected_count += int(selected)
        strict_eligible = bool(
            mode == "strict"
            and supported_site_qualified
            and selected
            and pose_qualified
            and record.get("assembly_eligible") is True
            and record.get("parent_redocking_supported") is True
        )
        row = {
            "analog_key": f"{export['parent_id']}/{export['job_id']}/{analog_id}",
            "analog_id": analog_id,
            "parent_id": export["parent_id"],
            "job_id": export["job_id"],
            "evidence_source_id": source_id,
            "canonical_smiles": exact_smiles,
            "transformation_class": record.get("transformation_class"),
            "broad_family": family,
            "site_maps": site_maps,
            "cheap_filter_valid": cheap_valid,
            "explicitly_qualified_for_counts": explicitly_qualified,
            "supported_site": supported,
            "actual_generated": actual_generated,
            "supported_site_qualified": supported_site_qualified,
            "selected": selected,
            "pose_qualified": pose_qualified,
            "strict_eligible": strict_eligible,
            "docking_metrics": docking_metrics,
            "docking_failure": docking_failure,
        }
        analog_rows.append(row)
        analog_by_id[analog_id] = row

    candidate_rows: list[dict[str, Any]] = []
    candidate_seen: dict[str, bytes] = {}
    branch_counts: dict[str, int] = {}
    for record in candidates:
        candidate_id = record.get("candidate_id") or record.get("id")
        _require_identifier(candidate_id, "candidate_id")
        if candidate_id in candidate_seen:
            raise EvidenceValidationError(
                f"duplicate candidate record for {export['job_id']}/{candidate_id}"
            )
        candidate_seen[candidate_id] = _canonical_bytes(record)
        exact_smiles = _record_smiles(record)
        exact_graph = _constitutional_smiles(exact_smiles) if exact_smiles else None
        analog_id = record.get("warhead_analog_id")
        analog = analog_by_id.get(analog_id) if isinstance(analog_id, str) else None
        e3_type = record.get("e3_type")
        status_values = (record.get("status"), record.get("assembly_status"))
        assembly_failed = (
            record.get("preview") is True
            or record.get("error") is not None
            or record.get("failure") is not None
            or any(
                isinstance(value, str)
                and any(token in value.lower() for token in ("failed", "error", "preview"))
                for value in status_values
            )
        )
        assembled = bool(
            analog
            and exact_graph
            and record.get("assembly_mode") == "pose_supported_hypothesis"
            and not assembly_failed
        )
        if assembled and isinstance(e3_type, str):
            branch_counts[e3_type] = branch_counts.get(e3_type, 0) + 1
        selected = bool(analog and analog["selected"])
        pose_qualified = bool(analog and analog["pose_qualified"])
        strict_eligible = bool(mode == "strict" and analog and analog["strict_eligible"] and assembled)
        candidate_rows.append(
            {
                "candidate_key": f"{export['parent_id']}/{export['job_id']}/{candidate_id}",
                "candidate_id": candidate_id,
                "parent_id": export["parent_id"],
                "job_id": export["job_id"],
                "evidence_source_id": source_id,
                "canonical_smiles": exact_smiles,
                "e3_type": e3_type,
                "warhead_analog_id": analog_id,
                "linker_id": record.get("linker_id"),
                "orientation": record.get("orientation"),
                "attachment_metadata": {
                    key: copy.deepcopy(record["attachment_metadata"][key])
                    for key in (
                        "selected_attachment_map",
                        "warhead_linker_bond",
                        "recruiter_linker_bond",
                    )
                    if isinstance(record.get("attachment_metadata"), dict)
                    and key in record["attachment_metadata"]
                },
                "risk_flags": copy.deepcopy(record.get("risk_flags")),
                "selected": selected,
                "pose_qualified": pose_qualified,
                "strict_eligible": strict_eligible,
                "assembled": assembled,
                "assembly_failure": record.get("failure") or record.get("assembly_failure"),
                "selected_analog_docking_metrics": analog["docking_metrics"] if analog and selected else {},
                "selected_analog_docking_failure": analog["docking_failure"] if analog and selected else None,
            }
        )

    supplied_summaries = []
    for document in documents:
        if isinstance(document, dict):
            result = document.get("result")
            if isinstance(result, dict) and isinstance(result.get("summary"), dict):
                supplied_summaries.append(_compact_diagnostic(result["summary"]))
            elif isinstance(document.get("summary"), dict):
                supplied_summaries.append(_compact_diagnostic(document["summary"]))

    parent = {
        "parent_id": export["parent_id"],
        "job_id": export["job_id"],
        "evidence_source_id": source_id,
        "mode": mode,
        "runtime_freshness_historical_as_supplied": _extract_runtime_freshness(receipt),
        "runtime_freshness_currently_verified": False,
        "approval": False,
        "counts_recomputed_from_records": {
            "analog_records": len(analog_rows),
            "distinct_graphs_by_supported_site": {
                str(key): len(value) for key, value in sorted(graph_sets.items())
            },
            "actual_generated_family_count": len(actual_generated_families),
            "actual_generated_families": sorted(actual_generated_families),
            "supported_site_family_count": len(supported_site_families),
            "supported_site_families": sorted(supported_site_families),
            "docking_attempts": attempted,
            "docking_completed": completed,
            "failed_docking": failed_docking,
            "computed_but_pose_filter_rejected": failed_posefilter,
            "docking_skipped": skipped,
            "docking_unknown_statuses": dict(sorted(unknown_statuses.items())),
            "selected_analogs": selected_count,
            "assemblies": sum(row["assembled"] for row in candidate_rows),
            "e3_branch_assemblies": dict(sorted(branch_counts.items())),
        },
        "supplied_summary_diagnostics_untrusted_for_counts": supplied_summaries,
    }
    analog_rows.sort(key=lambda row: row["analog_key"])
    candidate_rows.sort(key=lambda row: row["candidate_key"])
    return parent, analog_rows, candidate_rows


def _strict_source(path: Path) -> dict[str, Any]:
    raw = _read_regular(path, "caller-provided strict receipt", _MAX_JSON_BYTES)
    value = _loads_json(raw, "caller-provided strict receipt")
    if not isinstance(value, dict):
        raise EvidenceValidationError("caller-provided strict receipt must contain an object")
    return {
        "raw": raw,
        "value": value,
        "sha256": _sha256(raw),
        "leaf": path.name,
    }


def build_campaign_evidence(
    assessment_path: Path,
    parent_export_dirs: list[Path],
    strict_receipt_path: Path | None = None,
) -> dict[str, Any]:
    """Build a deterministic, compact evidence bundle from verified exports.

    Parent exports are trusted only for byte integrity and closure completeness.
    Scientific acceptance, current runtime freshness, and human approval are never
    inferred. The optional strict receipt has no manifest and is therefore recorded
    only as caller-provided context; it never contributes scientific counts.
    """
    assessment_path = Path(assessment_path)
    assessment_raw = _read_regular(assessment_path, "assessment", _MAX_JSON_BYTES)
    assessment = _loads_json(assessment_raw, "assessment")
    if not isinstance(assessment, dict):
        raise EvidenceValidationError("assessment must contain a JSON object")

    exports = [_verify_export(Path(path)) for path in parent_export_dirs]
    exports.sort(key=lambda item: (item["parent_id"], item["job_id"], item["receipt_sha256"]))
    seen_jobs: set[str] = set()
    seen_parents: set[str] = set()
    for export in exports:
        if export["job_id"] in seen_jobs:
            raise EvidenceValidationError(f"duplicate job_id across exports: {export['job_id']}")
        if export["parent_id"] in seen_parents:
            raise EvidenceValidationError(f"duplicate parent_id across exports: {export['parent_id']}")
        seen_jobs.add(export["job_id"])
        seen_parents.add(export["parent_id"])

    strict = _strict_source(Path(strict_receipt_path)) if strict_receipt_path is not None else None
    evidence_sources = [
        {
            "source_id": "E001",
            "kind": "assessment_snapshot",
            "path_leaf": assessment_path.name,
            "sha256": _sha256(assessment_raw),
            "closure_verified": False,
            "scientific_trust": False,
            "approval": False,
            "provenance_note": "Imported snapshot, not live currentness or authenticated scientific authority.",
        }
    ]
    export_source_ids: dict[tuple[str, str], str] = {}
    for export in exports:
        source_id = f"E{len(evidence_sources) + 1:03d}"
        export_source_ids[(export["parent_id"], export["job_id"])] = source_id
        evidence_sources.append(
            {
                "source_id": source_id,
                "kind": "verified_parent_export",
                "path_leaf": export["root"].name,
                "manifest_sha256": export["manifest_sha256"],
                "raw_receipt_sha256": export["receipt_sha256"],
                "job_id": export["job_id"],
                "parent_id": export["parent_id"],
                "closure_verified": True,
                "artifact_count": len(export["artifacts"]),
                "scientific_trust": False,
                "approval": False,
                "provenance_note": "Byte integrity and artifact closure verified; scientific claims and current runtime freshness are not verified.",
            }
        )
    if strict is not None:
        evidence_sources.append(
            {
                "source_id": f"E{len(evidence_sources) + 1:03d}",
                "kind": "caller_provided_strict_receipt",
                "path_leaf": strict["leaf"],
                "sha256": strict["sha256"],
                "mode": "strict",
                "closure_verified": False,
                "closure_status": "unverified_no_manifest",
                "scientific_trust": False,
                "counts_toward_scientific_evidence": False,
                "approval": False,
                "provenance_note": "Caller-provided receipt without a manifest; retained only as unverified context.",
            }
        )

    parents: list[dict[str, Any]] = []
    analog_rows: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, Any]] = []
    for export in exports:
        parent, analogs, candidates = _build_parent(
            export, export_source_ids[(export["parent_id"], export["job_id"])]
        )
        parents.append(parent)
        analog_rows.extend(analogs)
        candidate_rows.extend(candidates)

    criteria = _criteria_snapshot(assessment)
    generated_family_union = sorted(
        {
            family
            for parent in parents
            for family in parent["counts_recomputed_from_records"]["actual_generated_families"]
        }
    )
    supported_family_union = sorted(
        {
            family
            for parent in parents
            for family in parent["counts_recomputed_from_records"]["supported_site_families"]
        }
    )
    campaign_counts = {
        "actual_generated_family_count": len(generated_family_union),
        "actual_generated_families": generated_family_union,
        "supported_site_family_count": len(supported_family_union),
        "supported_site_families": supported_family_union,
        "docking_attempts": sum(
            parent["counts_recomputed_from_records"]["docking_attempts"] for parent in parents
        ),
        "docking_completed": sum(
            parent["counts_recomputed_from_records"]["docking_completed"] for parent in parents
        ),
        "failed_docking": sum(
            parent["counts_recomputed_from_records"]["failed_docking"] for parent in parents
        ),
        "computed_but_pose_filter_rejected": sum(
            parent["counts_recomputed_from_records"]["computed_but_pose_filter_rejected"]
            for parent in parents
        ),
        "selected_analogs": sum(row["selected"] for row in analog_rows),
        "assemblies": sum(row["assembled"] for row in candidate_rows),
    }
    bundle: dict[str, Any] = {
        "format": FORMAT,
        "generated_at": _datetime.datetime.now(_datetime.timezone.utc).isoformat(),
        "approval": False,
        "scientific_accepted": False,
        "criteria_relaxed": False,
        "assessment_snapshot": {
            "assessment_id": assessment.get("id"),
            "job_id": assessment.get("job_id"),
            "revision": assessment.get("revision"),
            "original_scientific_accepted": assessment.get("scientific_accepted"),
            "source_runtime_current_as_supplied": assessment.get("source_runtime_current"),
            "currentness_verified_now": False,
            "criteria14": criteria,
            "criteria_total_in_snapshot": len(assessment.get("criteria", [])),
            "status_policy": "Original criterion statuses are preserved; recomputed campaign aggregates do not create policy passes.",
        },
        "evidence_sources": evidence_sources,
        "parents": parents,
        "campaign_counts_recomputed_from_records": campaign_counts,
        "analogs": analog_rows,
        "candidates": candidate_rows,
        "retrieval_indexes": {
            "candidate_keys": [row["candidate_key"] for row in candidate_rows],
            "analog_keys": [row["analog_key"] for row in analog_rows],
            "candidates_by_job": {
                job_id: [row["candidate_key"] for row in candidate_rows if row["job_id"] == job_id]
                for job_id in sorted(seen_jobs)
            },
            "candidates_by_parent": {
                parent_id: [row["candidate_key"] for row in candidate_rows if row["parent_id"] == parent_id]
                for parent_id in sorted(seen_parents)
            },
            "evidence_by_job": {
                export["job_id"]: export_source_ids[(export["parent_id"], export["job_id"])]
                for export in exports
            },
        },
        "limitations": [
            "Integrity verification is not scientific validation or human approval.",
            "Runtime freshness is historical as supplied and is not claimed current.",
            "Family and graph labels are deterministic record-derived diagnostics, not human validation.",
            "CRBN and VHL branches remain separate; no cross-E3 composite score is produced.",
            "Caller-provided strict receipts without manifests never contribute scientific trust or counts.",
        ],
    }
    digest_input = copy.deepcopy(bundle)
    digest_input.pop("generated_at", None)
    digest_input.pop("digest", None)
    bundle["digest"] = _sha256(_canonical_bytes(digest_input))
    return bundle
