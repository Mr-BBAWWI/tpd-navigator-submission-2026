"""Build and verify a read-only, source-bound nine-parent funnel dossier.

The dossier is a computed diagnostic for later expert review.  It never performs
human selection, scientific approval, policy acceptance, or M2 registration.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shlex
import shutil
import stat
import zipfile
from collections import Counter
from pathlib import Path
from types import MappingProxyType
from typing import Any, Iterable, Mapping
from xml.etree import ElementTree

from rdkit import Chem

from packages.science.parent_sar_probe import (
    EvidenceError,
    _summarise_result,
    sha256_file,
    verify_parent_export,
)
from packages.science.reference_parents import load_reference_parent

FORMAT = "parent-funnel-dossier/1.0"
MANIFEST_FORMAT = "parent-funnel-bundle-manifest/1.0"
EXPECTED_PARENT_COUNT = 9
TRUSTED_EXPERT_SOURCE_SHA256 = "4132ea99e82d2711ac49e3921c9d095e2cef631ef1a0487a46a7cde45380570f"
_SHA256 = re.compile(r"[a-f0-9]{64}\Z")
_DOI = re.compile(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+\Z", re.IGNORECASE)
_FORBIDDEN_TRUE_KEYS = {
    "scientific_approval",
    "strict_approval",
    "formal_accept",
    "accepted",
    "human_selected",
    "m2_registered",
}
_REASON_CLASSES = (
    "NO_APPLICABLE_SITE",
    "PROTECTED_TOUCHED",
    "UNKNOWN_REQUIRES_EXPLORATORY",
    "quota",
    "chemical_filter",
)


class FunnelError(EvidenceError):
    """Raised when a funnel source or generated bundle is unsafe or inconsistent."""


def _reject_constant(value: str) -> None:
    raise FunnelError(f"NONFINITE_JSON_NUMBER:{value}")


def _no_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise FunnelError(f"DUPLICATE_JSON_KEY:{key}")
        result[key] = value
    return result


def _check_finite(value: Any) -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise FunnelError("NONFINITE_JSON_NUMBER")
    if isinstance(value, dict):
        for child in value.values():
            _check_finite(child)
    elif isinstance(value, list):
        for child in value:
            _check_finite(child)


def _strict_json_bytes(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=_no_duplicate_pairs,
            parse_constant=_reject_constant,
        )
    except FunnelError:
        raise
    except Exception as exc:
        raise FunnelError(f"INVALID_JSON:{label}") from exc
    if not isinstance(value, dict):
        raise FunnelError(f"JSON_OBJECT_REQUIRED:{label}")
    _check_finite(value)
    return value


def _strict_json(path: Path) -> dict[str, Any]:
    try:
        return _strict_json_bytes(path.read_bytes(), str(path))
    except OSError as exc:
        raise FunnelError(f"JSON_READ_FAILED:{path}") from exc


def _json_bytes(value: Any) -> bytes:
    _check_finite(value)
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")


def _write_json(path: Path, value: Any) -> None:
    path.write_bytes(_json_bytes(value))


def _canonical_hash(value: Any) -> str:
    return hashlib.sha256(_json_bytes(value)).hexdigest()


def _walk(value: Any, path: tuple[str, ...] = ()) -> Iterable[tuple[tuple[str, ...], Any]]:
    yield path, value
    if isinstance(value, dict):
        for key, child in value.items():
            yield from _walk(child, path + (str(key),))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from _walk(child, path + (str(index),))


def _family(row: Mapping[str, Any]) -> str | None:
    value = row.get("broad_family") or row.get("transformation_class") or row.get("family")
    if value in {"ring_expansion", "ring_contraction"}:
        return "ring_modification"
    return value if isinstance(value, str) and value else None


def _candidate_id(row: Mapping[str, Any]) -> str | None:
    for key in ("candidate_id", "warhead_id", "analog_id", "compound_id", "id", "title"):
        value = row.get(key)
        if isinstance(value, str) and value:
            return value
    return None


def _safe_component(value: str) -> str:
    if not value or value in {".", ".."} or "/" in value or "\\" in value:
        raise FunnelError(f"UNSAFE_PATH_COMPONENT:{value}")
    if not re.fullmatch(r"[A-Za-z0-9._-]+", value):
        raise FunnelError(f"UNSAFE_PATH_COMPONENT:{value}")
    return value


def _is_link_or_reparse(path: Path) -> bool:
    try:
        info = path.lstat()
    except OSError as exc:
        raise FunnelError(f"PATH_STAT_FAILED:{path}") from exc
    if stat.S_ISLNK(info.st_mode):
        return True
    attributes = getattr(info, "st_file_attributes", 0)
    reparse_flag = getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0)
    if not reparse_flag or not attributes & reparse_flag:
        return False
    # Cloud placeholders (including ordinary OneDrive files) are not links.
    # Reject only name-surrogate reparses and mount-point junctions.
    tag = getattr(info, "st_reparse_tag", 0)
    return bool(tag & 0x20000000) or tag == 0xA0000003


def _checked_chain(path: Path, root: Path, label: str) -> tuple[Path, Path]:
    path_absolute = path.absolute()
    root_absolute = root.absolute()
    try:
        relative = path_absolute.relative_to(root_absolute)
    except ValueError as exc:
        raise FunnelError(f"{label}_OUTSIDE_ROOT:{path}") from exc
    current = root_absolute
    for part in (Path(), *relative.parts):
        if part != Path():
            current = current / part
        if _is_link_or_reparse(current):
            raise FunnelError(f"{label}_LINK_FORBIDDEN:{current}")
    try:
        resolved_root = root_absolute.resolve(strict=True)
        resolved = path_absolute.resolve(strict=True)
        resolved.relative_to(resolved_root)
    except (OSError, ValueError) as exc:
        raise FunnelError(f"{label}_OUTSIDE_ROOT:{path}") from exc
    return resolved, resolved_root


def _regular_file(path: Path, root: Path, label: str) -> Path:
    resolved, _ = _checked_chain(path, root, label)
    try:
        mode = resolved.stat().st_mode
    except OSError as exc:
        raise FunnelError(f"{label}_STAT_FAILED:{path}") from exc
    if not stat.S_ISREG(mode):
        raise FunnelError(f"{label}_NOT_REGULAR_FILE:{path}")
    return resolved


def _safe_directory(path: Path, label: str) -> Path:
    resolved, _ = _checked_chain(path, path, label)
    if not resolved.is_dir():
        raise FunnelError(f"{label}_NOT_DIRECTORY:{path}")
    return resolved


def _bundle_parts(relative: str) -> tuple[str, ...]:
    if not isinstance(relative, str) or not relative or "\\" in relative or ":" in relative:
        raise FunnelError("INVALID_BUNDLE_RELATIVE_PATH")
    if relative.startswith("/") or Path(relative).is_absolute():
        raise FunnelError("INVALID_BUNDLE_RELATIVE_PATH")
    parts = tuple(relative.split("/"))
    if any(part in {"", ".", ".."} for part in parts):
        raise FunnelError("BUNDLE_PATH_TRAVERSAL")
    return parts


def _copy_verified_file(source: Path, target: Path, root: Path, label: str, expected_sha256: str, expected_bytes: int | None = None) -> None:
    source = _regular_file(source, root, label)
    before_hash = sha256_file(source)
    before_size = source.stat().st_size
    if before_hash != expected_sha256 or (expected_bytes is not None and before_size != expected_bytes):
        raise FunnelError(f"{label}_PRECOPY_MISMATCH:{source}")
    shutil.copyfile(source, target)
    if sha256_file(source) != before_hash or source.stat().st_size != before_size:
        raise FunnelError(f"{label}_CHANGED_DURING_COPY:{source}")
    if sha256_file(target) != before_hash or target.stat().st_size != before_size:
        raise FunnelError(f"{label}_COPY_MISMATCH:{target}")


def _docx_paragraphs(path: Path) -> list[str]:
    if path.suffix.lower() != ".docx":
        raise FunnelError("EXPERT_SOURCE_MUST_BE_DOCX")
    if _is_link_or_reparse(path) or not path.is_file():
        raise FunnelError("EXPERT_SOURCE_NOT_REGULAR_FILE")
    try:
        with zipfile.ZipFile(path, "r") as archive:
            names = archive.namelist()
            if len(names) != len(set(names)) or "word/document.xml" not in names:
                raise FunnelError("INVALID_EXPERT_DOCX_CLOSURE")
            xml = archive.read("word/document.xml")
    except FunnelError:
        raise
    except Exception as exc:
        raise FunnelError("INVALID_EXPERT_DOCX") from exc
    try:
        root = ElementTree.fromstring(xml)
    except ElementTree.ParseError as exc:
        raise FunnelError("INVALID_EXPERT_DOCUMENT_XML") from exc
    namespace = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    paragraphs: list[str] = []
    for paragraph in root.iter(namespace + "p"):
        text = "".join(node.text or "" for node in paragraph.iter(namespace + "t")).strip()
        if text:
            paragraphs.append(text)
    if not paragraphs:
        raise FunnelError("EXPERT_DOCX_HAS_NO_READABLE_PARAGRAPHS")
    return paragraphs


def _expert_scope(paragraphs: list[str], source_sha256: str) -> dict[str, Any]:
    if source_sha256 != TRUSTED_EXPERT_SOURCE_SHA256:
        raise FunnelError("UNTRUSTED_EXPERT_SOURCE_SHA256")
    joined = "\n".join(paragraphs)
    ids = sorted(set(re.findall(r"SMARCA2-[A-Za-z0-9-]+", joined)))
    return {
        "document_paragraphs": paragraphs,
        "summary": "전문가 회신은 ready reference parent 9개를 exploratory pool로 허용하고 최소 5개 비교를 요구하지만, strict parent 승인이나 formal acceptance를 부여하지 않는다.",
        "exploratory_pool_permitted": "9개" in joined or len(ids) >= EXPECTED_PARENT_COUNT,
        "minimum_comparison_parent_count": 5 if "최소 5개" in joined else None,
        "strict_approval_granted": False,
        "formal_accept_submitted": False,
        "parsed_parent_ids": ids,
    }


def _citation_records(*values: Any) -> list[dict[str, Any]]:
    records: dict[bytes, dict[str, Any]] = {}
    relevant = {
        "doi", "filename", "relative_path", "path", "source_row_locators",
        "record_number_1_based_including_header", "assay_source_locator",
        "source_column_label", "compound_id_column", "source_record_id",
    }
    for value in values:
        for path, node in _walk(value):
            if not isinstance(node, dict) or not relevant.intersection(node):
                continue
            row = {key: node[key] for key in sorted(relevant.intersection(node))}
            doi = row.get("doi")
            if doi is not None and (not isinstance(doi, str) or _DOI.fullmatch(doi.strip()) is None):
                row.pop("doi", None)
            if not row:
                continue
            row["evidence_path"] = "/".join(path)
            key = _json_bytes(row)
            records[key] = row
    return [records[key] for key in sorted(records)]


def _reason_class(row: Mapping[str, Any]) -> str:
    text = " ".join(str(value) for _, value in _walk(row) if isinstance(value, (str, int))).upper()
    if "NO_APPLICABLE_SITE" in text:
        return "NO_APPLICABLE_SITE"
    if "PROTECTED_TOUCHED" in text:
        return "PROTECTED_TOUCHED"
    if "UNKNOWN_REQUIRES_EXPLORATORY" in text:
        return "UNKNOWN_REQUIRES_EXPLORATORY"
    if "QUOTA" in text or "LIMIT" in text:
        return "quota"
    if "CHEAP_FILTER" in text or "CHEMICAL_FILTER" in text or "INVALID" in text:
        return "chemical_filter"
    return "other"


def _rule_id(row: Mapping[str, Any]) -> str:
    for key in ("rule_id", "transformation_rule_id", "rule", "generator_id"):
        value = row.get(key)
        if isinstance(value, str) and value:
            return value
    return "unidentified_rule"


def _named_values(row: Mapping[str, Any], names: set[str]) -> list[dict[str, Any]]:
    values: list[dict[str, Any]] = []
    for path, value in _walk(row):
        if path and path[-1] in names:
            values.append({"path": "/".join(path), "value": value})
    return values


def _rule_matrix(result: Mapping[str, Any]) -> list[dict[str, Any]]:
    analogs = [row for row in result.get("analogs", []) if isinstance(row, dict)] if isinstance(result.get("analogs"), list) else []
    rejections = [row for row in result.get("rejections", []) if isinstance(row, dict)] if isinstance(result.get("rejections"), list) else []
    catalog = result.get("rule_catalog")
    if isinstance(catalog, dict):
        catalog_rows = []
        for key, value in catalog.items():
            if isinstance(value, dict):
                row = dict(value)
                if _rule_id(row) == "unidentified_rule":
                    row["rule_id"] = str(key)
                catalog_rows.append(row)
    else:
        catalog_rows = [row for row in catalog if isinstance(row, dict)] if isinstance(catalog, list) else []
    rule_ids = {_rule_id(row) for row in analogs + rejections + catalog_rows}
    matrix: list[dict[str, Any]] = []
    for rule_id in sorted(rule_ids):
        emitted = [row for row in analogs if _rule_id(row) == rule_id]
        rejected = [row for row in rejections if _rule_id(row) == rule_id]
        catalog_matches = [row for row in catalog_rows if _rule_id(row) == rule_id]
        reasons = Counter(_reason_class(row) for row in rejected)
        families = sorted({family for family in (_family(row) for row in emitted + rejected) if family})
        matrix.append({
            "rule_id": rule_id,
            "catalog_scope_records": catalog_matches,
            "actual_families": families,
            "emitted_analog_count": len(emitted),
            "emitted_candidate_ids": sorted(filter(None, (_candidate_id(row) for row in emitted))),
            "rejection_count": len(rejected),
            "rejection_reason_counts": {key: reasons.get(key, 0) for key in (*_REASON_CLASSES, "other")},
            "rejection_records": rejected,
            "rejection_touched_maps": [_named_values(row, {"touched_maps", "touched_atom_maps", "affected_atom_maps"}) for row in rejected],
            "rejection_protected_maps": [_named_values(row, {"protected_maps", "protected_atom_maps", "protected_touched_maps"}) for row in rejected],
            "rejection_unknown_maps": [_named_values(row, {"unknown_maps", "unknown_atom_maps", "unresolved_atom_maps"}) for row in rejected],
            "cheap_filter_reasons": [_named_values(row, {"cheap_filter_reason", "cheap_filter_reasons", "chemical_filter_reason", "filter_reason"}) for row in rejected],
        })
    return matrix


def _redock(result: Mapping[str, Any]) -> dict[str, Any]:
    summary = result.get("summary") if isinstance(result.get("summary"), dict) else {}
    stages = summary.get("stage_counts") if isinstance(summary.get("stage_counts"), dict) else {}
    stage = stages.get("parent_redock_and_pose_filter") if isinstance(stages.get("parent_redock_and_pose_filter"), dict) else {}
    value = stage.get("parent_pass")
    return {
        "required": (result.get("parent_scope") or {}).get("parent_redock_required_for_qualification") is True,
        "actual_pass": value if type(value) is bool else None,
        "actual_status": "pass" if value is True else "fail" if value is False else "not_recorded",
        "source_stage_record": stage,
    }


def _evidence_tiers(warhead: Mapping[str, Any], linkage: Mapping[str, Any]) -> dict[str, bool]:
    locators = linkage.get("source_sar") if isinstance(linkage.get("source_sar"), list) else []
    supported_pair = any(
        isinstance(row, dict)
        and bool(row.get("source_row_locators"))
        and bool(row.get("selected_parent_affected_atom_maps"))
        and row.get("source_site_consensus") is True
        for row in locators
    )
    exact = linkage.get("status") == "exact_measured_parent_join"
    technical = warhead.get("evidence_tier") == "co_crystal" and bool(warhead.get("pdb")) and bool(warhead.get("ccd"))
    return {
        "co_crystal_technical_ready": technical,
        "exact_source_join": exact,
        "source_locator_atommap_supported_pair": supported_pair,
        "missing_evidence": not (technical and exact and supported_pair),
    }


def _map_free_canonical_isomeric_smiles(mol: Chem.Mol) -> str:
    copied = Chem.Mol(mol)
    for atom in copied.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.AssignStereochemistry(copied, cleanIt=True, force=True)
    return Chem.MolToSmiles(copied, canonical=True, isomericSmiles=True)


def _reference_check(parent_id: str, result: Mapping[str, Any], reference_root: Path | None) -> dict[str, Any]:
    if reference_root is None:
        return {"status": "not_requested", "independently_reverified_from_bundle": False, "exact_identity_consistent": None}
    mol, metadata = load_reference_parent(parent_id, reference_root)
    warhead = result.get("warhead") if isinstance(result.get("warhead"), dict) else {}
    stored = warhead.get("canonical_isomeric_smiles") or (warhead.get("identity") or {}).get("canonical_isomeric_smiles")
    metadata_smiles = metadata.get("canonical_isomeric_smiles")
    mappings = warhead.get("atom_mapping")
    metadata_mappings = metadata.get("atom_mapping")
    if not isinstance(stored, str) or not stored or not isinstance(metadata_smiles, str) or not metadata_smiles:
        raise FunnelError(f"REFERENCE_PARENT_CANONICAL_SMILES_MISSING:{parent_id}")
    receipt_mol = Chem.MolFromSmiles(stored)
    metadata_mol = Chem.MolFromSmiles(metadata_smiles)
    if receipt_mol is None or metadata_mol is None:
        raise FunnelError(f"INVALID_REFERENCE_OR_RECEIPT_WARHEAD_SMILES:{parent_id}")
    ref_smiles = _map_free_canonical_isomeric_smiles(mol)
    receipt_smiles = _map_free_canonical_isomeric_smiles(receipt_mol)
    metadata_canonical = _map_free_canonical_isomeric_smiles(metadata_mol)
    map_shape = isinstance(mappings, list) and len(mappings) == receipt_mol.GetNumAtoms()
    if map_shape:
        maps = [row.get("atom_map") for row in mappings if isinstance(row, dict)]
        indices = [row.get("rdkit_index_zero_based") for row in mappings if isinstance(row, dict)]
        map_shape = (
            len(maps) == len(mappings)
            and all(type(value) is int and value > 0 for value in maps)
            and len(set(maps)) == len(maps)
            and sorted(indices) == list(range(receipt_mol.GetNumAtoms()))
        )
    exact = receipt_smiles == ref_smiles == metadata_canonical and mappings == metadata_mappings and map_shape
    if not exact:
        raise FunnelError(f"REFERENCE_PARENT_IDENTITY_OR_MAP_MISMATCH:{parent_id}")
    return {
        "status": "verified_at_build_time_external_catalog",
        "independently_reverified_from_bundle": False,
        "exact_identity_consistent": True,
        "receipt_map_free_canonical_isomeric_smiles": receipt_smiles,
        "reference_map_free_canonical_isomeric_smiles": ref_smiles,
        "metadata_map_free_canonical_isomeric_smiles": metadata_canonical,
        "atom_mapping_exactly_matches_loaded_metadata": True,
        "atom_mapping_key_shape_consistent": True,
        "mapped_atom_count": len(mappings),
    }


def _copy_export(verified: Mapping[str, Any], source_root: Path, output_root: Path) -> dict[str, Any]:
    parent_id = _safe_component(str(verified["parent_id"]))
    destination = output_root / "exports" / parent_id
    destination.mkdir(parents=True, exist_ok=False)
    source_manifest = Path(str(verified["manifest_path"]))
    source_receipt = Path(str(verified["receipt_path"]))
    _copy_verified_file(
        source_manifest,
        destination / "manifest.json",
        source_root,
        "SOURCE_MANIFEST",
        str(verified["manifest_sha256"]),
    )
    _copy_verified_file(
        source_receipt,
        destination / "job-receipt.json",
        source_root,
        "SOURCE_RECEIPT",
        str(verified["raw_receipt_sha256"]),
    )
    blob_dir = destination / "blobs"
    blob_dir.mkdir()
    proof_blobs: list[dict[str, Any]] = []
    for blob in verified["artifact_blob_closure"]:
        artifact_id = _safe_component(str(blob["artifact_id"]))
        source = Path(str(blob["blob_path"]))
        target = blob_dir / artifact_id
        _copy_verified_file(
            source,
            target,
            source_manifest.parent,
            "SOURCE_BLOB",
            str(blob["sha256"]),
            blob["bytes"] if type(blob.get("bytes")) is int else None,
        )
        proof_blobs.append({
            "artifact_id": artifact_id,
            "relative_path": target.relative_to(output_root).as_posix(),
            "sha256": blob["sha256"],
            "bytes": blob["bytes"],
            "version": blob["version"],
        })
    return {
        "export_directory": destination.relative_to(output_root).as_posix(),
        "manifest_path": (destination / "manifest.json").relative_to(output_root).as_posix(),
        "receipt_path": (destination / "job-receipt.json").relative_to(output_root).as_posix(),
        "artifact_blobs": proof_blobs,
    }


def _parent_record(verified: Mapping[str, Any], copy_proof: Mapping[str, Any], reference_root: Path | None) -> dict[str, Any]:
    result = verified["result"]
    summary = _summarise_result(result)
    warhead = result.get("warhead") if isinstance(result.get("warhead"), dict) else {}
    scope = result.get("parent_scope") if isinstance(result.get("parent_scope"), dict) else {}
    linkage = result.get("selected_parent_measured_evidence") if isinstance(result.get("selected_parent_measured_evidence"), dict) else {}
    input_binding = result.get("input_binding")
    records = result.get("analogs") if isinstance(result.get("analogs"), list) else []
    selected_count = sum(1 for row in records if isinstance(row, dict) and row.get("selected") is True)
    return {
        "parent_id": verified["parent_id"],
        "source_proof": {
            "manifest_sha256": verified["manifest_sha256"],
            "raw_receipt_sha256": verified["raw_receipt_sha256"],
            "job_id": verified["job_id"],
            "result_canonical_json_sha256": _canonical_hash(result),
            "input_binding": input_binding,
            "bundle_paths": copy_proof,
        },
        "parent_scope": scope,
        "warhead": {
            "id": warhead.get("id"),
            "pdb": warhead.get("pdb"),
            "ccd": warhead.get("ccd"),
            "canonical_isomeric_smiles": warhead.get("canonical_isomeric_smiles"),
            "atom_mapping": warhead.get("atom_mapping"),
            "source_input_hashes": warhead.get("source_input_hashes"),
            "measured_evidence": warhead.get("measured_evidence"),
            "evidence_tier": warhead.get("evidence_tier"),
            "limitations": warhead.get("limitations"),
            "suitability": warhead.get("suitability"),
        },
        "selected_parent_measured_evidence": linkage,
        "citations": _citation_records(warhead, linkage),
        "evidence_tiers": _evidence_tiers(warhead, linkage),
        "parent_redocking": _redock(result),
        "actual_counts": {
            "cheap_passed": summary["cheap_passed_count"],
            "pose_qualified": len(summary["pose_qualified_candidate_ids"]),
            "selected": selected_count,
            "strict_qualified": len(summary["strict_qualified_candidate_ids"]),
        },
        "actual_families": {
            "cheap_passed": summary["cheap_passed_families"],
            "pose_qualified": summary["pose_qualified_families"],
            "strict_qualified": summary["strict_qualified_families"],
            "ring_expansion_and_contraction_merged_as": "ring_modification",
        },
        "stored_summary_count_mismatches": summary["stored_summary_count_mismatches"],
        "family_rule_reason_matrix": _rule_matrix(result),
        "reference_parent_verification": _reference_check(str(verified["parent_id"]), result, reference_root),
        "human_selection_pending": True,
        "atom_specific_missing_primary_evidence": not _evidence_tiers(warhead, linkage)["source_locator_atommap_supported_pair"],
        "actual_redock_failed": _redock(result)["actual_pass"] is False,
        "strict_trusted_selected_count": 0,
        "scientific_approval": False,
    }


def _compact_projection(parents: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "format": FORMAT,
        "status": "computed_diagnostic_for_expert_review",
        "parent_count": EXPECTED_PARENT_COUNT,
        "parents": [{
            "parent_id": row["parent_id"],
            "evidence_tiers": row["evidence_tiers"],
            "parent_redocking": row["parent_redocking"]["actual_status"],
            "actual_counts": row["actual_counts"],
            "actual_families": row["actual_families"],
        } for row in parents],
        "human_selection_count": 0,
        "strict_trusted_selected_count": 0,
        "scientific_approval": False,
    }


def _manifest_files(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    aliases: set[str] = set()
    for path in sorted(root.rglob("*")):
        if path == root / "manifest.json":
            continue
        if _is_link_or_reparse(path):
            raise FunnelError(f"OUTPUT_LINK_FORBIDDEN:{path}")
        if path.is_dir():
            continue
        relative = path.relative_to(root).as_posix()
        alias = relative.casefold()
        if alias in aliases:
            raise FunnelError(f"OUTPUT_PATH_ALIAS:{relative}")
        aliases.add(alias)
        rows.append({"path": relative, "sha256": sha256_file(path), "bytes": path.stat().st_size})
    return rows


def build_parent_funnel(export_root: str | Path, expert_source: str | Path, output: str | Path, reference_root: str | Path | None = "cases/reference_parents") -> dict[str, Any]:
    """Build a fresh nine-parent diagnostic bundle and return its manifest summary."""
    source_root = _safe_directory(Path(export_root).absolute(), "EXPORT_ROOT")
    expert_argument = Path(expert_source).absolute()
    expert_path = _regular_file(expert_argument, expert_argument.parent, "EXPERT_SOURCE")
    output_root = Path(output).absolute()
    if output_root.exists():
        raise FunnelError("OUTPUT_MUST_NOT_EXIST")
    output_root.mkdir(parents=True, exist_ok=False)
    manifests = sorted(source_root.glob("*/manifest.json"))
    if len(manifests) != EXPECTED_PARENT_COUNT:
        raise FunnelError(f"EXPECTED_9_PARENT_EXPORTS_FOUND_{len(manifests)}")
    ref = _safe_directory(Path(reference_root).absolute(), "REFERENCE_ROOT") if reference_root is not None else None
    expert_sha256 = sha256_file(expert_path)
    paragraphs = _docx_paragraphs(expert_path)
    expert = {
        "relative_path": "sources/expert-reply.docx",
        "sha256": expert_sha256,
        **_expert_scope(paragraphs, expert_sha256),
    }
    expert_destination = output_root / "sources" / "expert-reply.docx"
    expert_destination.parent.mkdir()
    _copy_verified_file(expert_path, expert_destination, expert_path.parent, "EXPERT_SOURCE", expert_sha256)

    parents: list[dict[str, Any]] = []
    parent_ids: set[str] = set()
    for manifest_path in manifests:
        _strict_json(manifest_path)
        _strict_json(manifest_path.parent / "job-receipt.json")
        verified = verify_parent_export(manifest_path, source_root)
        if verified["parent_id"] in parent_ids:
            raise FunnelError(f"DUPLICATE_PARENT_ID:{verified['parent_id']}")
        parent_ids.add(verified["parent_id"])
        proof = _copy_export(verified, source_root, output_root)
        copied = verify_parent_export(output_root / proof["manifest_path"], output_root / "exports")
        if copied["manifest_sha256"] != verified["manifest_sha256"] or copied["raw_receipt_sha256"] != verified["raw_receipt_sha256"]:
            raise FunnelError("COPIED_EXPORT_PROOF_MISMATCH")
        parents.append(_parent_record(verified, proof, ref))

    parents.sort(key=lambda row: row["parent_id"])
    funnel = {
        "format": FORMAT,
        "status": "computed_diagnostic_for_expert_review",
        "parent_count": EXPECTED_PARENT_COUNT,
        "parents": parents,
        "expert_source": expert,
        "aggregate_gate_substitution": False,
        "strict_trusted_selected_count": 0,
        "human_selection_count": 0,
        "human_selection_pending": True,
        "scientific_approval": False,
        "strict_approval": False,
        "formal_accept": False,
        "m2_registered": False,
    }
    compact = _compact_projection(parents)
    review = {
        "format": "parent-funnel-review-request/1.0",
        "status": "not_submitted",
        "parent_count": EXPECTED_PARENT_COUNT,
        "human_selection_count": 0,
        "expert_selected": None,
        "proposed_parents": [{
            "parent_id": row["parent_id"],
            "draft_rationale": "Verified exploratory comparison candidate; not a human selection or strict promotion.",
            "needed_source_references": [] if row["evidence_tiers"]["source_locator_atommap_supported_pair"] else ["atom-specific source + locator + structure atom-map evidence"],
            "fields_still_required": [
                *([] if row["evidence_tiers"]["source_locator_atommap_supported_pair"] else ["primary atom-specific evidence"]),
                *([] if row["parent_redocking"]["actual_pass"] is True else ["successful actual parent redocking"]),
                "expert_selected",
                "expert rationale",
            ],
            "expert_selected": None,
        } for row in parents],
        "valid_accept_request": False,
        "scientific_approval": False,
        "m2_write_requested": False,
    }
    _write_json(output_root / "funnel.json", funnel)
    _write_json(output_root / "compact-summary.json", compact)
    _write_json(output_root / "review-request.json", review)
    report = [
        "# 9-parent funnel 계산 진단 보고서",
        "",
        "## 추가 개발 필요 — 검토 상태",
        "본 번들은 검증된 export에서 계산한 expert-review 준비 자료이다. 사람의 선택, strict 승격, 과학적 승인 또는 formal acceptance가 아니다.",
        "전문가 문서는 9개 exploratory parent pool 사용과 최소 5개 비교를 허용하지만 strict approval을 부여하지 않는다.",
        "",
        "## 추가 개발 필요 — parent별 불완전 항목",
    ]
    for row in parents:
        gaps = []
        if row["atom_specific_missing_primary_evidence"]:
            gaps.append("atom-specific primary source/locator/map evidence 부족")
        if row["actual_redock_failed"]:
            gaps.append("actual parent redocking 실패")
        if not gaps:
            gaps.append("human expert selection 및 strict 근거 검토 대기")
        report.append(f"- `{row['parent_id']}`: " + "; ".join(gaps) + ". 이는 화학적 불가능성을 뜻하지 않는다.")
    report.extend([
        "",
        "## 재빌드 명령",
        "`python scripts/build_parent_funnel.py "
        f"--export-root {shlex.quote(str(source_root))} "
        f"--expert-source {shlex.quote(str(expert_path))} "
        "--output '<fresh-output-directory>' "
        f"--reference-root {shlex.quote(str(ref))}`",
        "",
        "## 통합 경계",
        "M2 root가 등록을 소유한다. 통합 시 `verify_funnel_bundle(root, expected_manifest_sha256)` 반환값을 computed diagnostic으로만 등록한다.",
    ])
    (output_root / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    manifest = {
        "format": MANIFEST_FORMAT,
        "status": "complete",
        "parent_count": EXPECTED_PARENT_COUNT,
        "files": _manifest_files(output_root),
        "claims": {
            "scientific_approval": False,
            "strict_approval": False,
            "formal_accept": False,
            "human_selection": False,
            "m2_registered": False,
        },
    }
    _write_json(output_root / "manifest.json", manifest)
    digest = sha256_file(output_root / "manifest.json")
    verify_funnel_bundle(output_root, digest)
    return {"manifest_sha256": digest, "parent_count": EXPECTED_PARENT_COUNT, "output": str(output_root)}


def _immutable(value: Any) -> Any:
    if isinstance(value, dict):
        return MappingProxyType({key: _immutable(child) for key, child in value.items()})
    if isinstance(value, list):
        return tuple(_immutable(child) for child in value)
    return value


def _assert_false_claims(value: Any) -> None:
    for path, node in _walk(value):
        if path and path[-1] in _FORBIDDEN_TRUE_KEYS and node is True:
            raise FunnelError(f"FORBIDDEN_TRUE_CLAIM:{'/'.join(path)}")


def _exact_int(value: Any, expected: int, label: str) -> None:
    if type(value) is not int or value != expected:
        raise FunnelError(f"{label}_INVALID")


def _expected_copy_proof(copied: Mapping[str, Any], parent_id: str) -> dict[str, Any]:
    base = f"exports/{parent_id}"
    blobs = []
    for blob in copied["artifact_blob_closure"]:
        artifact_id = _safe_component(str(blob["artifact_id"]))
        blobs.append({
            "artifact_id": artifact_id,
            "relative_path": f"{base}/blobs/{artifact_id}",
            "sha256": blob["sha256"],
            "bytes": blob["bytes"],
            "version": blob["version"],
        })
    return {
        "export_directory": base,
        "manifest_path": f"{base}/manifest.json",
        "receipt_path": f"{base}/job-receipt.json",
        "artifact_blobs": blobs,
    }


def _assert_exact_parent_record(stored: Mapping[str, Any], derived: Mapping[str, Any]) -> None:
    stored_copy = dict(stored)
    derived_copy = dict(derived)
    reference = stored_copy.pop("reference_parent_verification", None)
    derived_copy.pop("reference_parent_verification", None)
    if stored_copy != derived_copy:
        raise FunnelError("PARENT_SCIENTIFIC_DERIVATION_MISMATCH")
    if not isinstance(reference, dict):
        raise FunnelError("REFERENCE_BUILD_CHECK_MISSING")
    status = reference.get("status")
    if status == "verified_at_build_time_external_catalog":
        if reference.get("exact_identity_consistent") is not True or reference.get("independently_reverified_from_bundle") is not False:
            raise FunnelError("REFERENCE_BUILD_CHECK_INVALID")
    elif status == "not_requested":
        if reference != {"status": "not_requested", "independently_reverified_from_bundle": False, "exact_identity_consistent": None}:
            raise FunnelError("REFERENCE_BUILD_CHECK_INVALID")
    else:
        raise FunnelError("REFERENCE_BUILD_CHECK_INVALID")


def verify_funnel_bundle(root: str | Path, expected_manifest_sha256: str) -> Mapping[str, Any]:
    """Verify full closure and return an immutable computed-diagnostic value.

    The caller, not this helper, owns any M2 root registration operation.
    """
    if not isinstance(expected_manifest_sha256, str) or _SHA256.fullmatch(expected_manifest_sha256) is None:
        raise FunnelError("EXPECTED_MANIFEST_SHA256_INVALID")
    bundle = _safe_directory(Path(root).absolute(), "BUNDLE_ROOT")
    manifest_path = _regular_file(bundle / "manifest.json", bundle, "BUNDLE_MANIFEST")
    if sha256_file(manifest_path) != expected_manifest_sha256:
        raise FunnelError("BUNDLE_MANIFEST_HASH_MISMATCH")
    manifest = _strict_json(manifest_path)
    if manifest.get("format") != MANIFEST_FORMAT:
        raise FunnelError("INVALID_BUNDLE_MANIFEST_SCHEMA")
    _exact_int(manifest.get("parent_count"), EXPECTED_PARENT_COUNT, "BUNDLE_MANIFEST_PARENT_COUNT")
    claims = manifest.get("claims")
    if not isinstance(claims, dict) or any(claims.get(key) is not False for key in ("scientific_approval", "strict_approval", "formal_accept", "human_selection", "m2_registered")):
        raise FunnelError("BUNDLE_CLAIMS_MUST_BE_FALSE")
    declared: dict[str, dict[str, Any]] = {}
    aliases: set[str] = set()
    for row in manifest.get("files", []):
        if not isinstance(row, dict) or set(row) != {"path", "sha256", "bytes"}:
            raise FunnelError("INVALID_FILE_MANIFEST_RECORD")
        relative, digest, size = row["path"], row["sha256"], row["bytes"]
        parts = _bundle_parts(relative)
        alias = relative.casefold()
        if relative in declared or alias in aliases:
            raise FunnelError(f"DUPLICATE_OR_ALIAS_PATH:{relative}")
        aliases.add(alias)
        if not isinstance(digest, str) or _SHA256.fullmatch(digest) is None or type(size) is not int or size < 0:
            raise FunnelError("INVALID_FILE_HASH_OR_SIZE")
        path = _regular_file(bundle.joinpath(*parts), bundle, "BUNDLE_FILE")
        if path.stat().st_size != size or sha256_file(path) != digest:
            raise FunnelError(f"BUNDLE_FILE_MISMATCH:{relative}")
        declared[relative] = row
    actual: set[str] = set()
    for path in bundle.rglob("*"):
        if _is_link_or_reparse(path):
            raise FunnelError(f"BUNDLE_LINK_FORBIDDEN:{path}")
        if path.is_file() and path != manifest_path:
            _regular_file(path, bundle, "BUNDLE_CLOSURE_FILE")
            actual.add(path.relative_to(bundle).as_posix())
        elif not path.is_dir() and path != manifest_path:
            raise FunnelError(f"BUNDLE_NONREGULAR_ENTRY:{path}")
    if actual != set(declared):
        raise FunnelError(f"BUNDLE_CLOSURE_MISMATCH:unexpected={sorted(actual-set(declared))}:missing={sorted(set(declared)-actual)}")

    funnel = _strict_json(bundle / "funnel.json")
    compact = _strict_json(bundle / "compact-summary.json")
    review = _strict_json(bundle / "review-request.json")
    for value in (manifest, funnel, compact, review):
        _assert_false_claims(value)
    parents = funnel.get("parents")
    _exact_int(funnel.get("parent_count"), EXPECTED_PARENT_COUNT, "FUNNEL_PARENT_COUNT")
    if not isinstance(parents, list) or len(parents) != EXPECTED_PARENT_COUNT:
        raise FunnelError("FUNNEL_PARENT_COUNT_INVALID")
    _exact_int(funnel.get("human_selection_count"), 0, "FUNNEL_HUMAN_SELECTION_COUNT")
    _exact_int(funnel.get("strict_trusted_selected_count"), 0, "FUNNEL_STRICT_SELECTED_COUNT")
    if review.get("status") != "not_submitted" or review.get("expert_selected") is not None:
        raise FunnelError("REVIEW_REQUEST_MUST_REMAIN_UNSUBMITTED")
    _exact_int(review.get("parent_count"), EXPECTED_PARENT_COUNT, "REVIEW_PARENT_COUNT")
    _exact_int(review.get("human_selection_count"), 0, "REVIEW_SELECTION_COUNT")

    expert_relative = "sources/expert-reply.docx"
    expert_file = _regular_file(bundle / expert_relative, bundle, "EXPERT_BUNDLE_SOURCE")
    expert_sha256 = sha256_file(expert_file)
    paragraphs = _docx_paragraphs(expert_file)
    expected_expert = {
        "relative_path": expert_relative,
        "sha256": expert_sha256,
        **_expert_scope(paragraphs, expert_sha256),
    }
    if funnel.get("expert_source") != expected_expert:
        raise FunnelError("EXPERT_SOURCE_DERIVATION_MISMATCH")

    export_root = bundle / "exports"
    seen: set[str] = set()
    for parent in parents:
        if not isinstance(parent, dict):
            raise FunnelError("INVALID_FUNNEL_PARENT_RECORD")
        parent_id = _safe_component(str(parent.get("parent_id", "")))
        expected_manifest = f"exports/{parent_id}/manifest.json"
        expected_receipt = f"exports/{parent_id}/job-receipt.json"
        proof = (parent.get("source_proof") or {}).get("bundle_paths")
        if not isinstance(proof, dict):
            raise FunnelError("PARENT_EXPORT_PROOF_MISSING")
        if proof.get("manifest_path") != expected_manifest or proof.get("receipt_path") != expected_receipt:
            raise FunnelError("PARENT_EXPORT_PROOF_PATH_INVALID")
        _bundle_parts(expected_manifest)
        _bundle_parts(expected_receipt)
        copied = verify_parent_export(bundle.joinpath(*_bundle_parts(expected_manifest)), export_root)
        if copied["parent_id"] != parent_id:
            raise FunnelError("COPIED_PARENT_ID_MISMATCH")
        expected_proof = _expected_copy_proof(copied, parent_id)
        if proof != expected_proof:
            raise FunnelError("PARENT_EXPORT_PROOF_CLOSURE_MISMATCH")
        derived = _parent_record(copied, expected_proof, None)
        _assert_exact_parent_record(parent, derived)
        seen.add(parent_id)
    if len(seen) != EXPECTED_PARENT_COUNT:
        raise FunnelError("COPIED_EXPORT_PARENT_COUNT_INVALID")
    expected_compact = _compact_projection(parents)
    if compact != expected_compact:
        raise FunnelError("COMPACT_SUMMARY_DERIVATION_MISMATCH")
    return _immutable({
        "format": FORMAT,
        "status": "verified_computed_diagnostic",
        "manifest_sha256": expected_manifest_sha256,
        "parent_count": EXPECTED_PARENT_COUNT,
        "parent_ids": sorted(seen),
        "parents": parents,
        "compact_summary": compact,
        "reference_verification": "build_time_external_check_only_not_independently_reverified",
        "human_selection_count": 0,
        "strict_trusted_selected_count": 0,
        "scientific_approval": False,
        "m2_registered": False,
    })


__all__ = [
    "EXPECTED_PARENT_COUNT",
    "FORMAT",
    "FunnelError",
    "TRUSTED_EXPERT_SOURCE_SHA256",
    "build_parent_funnel",
    "verify_funnel_bundle",
]
