"""Import a verified local computed ligand-pKa evidence pack."""
from __future__ import annotations

import hashlib
import io
import json
import os
import stat
import zipfile
from pathlib import Path
from typing import Any

from packages.contracts import ContractError, encoded
from packages.platform.scientific_acceptance import (
    MAX_MICROSTATE_SOURCE_BYTES,
    ScientificAcceptanceService,
)

FORMAT = "ligand-pka-evidence/20261001.2"
MANIFEST_NAME = "manifest.json"
MAX_PACK_BYTES = 64 * 1024 * 1024
ALLOWED_PATHS = frozenset({
    "compact-summary.json",
    "core-manifest.json",
    "inputs/source-pose0.sdf",
    "LICENSE.md",
    "model-provenance.json",
    "quantitative-source.json",
    "raw-predictions.json",
    "README.md",
    "upstream-manifest.json",
    "worker-receipt.json",
})


def _require(condition: Any, code: str) -> None:
    if not condition:
        raise ContractError(code)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _canonical_sha(value: Any) -> str:
    return _sha(encoded(value))


def _is_link(path: Path) -> bool:
    junction = getattr(path, "is_junction", None)
    return path.is_symlink() or bool(junction and junction())


def _require_unlinked_chain(path: Path) -> None:
    current = path
    while True:
        _require(not _is_link(current), "COMPUTED_LIGAND_PKA_LINK")
        if current.parent == current:
            return
        current = current.parent


def _strict_json(raw: bytes):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError("duplicate JSON key")
            result[key] = value
        return result
    def constant(value):
        raise ValueError(f"invalid JSON constant: {value}")
    return json.loads(raw.decode("utf-8"), object_pairs_hook=pairs,
                      parse_constant=constant)


def _safe_regular_file(root: Path, relative_name: str) -> Path:
    relative = Path(relative_name)
    _require(
        not relative.is_absolute()
        and relative.as_posix() == relative_name
        and relative_name not in {"", ".", ".."}
        and all(part not in {"", ".", ".."} for part in relative.parts)
        and "\\" not in relative_name,
        "COMPUTED_LIGAND_PKA_PATH",
    )
    candidate = root / relative
    try:
        _require_unlinked_chain(root)
        root_resolved = root.resolve(strict=True)
        current = root
        for part in relative.parts:
            current = current / part
            info = current.lstat()
            _require(not stat.S_ISLNK(info.st_mode) and not _is_link(current),
                     "COMPUTED_LIGAND_PKA_LINK")
        resolved = candidate.resolve(strict=True)
        _require(resolved.parent == root_resolved or root_resolved in resolved.parents,
                 "COMPUTED_LIGAND_PKA_PATH")
        _require(candidate.is_file(), "COMPUTED_LIGAND_PKA_FILE")
    except ContractError:
        raise
    except OSError as error:
        raise ContractError("COMPUTED_LIGAND_PKA_FILE") from error
    return candidate


def _read_pack(pack_root: Path, expected_manifest_sha256: str):
    _require(isinstance(expected_manifest_sha256, str)
             and len(expected_manifest_sha256) == 64
             and all(c in "0123456789abcdef" for c in expected_manifest_sha256),
             "COMPUTED_LIGAND_PKA_MANIFEST_SHA256")
    try:
        _require_unlinked_chain(pack_root)
        _require(pack_root.is_dir() and not _is_link(pack_root),
                 "COMPUTED_LIGAND_PKA_PACK")
        strict_root = pack_root.resolve(strict=True)
        _require(strict_root == pack_root.resolve(), "COMPUTED_LIGAND_PKA_PACK")
        names = []
        for directory, directories, files in os.walk(pack_root, followlinks=False):
            base = Path(directory)
            for name in directories:
                child = base / name
                _require(not _is_link(child), "COMPUTED_LIGAND_PKA_LINK")
            for name in files:
                child = base / name
                _require(not _is_link(child), "COMPUTED_LIGAND_PKA_LINK")
                relative = child.relative_to(pack_root).as_posix()
                _safe_regular_file(pack_root, relative)
                names.append(relative)
    except ContractError:
        raise
    except OSError as error:
        raise ContractError("COMPUTED_LIGAND_PKA_PACK") from error

    expected_names = set(ALLOWED_PATHS) | {MANIFEST_NAME}
    _require(len(names) == 11 and len(set(names)) == 11
             and set(names) == expected_names,
             "COMPUTED_LIGAND_PKA_PACK_CLOSURE")
    raw_files = {}
    total = 0
    for name in sorted(names):
        try:
            path = _safe_regular_file(pack_root, name)
            with path.open("rb") as handle:
                raw = handle.read(MAX_PACK_BYTES - total + 1)
        except ContractError:
            raise
        except OSError as error:
            raise ContractError("COMPUTED_LIGAND_PKA_FILE") from error
        _require(len(raw) <= MAX_PACK_BYTES - total,
                 "COMPUTED_LIGAND_PKA_PACK_SIZE")
        total += len(raw)
        raw_files[name] = raw

    manifest_raw = raw_files[MANIFEST_NAME]
    _require(_sha(manifest_raw) == expected_manifest_sha256,
             "COMPUTED_LIGAND_PKA_MANIFEST_HASH")
    try:
        manifest = _strict_json(manifest_raw)
    except (UnicodeError, ValueError, TypeError) as error:
        raise ContractError("COMPUTED_LIGAND_PKA_MANIFEST_JSON") from error
    _require(isinstance(manifest, dict), "COMPUTED_LIGAND_PKA_MANIFEST")
    _require(set(manifest) == {
        "format", "files", "final_state_selection_performed",
        "population_prediction_performed", "scientific_approved",
    }, "COMPUTED_LIGAND_PKA_MANIFEST")
    _require(manifest["format"] == FORMAT
             and manifest["final_state_selection_performed"] is False
             and manifest["population_prediction_performed"] is False
             and manifest["scientific_approved"] is False,
             "COMPUTED_LIGAND_PKA_MANIFEST")
    entries = manifest["files"]
    _require(isinstance(entries, list) and len(entries) == 10,
             "COMPUTED_LIGAND_PKA_MANIFEST_FILES")
    seen = set()
    for entry in entries:
        _require(isinstance(entry, dict) and set(entry) == {"path", "sha256", "bytes"},
                 "COMPUTED_LIGAND_PKA_MANIFEST_FILE")
        name = entry["path"]
        _require(isinstance(name, str) and name in ALLOWED_PATHS and name not in seen,
                 "COMPUTED_LIGAND_PKA_MANIFEST_FILE")
        seen.add(name)
        raw = raw_files[name]
        _require(type(entry["bytes"]) is int and entry["bytes"] == len(raw)
                 and isinstance(entry["sha256"], str)
                 and entry["sha256"] == _sha(raw),
                 "COMPUTED_LIGAND_PKA_FILE_INTEGRITY")
    _require(seen == ALLOWED_PATHS, "COMPUTED_LIGAND_PKA_MANIFEST_FILES")
    return raw_files


def _archive(raw_files: dict[str, bytes]) -> bytes:
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_STORED) as archive:
        for name in sorted(raw_files):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_STORED
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            archive.writestr(info, raw_files[name])
    return output.getvalue()


def _validate_computed_document(document: dict) -> None:
    _require(document.get("scientific_approved") is False
             and document.get("state_selection_performed") is False,
             "COMPUTED_LIGAND_PKA_DOCUMENT_FLAGS")
    rows = document.get("records")
    _require(isinstance(rows, list), "COMPUTED_LIGAND_PKA_RECORDS")
    for row in rows:
        _require(isinstance(row, dict)
                 and row.get("source_kind") == "computed"
                 and "pKa" in row and "state_population" not in row,
                 "COMPUTED_LIGAND_PKA_RECORD")


def _collision() -> None:
    raise ContractError("COMPUTED_LIGAND_PKA_IDENTITY_COLLISION")


def import_computed_ligand_pka(store, project, job_id, pack_root, *,
                               expected_manifest_sha256) -> dict:
    service = ScientificAcceptanceService(store, project)
    raw_files = _read_pack(Path(pack_root), expected_manifest_sha256)
    archive_raw = _archive(raw_files)
    source_raw = raw_files["quantitative-source.json"]
    _require(len(source_raw) <= MAX_MICROSTATE_SOURCE_BYTES,
             "SCIENTIFIC_MICROSTATE_SOURCE_SIZE")
    document = service._strict_quantitative_document(source_raw)
    _validate_computed_document(document)
    canonical_sha = _canonical_sha(document)
    upload_sha = _sha(source_raw)
    paths = []
    try:
        with store.db() as db:
            db.execute("BEGIN IMMEDIATE")
            _, result, binding, source = service._verified_result(
                db, job_id, require_current=False)
            _require(source["source_inputs_current"] is True,
                     "SCIENTIFIC_DESIGN_SOURCE_CHANGED")
            service._validate_quantitative_source_document(
                db, job_id, result, document)
            rows = db.execute("""SELECT id,binding,ref,superseded
                FROM scientific_evidence WHERE project=? AND job_id=?
                AND kind='microstate_population_source' AND identity_key=?
                ORDER BY rowid""", (project, job_id, canonical_sha)).fetchall()
            if rows:
                if len(rows) != 1 or rows[0]["superseded"] != 0:
                    _collision()
                try:
                    old_binding = json.loads(rows[0]["binding"])
                    old_ref = json.loads(rows[0]["ref"])
                    registered = service._registered_population_source(
                        db, job_id, result, old_ref)
                    proof_ref = old_binding["proof_archive_ref"]
                    proof_raw = service._read_registered_ref(db, proof_ref)
                except (ContractError, KeyError, TypeError, ValueError,
                        json.JSONDecodeError):
                    _collision()
                expected = {
                    **binding,
                    "actor_id": None,
                    "import_mode": "explicit_local_computed_pack",
                    "source_kind": "computed",
                    "source_upload_sha256": upload_sha,
                    "source_document_sha256": canonical_sha,
                    "proof_pack_manifest_sha256": expected_manifest_sha256,
                    "source_inputs_current": source["source_inputs_current"],
                    "source_runtime_current": source["source_runtime_current"],
                    "source_state": "unverified_quantity_pending_expert_review",
                    "scientific_approval": False,
                    "formal_acceptance": False,
                }
                if (registered["id"] != rows[0]["id"]
                        or old_ref != registered["ref"]
                        or old_binding.get("proof_archive_ref") != proof_ref
                        or proof_raw != archive_raw
                        or any(old_binding.get(key) != value
                               for key, value in expected.items())):
                    _collision()
                made = registered
                created = False
            else:
                proof_ref = service._write(
                    db, paths, archive_raw, media="application/zip",
                    provenance="computed")
                evidence_binding = {
                    **binding, "actor_id": None,
                    "import_mode": "explicit_local_computed_pack",
                    "source_kind": "computed",
                    "source_upload_sha256": upload_sha,
                    "source_document_sha256": canonical_sha,
                    "proof_pack_manifest_sha256": expected_manifest_sha256,
                    "proof_archive_ref": proof_ref,
                    "source_inputs_current": source["source_inputs_current"],
                    "source_runtime_current": source["source_runtime_current"],
                    "source_state": "unverified_quantity_pending_expert_review",
                    "scientific_approval": False, "formal_acceptance": False,
                }
                made = service._register_evidence(
                    db, paths, job_id, "microstate_population_source",
                    canonical_sha, evidence_binding, document,
                    provenance="computed")
                created = made["created"]
            result_binding = made["binding"]
            return {
                "id": made["id"], "ref": made["ref"],
                "binding": result_binding, "record_count": len(document["records"]),
                "created": created, "scientific_approved": False,
            }
    except BaseException:
        for path in paths:
            path.unlink(missing_ok=True)
        raise
