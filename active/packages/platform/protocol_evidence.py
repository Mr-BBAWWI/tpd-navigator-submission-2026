"""Immutable registration of verified protocol diagnostic evidence packs."""
from __future__ import annotations

import hashlib
import io
import json
import math
import os
import sqlite3
import uuid
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any

from jsonschema import Draft202012Validator

from packages.contracts import ContractError, encoded, now, parse_json, validate_payload
from packages.platform.scientific_acceptance import ScientificAcceptanceService
from packages.science import novel_ternary
from packages.science.protocol_evidence import ProtocolEvidenceError, verify_pack
from scripts import run_novel_ternary

ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = ROOT / "contracts" / "protocol-evidence" / "v0.1.0" / "schema.json"
SCHEMA_ID = "urn:tpd-navigator:raw:1"
STATUS = "computed_diagnostic_unreviewed"
OPERATOR_KIND = "explicit_local_cli"
TABLE_SQL = """
CREATE TABLE IF NOT EXISTS protocol_diagnostics(
    id TEXT PRIMARY KEY,
    project TEXT NOT NULL,
    job_id TEXT NOT NULL,
    protocol_id TEXT NOT NULL,
    kind TEXT NOT NULL,
    pack_manifest_sha256 TEXT NOT NULL,
    binding TEXT NOT NULL,
    record_ref TEXT NOT NULL,
    archive_ref TEXT NOT NULL,
    created_at TEXT NOT NULL,
    UNIQUE(project,job_id,protocol_id));
"""


def _require(condition: bool, code: str) -> None:
    if not condition:
        raise ContractError(code)


def _strict_loads(raw: str | bytes) -> Any:
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ContractError("PROTOCOL_EVIDENCE_DUPLICATE_JSON_KEY")
            result[key] = value
        return result

    def constant(_value):
        raise ContractError("PROTOCOL_EVIDENCE_NONFINITE_JSON")

    try:
        value = json.loads(raw, object_pairs_hook=pairs, parse_constant=constant)
    except ContractError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError) as error:
        raise ContractError("PROTOCOL_EVIDENCE_INVALID_JSON") from error
    _require(_finite(value), "PROTOCOL_EVIDENCE_NONFINITE_JSON")
    return value


def _finite(value: Any) -> bool:
    if isinstance(value, float):
        return math.isfinite(value)
    if isinstance(value, dict):
        return all(isinstance(key, str) and _finite(item) for key, item in value.items())
    if isinstance(value, list):
        return all(_finite(item) for item in value)
    return True


def _canonical(value: Any) -> bytes:
    _require(_finite(value), "PROTOCOL_EVIDENCE_NONFINITE_JSON")
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
                      allow_nan=False).encode("utf-8")


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _digest(value: Any) -> str:
    return _sha(_canonical(value))


def _valid_sha(value: Any) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def _safe_relative(value: Any) -> str:
    _require(isinstance(value, str) and value != "", "PROTOCOL_EVIDENCE_FILE_PATH")
    _require("\\" not in value and "\x00" not in value and ":" not in value,
             "PROTOCOL_EVIDENCE_FILE_PATH")
    path = PurePosixPath(value)
    _require(not path.is_absolute() and value == path.as_posix(), "PROTOCOL_EVIDENCE_FILE_PATH")
    _require(all(part not in ("", ".", "..") for part in path.parts),
             "PROTOCOL_EVIDENCE_FILE_PATH")
    return value


def _artifact_ref(identifier: str, raw: bytes, media: str, provenance: str) -> dict:
    ref = {
        "artifact_id": identifier,
        "version": 1,
        "sha256": _sha(raw),
        "media_type": media,
        "schema_id": SCHEMA_ID,
        "provenance": "computed",
    }
    validate_payload("ArtifactRef", ref)
    return ref


class ProtocolEvidenceService:
    """Register diagnostic-only packs without changing scientific acceptance state."""

    def __init__(self, store, project):
        self.store = store
        self.project = project
        schema = _strict_loads(SCHEMA_PATH.read_bytes())
        self.validator = Draft202012Validator(schema)
        # Construct this service before this component opens its own write transaction;
        # its constructor may perform schema initialization.
        self.acceptance = ScientificAcceptanceService(self.store, self.project)
        with self.store.db() as db:
            db.execute(TABLE_SQL)
        sources = {
            "packages/platform/protocol_evidence.py": Path(__file__).resolve(),
            "packages/science/protocol_evidence.py":
                ROOT / "packages" / "science" / "protocol_evidence.py",
            "contracts/protocol-evidence/v0.1.0/schema.json": SCHEMA_PATH,
        }
        self.verification_engine = {
            "version": 1,
            "source_hashes": {name: _sha(path.read_bytes()) for name, path in sources.items()},
        }

    def _validate_record(self, record: dict) -> None:
        errors = sorted(self.validator.iter_errors(record), key=lambda error: list(error.path))
        if errors:
            path = ".".join(str(part) for part in errors[0].absolute_path)
            raise ContractError("PROTOCOL_EVIDENCE_RECORD_SCHEMA" + (":" + path if path else ""))
        _require(_finite(record), "PROTOCOL_EVIDENCE_NONFINITE_JSON")

    def _candidate_graphs(self, result: dict) -> dict:
        candidates = result.get("protac_candidates")
        _require(isinstance(candidates, list), "PROTOCOL_EVIDENCE_CANDIDATES")
        graphs = {}
        for item in candidates:
            _require(isinstance(item, dict) and isinstance(item.get("candidate_id"), str),
                     "PROTOCOL_EVIDENCE_CANDIDATE")
            candidate_id = item["candidate_id"]
            _require(candidate_id not in graphs, "PROTOCOL_EVIDENCE_CANDIDATE_DUPLICATE")
            try:
                candidate = run_novel_ternary._candidate(result, candidate_id)
                graphs[candidate_id] = novel_ternary._mapped_graph(candidate)
            except Exception as error:
                raise ContractError("PROTOCOL_EVIDENCE_CANDIDATE_GRAPH") from error
        return graphs

    def _assessment_link(self, db, assessment_id: str | None, job_id: str,
                         binding: dict) -> dict | None:
        if assessment_id is None:
            return None
        row = db.execute(
            "SELECT id,project,job_id,revision,ref FROM scientific_assessments "
            "WHERE id=? AND project=? AND job_id=?",
            (assessment_id, self.project, job_id),
        ).fetchone()
        _require(row is not None, "PROTOCOL_EVIDENCE_ASSESSMENT_NOT_FOUND")
        try:
            ref = _strict_loads(row["ref"])
            document = _strict_loads(self.store.read(self.project, ref))
        except ContractError:
            raise
        except Exception as error:
            raise ContractError("PROTOCOL_EVIDENCE_ASSESSMENT_CORRUPT") from error
        assessment_binding = document.get("source_binding") or document.get("binding")
        _require(isinstance(assessment_binding, dict), "PROTOCOL_EVIDENCE_ASSESSMENT_BINDING")
        for key in ("project", "job_id", "input_sha256", "result_sha256"):
            _require(assessment_binding.get(key) == binding.get(key),
                     "PROTOCOL_EVIDENCE_ASSESSMENT_BINDING")
        status = document.get("status")
        if not isinstance(status, str):
            status = document.get("summary", {}).get("status")
        _require(status is None or isinstance(status, str), "PROTOCOL_EVIDENCE_ASSESSMENT_STATUS")
        return {"assessment_id": row["id"], "revision": row["revision"],
                "original_status_summary": status}

    def _verified_source(self, db, job_id: str, assessment_id: str | None):
        job, result, binding, _source = self.acceptance._verified_result(
            db, job_id, require_current=False
        )
        link = self._assessment_link(db, assessment_id, job_id, binding)
        return job, result, binding, self._candidate_graphs(result), link

    def _normalize_bundle(self, bundle: Any, expected_sha: str) -> tuple[dict, list[dict]]:
        _require(isinstance(bundle, dict), "PROTOCOL_EVIDENCE_VERIFIER_RESULT")
        _require(bundle.get("pack_manifest_sha256") == expected_sha,
                 "PROTOCOL_EVIDENCE_MANIFEST_MISMATCH")
        _require(_valid_sha(bundle.get("compact_sha256")), "PROTOCOL_EVIDENCE_COMPACT_HASH")
        files = bundle.get("files")
        protocols = bundle.get("protocols")
        _require(isinstance(files, list) and isinstance(protocols, list) and len(protocols) == 2,
                 "PROTOCOL_EVIDENCE_BUNDLE_SHAPE")
        normalized_files = []
        seen_paths = set()
        for item in files:
            _require(isinstance(item, dict) and set(item) >= {"path", "sha256", "bytes"},
                     "PROTOCOL_EVIDENCE_FILE_ENTRY")
            relative = _safe_relative(item["path"])
            folded = relative.casefold()
            _require(folded not in seen_paths, "PROTOCOL_EVIDENCE_FILE_DUPLICATE")
            _require(_valid_sha(item["sha256"]) and type(item["bytes"]) is int and item["bytes"] >= 0,
                     "PROTOCOL_EVIDENCE_FILE_ENTRY")
            seen_paths.add(folded)
            normalized_files.append({"path": relative, "sha256": item["sha256"],
                                     "bytes": item["bytes"]})
        _require("output-manifest.json".casefold() in seen_paths,
                 "PROTOCOL_EVIDENCE_MANIFEST_NOT_LISTED")
        normalized_protocols = []
        seen_protocols = set()
        for protocol in protocols:
            _require(isinstance(protocol, dict), "PROTOCOL_EVIDENCE_PROTOCOL")
            _require(set(protocol) >= {"protocol_id", "kind", "summary", "measurement_verification",
                                      "scientific_approved", "gates_affected"},
                     "PROTOCOL_EVIDENCE_PROTOCOL")
            _require(isinstance(protocol["protocol_id"], str) and protocol["protocol_id"] and
                     isinstance(protocol["kind"], str) and protocol["kind"],
                     "PROTOCOL_EVIDENCE_PROTOCOL")
            _require(protocol["protocol_id"] not in seen_protocols,
                     "PROTOCOL_EVIDENCE_PROTOCOL_DUPLICATE")
            _require(protocol["scientific_approved"] is False and
                     protocol["gates_affected"] is False,
                     "PROTOCOL_EVIDENCE_NEGATIVE_FLAGS")
            _require(isinstance(protocol["summary"], dict) and
                     isinstance(protocol["measurement_verification"], dict),
                     "PROTOCOL_EVIDENCE_PROTOCOL")
            _require(_finite(protocol), "PROTOCOL_EVIDENCE_NONFINITE_JSON")
            seen_protocols.add(protocol["protocol_id"])
            normalized_protocols.append(protocol)
        return {"pack_manifest_sha256": expected_sha,
                "compact_sha256": bundle["compact_sha256"], "files": normalized_files}, normalized_protocols

    def _reread_and_archive(self, root: Path, files: list[dict], expected_sha: str) -> bytes:
        _require(root.is_dir() and not root.is_symlink(), "PROTOCOL_EVIDENCE_PACK_ROOT")
        strict_root = root.resolve(strict=True)
        expected = {item["path"] for item in files}
        actual = set()
        for current, directories, names in os.walk(root, followlinks=False):
            base = Path(current)
            for name in directories:
                child = base / name
                _require(not child.is_symlink(), "PROTOCOL_EVIDENCE_FILE_SYMLINK")
                resolved = child.resolve(strict=True)
                _require(strict_root in resolved.parents,
                         "PROTOCOL_EVIDENCE_FILE_PATH")
            for name in names:
                child = base / name
                _require(not child.is_symlink(), "PROTOCOL_EVIDENCE_FILE_SYMLINK")
                resolved = child.resolve(strict=True)
                _require(resolved.is_file() and strict_root in resolved.parents,
                         "PROTOCOL_EVIDENCE_FILE_PATH")
                actual.add(child.relative_to(root).as_posix())
        _require(actual == expected, "PROTOCOL_EVIDENCE_FILE_SET_CHANGED")
        _require(len({name.casefold() for name in actual}) == len(actual),
                 "PROTOCOL_EVIDENCE_FILE_DUPLICATE")
        payloads = []
        for item in files:
            candidate = root.joinpath(*PurePosixPath(item["path"]).parts)
            try:
                _require(not candidate.is_symlink(), "PROTOCOL_EVIDENCE_FILE_SYMLINK")
                resolved = candidate.resolve(strict=True)
                _require(resolved.is_file() and not resolved.is_symlink() and
                         (resolved == strict_root or strict_root in resolved.parents),
                         "PROTOCOL_EVIDENCE_FILE_PATH")
                raw = resolved.read_bytes()
            except ContractError:
                raise
            except OSError as error:
                raise ContractError("PROTOCOL_EVIDENCE_FILE_READ") from error
            _require(len(raw) == item["bytes"] and _sha(raw) == item["sha256"],
                     "PROTOCOL_EVIDENCE_FILE_CHANGED")
            payloads.append((item["path"], raw))
        manifest = next(raw for path, raw in payloads if path == "output-manifest.json")
        _require(_sha(manifest) == expected_sha, "PROTOCOL_EVIDENCE_MANIFEST_MISMATCH")
        stream = io.BytesIO()
        with zipfile.ZipFile(stream, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for path, raw in sorted(payloads):
                info = zipfile.ZipInfo(path, (1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o100600 << 16
                archive.writestr(info, raw)
        return stream.getvalue()

    def _write_artifact(self, db, paths: list[Path], raw: bytes, media: str,
                        provenance: str) -> dict:
        identifier = "a-" + uuid.uuid4().hex
        ref = _artifact_ref(identifier, raw, media, provenance)
        blob_root = self.store.root / "blobs"
        _require(blob_root.is_dir() and not blob_root.is_symlink(),
                 "PROTOCOL_EVIDENCE_BLOB_ROOT")
        path = blob_root / identifier
        with path.open("xb") as handle:
            # Once exclusive creation succeeds, every subsequent failure must clean it up.
            paths.append(path)
            handle.write(raw)
            handle.flush()
            os.fsync(handle.fileno())
        db.execute("INSERT INTO artifacts(id,project,metadata) VALUES(?,?,?)",
                   (identifier, self.project,
                    json.dumps(ref, ensure_ascii=False, sort_keys=True, allow_nan=False)))
        return ref

    def _public(self, record: dict) -> dict:
        return dict(record)

    def register_pack(self, job_id, pack_root, *, expected_manifest_sha256,
                      expected_assessment_id=None):
        _require(_valid_sha(expected_manifest_sha256), "PROTOCOL_EVIDENCE_EXPECTED_SHA256")
        root = Path(pack_root)
        created_paths: list[Path] = []
        try:
            with self.store.db() as db:
                db.execute("BEGIN IMMEDIATE")
                _job, result, binding, graphs, assessment_link = self._verified_source(
                    db, job_id, expected_assessment_id
                )
                try:
                    verified = verify_pack(
                        root,
                        expected_manifest_sha256=expected_manifest_sha256,
                        expected_binding=binding,
                        expected_candidate_graphs=graphs,
                    )
                except ProtocolEvidenceError as error:
                    raise ContractError("PROTOCOL_EVIDENCE_PACK_INVALID") from error
                bundle, protocols = self._normalize_bundle(verified, expected_manifest_sha256)
                archive_raw = self._reread_and_archive(root, bundle["files"], expected_manifest_sha256)

                # Recheck the archived design binding while holding the write lock.
                _job2, _result2, binding2, _graphs2, link2 = self._verified_source(
                    db, job_id, expected_assessment_id
                )
                _require(binding2 == binding and link2 == assessment_link,
                         "PROTOCOL_EVIDENCE_SOURCE_CHANGED")

                existing = {}
                missing = []
                for protocol in protocols:
                    row = db.execute(
                        "SELECT * FROM protocol_diagnostics WHERE project=? AND job_id=? AND protocol_id=?",
                        (self.project, job_id, protocol["protocol_id"]),
                    ).fetchone()
                    if row is None:
                        missing.append(protocol["protocol_id"])
                        continue
                    stored = self._read_row(row)
                    _require(
                        stored["kind"] == protocol["kind"] and
                        stored["source_binding"] == binding and
                        stored.get("assessment_link") == assessment_link and
                        stored["summary"] == protocol["summary"] and
                        stored["measurement_verification"] == protocol["measurement_verification"],
                        "PROTOCOL_EVIDENCE_PROTOCOL_CONFLICT",
                    )
                    existing[protocol["protocol_id"]] = stored

                archive_ref = None
                created_at = None
                if missing:
                    archive_ref = self._write_artifact(
                        db, created_paths, archive_raw, "application/zip",
                        "verified_protocol_evidence_bundle"
                    )
                    created_at = now()

                records = []
                for protocol in protocols:
                    stored = existing.get(protocol["protocol_id"])
                    if stored is not None:
                        records.append(stored)
                        continue
                    record = {
                        "id": "pdiag-" + uuid.uuid4().hex,
                        "project_id": self.project,
                        "job_id": job_id,
                        "protocol_id": protocol["protocol_id"],
                        "kind": protocol["kind"],
                        "created_at": created_at,
                        "source_binding": binding,
                        "assessment_link": assessment_link,
                        "summary": protocol["summary"],
                        "measurement_verification": protocol["measurement_verification"],
                        "pack_manifest_sha256": expected_manifest_sha256,
                        "compact_sha256": bundle["compact_sha256"],
                        "archive_ref": archive_ref,
                        "verification_engine": self.verification_engine,
                        "scientific_approved": False,
                        "formal_acceptance": False,
                        "gates_affected": False,
                        "status": STATUS,
                        "operator_kind": OPERATOR_KIND,
                    }
                    self._validate_record(record)
                    record_ref = self._write_artifact(
                        db, created_paths, _canonical(record), "application/json",
                        "protocol_diagnostic_record"
                    )
                    db.execute(
                        "INSERT INTO protocol_diagnostics "
                        "(id,project,job_id,protocol_id,kind,pack_manifest_sha256,binding,record_ref,archive_ref,created_at) "
                        "VALUES(?,?,?,?,?,?,?,?,?,?)",
                        (record["id"], self.project, job_id, protocol["protocol_id"],
                         protocol["kind"], expected_manifest_sha256,
                         json.dumps(binding, ensure_ascii=False, sort_keys=True, allow_nan=False),
                         json.dumps(record_ref, ensure_ascii=False, sort_keys=True, allow_nan=False),
                         json.dumps(archive_ref, ensure_ascii=False, sort_keys=True, allow_nan=False),
                         created_at),
                    )
                    records.append(record)
                return {"created": bool(missing), "records": records}
        except Exception:
            for path in created_paths:
                try:
                    path.unlink()
                except FileNotFoundError:
                    pass
            raise

    def _read_row(self, row) -> dict:
        try:
            binding = _strict_loads(row["binding"])
            record_ref = _strict_loads(row["record_ref"])
            archive_ref = _strict_loads(row["archive_ref"])
            raw = self.store.read(self.project, record_ref)
            record = _strict_loads(raw)
            archive_raw = self.store.read(self.project, archive_ref)
        except ContractError:
            raise
        except Exception as error:
            raise ContractError("PROTOCOL_EVIDENCE_RECORD_CORRUPT") from error
        self._validate_record(record)
        _require(_canonical(record) == raw, "PROTOCOL_EVIDENCE_RECORD_ENCODING")
        _require(record["id"] == row["id"] and record["project_id"] == self.project and
                 record["job_id"] == row["job_id"] and
                 record["protocol_id"] == row["protocol_id"] and
                 record["kind"] == row["kind"] and
                 record["created_at"] == row["created_at"] and
                 record["pack_manifest_sha256"] == row["pack_manifest_sha256"] and
                 record["source_binding"] == binding and record["archive_ref"] == archive_ref,
                 "PROTOCOL_EVIDENCE_ROW_MISMATCH")
        _require(_sha(archive_raw) == archive_ref["sha256"], "PROTOCOL_EVIDENCE_ARCHIVE_HASH")
        return self._public(record)

    def list(self, job_id):
        with self.store.db() as db:
            jobs_table = db.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='workbench_jobs'"
            ).fetchone()
            if jobs_table is None:
                raise KeyError(job_id)
            job = db.execute(
                "SELECT 1 FROM workbench_jobs WHERE id=? AND project=?",
                (job_id, self.project),
            ).fetchone()
            if job is None:
                raise KeyError(job_id)
            table = db.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='protocol_diagnostics'"
            ).fetchone()
            if table is None:
                return []
            rows = db.execute(
                "SELECT * FROM protocol_diagnostics WHERE project=? AND job_id=? ORDER BY created_at,id",
                (self.project, job_id),
            ).fetchall()
        return [self._read_row(row) for row in rows]

    def view(self, diagnostic_id):
        with self.store.db() as db:
            table = db.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='protocol_diagnostics'"
            ).fetchone()
            if table is None:
                raise KeyError(diagnostic_id)
            row = db.execute(
                "SELECT * FROM protocol_diagnostics WHERE project=? AND id=?",
                (self.project, diagnostic_id),
            ).fetchone()
        if row is None:
            raise KeyError(diagnostic_id)
        return self._read_row(row)
