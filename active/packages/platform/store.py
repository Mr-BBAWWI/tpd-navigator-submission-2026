"""Immutable local artifacts and transactional run records for the I2 loopback app."""
import hashlib
import json
import sqlite3
import uuid
from contextlib import contextmanager
from pathlib import Path

from packages.contracts import DOMAIN_URI, ContractError, parse_json, validate_payload, encoded, now


class Store:
    def __init__(self, root):
        self.root = Path(root).resolve()
        (self.root / "blobs").mkdir(parents=True, exist_ok=True)
        self.db_path = self.root / "index.sqlite3"
        with self.db() as db:
            db.executescript("""
              CREATE TABLE IF NOT EXISTS artifacts(id TEXT PRIMARY KEY, project TEXT NOT NULL, metadata TEXT NOT NULL);
              CREATE TABLE IF NOT EXISTS runs(id TEXT PRIMARY KEY, project TEXT NOT NULL, request_key TEXT NOT NULL,
                request_hash TEXT NOT NULL, revision INTEGER NOT NULL, state TEXT NOT NULL, body TEXT NOT NULL,
                UNIQUE(project,request_key));
            """)

    @contextmanager
    def db(self):
        connection = sqlite3.connect(self.db_path, timeout=20)
        connection.row_factory = sqlite3.Row
        try:
            with connection:
                yield connection
        finally:
            connection.close()

    def scope(self, project_id):
        return ArtifactPort(self, project_id)

    def put(self, project, data, media, schema, provenance):
        identifier = "a-" + uuid.uuid4().hex
        digest = hashlib.sha256(data).hexdigest()
        reference = {"artifact_id": identifier, "version": 1, "sha256": digest, "media_type": media,
                     "schema_id": schema, "provenance": provenance}
        validate_payload("ArtifactRef", reference)
        path = self.root / "blobs" / identifier
        with path.open("xb") as handle:
            handle.write(data)
        with self.db() as db:
            db.execute("INSERT INTO artifacts VALUES(?,?,?)", (identifier, project, json.dumps(reference)))
        return reference

    def reference(self, project, artifact_id):
        with self.db() as db:
            row = db.execute("SELECT metadata FROM artifacts WHERE id=? AND project=?", (artifact_id, project)).fetchone()
        if row is None:
            raise KeyError("Artifact not registered in this project")
        return json.loads(row["metadata"])

    def read(self, project, reference):
        registered = self.reference(project, reference["artifact_id"])
        if reference != registered:
            raise ContractError("REGISTERED_ARTIFACT_MISMATCH")
        data = (self.root / "blobs" / registered["artifact_id"]).read_bytes()
        if hashlib.sha256(data).hexdigest() != registered["sha256"]:
            raise ContractError("STORED_ARTIFACT_CORRUPTED")
        return data

    def create_run(self, project, request_key, request):
        fingerprint = hashlib.sha256(encoded(request)).hexdigest()
        with self.db() as db:
            db.execute("BEGIN IMMEDIATE")
            existing = db.execute("SELECT request_hash,body FROM runs WHERE project=? AND request_key=?", (project, request_key)).fetchone()
            if existing:
                if existing["request_hash"] != fingerprint:
                    raise ValueError("같은 요청 키에 다른 입력을 보낼 수 없습니다.")
                return json.loads(existing["body"]), False
            pending = db.execute("SELECT COUNT(*) FROM runs WHERE project=? AND state IN ('queued','resolving','collecting')", (project,)).fetchone()[0]
            if pending >= 3:
                raise ValueError("처리 중인 작업이 많습니다. 기존 작업이 끝난 뒤 다시 실행해 주세요.")
            identifier = "r-" + uuid.uuid4().hex
            run = {"id": identifier, "project_id": project, "revision": 1, "state": "queued", "created_at": now(),
                   "updated_at": now(), "request": request, "target_ref": None, "bundle_ref": None,
                   "target_query_ref": None, "jobs": [], "issues": [], "selection": None,
                   "message": "표적 식별 대기 중", "data_mode": "real"}
            db.execute("INSERT INTO runs VALUES(?,?,?,?,?,?,?)", (identifier, project, request_key, fingerprint, 1, "queued", json.dumps(run)))
        return run, True

    def get_run(self, project, identifier):
        with self.db() as db:
            row = db.execute("SELECT body FROM runs WHERE project=? AND id=?", (project, identifier)).fetchone()
        if row is None:
            raise KeyError("Run not found in this project")
        return json.loads(row["body"])

    def update(self, project, identifier, updates, expected_revision=None, expected_state=None):
        with self.db() as db:
            db.execute("BEGIN IMMEDIATE")
            row = db.execute("SELECT body FROM runs WHERE project=? AND id=?", (project, identifier)).fetchone()
            if row is None:
                raise KeyError(identifier)
            run = json.loads(row["body"])
            if expected_revision is not None and run["revision"] != expected_revision:
                raise ValueError("입력 버전이 바뀌었습니다. 화면을 새로 읽어 주세요.")
            if expected_state is not None and run["state"] != expected_state:
                raise ValueError("이미 처리되었거나 현재 상태에서 선택할 수 없습니다.")
            run.update(updates)
            run["updated_at"] = now()
            db.execute("UPDATE runs SET revision=?,state=?,body=? WHERE id=?", (run["revision"], run["state"], json.dumps(run), identifier))
            return run

    def list_runs(self, project):
        with self.db() as db:
            rows = db.execute("SELECT body FROM runs WHERE project=? ORDER BY rowid DESC LIMIT 30", (project,)).fetchall()
        return [json.loads(r["body"]) for r in rows]

    def interrupt_pending(self, project):
        with self.db() as db:
            rows = db.execute("SELECT id FROM runs WHERE project=? AND state IN ('queued','resolving','collecting')", (project,)).fetchall()
        for row in rows:
            self.update(project, row["id"], {"state": "interrupted", "message": "서버가 재시작되어 중단되었습니다. 저장된 결과는 유지됩니다."})


class ArtifactPort:
    """A scoped port injected into M3/M4, without run-state or SQL access."""
    def __init__(self, store, project):
        self._store, self.project = store, project

    def read(self, reference):
        return self._store.read(self.project, reference)

    def put_raw(self, data, media="application/json", provenance="source", schema="urn:tpd-navigator:raw:1"):
        return self._store.put(self.project, data, media, schema, provenance)

    def put_json(self, value, kind=None, provenance="computed"):
        if kind:
            validate_payload(kind, value)
        return self.put_raw(encoded(value), "application/json", provenance, DOMAIN_URI + "#/$defs/" + kind if kind else "urn:tpd-navigator:raw:1")

    def json(self, reference):
        return parse_json(self.read(reference))
