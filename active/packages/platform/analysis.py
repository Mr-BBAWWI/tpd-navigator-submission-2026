"""M2-owned durable analysis jobs; M3 receives only an artifact port."""
import hashlib
import json
import uuid
import os
from contextlib import contextmanager

from packages.contracts import encoded, now, ContractError
from packages.agents.provider import AgentError
from packages.literature.assessment import LiteratureAssessment
from packages.review_contracts import ReviewReader, save


@contextmanager
def worker_lock(root):
    """OS-released lock: an app startup cannot interrupt a live CLI/API worker."""
    with (root / "analysis-worker.lock").open("a+b") as handle:
        handle.seek(0, 2)
        if handle.tell() == 0:
            handle.write(b"0")
            handle.flush()
        handle.seek(0)
        try:
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            yield False
            return
        try:
            yield True
        finally:
            handle.seek(0)
            if os.name == "nt":
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


class AnalysisService:
    def __init__(self, store, project, provider_factory=None):
        self.store, self.project, self.provider_factory = store, project, provider_factory
        self.port = store.scope(project)
        with store.db() as db:
            db.execute("""CREATE TABLE IF NOT EXISTS analyses(
                id TEXT PRIMARY KEY, project TEXT NOT NULL, run_id TEXT NOT NULL, request_key TEXT NOT NULL,
                fingerprint TEXT NOT NULL, state TEXT NOT NULL, body TEXT NOT NULL,
                UNIQUE(project,run_id,request_key))""")

    def create(self, run_id, revision, segment_ids, request_key):
        if self.provider_factory is None:
            raise ValueError("AI 분석이 비활성 상태입니다. --enable-llm으로 로컬 서버를 시작해 주세요.")
        if not 1 <= len(segment_ids) <= 12 or len(set(segment_ids)) != len(segment_ids):
            raise ValueError("서로 다른 원문 구간 1~12개를 선택해 주세요.")
        fingerprint = hashlib.sha256(encoded({"revision": revision, "segments": sorted(segment_ids)})).hexdigest()
        # Validate source graph before queueing; no external calls here.
        run = self.store.get_run(self.project, run_id)
        if run.get('replay_only'):
            raise ValueError('가져온 재생 자료에서는 새 AI 분석을 실행할 수 없습니다.')
        if run["revision"] != revision or run["state"] not in {"completed", "partial"} or not run["bundle_ref"]:
            raise ValueError("문헌 수집 완료 상태와 현재 입력 버전을 확인해 주세요.")
        bundle = self.port.json(run["bundle_ref"])
        segments = {s["segment_id"]: s for s in bundle["segments"]}
        if not set(segment_ids) <= set(segments):
            raise ValueError("현재 문헌에 없는 구간이 포함되어 있습니다.")
        if sum(len(self.port.read(segments[s]["text_artifact_ref"])) for s in segment_ids) > 24000:
            raise ValueError("선택한 원문이 너무 깁니다. 구간 수를 줄여 주세요(합계 24 KB).")
        context = {"project_id": self.project, "run_id": run_id, "data_mode": run["data_mode"]}
        policy_ref = save(self.port, "PolicySnapshot", context, target_resolution_ref=run["target_ref"], fields=[],
                          limitations=["문헌 분석 전용. 과학 정책 및 G1 실행 조건은 아직 구성하지 않음."])
        request_ref = save(self.port, "LiteratureAssessmentRequest", context, input_revision=revision,
            evidence_bundle_ref=run["bundle_ref"], target_resolution_ref=run["target_ref"], base_policy_ref=policy_ref,
            allowed_field_ids=[], requested_segment_ids=segment_ids, research_context_ref=run["target_query_ref"],
            limits={"max_delivered_segments": 12, "max_claims": 8, "max_proposals": 0})
        ReviewReader(self.port.read).verify(request_ref)
        with self.store.db() as db:
            db.execute("BEGIN IMMEDIATE")
            existing = db.execute("SELECT body,fingerprint FROM analyses WHERE project=? AND run_id=? AND request_key=?",
                                  (self.project, run_id, request_key)).fetchone()
            if existing:
                if existing["fingerprint"] != fingerprint:
                    raise ValueError("같은 요청 키의 분석 입력을 바꿀 수 없습니다.")
                return json.loads(existing["body"]), False
            current = json.loads(db.execute("SELECT body FROM runs WHERE id=? AND project=?", (run_id, self.project)).fetchone()[0])
            if not self.matches(current, run):
                raise ValueError("입력이 바뀌었습니다. 새로 읽은 뒤 분석해 주세요.")
            pending = db.execute("SELECT COUNT(*) FROM analyses WHERE project=? AND state IN ('queued','running')", (self.project,)).fetchone()[0]
            if pending:
                raise ValueError("진행 중인 AI 분석이 있습니다. 완료 후 다시 요청해 주세요.")
            row = {"id": "analysis-" + uuid.uuid4().hex, "run_id": run_id, "project_id": self.project,
                   "input_revision": revision, "bundle_ref": run["bundle_ref"], "request_ref": request_ref,
                   "target_ref": run["target_ref"], "target_query_ref": run["target_query_ref"],
                   "data_mode": run["data_mode"], "state": "queued", "stage": "queued", "created_at": now(), "updated_at": now(),
                   "assessment_ref": None, "review_ref": None, "code_check_ref": None, "receipt_ref": None,
                   "trace": [], "calls": 0, "total_tokens": 0, "usage_status": "not_started", "message": "AI 분석 대기 중"}
            db.execute("INSERT INTO analyses VALUES(?,?,?,?,?,?,?)", (row["id"], self.project, run_id, request_key, fingerprint, row["state"], json.dumps(row)))
        return row, True

    @staticmethod
    def matches(current, snapshot):
        return (current["revision"] == snapshot.get("input_revision", snapshot.get("revision"))
                and current["bundle_ref"] == snapshot["bundle_ref"]
                and all(current[k] == snapshot[k] for k in ('target_ref','target_query_ref','data_mode') if k in snapshot))

    def get(self, identifier):
        with self.store.db() as db:
            row = db.execute("SELECT body FROM analyses WHERE id=? AND project=?", (identifier, self.project)).fetchone()
        if row is None:
            raise KeyError(identifier)
        return json.loads(row["body"])

    def update(self, identifier, fields, expected=None, finish=False):
        with self.store.db() as db:
            db.execute("BEGIN IMMEDIATE")
            row = db.execute("SELECT body FROM analyses WHERE id=? AND project=?", (identifier, self.project)).fetchone()
            if row is None:
                raise KeyError(identifier)
            value = json.loads(row["body"])
            if expected and value["state"] != expected:
                raise ValueError("ANALYSIS_STATE_CHANGED")
            value.update(fields, updated_at=now())
            if finish:
                current = json.loads(db.execute("SELECT body FROM runs WHERE id=? AND project=?", (value["run_id"], self.project)).fetchone()[0])
                if not self.matches(current, value):
                    value.update(state="stale", message="분석 중 입력이 바뀌어 현재 결과로 사용할 수 없습니다.")
            db.execute("UPDATE analyses SET state=?,body=? WHERE id=? AND project=?", (value["state"], json.dumps(value), identifier, self.project))
            return value

    def execute(self, identifier):
        with worker_lock(self.store.root) as acquired:
            if acquired:
                self._execute(identifier)

    def _execute(self, identifier):
        try:
            job = self.update(identifier, {"state": "running", "message": "원문 근거 분석 중"}, expected="queued")
        except ValueError:
            return  # Duplicate worker does not make another model call.
        provider = None
        try:
            provider = self.provider_factory()
            runtime = LiteratureAssessment(self.port, provider)
            result = runtime.run(job["request_ref"],
                checkpoint=lambda **fields: self.update(identifier, fields),
                is_current=lambda: self.get(identifier)["state"] == "running" and self.matches(self.store.get_run(self.project, job["run_id"]), job))
            self.update(identifier, result, expected="running", finish=True)
        except Exception as exc:
            # Never propagate source/model text or provider credentials in error messages.
            code = str(exc) if isinstance(exc, AgentError) else "LITERATURE_CONTRACT_REJECTED" if isinstance(exc, ContractError) else "LITERATURE_EXECUTION_FAILED"
            # AgentError originates in our adapters, not arbitrary exception text from documents.
            if not code.replace("_", "").isalnum() or len(code) > 100:
                code = "LITERATURE_EXECUTION_FAILED"
            if self.get(identifier)["state"] == "running":
                self.update(identifier, {"state": "failed", "error_code": code,
                    "message": "분석을 보류했습니다. 저장된 실행 기록을 확인해 주세요.",
                    "usage_status": "unreconciled" if code.startswith("API_") else self.get(identifier)["usage_status"]}, expected="running", finish=True)
        finally:
            if provider:
                provider.close()

    def interrupt_pending(self):
        with worker_lock(self.store.root) as acquired:
            if acquired:
                self._interrupt_pending()

    def _interrupt_pending(self):
        with self.store.db() as db:
            db.execute("BEGIN IMMEDIATE")
            rows = db.execute("SELECT body FROM analyses WHERE project=? AND state IN ('queued','running')", (self.project,)).fetchall()
            for row in rows:
                value = json.loads(row[0])
                value.update(state="interrupted", updated_at=now(), message="서버 재시작으로 중단. 자동 재호출하지 않습니다.")
                if value.get("pending_prompt_ref"):
                    value["usage_status"] = "unreconciled"
                db.execute("UPDATE analyses SET state=?,body=? WHERE id=?", (value["state"], json.dumps(value), value["id"]))

    def list(self, run_id):
        self.store.get_run(self.project, run_id)
        with self.store.db() as db:
            rows = db.execute("SELECT body FROM analyses WHERE project=? AND run_id=? ORDER BY rowid DESC LIMIT 20", (self.project, run_id)).fetchall()
        return [json.loads(r[0]) for r in rows]

    def view(self, identifier):
        value = self.get(identifier)
        # Old saved rows predate explicit target/context snapshots; their immutable
        # request still supplies those bindings without rewriting historical results.
        request = self.port.json(value['request_ref'])
        value['target_ref'] = request['target_resolution_ref']
        value['target_query_ref'] = request['research_context_ref']
        current = self.store.get_run(self.project, value["run_id"])
        value["current_input"] = self.matches(current, value)
        for name in ("assessment", "review", "code_check", "receipt"):
            ref = value.get(name + "_ref")
            if ref:
                ReviewReader(self.port.read).verify(ref)
                value[name] = self.port.json(ref)
        if value.get("assessment"):
            value["annotated_bundle"] = self.port.json(value["assessment"]["annotated_evidence_ref"])
        return value

    def g1_status(self, run_id):
        self.store.get_run(self.project, run_id)
        rows = self.list(run_id)
        latest = self.view(rows[0]["id"]) if rows else None
        reasons = ["B_PHASE1_AND_CPU_PREVIEW_NOT_CONNECTED", "COST_SCOPE_NOT_CONFIGURED", "HUMAN_REVIEW_NOT_RECORDED"]
        if latest is None or latest["state"] != "review_ready" or not latest["current_input"]:
            reasons.insert(0, "CURRENT_LITERATURE_REVIEW_NOT_READY")
        return {"gate": "G1", "ready": False, "dispatch_authorized": False, "packet_ref": None,
                "latest_analysis_id": latest["id"] if latest else None, "blocking_reasons": reasons}
