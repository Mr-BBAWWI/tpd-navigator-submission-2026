"""M2 local CPU runner: two explicit capabilities, immutable contracts, saved views."""
import time
import uuid
import json
from pathlib import Path

from packages.contracts import ArtifactReader, ContractError, encoded, issue, missing, known, payload, validate_exchange
from packages.literature.collector import Collector, questions_for
from packages.science.target import resolve_target, select_target
from packages.transport import PublicTransport, FetchError

CAPABILITIES = {"resolve_target": ("TargetQuery", "TargetResolution", "science"),
                "collect_literature": ("LiteratureCollectionRequest", "EvidenceBundle", "literature")}


class Orchestrator:
    def __init__(self, store, project="local-research", transport_factory=PublicTransport):
        self.store, self.project = store, project
        self.transport_factory = transport_factory

    def record_job(self, run, operation, request_ref, inputs, policy_ref, work):
        if operation not in CAPABILITIES:
            raise ValueError("I2에서 허용되지 않은 작업입니다.")
        registry = json.loads((Path(__file__).resolve().parents[2] / "contracts/v0.1.0/operations.json").read_text(encoding="utf-8"))
        capability = next(item for item in registry["operations"] if item["operation_id"] == operation)
        if not capability["live_execution_ready"] or "local-single-user-i2" not in capability.get("execution_profiles", []):
            raise ValueError("현재 로컬 실행 프로필에서 비활성화된 작업입니다.")
        port = self.store.scope(self.project)
        _, output_kind, owner = CAPABILITIES[operation]
        parameters = port.put_json(payload("OperationParameters", profile_id="local-cpu-i2-v1"), "OperationParameters")
        manifest = payload("InputManifest", boundary_version="0.1.0", project_id=self.project, run_id=run["id"],
            job_id="j-" + uuid.uuid4().hex, input_revision=run["revision"], operation_id=operation,
            execution_mode="development_live", data_mode="real", payload_ref=request_ref, input_artifacts=inputs,
            parameters_ref=parameters, policy_ref=policy_ref)
        manifest_ref = port.put_json(manifest, "InputManifest")
        job = {"contract_version": "0.1.0", **{k: v for k, v in manifest.items() if k not in {"payload_type", "payload_version", "boundary_version", "payload_ref"}},
               "attempt": 1, "input_manifest": manifest_ref, "input_digest": manifest_ref["sha256"], "approval_ref": None, "reservation_id": None,
               "limits": {"wall_time_seconds": 240, "max_llm_requests": 0, "max_llm_tokens": 0, "max_gpu_seconds": 0, "max_output_bytes": 30_000_000}}
        reader = ArtifactReader(port.read)
        validate_exchange(job, reader)
        job_ref = port.put_raw(encoded(job), "application/json", "computed", "urn:tpd-navigator:boundary:0.1.0#/$defs/JobSpec")
        start = time.monotonic()
        output_ref, output = None, None
        try:
            output = work()
            output_ref = port.put_json(output, output_kind)
            notices = output.get("issues", [])
            status = "partial" if notices else "succeeded"
        except Exception as exc:
            code = exc.code if isinstance(exc, FetchError) else "OUTPUT_CONTRACT_REJECTED" if isinstance(exc, ContractError) else "MODULE_EXECUTION_FAILED"
            notices = [issue(code, str(exc) if isinstance(exc, (FetchError, ValueError)) else "자료 처리 중 오류가 발생했습니다.", "transient", retry="same_input_may_retry")]
            status = "failed"
        result = {k: job[k] for k in ("contract_version", "project_id", "run_id", "job_id", "attempt", "input_revision", "input_digest", "operation_id", "data_mode")}
        result.update(status=status, output_artifacts=[output_ref] if output_ref else [], issues=notices,
            actual_usage={"wall_time_ms": int((time.monotonic() - start) * 1000), "llm_requests": 0, "llm_tokens": 0, "gpu_seconds": 0, "usage_complete": True},
            producer={"module": owner, "code_version": "i2-local-0.1.0", "tool_versions_ref": None})
        try:
            validate_exchange(job, reader, result)
        except ContractError:
            # Reject the output, but preserve the failed attempt in the run journal.
            output_ref, output, status = None, None, "failed"
            notices = [issue("OUTPUT_CONTRACT_REJECTED", "모듈 출력이 계약/참조 검사를 통과하지 못해 사용하지 않았습니다.")]
            result.update(status=status, output_artifacts=[], issues=notices)
            validate_exchange(job, reader, result)
        result_ref = port.put_raw(encoded(result), "application/json", "computed", "urn:tpd-navigator:boundary:0.1.0#/$defs/ModuleResult")
        current = self.store.get_run(self.project, run["id"])
        self.store.update(self.project, run["id"], {"jobs": current["jobs"] + [{"operation_id": operation, "job_ref": job_ref, "result_ref": result_ref, "status": status}]})
        if status == "failed":
            raise FetchError(notices[0]["code"], notices[0]["message"])
        return output_ref, output

    def execute(self, identifier):
        run = self.store.get_run(self.project, identifier)
        port = self.store.scope(self.project)
        transport = self.transport_factory()
        try:
            self.store.update(self.project, identifier, {"state": "resolving", "message": "UniProt에서 표적을 식별하는 중"}, expected_revision=run["revision"], expected_state="queued")
        except ValueError:
            return  # A duplicate worker cannot overwrite the active/completed run.
        try:
            if run["selection"]:
                previous = port.json(run["target_ref"])
                query_ref = run["target_query_ref"]
                annotation_ref = run["selection"]["annotation_ref"]
                inputs = [query_ref, run["target_ref"], annotation_ref]
                work = lambda: select_target(previous, run["selection"]["candidate_id"], annotation_ref, port)
            else:
                request = run["request"]
                context = {"objective": known(request["objective"]) if request["objective"] else missing("연구 목적 미입력", "not_provided"),
                           "cell_line": known(request["cell_line"]) if request["cell_line"] else missing("세포주 미입력", "not_provided"), "tissue": missing("조직 미입력", "not_provided")}
                query = payload("TargetQuery", query=request["query"], query_kind=request["query_kind"], taxon_id=9606,
                    isoform=missing("isoform 미입력", "not_provided"), mutation=missing("변이 미입력", "not_provided"),
                    construct=missing("construct 미입력", "not_provided"), research_context=context)
                query_ref = port.put_json(query, "TargetQuery", "human_annotation")
                inputs = [query_ref]
                work = lambda: resolve_target(query_ref, port, transport)
            policy_ref = port.put_json({"profile": "local-human-target-v1", "taxon_id": 9606, "exclude_annotated_transmembrane": True})
            target_ref, target = self.record_job(run, "resolve_target", query_ref, inputs, policy_ref, work)
            self.store.update(self.project, identifier, {"target_ref": target_ref, "target_query_ref": query_ref})
            if target["resolution_status"] != "resolved":
                state = "awaiting_selection" if target["resolution_status"] == "ambiguous" else "not_found"
                self.store.update(self.project, identifier, {"state": state, "message": "표적 후보를 선택해 주세요." if state == "awaiting_selection" else "일치하는 표적을 찾지 못했습니다."})
                return
            if target["scope_status"] == "outside_current_scope":
                self.store.update(self.project, identifier, {"state": "unsupported", "message": "현재 제품의 표적 범위에 포함되지 않습니다."})
                return
            self.store.update(self.project, identifier, {"state": "collecting", "message": "문헌 수집 준비 중"})
            candidate = next(c for c in target["candidates"] if c["candidate_id"] == target["selected_candidate_id"])
            uploads = [self.store.reference(self.project, identifier) for identifier in run["request"]["pdf_artifact_ids"]]
            policy = payload("CollectionPolicy", provider="europe_pmc", max_records=run["request"]["max_documents"],
                max_documents=run["request"]["max_documents"], max_segments=180, max_source_bytes=12_000_000)
            policy_ref = port.put_json(policy, "CollectionPolicy")
            request = payload("LiteratureCollectionRequest", target_resolution_ref=target_ref,
                research_context=port.json(query_ref)["research_context"], questions=questions_for(candidate["identity"]),
                attached_source_refs=uploads, collection_policy_ref=policy_ref)
            request_ref = port.put_json(request, "LiteratureCollectionRequest")
            progress = lambda text: self.store.update(self.project, identifier, {"message": text})
            bundle_ref, bundle = self.record_job(run, "collect_literature", request_ref, [request_ref, target_ref] + uploads,
                policy_ref, lambda: Collector(port, transport).collect(request_ref, progress))
            state = "partial" if bundle["issues"] else "completed"
            self.store.update(self.project, identifier, {"bundle_ref": bundle_ref, "state": state, "issues": bundle["issues"],
                "message": "원문 수집과 텍스트 추출을 마쳤습니다. AI 분석은 아직 실행하지 않았습니다."})
        except Exception as exc:
            self.store.update(self.project, identifier, {"state": "failed", "message": "작업이 완료되지 않았습니다. 저장된 상태와 오류를 확인해 주세요.",
                "issues": [issue(exc.code if isinstance(exc, FetchError) else "PIPELINE_FAILED", str(exc), "internal_error", retry="same_input_may_retry")]})

    def choose(self, identifier, candidate_id, revision):
        run = self.store.get_run(self.project, identifier)
        if run["state"] != "awaiting_selection" or run["revision"] != revision:
            raise ValueError("선택 요청이 이미 처리되었거나 입력 버전이 바뀌었습니다.")
        port = self.store.scope(self.project)
        target = port.json(run["target_ref"])
        if candidate_id not in {c["candidate_id"] for c in target["candidates"]}:
            raise ValueError("표시된 후보에서 선택해 주세요.")
        annotation = port.put_json({"kind": "target_identity_selection", "candidate_id": candidate_id, "target_resolution_ref": run["target_ref"]}, provenance="human_annotation")
        return self.store.update(self.project, identifier, {"state": "queued", "revision": revision + 1,
            "selection": {"candidate_id": candidate_id, "annotation_ref": annotation}, "message": "선택한 표적으로 수집을 준비합니다."},
            expected_revision=revision, expected_state="awaiting_selection")

    def view(self, identifier):
        run = self.store.get_run(self.project, identifier)
        port = self.store.scope(self.project)
        return {**run, "target": port.json(run["target_ref"]) if run["target_ref"] else None,
                "bundle": port.json(run["bundle_ref"]) if run["bundle_ref"] else None}
