"""Loopback research app with opt-in literature LLM analysis; no GPU dispatch."""
import os
import re
from concurrent.futures import ThreadPoolExecutor
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Annotated, Literal, Union
from urllib.parse import urlsplit

from fastapi import FastAPI, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, ConfigDict, Field, field_validator

from packages.contracts import ContractError
from packages.platform.store import Store
from packages.platform.orchestrator import Orchestrator
from packages.platform.analysis import AnalysisService
from packages.platform.dossiers import DossierService
from packages.platform.workflows import WorkflowService
from packages.platform.replay import export_bytes
from packages.platform.saved_results import SavedResultsService
from packages.platform.evidence_reports import EvidenceReportService
from packages.platform.report_reviews import ReportReviewService
from packages.platform.review_identity import ReviewAuthError
from packages.platform.scientific_acceptance import ScientificAcceptanceService
from packages.platform.protocol_evidence import ProtocolEvidenceService
from packages.science.protocol_evidence import ProtocolEvidenceError
from packages.platform.workbench import WorkbenchService
from packages.science.handoff import HandoffError
from packages.transport import PublicTransport

ROOT = Path(__file__).resolve().parents[2]
PROJECT = "local-research"


class WorkbenchInput(BaseModel):
    model_config = ConfigDict(extra='forbid',strict=True)
    result_id: str = Field(min_length=1,max_length=80)
    operation: Literal['verify','report','review','candidate','chemistry','design_panel']
    request_key: str = Field(min_length=8,max_length=120)
    parameters: dict = Field(default_factory=dict)
    retry_of: str | None = Field(default=None,max_length=80)


class RunInput(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)
    query: str = Field(min_length=1, max_length=180)
    query_kind: Literal["name", "gene_symbol", "uniprot_accession"]
    objective: str = Field(default="", max_length=2000)
    cell_line: str = Field(default="", max_length=160)
    pdf_artifact_ids: list[str] = Field(default_factory=list, max_length=3)
    max_documents: int = Field(default=6, ge=2, le=12)
    request_key: str = Field(min_length=8, max_length=100, pattern=r"^[A-Za-z0-9_-]+$")


class SelectionInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    candidate_id: str = Field(min_length=1, max_length=128)
    revision: int = Field(ge=1)


class AnalysisInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    revision: int = Field(ge=1)
    segment_ids: list[str] = Field(min_length=1, max_length=12)
    request_key: str = Field(min_length=8, max_length=100, pattern=r"^[A-Za-z0-9_-]+$")


class WorkflowInput(BaseModel):
    model_config = ConfigDict(extra='forbid')
    b_input_id: str = Field(pattern=r'^binput-[a-f0-9]{32}$')
    request_key: str = Field(min_length=8,max_length=100,pattern=r'^[A-Za-z0-9_-]+$')


class EvidenceReportInput(BaseModel):
    model_config = ConfigDict(extra='forbid')
    expected_digest: str = Field(pattern=r'^[a-f0-9]{64}$')


class ReviewLoginInput(BaseModel):
    model_config = ConfigDict(extra='forbid')
    reviewer_id: str = Field(min_length=3,max_length=64)
    access_key: str = Field(min_length=1,max_length=128)


class HumanReviewRequestInput(BaseModel):
    model_config = ConfigDict(extra='forbid')
    expected_report_sha256: str = Field(pattern=r'^[a-f0-9]{64}$')
    scopes: list[Literal['chemical_state','structural_comparison','literature_conditions','synthesis_evidence']] = Field(min_length=1,max_length=4)
    supersedes: str | None = Field(default=None,pattern=r'^hreview-[a-f0-9]{32}$')


class ScopeJudgmentInput(BaseModel):
    model_config = ConfigDict(extra='forbid')
    scope_id: Literal['chemical_state','structural_comparison','literature_conditions','synthesis_evidence']
    verdict: Literal['accept_scope','request_changes','defer']
    reason: str = Field(min_length=1,max_length=5000)


class HumanReviewDecisionInput(BaseModel):
    model_config = ConfigDict(extra='forbid')
    expected_request_sha256: str = Field(pattern=r'^[a-f0-9]{64}$')
    expected_revision: int = Field(ge=0)
    idempotency_key: str = Field(pattern=r'^[A-Za-z0-9_-]{8,100}$')
    action: Literal['record','withdraw']
    judgments: list[ScopeJudgmentInput] = Field(max_length=4)
    acknowledgement_ids: list[Literal['scope_only','no_efficacy','no_execution']] = Field(max_length=3)
    withdrawal_reason: str = Field(max_length=5000)


class ScientificAssessmentInput(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    expected_result_sha256: str | None = Field(default=None, pattern=r'^[a-f0-9]{64}$')
    protocol_ids: list[str] | None = Field(default=None, max_length=8)

    @field_validator('protocol_ids')
    @classmethod
    def validate_protocol_ids(cls, value):
        if value is None:
            return value
        if any(not item.strip() or len(item) > 256 for item in value):
            raise ValueError('protocol_ids must contain non-empty strings of at most 256 characters')
        if len(value) != len(set(value)):
            raise ValueError('protocol_ids must not contain duplicates')
        return value

class ScientificComputeInput(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    candidate_ids: list[str] | None = Field(default=None, max_length=256)
    analog_ids: list[str] | None = Field(default=None, max_length=256)

class ScientificStatementInput(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    text: str = Field(min_length=1, max_length=50000)
    locator: str = Field(min_length=1, max_length=2000)

class MicrostatePopulationSourceInput(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    text: str = Field(min_length=1, max_length=1024 * 1024)

class MicrostatePopulationEvidenceInput(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    analog_id: str = Field(min_length=1, max_length=256)
    report_id: str = Field(min_length=1, max_length=256)
    state_index: int = Field(ge=0)
    pH: float
    evidence_ref: dict
    locator: str = Field(min_length=1, max_length=2000)
    review_ref: dict
    interpretation: str = Field(min_length=1, max_length=10000)

class ParentChoice(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    parent_id: str
    expert_selected: Literal[True]
    rationale: str
    source_ref: dict
class ParentData(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    parents: list[ParentChoice] = Field(min_length=1, max_length=10)
class SiteChoice(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    atom_map: int = Field(ge=1)
    state: Literal['MODIFIABLE','PROTECTED','UNKNOWN']
    allowed_rule_ids: list[str] = Field(min_length=1)
    rationale: str
    source_ref: dict
    allow_release_protected: bool
    release_justification: str | None = None
    release_source_ref: dict | None = None
class SiteData(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    parent_id: str
    sites: list[SiteChoice] = Field(min_length=1)
class Requirement(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    id: str
    kind: Literal['directional_hbond','salt_bridge','aromatic','hydrophobic','proximity']
    ligand_maps: list[int] = Field(min_length=1)
    protein_atom_ids: list[str] = Field(min_length=1)
    required: bool
class InteractionData(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    parent_id: str
    analog_id: str
    requirements: list[Requirement] = Field(min_length=1)
class PHConditions(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    pH: float = Field(ge=0, le=14)
    description: str
class MicrostateData(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    state_choices: dict[str,int]
    pH_conditions: PHConditions
    population_evidence: list[MicrostatePopulationEvidenceInput] = Field(
        default_factory=list, max_length=256)
class CalibrationData(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    metric: str
    threshold: float
    comparison: Literal['max_lte','min_gte']
    expected_seeds: list[int] = Field(min_length=1)
    scope: str
class GeometryData(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    disposition: Literal['accept','reject']
    metric: str
    comparison: Literal['lte','gte','max_lte','min_gte']
    threshold: float
    reference_free: Literal[True]
    rationale: str
    source_ref: dict
class TernaryData(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    chosen_per_e3: dict[Literal['CRBN','VHL'],str]
    expected_unique_seeds: list[int] = Field(min_length=2, max_length=32)
    geometry: GeometryData
class SynthesisData(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    chosen_candidate_ids: list[str] = Field(min_length=1)
    candidate_policies: list[dict] = Field(min_length=1)
    exact_route_records: list[dict]
class EmptyData(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
class ScientificDecisionBase(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    action: Literal['accept']
    intent: str
    reason: str
    expected_policy_revision: int = Field(ge=0)
class ParentDecision(ScientificDecisionBase):
    kind: Literal['parent_funnel']
    data: ParentData
class SiteDecision(ScientificDecisionBase):
    kind: Literal['site_policy']
    data: SiteData
class InteractionDecision(ScientificDecisionBase):
    kind: Literal['interaction_requirements']
    data: InteractionData
class MicrostateDecision(ScientificDecisionBase):
    kind: Literal['microstate_decision']
    data: MicrostateData
class CalibrationDecision(ScientificDecisionBase):
    kind: Literal['calibration_criterion']
    data: CalibrationData
class TernaryDecision(ScientificDecisionBase):
    kind: Literal['ternary_criterion']
    data: TernaryData
class SynthesisDecision(ScientificDecisionBase):
    kind: Literal['synthesis_policy']
    data: SynthesisData
class FormalDecision(BaseModel):
    model_config = ConfigDict(extra='forbid', strict=True)
    kind: Literal['formal_decision']
    action: Literal['accept','reject','defer']
    intent: str
    reason: str
    expected_policy_revision: int = Field(ge=0)
    data: EmptyData
ScientificDecisionInput = Annotated[Union[ParentDecision,SiteDecision,InteractionDecision,MicrostateDecision,CalibrationDecision,TernaryDecision,SynthesisDecision,FormalDecision],Field(discriminator='kind')]


def create_app(data_root=None, transport_factory=PublicTransport, enable_worker=True, provider_factory=None, read_only=False, review_only=False):
    store = Store(data_root or os.environ.get("TPD_DATA_DIR", ROOT / ".localdata"))
    orchestrator = Orchestrator(store, PROJECT, transport_factory)
    analysis = AnalysisService(store, PROJECT, provider_factory)
    dossiers = DossierService(store, PROJECT)
    workflows = WorkflowService(store, PROJECT, provider_factory)
    saved_results = SavedResultsService(store, PROJECT)
    evidence_reports = EvidenceReportService(store, PROJECT)
    report_reviews = ReportReviewService(store, PROJECT)
    scientific_acceptance = ScientificAcceptanceService(store, PROJECT)
    protocol_evidence = ProtocolEvidenceService(store, PROJECT)
    workbench = WorkbenchService(store, PROJECT, provider_factory)
    if review_only: enable_worker=False
    executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="tpd-i2")

    @asynccontextmanager
    async def lifespan(app):
        if not read_only and not review_only:
            store.interrupt_pending(PROJECT)
            analysis.interrupt_pending()
            workflows.interrupt_pending()
            workbench.recover()
        yield
        executor.shutdown(wait=False, cancel_futures=True)

    app = FastAPI(title="TPD Navigator local I2", version="0.1.0", lifespan=lifespan, docs_url=None, redoc_url=None)
    app.state.store = store
    app.state.orchestrator = orchestrator
    app.state.analysis = analysis
    app.state.dossiers = dossiers
    app.state.workflows = workflows
    app.state.saved_results = saved_results
    app.state.evidence_reports = evidence_reports
    app.state.report_reviews = report_reviews
    app.state.scientific_acceptance = scientific_acceptance
    app.state.protocol_evidence = protocol_evidence
    app.state.workbench = workbench

    @app.middleware("http")
    async def local_boundary(request: Request, call_next):
        if request.url.hostname not in {"127.0.0.1", "localhost", "::1", "testserver"}:
            return JSONResponse({"detail": "이 실행기는 로컬 개발용입니다."}, status_code=403)
        if request.method not in {"GET", "HEAD", "OPTIONS"}:
            if read_only:
                return JSONResponse({"detail": "저장 자료 열람 전용입니다."}, status_code=403)
            review_mutation = request.url.path in {'/api/review-auth/login','/api/review-auth/logout','/api/scientific-statements'} or bool(re.fullmatch(r'/api/(?:report-reviews/hreview-[a-f0-9]{32}/decisions|evidence-reports/report-[a-f0-9]{32}/review-requests|lab/jobs/[^/]+/(?:scientific-assessments|microstate-population-sources)|scientific-assessments/[^/]+/(?:decisions|export)|scientific-decisions/[^/]+/withdraw)',request.url.path))
            if review_only and not review_mutation:
                return JSONResponse({'detail':'검토 기록 전용 환경입니다.'},status_code=403)
            origin = request.headers.get("origin")
            if request.headers.get("x-tpd-local") != "1" or (origin and urlsplit(origin).netloc != request.headers.get("host")):
                return JSONResponse({"detail": "같은 로컬 화면에서 요청해 주세요."}, status_code=403)
        response = await call_next(request)
        response.headers["Cache-Control"] = "no-store"
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["Content-Security-Policy"] = "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; object-src 'none'; frame-ancestors 'none'"
        return response

    @app.exception_handler(KeyError)
    async def missing(request, exc):
        return JSONResponse({"detail": "이 프로젝트에서 자료를 찾지 못했습니다."}, status_code=404)

    @app.exception_handler(ReviewAuthError)
    async def review_auth_error(request, exc):
        return JSONResponse({'detail':'검수자 로그인이 필요하거나 세션이 만료되었습니다.'},status_code=401)

    @app.exception_handler(ContractError)
    @app.exception_handler(HandoffError)
    @app.exception_handler(ProtocolEvidenceError)
    async def contract_error(request, exc):
        return JSONResponse({"detail": "저장된 자료의 계약/무결성 검사에 실패했습니다."}, status_code=409)

    @app.get("/api/health")
    def health():
        return {"workbench_version":"20260930.design.1", "profile": "local-literature-review" if provider_factory else "local-single-user-i2",
                "llm_enabled": provider_factory is not None, "gpu_enabled": False,
                "server_pid":os.getpid(),
                "read_only": read_only, "review_only": review_only,
                "public_deployment_ready": False, "contracts": {"boundary": "0.1.0", "payload": "0.1.0", "review": "0.1.0"}}

    @app.post("/api/uploads")
    async def upload(request: Request):
        # Stream admission cap applies before multipart parsing/spooling.
        chunks, size = [], 0
        async for chunk in request.stream():
            size += len(chunk)
            if size > 6_100_000:
                raise HTTPException(413, "PDF 파일은 한 개당 6 MB 이하로 올려 주세요.")
            chunks.append(chunk)
        request._body = b"".join(chunks)
        form = await request.form(max_files=1, max_fields=1, max_part_size=6_000_000)
        file = form.get("file")
        if file is None or not hasattr(file, "read"):
            raise HTTPException(422, "PDF 파일을 선택해 주세요.")
        raw = await file.read(6_000_001)
        await file.close()
        if len(raw) > 6_000_000 or not raw.startswith(b"%PDF-"):
            raise HTTPException(422, "PDF 형식과 크기를 확인해 주세요.")
        reference = store.scope(PROJECT).put_raw(raw, "application/pdf", "source")
        return {"artifact": reference, "filename": Path(file.filename or "attachment.pdf").name, "bytes": len(raw)}

    @app.post("/api/runs", status_code=202)
    def create_run(body: RunInput):
        values = body.model_dump()
        request_key = values.pop("request_key")
        if len(values["pdf_artifact_ids"]) >= values["max_documents"]:
            raise HTTPException(422, "기본 문헌 검색을 위해 첨부 수보다 큰 문서 한도를 선택해 주세요.")
        if len(values["pdf_artifact_ids"]) != len(set(values["pdf_artifact_ids"])):
            raise HTTPException(422, "같은 첨부를 중복 선택할 수 없습니다.")
        upload_bytes = 0
        for identifier in values["pdf_artifact_ids"]:
            reference = store.reference(PROJECT, identifier)
            if reference["media_type"] != "application/pdf" or reference["provenance"] == "test_fixture":
                raise HTTPException(422, "실제 PDF 첨부 자료를 선택해 주세요.")
            upload_bytes += len(store.read(PROJECT, reference))
        if upload_bytes > 10_000_000:
            raise HTTPException(413, "검색 자료를 받을 공간을 남기기 위해 PDF 합계는 10 MB 이하로 제한합니다.")
        try:
            run, fresh = store.create_run(PROJECT, request_key, values)
        except ValueError as exc:
            raise HTTPException(409, str(exc))
        if fresh and enable_worker:
            executor.submit(orchestrator.execute, run["id"])
        return {"run_id": run["id"], "created": fresh}

    @app.post("/api/runs/{run_id}/select", status_code=202)
    def choose(run_id: str, body: SelectionInput):
        try:
            run = orchestrator.choose(run_id, body.candidate_id, body.revision)
        except ValueError as exc:
            raise HTTPException(409, str(exc))
        if enable_worker:
            executor.submit(orchestrator.execute, run["id"])
        return {"run_id": run["id"]}

    @app.get("/api/runs")
    def runs():
        return [{"id": r["id"], "query": r["request"]["query"], "state": r["state"], "updated_at": r["updated_at"]} for r in store.list_runs(PROJECT)]

    @app.get("/api/runs/{run_id}")
    def run(run_id: str):
        value = orchestrator.view(run_id)
        value["analyses"] = analysis.list(run_id)
        value["llm_enabled"] = provider_factory is not None and not value.get('replay_only')
        value['workflows'] = workflows.list(run_id)
        value["g1"] = analysis.g1_status(run_id)
        return value

    @app.post("/api/runs/{run_id}/analyses", status_code=202)
    def analyze(run_id: str, body: AnalysisInput):
        try:
            job, fresh = analysis.create(run_id, body.revision, body.segment_ids, body.request_key)
        except ValueError as exc:
            raise HTTPException(409, str(exc))
        if fresh and enable_worker:
            executor.submit(analysis.execute, job["id"])
        return {"analysis_id": job["id"], "created": fresh}

    @app.get("/api/analyses/{analysis_id}")
    def analysis_result(analysis_id: str):
        return analysis.view(analysis_id)

    @app.get("/api/runs/{run_id}/g1")
    def gate_readiness(run_id: str):
        return analysis.g1_status(run_id)

    @app.get("/api/runs/{run_id}/dossiers")
    def candidate_dossiers(run_id: str):
        return dossiers.list(run_id)

    @app.get("/api/dossiers/{dossier_id}")
    def candidate_dossier(dossier_id: str):
        return dossiers.view(dossier_id)

    @app.get('/api/dossiers/{dossier_id}/b-inputs')
    def b_inputs(dossier_id: str):
        return workflows.handoffs.list(dossier_id)

    @app.get('/api/dossiers/{dossier_id}/saved-results')
    def saved_result_history(dossier_id: str):
        return saved_results.list(dossier_id)

    @app.get('/api/saved-results/{identifier}')
    def saved_result(identifier: str):
        return saved_results.view(identifier)

    @app.get('/api/saved-results/{identifier}/files')
    def saved_files(identifier: str):
        return saved_results.files(identifier)

    @app.get('/api/saved-results/{identifier}/report-preparation')
    def prepare_evidence_report(identifier: str):
        return evidence_reports.prepare(identifier)

    @app.get('/api/saved-results/{identifier}/reports')
    def evidence_report_history(identifier: str):
        return evidence_reports.list(identifier)

    @app.post('/api/saved-results/{identifier}/reports')
    def create_evidence_report(identifier: str, body: EvidenceReportInput):
        report, fresh = evidence_reports.create(identifier,body.expected_digest)
        return {'report_id':report['id'],'created':fresh,'freshness':report['freshness']}

    @app.get('/api/evidence-reports/{identifier}')
    def evidence_report(identifier: str):
        return evidence_reports.view(identifier)

    @app.get('/api/review-auth/session')
    def review_session(request: Request):
        try: actor=report_reviews.identities.session(request.cookies.get('tpd_review_session'))
        except ReviewAuthError: actor=None
        return {'actor':actor,'reviewers':report_reviews.identities.actors()}

    @app.post('/api/review-auth/login')
    def review_login(body: ReviewLoginInput, request: Request):
        token=report_reviews.identities.login(body.reviewer_id,body.access_key)
        report_reviews.identities.logout(request.cookies.get('tpd_review_session'))
        response=JSONResponse({'actor':report_reviews.identities.session(token)})
        response.set_cookie('tpd_review_session',token,httponly=True,samesite='strict',secure=request.url.scheme=='https',
                            max_age=report_reviews.identities.SESSION_SECONDS,path='/api')
        return response

    @app.post('/api/review-auth/logout')
    def review_logout(request: Request):
        report_reviews.identities.logout(request.cookies.get('tpd_review_session'))
        response=JSONResponse({'logged_out':True});response.delete_cookie('tpd_review_session',path='/api')
        return response

    @app.get('/api/evidence-reports/{identifier}/review-requests')
    def report_review_history(identifier: str):
        return report_reviews.list(identifier)

    @app.get('/api/evidence-reports/{identifier}/review-predecessors')
    def report_review_predecessors(identifier: str):
        return report_reviews.predecessors(identifier)

    @app.post('/api/evidence-reports/{identifier}/review-requests')
    def create_report_review(identifier: str, body: HumanReviewRequestInput, request: Request):
        actor=report_reviews.identities.session(request.cookies.get('tpd_review_session'))
        row,fresh=report_reviews.create(identifier,body.expected_report_sha256,body.scopes,actor['id'],body.supersedes,
                                       request.cookies.get('tpd_review_session'))
        return {'request_id':row['id'],'created':fresh}

    @app.get('/api/report-reviews/{identifier}')
    def report_review(identifier: str):
        return report_reviews.view(identifier)

    @app.post('/api/report-reviews/{identifier}/decisions')
    def report_review_decision(identifier: str, body: HumanReviewDecisionInput, request: Request):
        return report_reviews.decide(identifier,body.model_dump(),request.cookies.get('tpd_review_session'))

    @app.get('/api/b-inputs/{input_id}/review-preparation')
    def review_preparation(input_id: str):
        return workflows.handoffs.packet(input_id)

    @app.post('/api/workflows',status_code=202)
    def start_workflow(body: WorkflowInput):
        try:
            row,fresh=workflows.create(body.b_input_id,body.request_key)
        except ValueError as exc:
            raise HTTPException(409,str(exc))
        if fresh and enable_worker: executor.submit(workflows.execute,row['id'])
        return {'workflow_id':row['id'],'created':fresh}

    @app.get('/api/workflows/{identifier}')
    def workflow(identifier: str):
        return workflows.view(identifier)

    @app.post('/api/workflows/{identifier}/cancel')
    def cancel_workflow(identifier: str):
        try:return workflows.cancel(identifier)
        except ValueError as exc:raise HTTPException(409,str(exc))

    @app.post('/api/workflows/{identifier}/resume',status_code=202)
    def resume_workflow(identifier: str):
        try:row=workflows.resume(identifier)
        except ValueError as exc:raise HTTPException(409,str(exc))
        if enable_worker:executor.submit(workflows.execute,identifier)
        return {'workflow_id':row['id']}

    @app.get('/api/runs/{run_id}/replay-package')
    def replay_package(run_id: str):
        raw,_=export_bytes(store,PROJECT,run_id)
        return Response(raw,media_type='application/zip',headers={'Content-Disposition':'attachment; filename="tpd-internal-replay.zip"'})

    @app.get('/api/lab/cases')
    def lab_cases():
        with store.db() as db:
            rows=db.execute('SELECT id,ref FROM saved_results WHERE project=? ORDER BY rowid DESC',(PROJECT,)).fetchall()
        result=[]
        for row in rows:
            import json
            record=store.scope(PROJECT).json(json.loads(row['ref']))
            data=store.scope(PROJECT).json(record['projection_ref'])
            run=store.get_run(PROJECT,record['run_id'])
            title=run['request']['query']
            if run.get('target_ref'):
                target=store.scope(PROJECT).json(run['target_ref'])
                selected=next((c for c in target.get('candidates',[]) if c['candidate_id']==target.get('selected_candidate_id')),None)
                if selected:
                    symbol=selected['identity'].get('gene_symbol',{}).get('value')
                    if symbol:title=symbol+' · '+selected['identity']['uniprot_accession']
            result.append({'id':row['id'],'dossier_id':record['dossier_id'],'run_id':record['run_id'],
                           'title':title,'summary':data['summary'],'replay_only':bool(run.get('replay_only'))})
        return result

    @app.get('/api/lab/linker-library')
    def lab_linker_library():
        from packages.science.linker_design import library
        return library()

    @app.get('/api/design/catalog')
    def design_catalog():
        from packages.science.dual_e3 import catalog
        return catalog()

    @app.get('/design')
    def design_page():
        return FileResponse(ROOT/'apps/web/design-lab.html')

    @app.get('/api/lab/cases/{identifier}')
    def lab_case(identifier:str):
        data=saved_results.view(identifier)
        data['reports']=evidence_reports.list(identifier)
        data['replay_only']=bool(store.get_run(PROJECT,data['run_id']).get('replay_only'))
        return data

    @app.get('/api/lab/jobs')
    def lab_jobs(result_id:str|None=None):
        return workbench.list(result_id)

    @app.get('/api/lab/cases/{identifier}/replay')
    def lab_replay(identifier:str):
        from packages.platform.workbench_replay import export_bytes as export_workbench
        raw,_=export_workbench(store,PROJECT,saved_results.view(identifier)['run_id'])
        return Response(raw,media_type='application/zip',headers={'Content-Disposition':'attachment; filename="tpd-workbench-replay.zip"'})

    @app.get('/api/lab/jobs/{identifier}')
    def lab_job(identifier:str):
        return workbench.view(identifier)

    def scientific_session(request:Request):
        return report_reviews.identities.session(request.cookies.get('tpd_review_session'))
    def clean_scientific(value):
        blocked={'access_key','secret_key','password','traceback','raw_traceback','local_path','filesystem_path'}
        if isinstance(value,dict): return {k:clean_scientific(v) for k,v in value.items() if k not in blocked and not k.endswith('_local_path')}
        if isinstance(value,(list,tuple)): return [clean_scientific(v) for v in value]
        return value
    def created(value): return value if isinstance(value,tuple) else (value,True)

    @app.get('/api/lab/jobs/{identifier}/protocol-diagnostics')
    def list_protocol_diagnostics(identifier:str):
        return clean_scientific(protocol_evidence.list(identifier))
    @app.get('/api/protocol-diagnostics/{identifier}')
    def view_protocol_diagnostic(identifier:str):
        return clean_scientific(protocol_evidence.view(identifier))

    @app.get('/api/lab/jobs/{identifier}/scientific-assessments')
    def list_scientific(identifier:str,compact:bool=False):
        values=(scientific_acceptance.list(identifier,compact=True) if compact
                else scientific_acceptance.list(identifier))
        return clean_scientific(values)
    @app.get('/api/lab/jobs/{identifier}/scientific-policy')
    def get_scientific_policy(identifier:str): return clean_scientific(scientific_acceptance.policy(identifier))
    @app.get('/api/lab/jobs/{identifier}/microstate-population-sources')
    def list_microstate_population_sources(identifier:str,request:Request):
        scientific_session(request)
        return clean_scientific(scientific_acceptance.list_quantitative_sources(
            identifier,request.cookies.get('tpd_review_session')))
    @app.post('/api/lab/jobs/{identifier}/microstate-population-sources')
    def register_microstate_population_source(identifier:str,body:MicrostatePopulationSourceInput,request:Request):
        scientific_session(request)
        return clean_scientific(scientific_acceptance.register_quantitative_source(
            identifier,body.text,request.cookies.get('tpd_review_session')))
    @app.post('/api/lab/jobs/{identifier}/scientific-assessments')
    def create_scientific(identifier:str,body:ScientificAssessmentInput,request:Request):
        actor=scientific_session(request)
        create_options={'reviewer_id':actor['id']}
        if body.protocol_ids:
            create_options['protocol_ids']=body.protocol_ids
        row,is_new=created(scientific_acceptance.create(
            identifier,body.expected_result_sha256,**create_options))
        return {'assessment_id':row['id'],'created':is_new}
    @app.get('/api/scientific-assessments/{identifier}')
    def view_scientific(identifier:str): return clean_scientific(scientific_acceptance.view(identifier))
    @app.post('/api/scientific-statements')
    def create_scientific_statement(body:ScientificStatementInput,request:Request):
        scientific_session(request)
        return clean_scientific(scientific_acceptance.create_statement(body.text,body.locator,request.cookies.get('tpd_review_session')))
    @app.post('/api/scientific-assessments/{identifier}/decisions')
    def record_scientific(identifier:str,body:ScientificDecisionInput,request:Request):
        scientific_session(request)
        row,is_new=created(scientific_acceptance.record_decision(identifier,body.model_dump(),request.cookies.get('tpd_review_session')))
        return {'decision':clean_scientific(row),'created':is_new}
    @app.post('/api/scientific-decisions/{identifier}/withdraw')
    def withdraw_scientific(identifier:str,request:Request):
        scientific_session(request)
        return clean_scientific(scientific_acceptance.withdraw_decision(identifier,request.cookies.get('tpd_review_session')))
    @app.post('/api/scientific-assessments/{identifier}/compute')
    def compute_scientific(identifier:str,body:ScientificComputeInput,request:Request):
        scientific_session(request)
        assessment=scientific_acceptance.view(identifier)
        return clean_scientific(scientific_acceptance.compute(assessment['job_id'],request.cookies.get('tpd_review_session'),body.candidate_ids,body.analog_ids))
    @app.post('/api/scientific-assessments/{identifier}/export')
    def export_scientific(identifier:str,request:Request):
        scientific_session(request)
        return clean_scientific(scientific_acceptance.export(identifier,request.cookies.get('tpd_review_session')))

    @app.post('/api/lab/jobs',status_code=202)
    def lab_start(body:WorkbenchInput):
        job,fresh=workbench.create(body.result_id,body.operation,body.request_key,body.parameters,body.retry_of)
        if fresh and enable_worker:executor.submit(workbench.execute,job['id'])
        return {'job_id':job['id'],'created':fresh}

    @app.post('/api/lab/jobs/{identifier}/cancel')
    def lab_cancel(identifier:str):
        return workbench.cancel(identifier)

    @app.get('/api/lab/structures/{identifier}/{candidate_id}/{sample_index}')
    def lab_structure(identifier:str,candidate_id:str,sample_index:int):
        from packages.platform.structure_views import structure_view
        return structure_view(saved_results,identifier,candidate_id,sample_index)

    @app.get('/api/lab/artifacts/{identifier}/preview')
    def lab_preview(identifier:str):
        ref=store.reference(PROJECT,identifier)
        if ref['media_type'] not in {'image/png','text/markdown','text/plain','application/json'}:
            raise HTTPException(422,'미리보기를 지원하지 않는 자료입니다.')
        return Response(store.read(PROJECT,ref),media_type=ref['media_type'])

    @app.get("/api/artifacts/{artifact_id}/text")
    def text(artifact_id: str):
        reference = store.reference(PROJECT, artifact_id)
        if reference["media_type"] != "text/plain":
            raise HTTPException(422, "텍스트 구간 자료가 아닙니다.")
        return {"text": store.read(PROJECT, reference).decode("utf-8"), "artifact": reference}

    @app.get("/api/artifacts/{artifact_id}/download")
    def download(artifact_id: str):
        reference = store.reference(PROJECT, artifact_id)
        suffix = {"application/pdf": ".pdf", "application/xml": ".xml", "application/json": ".json", "text/plain": ".txt", "text/markdown": ".md",
                  "image/png": ".png", "chemical/x-mdl-sdfile": ".sdf", "chemical/x-mmcif": ".cif",
                  "application/vnd.openxmlformats-officedocument.wordprocessingml.document": ".docx", "application/zip": ".zip"}.get(reference["media_type"], ".bin")
        return Response(store.read(PROJECT, reference), media_type="application/octet-stream",
            headers={"Content-Disposition": 'attachment; filename="' + artifact_id + suffix + '"'})

    app.mount("/static", StaticFiles(directory=ROOT / "apps/web"), name="static")

    @app.get("/")
    def index():
        return FileResponse(ROOT / 'apps/web/workbench.html')

    @app.get('/literature')
    def literature_page():
        return FileResponse(ROOT / "apps/web/index.html")

    @app.get('/saved-review/{dossier_id}')
    def saved_review_page(dossier_id: str):
        return FileResponse(ROOT / 'apps/web/saved-review.html')

    @app.get('/evidence-reports/{result_id}')
    def evidence_report_page(result_id: str):
        return FileResponse(ROOT / 'apps/web/evidence-report.html')

    @app.get('/human-review/{report_id}')
    def human_review_page(report_id: str):
        return FileResponse(ROOT / 'apps/web/human-review.html')

    @app.get('/scientific-review/{job_id}')
    def scientific_review_page(job_id:str):
        return FileResponse(ROOT / 'apps/web/scientific-review.html')

    return app
