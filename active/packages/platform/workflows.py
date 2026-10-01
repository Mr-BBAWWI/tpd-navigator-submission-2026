"""M2-owned bounded review workflow with durable call receipts and explicit resume."""
import json
import time
import uuid
from pathlib import Path

from packages.contracts import encoded, now
from packages.agents.provider import AgentError, Completion
from packages.agents.runtime import ReviewRuntime, PROMPT_VERSION, CONTRACT
from packages.agents.workflow_bundle import build_service_bundle, VERSION as BUNDLE_VERSION
from packages.platform.analysis import worker_lock
from packages.platform.dossiers import digest
from packages.platform.handoffs import HandoffService

VERSION = 'review-workflow-20260924.2'
LIMITS = {'max_calls':14,'max_total_tokens':180000,'max_output_tokens':4000,'max_elapsed_seconds':600}


def runtime_identity():
    root=Path(__file__).resolve().parents[1]
    code={name:digest((root/name).read_bytes()) for name in ('platform/workflows.py','platform/handoffs.py','agents/workflow_bundle.py','agents/runtime.py')}
    return {'service':VERSION,'bundle':BUNDLE_VERSION,'prompt':PROMPT_VERSION,'schema':digest(encoded(CONTRACT)),'code':code}


class JournalProvider:
    """Only fully saved answers can be replayed; uncertain requests never auto-retry."""
    def __init__(self, service, identifier, provider):
        self.service, self.identifier, self.provider = service, identifier, provider
        self.mode, self.index = provider.mode, 0

    def complete(self, **request):
        s, identifier, index = self.service, self.identifier, self.index
        self.index += 1
        s.check_active(identifier)
        fingerprint = digest(encoded(request))
        with s.store.db() as db:
            previous = db.execute('SELECT body FROM workflow_calls WHERE workflow_id=? AND ordinal=?',(identifier,index)).fetchone()
        if previous:
            call = json.loads(previous[0])
            if call['fingerprint'] != fingerprint:
                raise AgentError('WORKFLOW_REPLAY_PROMPT_CHANGED')
            if not call.get('answer_ref') or call['usage_status']!='observed':
                raise AgentError('WORKFLOW_UNCERTAIN_CALL_NO_RETRY')
            answer = s.port.json(call['answer_ref'])
            return Completion(answer['text'], answer['metadata'])
        prompt_ref = s.port.put_json(request)
        call = {'ordinal':index,'fingerprint':fingerprint,'prompt_ref':prompt_ref,'answer_ref':None,
            'usage_status':'unreconciled','tokens':0,'requested_at':now()}
        with s.store.db() as db:
            db.execute('BEGIN IMMEDIATE')
            row = s._row(db,identifier)
            if row['state']!='running':
                raise AgentError('WORKFLOW_NOT_RUNNING')
            db.execute('INSERT INTO workflow_calls VALUES(?,?,?)',(identifier,index,json.dumps(call)))
            row.update(calls=row['calls']+1,usage_status='unreconciled',stage='model_call_'+str(index+1),
                       active_seconds=round(s._base_seconds+time.monotonic()-s._started,3))
            s._save(db,row)
        # The request is journaled before any network call, closing the duplicate billing window.
        answer = self.provider.complete(**request)
        answer_ref = s.port.put_json({'text':answer.text,'metadata':answer.metadata})
        tokens = answer.metadata.get('usage',{}).get('total_tokens')
        observed = type(tokens) is int and tokens>0
        call.update(answer_ref=answer_ref,usage_status='observed' if observed else 'unreconciled',
                    tokens=tokens if observed else 0,answered_at=now())
        with s.store.db() as db:
            db.execute('BEGIN IMMEDIATE')
            db.execute('UPDATE workflow_calls SET body=? WHERE workflow_id=? AND ordinal=?',(json.dumps(call),identifier,index))
            row = s._row(db,identifier)
            row.update(total_tokens=row['total_tokens']+call['tokens'],usage_status=call['usage_status'],
                       active_seconds=round(s._base_seconds+time.monotonic()-s._started,3))
            s._save(db,row)
        s.check_active(identifier)
        return answer


class WorkflowService:
    def __init__(self, store, project, provider_factory=None):
        self.store, self.project, self.provider_factory = store, project, provider_factory
        self.port = store.scope(project)
        self.handoffs = HandoffService(store, project)
        with store.db() as db:
            db.executescript('''CREATE TABLE IF NOT EXISTS workflows(
                id TEXT PRIMARY KEY, project TEXT NOT NULL, run_id TEXT NOT NULL,
                request_key TEXT NOT NULL, fingerprint TEXT NOT NULL, state TEXT NOT NULL, body TEXT NOT NULL,
                UNIQUE(project,run_id,request_key));
                CREATE TABLE IF NOT EXISTS workflow_calls(
                workflow_id TEXT NOT NULL, ordinal INTEGER NOT NULL, body TEXT NOT NULL,
                PRIMARY KEY(workflow_id,ordinal));''')

    def _row(self, db, identifier):
        row = db.execute('SELECT body FROM workflows WHERE project=? AND id=?',(self.project,identifier)).fetchone()
        if row is None:
            raise KeyError(identifier)
        return json.loads(row[0])

    def _save(self, db, row):
        row['updated_at'] = now()
        db.execute('UPDATE workflows SET state=?,body=? WHERE project=? AND id=?',(row['state'],json.dumps(row),self.project,row['id']))

    def get(self, identifier):
        with self.store.db() as db:
            return self._row(db,identifier)

    def create(self, input_id, request_key):
        if self.provider_factory is None:
            raise ValueError('AI 검토가 꺼져 있습니다.')
        value = self.handoffs.view(input_id)
        run = self.store.get_run(self.project,value['run_id'])
        if run.get('replay_only') or not value['current_input']:
            raise ValueError('현재 연구 입력에서만 새 검토를 실행할 수 있습니다.')
        bundle = build_service_bundle(self.store,self.project,input_id)
        fingerprint = digest(encoded({'input':bundle['input_digest'],'runtime':runtime_identity(),'limits':LIMITS}))
        bundle_ref = self.port.put_json(bundle)
        with self.store.db() as db:
            db.execute('BEGIN IMMEDIATE')
            old = db.execute('SELECT fingerprint,body FROM workflows WHERE project=? AND run_id=? AND request_key=?',
                (self.project,value['run_id'],request_key)).fetchone()
            if old:
                if old[0]!=fingerprint:
                    raise ValueError('같은 요청 키에 다른 검토 입력을 사용할 수 없습니다.')
                return json.loads(old[1]),False
            if db.execute("SELECT 1 FROM workflows WHERE project=? AND state IN ('queued','running')",(self.project,)).fetchone():
                raise ValueError('진행 중인 후보 검토를 먼저 완료해 주세요.')
            current = json.loads(db.execute('SELECT body FROM runs WHERE id=? AND project=?',(value['run_id'],self.project)).fetchone()[0])
            if self.handoffs.dossiers.snapshot(current)!=value['snapshot']:
                raise ValueError('입력이 변경됐습니다.')
            row = {'id':'workflow-'+uuid.uuid4().hex,'project_id':self.project,'run_id':value['run_id'],
                'b_input_id':input_id,'bundle_ref':bundle_ref,'input_digest':bundle['input_digest'],
                'runtime':runtime_identity(),'limits':dict(LIMITS),'model':'gpt-5.6-sol',
                'state':'queued','stage':'queued','created_at':now(),'updated_at':now(),'attempts':0,
                'calls':0,'total_tokens':0,'usage_status':'not_started','active_seconds':0.0,
                'result_ref':None,'report_ref':None,'packet_ref':None,'error_code':None,
                'human_approved':False,'dispatch_authorized':False}
            db.execute('INSERT INTO workflows VALUES(?,?,?,?,?,?,?)',(row['id'],self.project,row['run_id'],request_key,fingerprint,row['state'],json.dumps(row)))
        return row,True

    def check_active(self, identifier):
        row = self.get(identifier)
        if row['state']!='running':
            raise AgentError('WORKFLOW_NOT_RUNNING')
        if row['runtime']!=runtime_identity():
            raise AgentError('WORKFLOW_RUNTIME_CHANGED')
        if not self.handoffs.view(row['b_input_id'])['current_input']:
            raise AgentError('WORKFLOW_INPUT_STALE')
        if time.monotonic()-self._started+self._base_seconds > row['limits']['max_elapsed_seconds']:
            raise AgentError('WORKFLOW_TIME_LIMIT')
        return row['input_digest']

    def execute(self, identifier):
        # Shared OS lock serializes this workflow with literature analysis across processes.
        with worker_lock(self.store.root) as acquired:
            if not acquired:
                with self.store.db() as db:
                    db.execute('BEGIN IMMEDIATE');row=self._row(db,identifier)
                    if row['state']=='queued':
                        row.update(state='interrupted',stage='worker_busy',error_code='WORKFLOW_WORKER_BUSY');self._save(db,row)
                return
            with self.store.db() as db:
                db.execute('BEGIN IMMEDIATE'); row=self._row(db,identifier)
                if row['state']!='queued':
                    return
                row.update(state='running',attempts=row['attempts']+1); self._save(db,row)
            self._started=time.monotonic(); self._base_seconds=row['active_seconds']; provider=None
            try:
                self.check_active(identifier)
                provider=self.provider_factory()
                runtime=ReviewRuntime(JournalProvider(self,identifier,provider),
                    self.store.root/'workflow-attempts'/identifier/str(row['attempts']),model=row['model'],
                    **{k:v for k,v in row['limits'].items() if k!='max_elapsed_seconds'})
                bundle=self.port.json(row['bundle_ref'])
                result=runtime.run(bundle,revalidate=lambda:self.check_active(identifier))
                result_ref=self.port.put_json(result)
                report_path=runtime.output/'REVIEW.md'
                report_ref=self.port.put_raw(report_path.read_bytes(),'text/markdown','computed') if report_path.exists() else None
                # Keep failed attempts as artifacts too; journal call evidence remains accessible.
                trace_ref=self.port.put_raw((runtime.output/'trace.jsonl').read_bytes(),'text/plain','computed') if (runtime.output/'trace.jsonl').exists() else None
                with self.store.db() as db:
                    db.execute('BEGIN IMMEDIATE'); current=self._row(db,identifier)
                    state='review_ready' if result['status']=='draft_pending_human_review' else 'held'
                    if current['state']=='cancelled': state='cancelled'
                    run=json.loads(db.execute('SELECT body FROM runs WHERE project=? AND id=?',(self.project,row['run_id'])).fetchone()[0])
                    frozen=self.port.json(bundle['b_input_ref'])['snapshot']
                    if self.handoffs.dossiers.snapshot(run)!=frozen: state='stale'
                    current.update(stage='preparing_review_material',result_ref=result_ref,report_ref=report_ref,trace_ref=trace_ref,
                        error_code=result.get('error'),active_seconds=round(self._base_seconds+time.monotonic()-self._started,3))
                    self._save(db,current)
                # Keep the job active until the final packet is saved so the UI and exporter
                # cannot observe a completed workflow whose review material is still missing.
                packet=self.handoffs.packet(row['b_input_id'],{**self.summary(identifier),'state':state})
                packet_ref=self.port.put_json(packet)
                with self.store.db() as db:
                    db.execute('BEGIN IMMEDIATE'); current=self._row(db,identifier)
                    if current['state']=='cancelled': state='cancelled'
                    run=json.loads(db.execute('SELECT body FROM runs WHERE project=? AND id=?',(self.project,row['run_id'])).fetchone()[0])
                    if self.handoffs.dossiers.snapshot(run)!=frozen: state='stale'
                    current.update(state=state,stage='finished',packet_ref=packet_ref,
                                   active_seconds=round(self._base_seconds+time.monotonic()-self._started,3))
                    self._save(db,current)
            except Exception as exc:
                code=str(exc) if isinstance(exc,AgentError) else 'WORKFLOW_EXECUTION_FAILED'
                if not code.replace('_','').isalnum() or len(code)>100: code='WORKFLOW_EXECUTION_FAILED'
                with self.store.db() as db:
                    db.execute('BEGIN IMMEDIATE'); current=self._row(db,identifier)
                    if current['state']=='running':
                        current.update(state='held',stage='finished',error_code=code,active_seconds=round(self._base_seconds+time.monotonic()-self._started,3)); self._save(db,current)
            finally:
                if provider: provider.close()

    def summary(self, identifier):
        row=self.get(identifier)
        return {k:row[k] for k in ('id','b_input_id','state','calls','total_tokens','usage_status','result_ref','report_ref','input_digest','human_approved','dispatch_authorized')}

    def calls(self, identifier):
        self.get(identifier)
        with self.store.db() as db:
            return [json.loads(r[0]) for r in db.execute('SELECT body FROM workflow_calls WHERE workflow_id=? ORDER BY ordinal',(identifier,))]

    def view(self, identifier):
        row=self.get(identifier)
        row['current_input']=self.handoffs.view(row['b_input_id'])['current_input']
        row['calls_detail']=self.calls(identifier)
        for call in row['calls_detail']:
            self.port.read(call['prompt_ref'])
            if call['answer_ref']: self.port.read(call['answer_ref'])
        for name in ('bundle_ref','result_ref','report_ref','packet_ref'):
            if row.get(name): self.port.read(row[name])
        if row['result_ref']: row['result']=self.port.json(row['result_ref'])
        replay=self.store.get_run(self.project,row['run_id']).get('replay_only',False)
        row['can_resume']=not replay and row['state']=='interrupted' and row['current_input'] and row['runtime']==runtime_identity() and row['usage_status']!='unreconciled'
        return row

    def list(self, run_id):
        self.store.get_run(self.project,run_id)
        with self.store.db() as db:
            ids=[r[0] for r in db.execute('SELECT id FROM workflows WHERE project=? AND run_id=? ORDER BY rowid DESC',(self.project,run_id))]
        return [self.view(i) for i in ids]

    def cancel(self, identifier):
        with self.store.db() as db:
            db.execute('BEGIN IMMEDIATE'); row=self._row(db,identifier)
            run=json.loads(db.execute('SELECT body FROM runs WHERE project=? AND id=?',(self.project,row['run_id'])).fetchone()[0])
            if run.get('replay_only'):raise ValueError('가져온 재생 자료는 변경할 수 없습니다.')
            if row['state'] not in {'queued','running','interrupted'}:
                raise ValueError('현재 상태에서는 중단할 수 없습니다.')
            row.update(state='cancelled',stage='cancelled'); self._save(db,row)
        return row

    def resume(self, identifier):
        if self.provider_factory is None or not self.view(identifier)['can_resume']:
            raise ValueError('현재 입력과 호출 기록을 확인할 수 있는 중단 작업만 재개할 수 있습니다.')
        with self.store.db() as db:
            db.execute('BEGIN IMMEDIATE'); row=self._row(db,identifier)
            if row['state']!='interrupted' or db.execute("SELECT 1 FROM workflows WHERE project=? AND state IN ('queued','running')",(self.project,)).fetchone():
                raise ValueError('작업 상태가 바뀌었거나 진행 중인 검토가 있습니다.')
            row.update(state='queued',stage='resume_requested'); self._save(db,row)
        return row

    def interrupt_pending(self):
        with worker_lock(self.store.root) as acquired:
            if not acquired:return
            with self.store.db() as db:
                db.execute('BEGIN IMMEDIATE')
                ids=[r[0] for r in db.execute("SELECT id FROM workflows WHERE project=? AND state IN ('queued','running')",(self.project,))]
                for identifier in ids:
                    row=self._row(db,identifier);row.update(state='interrupted',stage='interrupted');self._save(db,row)
