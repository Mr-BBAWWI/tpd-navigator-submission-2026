"""Durable local tool/agent jobs. No arbitrary command execution or implicit approval."""
import importlib.metadata
import json
import time
import uuid
from pathlib import Path

from packages.contracts import encoded, now, ContractError
from packages.agents.provider import AgentError
from packages.agents.runtime import ReviewRuntime, PROMPT_VERSION
from packages.agents.saved_bundle import build_saved_bundle
from packages.platform.analysis import worker_lock
from packages.platform.dossiers import digest, require
from packages.platform.evidence_reports import EvidenceReportService
from packages.platform.saved_results import AUTHORITY


VERSION='workbench-jobs/20260930.design.4'
TOOLS={'verify':'저장 자료 검사','report':'통합 보고서 작성','review':'AI 근거 검토',
       'candidate':'C03 설계 가설 생성','chemistry':'리간드 수소 방향 시험',
       'design_panel':'Warhead SAR · CRBN/VHL 설계 Agent'}
LIMITS={'max_calls':14,'max_total_tokens':180000,'max_output_tokens':4000,'max_elapsed_seconds':900}
ROOT=Path(__file__).resolve().parents[2]


def identity():
    paths=['packages/platform/workbench.py','packages/agents/saved_bundle.py','packages/agents/runtime.py',
           'contracts/drafts/agent_review.schema.json','packages/platform/design_panel.py',
           'packages/science/warhead_sites.py','packages/science/analog_generation.py',
           'packages/science/analog_filters.py','packages/science/mapped_stereo.py',
           'packages/science/design_docking.py',
           'packages/science/dual_e3.py','packages/science/e3_benchmark.py',
           'packages/science/chemical_states.py','packages/science/linker_design.py',
           'packages/science/linker_assessment.py','packages/science/acceptance_evidence.py',
           'packages/science/reference_parents.py','packages/science/medchem_sar.py',
           'packages/science/molecules.py','packages/science/structures.py',
           'contracts/drafts/design_panel.schema.json','cases/design_sources/manifest.json']
    for name in ['packages/platform/cpu_tools.py','cases/linker_library_20260930.json']:
        if (ROOT/name).is_file():paths.append(name)
    dependencies={}
    for package in ('rdkit','meeko','gemmi','numpy','scipy'):
        try:dependencies[package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:dependencies[package]='not_installed'
    from packages.platform.design_panel import input_binding
    return {'version':VERSION,'prompt':PROMPT_VERSION,
            'files':{p:digest((ROOT/p).read_bytes()) for p in paths if (ROOT/p).is_file()},
            'dependencies':dependencies,'data_sources':input_binding()['source_files']}


def _observed_usage(metadata):
    metadata = metadata if type(metadata) is dict else {}
    usage = metadata.get('usage')
    if type(usage) is not dict:
        return None, 'unreconciled', 0
    total = usage.get('total_tokens')
    if type(total) is int and total >= 0:
        return usage, 'observed', total
    return usage, 'observed_partial', 0


class CallJournal:
    def __init__(self, service, job_id, provider):
        self.service,self.job_id,self.provider=service,job_id,provider
        self.mode=provider.mode

    def _record_answer(self, call, metadata, text):
        s=self.service
        usage,status,total=_observed_usage(metadata)
        ref=s.port.put_json({'text':text,'metadata':metadata})
        with s.store.db() as db:
            db.execute('BEGIN IMMEDIATE');row=s._row(db,self.job_id)
            update={'answer_ref':ref,'usage_status':status,'completed_at':now()}
            if usage is not None:update['usage']=usage
            call.update(update)
            db.execute('UPDATE workbench_calls SET body=? WHERE job_id=? AND ordinal=?',
                       (json.dumps(call),row['id'],call['ordinal']))
            row.update(total_tokens=row['total_tokens']+total,usage_status=status)
            s._save(db,row)

    def complete(self, **request):
        s=self.service;s.check_active(self.job_id)
        prompt=s.port.put_json(request)
        with s.store.db() as db:
            db.execute('BEGIN IMMEDIATE');row=s._row(db,self.job_id)
            require(row['state']=='running','JOB_NOT_RUNNING')
            ordinal=row['calls']+1
            require(ordinal<=row['limits']['max_calls'],'JOB_CALL_LIMIT')
            call={'ordinal':ordinal,'requested_at':now(),'prompt_ref':prompt,'answer_ref':None,'usage_status':'unreconciled'}
            db.execute('INSERT INTO workbench_calls VALUES(?,?,?)',(row['id'],ordinal,json.dumps(call)))
            row.update(calls=ordinal,usage_status='unreconciled',stage='agent_call_'+str(ordinal));s._save(db,row)
        # Cancellation after this point cannot recall the HTTP request; preserve any response or failure usage.
        try:
            answer=self.provider.complete(**request)
        except AgentError as error:
            metadata=getattr(error,'metadata',None)
            metadata=metadata if type(metadata) is dict else {}
            self._record_answer(call,metadata,'')
            raise
        metadata=answer.metadata if type(answer.metadata) is dict else {}
        self._record_answer(call,metadata,answer.text)
        s.check_active(self.job_id)
        return answer


class WorkbenchService:
    def __init__(self,store,project,provider_factory=None,scientific_acceptance_service=None):
        self.store,self.project,self.port=store,project,store.scope(project)
        self.provider_factory=provider_factory
        self.scientific_acceptance=scientific_acceptance_service
        self.reports=EvidenceReportService(store,project)
        self.saved=self.reports.saved
        with store.db() as db:
            db.executescript('''CREATE TABLE IF NOT EXISTS workbench_jobs(
                id TEXT PRIMARY KEY,project TEXT NOT NULL,run_id TEXT NOT NULL,result_id TEXT NOT NULL,
                request_key TEXT NOT NULL,fingerprint TEXT NOT NULL,state TEXT NOT NULL,body TEXT NOT NULL,
                UNIQUE(project,request_key));
                CREATE TABLE IF NOT EXISTS workbench_calls(job_id TEXT NOT NULL,ordinal INTEGER NOT NULL,
                body TEXT NOT NULL,PRIMARY KEY(job_id,ordinal));''')

    def _row(self,db,identifier):
        r=db.execute('SELECT body FROM workbench_jobs WHERE id=? AND project=?',(identifier,self.project)).fetchone()
        if r is None:raise KeyError(identifier)
        return json.loads(r[0])

    def _save(self,db,row):
        row['updated_at']=now()
        db.execute('UPDATE workbench_jobs SET state=?,body=? WHERE id=? AND project=?',
                   (row['state'],json.dumps(row),row['id'],self.project))

    def get(self,identifier):
        with self.store.db() as db:return self._row(db,identifier)

    @staticmethod
    def _public_row(row):
        result=dict(row)
        result.pop('_scientific_policy_snapshot',None)
        return result

    def list(self,result_id=None):
        with self.store.db() as db:
            if result_id:rows=db.execute('SELECT body FROM workbench_jobs WHERE project=? AND result_id=? ORDER BY rowid DESC',(self.project,result_id))
            else:rows=db.execute('SELECT body FROM workbench_jobs WHERE project=? ORDER BY rowid DESC LIMIT 100',(self.project,))
            return [self._public_row(json.loads(r[0])) for r in rows]

    def _policy_service(self):
        if self.scientific_acceptance is None:
            from packages.platform.scientific_acceptance import ScientificAcceptanceService
            self.scientific_acceptance=ScientificAcceptanceService(self.store,self.project)
        return self.scientific_acceptance

    def _resolve_design_policy(self,policy_id,parent_id,runtime=False):
        from packages.science.dual_e3 import validate_scientific_policy
        try:
            policy=self._policy_service().policy_for_design(policy_id)
            if type(policy) is not dict or policy.get('_validated_by_platform') is not True:
                raise ValueError('SCIENTIFIC_POLICY_NOT_PLATFORM_VALIDATED')
            policy=dict(policy)
            policy.pop('_validated_by_platform')
            return validate_scientific_policy(policy,project_id=self.project,
                policy_id=policy_id,parent_id=parent_id)
        except Exception:
            if runtime:raise AgentError('JOB_SCIENTIFIC_POLICY_INVALID')
            require(False,'JOB_SCIENTIFIC_POLICY_INVALID')

    def tool_snapshot(self,result_id):
        # Bounded selection: latest result for each explicit experiment setting.
        selected={}
        for row in self.list(result_id):
            if row['state']!='completed' or row['operation'] not in {'candidate','chemistry','design_panel'}:continue
            key=(row['operation'],json.dumps(row['parameters'],sort_keys=True))
            if key not in selected:selected[key]={'job_id':row['id'],'operation':row['operation'],
                'input_digest':row['input_digest'],'parameters':row['parameters'],'outputs':row['outputs'],'result':row['result']}
        return list(selected.values())[:32]

    def _prepare(self,result_id):
        if result_id=='design:SMARCA2':
            from packages.platform.design_panel import input_binding
            binding=input_binding()
            return {'can_create':True,'input_digest':binding['digest'],'binding':binding,
                    'dossier_id':None,'run_id':'design:SMARCA2'}
        return self.reports.prepare(result_id)

    def create(self,result_id,operation,request_key,parameters=None,retry_of=None):
        require(operation in TOOLS,'JOB_TOOL_NOT_ALLOWED')
        parameters={} if parameters is None else parameters
        require(type(parameters) is dict,'JOB_PARAMETERS')
        allowed=({'exploratory','dock','use_api','panel_size','linker_ids','target','parent_id','scientific_policy_id'} if operation=='design_panel' else
            {'parent','linker_extension','linker_id','orientation'} if operation=='candidate' else {'candidate_id','sample_index'} if operation=='chemistry' else set())
        require(set(parameters)<=allowed,'JOB_PARAMETERS')
        scientific_policy_id=None
        if operation=='design_panel':
            from packages.science.dual_e3 import validate
            scientific_policy_id=parameters.get('scientific_policy_id')
            require(scientific_policy_id is None or (type(scientific_policy_id) is str and scientific_policy_id),'JOB_DESIGN_PARAMETERS')
            design_parameters=dict(parameters);design_parameters.pop('scientific_policy_id',None)
            try:parameters=validate(design_parameters)
            except (ValueError,TypeError):require(False,'JOB_DESIGN_PARAMETERS')
            if scientific_policy_id is not None:parameters['scientific_policy_id']=scientific_policy_id
            if parameters['use_api']:require(self.provider_factory is not None,'JOB_API_NOT_CONFIGURED')
        if result_id=='design:SMARCA2':require(operation=='design_panel','JOB_DESIGN_SOURCE_TOOL')
        if operation=='candidate':
            from packages.science.linker_design import validate_parameters
            validate_parameters(parameters)
        if operation=='chemistry':
            require(parameters.get('candidate_id','C01') in {'C01','C02'} and type(parameters.get('sample_index',0)) is int
                    and 0<=parameters.get('sample_index',0)<10,'JOB_CHEMISTRY_PARAMETERS')
        if operation=='review':require(self.provider_factory is not None,'JOB_API_NOT_CONFIGURED')
        prepared=self._prepare(result_id)
        require(prepared['can_create'],'JOB_CURRENT_WRITABLE_INPUT_REQUIRED')
        runtime=identity()
        tool_inputs=self.tool_snapshot(result_id) if operation in {'review','report'} else []
        scientific_policy=None
        scientific_binding=None
        if scientific_policy_id is not None:
            scientific_policy=self._resolve_design_policy(scientific_policy_id,parameters.get('parent_id','SMARCA2-FX5'))
            scientific_binding={'id':scientific_policy_id,
                'revision':scientific_policy['scope']['policy_revision'],
                'digest':scientific_policy['scope']['policy_digest']}
        fingerprint=digest(encoded({'input':prepared['input_digest'],'tools':tool_inputs,'operation':operation,
            'parameters':parameters,'runtime':runtime,'retry_of':retry_of,'scientific_policy_binding':scientific_binding}))
        with self.store.db() as db:
            db.execute('BEGIN IMMEDIATE')
            prior=db.execute('SELECT fingerprint,body FROM workbench_jobs WHERE project=? AND request_key=?',(self.project,request_key)).fetchone()
            if prior:
                require(prior[0]==fingerprint,'JOB_REQUEST_KEY_CONFLICT');return json.loads(prior[1]),False
            if retry_of:
                old=self._row(db,retry_of)
                require(old['state'] in {'failed','held','cancelled','interrupted'} and old['result_id']==result_id
                        and old['operation']==operation and old['parameters']==parameters,'JOB_RETRY_CONTEXT')
                require(old['usage_status']!='unreconciled','JOB_UNKNOWN_USAGE_REQUIRES_RECONCILIATION')
            require(not db.execute("SELECT 1 FROM workbench_jobs WHERE project=? AND state IN ('queued','running')",(self.project,)).fetchone(),'JOB_ALREADY_ACTIVE')
            if result_id!='design:SMARCA2':
                require(self.reports._state(db,result_id,prepared['dossier_id'],prepared['run_id'])==prepared['binding']['state'],'JOB_INPUT_CHANGED')
            row={'id':'job-'+uuid.uuid4().hex,'project_id':self.project,'run_id':prepared['run_id'],'result_id':result_id,
                 'operation':operation,'title':TOOLS[operation],'parameters':parameters,'input_digest':prepared['input_digest'],
                 'binding':prepared['binding'],'runtime':runtime,'state':'queued','stage':'queued',
                 'tool_inputs':tool_inputs,
                 'created_at':now(),'updated_at':now(),'retry_of':retry_of,'attempt':old['attempt']+1 if retry_of else 1,
                 'calls':0,'total_tokens':0,'usage_status':'not_started','limits':dict(LIMITS),'outputs':[],
                 'error_code':None,'authority':AUTHORITY.copy(),'result':None,
                 'scientific_policy_binding':scientific_binding,
                 '_scientific_policy_snapshot':scientific_policy}
            db.execute('INSERT INTO workbench_jobs VALUES(?,?,?,?,?,?,?,?)',
                (row['id'],self.project,row['run_id'],result_id,request_key,fingerprint,row['state'],json.dumps(row)))
        return row,True

    def cancel(self,identifier):
        with self.store.db() as db:
            db.execute('BEGIN IMMEDIATE');row=self._row(db,identifier)
            if row['state'] in {'queued','running'}:row.update(state='cancelled',stage='cancel_requested');self._save(db,row)
        return row

    def recover(self):
        # Never interrupt a worker in another process holding the shared execution lock.
        with worker_lock(self.store.root) as acquired:
            if not acquired:return
            with self.store.db() as db:
                db.execute('BEGIN IMMEDIATE')
                for record in db.execute("SELECT body FROM workbench_jobs WHERE project=? AND state IN ('queued','running')",(self.project,)).fetchall():
                    row=json.loads(record[0]);row.update(state='interrupted',stage='server_restarted',error_code='JOB_INTERRUPTED');self._save(db,row)

    def check_active(self,identifier):
        row=self.get(identifier)
        if row['state']!='running':raise AgentError('JOB_CANCELLED_OR_NOT_RUNNING')
        if identity()!=row['runtime']:raise AgentError('JOB_CODE_CHANGED')
        if self._prepare(row['result_id'])['input_digest']!=row['input_digest']:raise AgentError('JOB_INPUT_CHANGED')
        if row['operation'] in {'report','review'} and row.get('tool_inputs',[])!=self.tool_snapshot(row['result_id']):raise AgentError('JOB_TOOL_INPUT_CHANGED')
        binding=row.get('scientific_policy_binding')
        if binding is not None:
            current=self._resolve_design_policy(binding['id'],row['parameters'].get('parent_id','SMARCA2-FX5'),runtime=True)
            actual={'id':binding['id'],'revision':current['scope']['policy_revision'],'digest':current['scope']['policy_digest']}
            if actual!=binding:raise AgentError('JOB_SCIENTIFIC_POLICY_CHANGED')
        if time.monotonic()-self._started>row['limits']['max_elapsed_seconds']:raise AgentError('JOB_TIME_LIMIT')
        return self._agent_digest

    def execute(self,identifier):
        with worker_lock(self.store.root) as acquired:
            if not acquired:
                with self.store.db() as db:
                    db.execute('BEGIN IMMEDIATE');row=self._row(db,identifier)
                    if row['state']=='queued':row.update(state='interrupted',stage='worker_busy',error_code='JOB_WORKER_BUSY');self._save(db,row)
                return
            with self.store.db() as db:
                db.execute('BEGIN IMMEDIATE');row=self._row(db,identifier)
                if row['state']!='queued':return
                row.update(state='running',stage='checking_inputs',started_at=now());self._save(db,row)
            self._started=time.monotonic();self._agent_digest=None;provider=None;outputs=[]
            try:
                self.check_active(identifier)
                def put(name,raw,media='application/json'):
                    ref=self.port.put_raw(raw,media,'computed');outputs.append({'name':name,'ref':ref});return ref
                if row['operation']=='verify':
                    result=self.saved.verify(row['result_id']);state='completed'
                elif row['operation']=='report':
                    report,_=self.reports.create(row['result_id'],row['input_digest'])
                    result={k:report[k] for k in ('id','counts','markdown_ref','content_ref')};state='completed'
                    outputs.extend([{'name':'통합보고서.md','ref':report['markdown_ref']},{'name':'근거.json','ref':report['content_ref']}])
                    if row.get('tool_inputs'):
                        supplement={'format':'tpd-tool-report/0.1.0','saved_report_id':report['id'],
                            'tools':row['tool_inputs'],'authority':AUTHORITY.copy()}
                        put('CPU-experiments.json',encoded(supplement))
                        text=['# CPU 실험 부록','원 계산 보고서와 별도인 설계 가설·리간드 국소 시험입니다.','']
                        for item in row['tool_inputs']:
                            r=item['result'];text.extend(['## '+r['candidate_id']+' · '+item['operation'],
                                '- 상태: '+r['status'],'- 작업: '+item['job_id'],
                                '- 효능·합성 가능성·정식 실행 승인: 미확립',''])
                            if item['operation']=='candidate':
                                from packages.science.linker_design import classification_label
                                text.extend(['- 분류: '+classification_label(r),'- 설계 근거: '+r['rationale'],
                                    '- linker 변형: '+r['transformation'].get('label','메틸렌 경계 연장'),
                                    '- 합성 평가: '+r['synthesis']['status'],
                                    '- 실제 결합 구조·복합체 접촉·효능: 미검증',
                                    '- 다음 검토: '+' / '.join(r['remaining']),''])
                            if item['operation']=='design_panel':
                                from packages.platform.design_panel import report_markdown
                                text.extend([report_markdown(r),''])
                        put('CPU-실험-부록.md','\n'.join(text).encode(),'text/markdown')
                elif row['operation']=='design_panel':
                    from packages.platform.design_panel import run_panel
                    if row['parameters']['use_api']:provider=self.provider_factory()
                    def stage(value):
                        self.check_active(identifier)
                        with self.store.db() as db:
                            db.execute('BEGIN IMMEDIATE');current=self._row(db,identifier)
                            require(current['state']=='running','JOB_NOT_RUNNING')
                            current['stage']=value;self._save(db,current)
                    panel_parameters=dict(row['parameters']);panel_parameters.pop('scientific_policy_id',None)
                    result=run_panel(panel_parameters,put,lambda:self.check_active(identifier),
                        CallJournal(self,identifier,provider) if provider else None,stage,
                        scientific_policy=row.get('_scientific_policy_snapshot'))
                    state='completed'
                elif row['operation']=='review':
                    data=self.saved.view(row['result_id']);dossier=self.saved.dossiers.view(data['dossier_id'])
                    bundle=build_saved_bundle(data,dossier,row.get('tool_inputs',[]));self._agent_digest=bundle['input_digest']
                    put('agent-input.json',encoded(bundle))
                    provider=self.provider_factory()
                    runtime=ReviewRuntime(CallJournal(self,identifier,provider),self.store.root/'workbench-attempts'/identifier,
                        model='gpt-5.6-sol',**{k:v for k,v in row['limits'].items() if k!='max_elapsed_seconds'})
                    result=runtime.run(bundle,revalidate=lambda:self.check_active(identifier))
                    for name,media in [('REVIEW.md','text/markdown'),('run.json','application/json')]:
                        p=runtime.output/name
                        if p.exists():put(name,p.read_bytes(),media)
                    state='review_ready' if result['status']=='draft_pending_human_review' else 'held'
                else:
                    from packages.platform.cpu_tools import run_tool
                    result=run_tool(row['operation'],self.saved,row['result_id'],row['parameters'],put,
                                    lambda:self.check_active(identifier));state='completed'
                put('result.json',encoded(result))
                self.check_active(identifier)
                with self.store.db() as db:
                    db.execute('BEGIN IMMEDIATE');current=self._row(db,identifier)
                    if current['state']=='running':
                        current.update(state=state,stage='finished',result=result,outputs=outputs,
                            error_code=result.get('error') if isinstance(result,dict) else None,finished_at=now())
                    else:current.update(outputs=outputs,result=result)
                    self._save(db,current)
            except Exception as exc:
                code=str(exc) if isinstance(exc,(AgentError,ContractError)) else 'JOB_EXECUTION_FAILED'
                if not re_safe_code(code):code='JOB_EXECUTION_FAILED'
                with self.store.db() as db:
                    db.execute('BEGIN IMMEDIATE');current=self._row(db,identifier)
                    if current['state']=='running':current.update(state='held' if isinstance(exc,AgentError) else 'failed',stage='finished',error_code=code)
                    current.update(outputs=outputs,finished_at=now());self._save(db,current)
            finally:
                if provider:provider.close()

    def view(self,identifier):
        row=self.get(identifier)
        with self.store.db() as db:row['calls_detail']=[json.loads(r[0]) for r in db.execute('SELECT body FROM workbench_calls WHERE job_id=? ORDER BY ordinal',(identifier,))]
        for output in row['outputs']:self.port.read(output['ref'])
        for call in row['calls_detail']:
            for key in ('prompt_ref','answer_ref'):
                if call.get(key):self.port.read(call[key])
        current=self._prepare(row['result_id'])
        row['freshness']={'current':current['input_digest']==row['input_digest'] and identity()==row['runtime']
            and (row['operation'] not in {'review','report'} or row.get('tool_inputs',[])==self.tool_snapshot(row['result_id']))}
        return self._public_row(row)


def re_safe_code(value):
    return isinstance(value,str) and len(value)<=100 and value.replace('_','').isalnum() and value.isascii()
