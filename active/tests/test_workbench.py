"""Synthetic job-control tests, never evidence of molecular efficacy."""
import copy
import io
import json
import time
import unittest
import zipfile
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch
from fastapi.testclient import TestClient
from apps.api.main import create_app,PROJECT
from packages.contracts import ContractError
from packages.agents.provider import Completion,AgentError
from packages.platform.workbench import WorkbenchService,CallJournal
from packages.platform.expert_responses import extract_answers
from packages.platform.review_identity import ReviewIdentityService
from packages.platform import workbench_replay as replay
import test_evidence_reports as fixture
import test_candidate_dossiers as dossier_fixture

class JobTests(unittest.TestCase):
    def setUp(self):
        fixture.EvidenceReportTests.setUp(self)
        self.jobs=WorkbenchService(self.store,PROJECT,lambda:None)
        self.jobs.reports=self.svc;self.jobs.saved=self.svc.saved

    def create(self,key='synthetic-key',operation='verify',parameters=None):
        return self.jobs.create('s',operation,key,parameters)

    def test_concurrent_idempotency_and_active_admission(self):
        with ThreadPoolExecutor(2) as pool:rows=list(pool.map(lambda _:self.create(),range(2)))
        self.assertEqual(sum(fresh for _,fresh in rows),1)
        self.assertEqual(rows[0][0]['id'],rows[1][0]['id'])
        with self.assertRaisesRegex(ContractError,'ALREADY_ACTIVE'):self.create('other-key')
        with self.assertRaisesRegex(ContractError,'KEY_CONFLICT'):self.create(operation='candidate')

    def test_invalid_parameters_and_nonallowlisted_tool(self):
        for op,p in [('shell',{}),('verify',{'command':'x'}),('candidate',{'parent':'C03'}),
                     ('candidate',{'linker_extension':True}),('chemistry',{'sample_index':-1})]:
            with self.assertRaises(ContractError):self.create(operation=op,parameters=p)

    def test_execution_hash_protected_output(self):
        j,_=self.create();self.jobs.execute(j['id']);v=self.jobs.view(j['id'])
        self.assertEqual(v['state'],'completed');self.assertEqual(v['result']['files_verified'],1)
        self.assertTrue(v['freshness']['current']);ref=v['outputs'][0]['ref']
        (self.store.root/'blobs'/ref['artifact_id']).write_bytes(b'changed')
        with self.assertRaisesRegex(ContractError,'CORRUPTED'):self.jobs.view(j['id'])

    def test_cancel_before_execution_and_explicit_retry(self):
        j,_=self.create();self.jobs.cancel(j['id']);self.jobs.execute(j['id'])
        self.assertEqual(self.jobs.get(j['id'])['calls'],0)
        retry,fresh=self.jobs.create('s','verify','new-key',{},j['id'])
        self.assertTrue(fresh);self.assertEqual(retry['attempt'],2)
        self.jobs.execute(retry['id']);self.assertEqual(self.jobs.get(retry['id'])['state'],'completed')

    def running(self):
        j,_=self.create(operation='review')
        with self.store.db() as db:j['state']='running';self.jobs._save(db,j)
        self.jobs._started=time.monotonic();self.jobs._agent_digest='synthetic'
        return j

    def test_cancel_during_call_preserves_observed_usage(self):
        j=self.running();jobs=self.jobs
        class Provider:
            mode='synthetic'
            def complete(self,**request):
                jobs.cancel(j['id']);return Completion('synthetic',{'usage':{'total_tokens':42}})
        with self.assertRaisesRegex(AgentError,'CANCELLED'):
            CallJournal(jobs,j['id'],Provider()).complete(model='synthetic',context={},instructions='synthetic',max_output_tokens=10)
        v=jobs.view(j['id']);self.assertEqual(v['state'],'cancelled');self.assertEqual(v['total_tokens'],42)
        self.assertEqual(v['usage_status'],'observed');self.assertIsNotNone(v['calls_detail'][0]['answer_ref'])

    def test_unknown_usage_survives_recovery_and_blocks_retry(self):
        j=self.running()
        class Provider:
            mode='synthetic'
            def complete(self,**request):raise AgentError('API_TRANSPORT_ERROR_USAGE_UNKNOWN')
        with self.assertRaises(AgentError):CallJournal(self.jobs,j['id'],Provider()).complete(context={},model='synthetic')
        self.jobs.recover();v=self.jobs.get(j['id']);self.assertEqual(v['state'],'interrupted')
        self.assertEqual(v['usage_status'],'unreconciled')
        with self.assertRaisesRegex(ContractError,'UNKNOWN_USAGE'):self.jobs.create('s','review','retry-key',{},j['id'])

    def test_stale_input_and_runtime_held_before_tool(self):
        for change in ['input','runtime']:
            j,_=self.create(key=change)
            if change=='input':self.run['revision']=2
            else:j['runtime']['version']='old'
            if change=='input':
                with self.store.db() as db:db.execute('UPDATE runs SET body=? WHERE id=?',(json.dumps(self.run),'r'))
            else:
                with self.store.db() as db:self.jobs._save(db,j)
            self.jobs.execute(j['id']);self.assertEqual(self.jobs.get(j['id'])['state'],'held')
            self.run['revision']=1
            with self.store.db() as db:db.execute('UPDATE runs SET body=? WHERE id=?',(json.dumps(self.run),'r'))

    def test_local_and_read_only_api_boundaries(self):
        body={'result_id':'s','operation':'verify','request_key':'test-12345678','parameters':{}}
        for readonly in (True,False):
            with TestClient(create_app(self.store.root,enable_worker=False,read_only=readonly)) as c:
                self.assertEqual(c.post('/api/lab/jobs',json=body).status_code,403)
                if readonly:self.assertEqual(c.post('/api/lab/jobs',json=body,headers={'X-TPD-Local':'1'}).status_code,403)
                self.assertEqual(c.post('/api/lab/jobs',json=body,headers={'X-TPD-Local':'1','Origin':'https://other.invalid'}).status_code,403)

class DocumentTests(unittest.TestCase):
    headings=['1. 원자형·수소결합 역할 차이','2. Protonation·tautomer·수소 방향',
              '3. RMSD·Boltz confidence·ternary 구조 비교','4. 문헌 관측값·조건·결측 표현','5. 합성 근거와 신규 C03 진행 기준']
    def doc(self,lines):
        out=io.BytesIO()
        with zipfile.ZipFile(out,'w') as z:z.writestr('word/document.xml','<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body>'+''.join('<w:p><w:r><w:t>'+x+'</w:t></w:r></w:p>' for x in lines)+'</w:body></w:document>')
        return out.getvalue()
    def test_five_explicit_sections_and_summary_boundary(self):
        lines=[v for i,h in enumerate(self.headings) for v in [h,'synthetic answer '+str(i)]]
        a=extract_answers(self.doc(lines+['개발팀 전달용 최종 요약','not Q05']))
        self.assertEqual(len(a),5);self.assertEqual(a['Q05'],'synthetic answer 4')
    def test_partial_and_duplicate_numbered_sections_rejected(self):
        with self.assertRaisesRegex(ContractError,'FIVE_SECTIONS'):extract_answers(self.doc(self.headings[:4]))
        with self.assertRaisesRegex(ContractError,'DUPLICATE'):extract_answers(self.doc(self.headings+self.headings[:1]))

class WorkbenchReplayTests(unittest.TestCase):
    def setUp(self):
        self.f=dossier_fixture.DossierTests();self.f.setUp();self.addCleanup(self.f.doCleanups);self.f.load()
        identity=ReviewIdentityService(self.f.store,PROJECT);self.secret=identity.register('synthetic-reviewer','SYNTHETIC')
        self.raw,self.manifest=replay.export_bytes(self.f.store,PROJECT,self.f.run['id'])
    def alter(self,fn):
        with zipfile.ZipFile(io.BytesIO(self.raw)) as z:files={n:z.read(n) for n in z.namelist()}
        fn(files);out=io.BytesIO()
        with zipfile.ZipFile(out,'w') as z:
            for n,b in files.items():z.writestr(n,b)
        return out.getvalue()
    def test_roundtrip_excludes_accounts_and_sessions(self):
        store,m=replay.import_archive(self.raw,self.f.root/'workbench-replay',PROJECT)
        with store.db() as db:self.assertEqual(db.execute('SELECT COUNT(*) FROM review_actors').fetchone()[0],0)
        self.assertTrue(store.get_run(PROJECT,self.f.run['id'])['replay_only'])
        self.assertNotIn(self.secret,encoded_json(m))
    def test_base_hash_and_unsafe_path_rejected_before_destination(self):
        for name in ['base.zip','../unsafe']:
            dest=self.f.root/'bad'
            with self.assertRaises(ContractError):replay.import_archive(self.alter(lambda f:f.update({name:b'altered'})),dest,PROJECT)
            self.assertFalse(dest.exists())
    def test_unknown_table_cannot_import_credentials(self):
        from packages.contracts import encoded
        from packages.platform.dossiers import digest
        def change(files):
            m=json.loads(files['manifest.json']);m['records']['review_actors']=[];m.pop('digest');m['digest']=digest(encoded(m));files['manifest.json']=encoded(m)
        with self.assertRaisesRegex(ContractError,'TABLES'):replay.inspect_archive(self.alter(change))

def encoded_json(value):return json.dumps(value)
