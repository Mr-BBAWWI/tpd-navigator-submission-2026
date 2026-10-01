"""Workflow failure/recovery tests use scripted responses, never a real LLM."""
import copy
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch
from fastapi.testclient import TestClient
from apps.api.main import create_app

from packages.agents.bundle import json_digest
from packages.agents.provider import AgentError
from packages.contracts import ContractError
from packages.platform.dossiers import DossierService, ROOT
from packages.platform.handoffs import HandoffService
from packages.platform.store import Store
from packages.platform.workflows import WorkflowService
import test_agents as agents_fixture


class Crash(BaseException):
    """Simulated process loss bypasses ordinary application exception handlers."""


class Provider(agents_fixture.ScriptedProvider):
    def __init__(self, action=None, **kwargs):
        super().__init__(**kwargs); self.count=0; self.action=action

    def complete(self, **kwargs):
        self.count+=1
        if self.action:self.action(self.count)
        return super().complete(**kwargs)

    def close(self):pass


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.store=Store(self.root/'store');self.project='local-research'
        self.run,_=self.store.create_run(self.project,'workflow-fixture',{'query':'fixture'})
        self.run=self.store.update(self.project,self.run['id'],{'state':'completed'})
        self.port=self.store.scope(self.project)
        self.dossier=json.loads((ROOT/'outputs/a_handoff_20260923/candidate_dossier.json').read_text('utf-8'))
        self.dossier.update(run_id=self.run['id'],snapshot=DossierService.snapshot(self.run),current_input=True,
                            analysis_ids=[])
        self.dossier['artifact_refs']={p.name:self.port.put_raw(p.read_bytes(),'application/json' if p.suffix=='.json' else 'application/octet-stream','computed')
            for p in (ROOT/'outputs/b_evidence_20260923/cpu').iterdir() if p.name!='handoff.json'}
        self.dossier['handoff_ref']=self.port.put_raw((ROOT/'outputs/b_evidence_20260923/cpu/handoff.json').read_bytes(),'application/json','computed')
        for c in self.dossier['candidates']:c['literature_evidence']=[]
        self.dossier['dossier_ref']=self.port.put_json({k:v for k,v in self.dossier.items() if k not in ('dossier_ref','current_input')})
        def view(service,identifier):
            if identifier!=self.dossier['id']:raise KeyError(identifier)
            d=copy.deepcopy(self.dossier)
            d['current_input']=DossierService.snapshot(self.store.get_run(self.project,self.run['id']))==d['snapshot']
            return d
        self.patcher=patch.object(DossierService,'view',view);self.patcher.start();self.addCleanup(self.patcher.stop)
        self.handoffs=HandoffService(self.store,self.project)
        self.input,_=self.handoffs.import_export(self.dossier['id'],ROOT/'outputs/b_evidence_20260923')
        b=agents_fixture.small_bundle();b['b_input_ref']=self.input['input_ref'];b['input_digest']=json_digest({k:v for k,v in b.items() if k!='input_digest'})
        self.bundle_patch=patch('packages.platform.workflows.build_service_bundle',return_value=b)
        self.bundle_patch.start();self.addCleanup(self.bundle_patch.stop)
        self.provider=Provider();self.service=WorkflowService(self.store,self.project,lambda:self.provider)

    def job(self):return self.service.create(self.input['id'],'test-workflow-key')[0]

    def test_import_idempotent_and_unfilled_values_have_no_authority(self):
        again,fresh=self.handoffs.import_export(self.dossier['id'],ROOT/'outputs/b_evidence_20260923')
        self.assertFalse(fresh);self.assertEqual(again['id'],self.input['id'])
        packet=self.handoffs.packet(again['id'])
        self.assertFalse(packet['ready_for_execution']);self.assertFalse(packet['dispatch_authorized'])
        self.assertTrue(all(x['value'] is None and x['status']=='pending' for x in packet['pending_inputs']))

    def test_corrupted_B_export_rejected(self):
        target=self.root/'bad';shutil.copytree(ROOT/'outputs/b_evidence_20260923',target)
        (target/'cpu/C01.smi').write_text('wrong')
        with self.assertRaises(ValueError):self.handoffs.import_export(self.dossier['id'],target)

    def test_wrong_molecule_even_with_matching_names_rejected(self):
        self.dossier['candidates'][0]['molecule_id']='wrong'
        with self.assertRaisesRegex(ContractError,'IDENTITY'):self.handoffs.import_export(self.dossier['id'],ROOT/'outputs/b_evidence_20260923')

    def test_service_full_run_idempotency_and_duplicate_worker(self):
        job=self.job();same,fresh=self.service.create(self.input['id'],'test-workflow-key')
        self.assertFalse(fresh);self.assertEqual(job['id'],same['id'])
        self.service.execute(job['id']);self.service.execute(job['id'])
        result=self.service.view(job['id'])
        self.assertEqual(result['state'],'review_ready');self.assertEqual(self.provider.count,7)
        self.assertEqual(result['calls'],7);self.assertEqual(result['total_tokens'],700)
        self.assertIsNotNone(result['report_ref']);self.assertIsNotNone(result['packet_ref'])
        self.assertFalse(result['human_approved']);self.assertFalse(result['dispatch_authorized'])

    def test_critic_one_repair_and_unresolved_hold(self):
        self.provider.revision=True;job=self.job();self.service.execute(job['id'])
        self.assertEqual(self.service.get(job['id'])['state'],'review_ready')
        self.assertEqual(self.provider.count,10)

    def test_completion_waits_until_review_material_is_saved(self):
        job=self.job(); original=self.service.handoffs.packet
        def observed(*args,**kwargs):
            self.assertEqual(self.service.get(job['id'])['state'],'running')
            return original(*args,**kwargs)
        with patch.object(self.service.handoffs,'packet',side_effect=observed):
            self.service.execute(job['id'])
        result=self.service.get(job['id'])
        self.assertEqual(result['state'],'review_ready');self.assertIsNotNone(result['packet_ref'])

    def test_review_material_failure_holds_and_preserves_completed_answers(self):
        job=self.job()
        with patch.object(self.service.handoffs,'packet',side_effect=ValueError('synthetic failure')):
            self.service.execute(job['id'])
        result=self.service.get(job['id'])
        self.assertEqual(result['state'],'held');self.assertEqual(result['error_code'],'WORKFLOW_EXECUTION_FAILED')
        self.assertIsNotNone(result['result_ref']);self.assertIsNone(result['packet_ref'])
        self.assertEqual(len(self.service.calls(job['id'])),7)

    def test_uncertain_network_failure_is_not_resumable(self):
        self.provider.action=lambda n:(_ for _ in ()).throw(AgentError('API_TIMEOUT'))
        job=self.job();self.service.execute(job['id'])
        view=self.service.view(job['id'])
        self.assertEqual(view['state'],'held');self.assertEqual(view['usage_status'],'unreconciled')
        with self.assertRaises(ValueError):self.service.resume(job['id'])
        self.assertEqual(self.provider.count,1)

    def test_process_loss_during_request_blocks_resume(self):
        self.provider.action=lambda n:(_ for _ in ()).throw(Crash())
        job=self.job()
        with self.assertRaises(Crash):self.service.execute(job['id'])
        self.service.interrupt_pending()
        self.assertEqual(self.service.get(job['id'])['state'],'interrupted')
        self.assertFalse(self.service.view(job['id'])['can_resume'])

    def test_resume_reuses_saved_answer_without_new_first_call(self):
        job=self.job();original=self.service.check_active
        def crash_after_saved(identifier):
            calls=self.service.calls(identifier)
            if calls and calls[-1]['answer_ref']:raise Crash()
            return original(identifier)
        with patch.object(self.service,'check_active',side_effect=crash_after_saved):
            with self.assertRaises(Crash):self.service.execute(job['id'])
        self.service.interrupt_pending();self.assertTrue(self.service.view(job['id'])['can_resume'])
        self.service.resume(job['id']);self.service.execute(job['id'])
        result=self.service.view(job['id'])
        self.assertEqual(result['state'],'review_ready');self.assertEqual(result['attempts'],2)
        self.assertEqual(self.provider.count,7);self.assertEqual(result['calls'],7)

    def test_cancel_during_call_retains_answer_and_no_more_calls(self):
        job=self.job();self.provider.action=lambda n:self.service.cancel(job['id'])
        self.service.execute(job['id']);result=self.service.view(job['id'])
        self.assertEqual(result['state'],'cancelled');self.assertEqual(self.provider.count,1)
        self.assertEqual(result['total_tokens'],100);self.assertIsNotNone(result['calls_detail'][0]['answer_ref'])

    def test_input_change_during_call_marks_stale(self):
        job=self.job()
        self.provider.action=lambda n:self.store.update(self.project,self.run['id'],{'revision':2})
        self.service.execute(job['id']);result=self.service.view(job['id'])
        self.assertEqual(result['state'],'stale');self.assertFalse(result['current_input'])
        self.assertEqual(self.provider.count,1)

    def test_second_project_cannot_read_or_cancel(self):
        job=self.job();other=WorkflowService(self.store,'other',lambda:self.provider)
        for method in (other.get,other.view,other.cancel):
            with self.assertRaises(KeyError):method(job['id'])

    def test_replay_collection_cannot_start(self):
        self.store.update(self.project,self.run['id'],{'replay_only':True})
        with self.assertRaises(ValueError):self.job()

    def test_runtime_change_blocks_resume(self):
        job=self.job();self.service.interrupt_pending()
        with patch('packages.platform.workflows.VERSION','changed'):
            with self.assertRaises(ValueError):self.service.resume(job['id'])

    def test_local_api_authority_rejection_and_read_never_calls_model(self):
        app=create_app(self.store.root,enable_worker=False,provider_factory=lambda:self.provider)
        with TestClient(app) as client:
            body={'b_input_id':self.input['id'],'request_key':'api-workflow-test'}
            self.assertEqual(client.post('/api/workflows',json=body).status_code,403)
            self.assertEqual(client.post('/api/workflows',json={**body,'human_approved':True},headers={'x-tpd-local':'1'}).status_code,422)
            response=client.post('/api/workflows',json=body,headers={'x-tpd-local':'1'})
            self.assertEqual(response.status_code,202)
            identifier=response.json()['workflow_id']
            self.assertEqual(client.get('/api/workflows/'+identifier).status_code,200)
            self.assertEqual(client.get('/api/b-inputs/'+self.input['id']+'/review-preparation').status_code,200)
            self.assertEqual(self.provider.count,0)
            self.assertEqual(client.post('/api/workflows/'+identifier+'/cancel',headers={'x-tpd-local':'1'}).status_code,200)

    def test_worker_busy_is_resumable_without_a_request(self):
        from packages.platform.analysis import worker_lock
        job=self.job()
        with worker_lock(self.store.root):self.service.execute(job['id'])
        self.assertTrue(self.service.view(job['id'])['can_resume']);self.assertEqual(self.provider.count,0)

    def test_review_packet_cannot_claim_human_approval(self):
        from packages.science.evidence_common import schema_check
        packet=self.handoffs.packet(self.input['id']);packet['human_approved']=True
        with self.assertRaises(ValueError):schema_check(packet,'a_review_handoff.schema.json')
