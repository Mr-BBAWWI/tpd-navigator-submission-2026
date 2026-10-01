"""Synthetic technical tests. These records are never used as scientific results."""
import copy
import json
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch

from fastapi.testclient import TestClient
from apps.api.main import PROJECT, create_app
from packages.contracts import ContractError, encoded
from packages.platform.store import Store
from packages.platform.evidence_reports import EvidenceReportService, implementation, bind_b_sources, references
from packages.platform.saved_results import AUTHORITY
from packages.reporting.saved_report import assemble, json_block, safe
from packages.science.evidence_common import schema_check


class EvidenceReportTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.store=Store(self.tmp.name);self.svc=EvidenceReportService(self.store,PROJECT)
        self.ref=self.store.scope(PROJECT).put_json({'synthetic':True})
        self.run={'revision':1,'bundle_ref':self.ref,'target_ref':self.ref,'target_query_ref':self.ref,'data_mode':'real'}
        self.dossier={'id':'d','run_id':'r','dossier_ref':self.ref,'anchors':{},'planned_candidates':[],
            'candidates':[{'id':'C01','paper_name':'synthetic 1','paper_compound_number':'2','pdb':'TEST','ccd':'TEST','molecule_id':'synthetic','attachment':{}}]}
        sample={'run_id':'synthetic-run','rank':0,'binding':{},'reported_execution':{'settings':{'seed':23,'use_potentials':False}},
            'summary':{'target_rmsd_A':0,'vhl_rmsd_after_target_alignment_A':None,'confidence_score':0},
            'bad_overlap_pairs':{'xray':0,'nuclear':None},'parts_rmsd_ranges_A':{},'reports':{'comparison':self.ref},'status':{'contacts':'not_provided'}}
        self.data={'id':'s','dossier_id':'d','run_id':'r','record_ref':self.ref,'dossier_ref':self.ref,'projection_ref':self.ref,
            'snapshot':self.run.copy(),'bundle_digest':'a'*64,'summary':{'candidate_count':1,'sample_count':1},
            'candidates':[{'compound_id':'C01','molecule_id':'synthetic','samples':[sample]}],
            'literature':{'records':[],'missing':[{'compound_id':'C01','value':None}]},'expert_reviews':[],
            'source_authority':AUTHORITY.copy()}
        with self.store.db() as db:
            db.execute('INSERT INTO runs VALUES(?,?,?,?,?,?,?)',('r',PROJECT,'key','hash',1,'completed',json.dumps(self.run)))
            db.execute('INSERT INTO dossiers VALUES(?,?,?,?,?)',('d',PROJECT,'r','hash',json.dumps(self.ref)))
            db.execute('INSERT INTO saved_results VALUES(?,?,?,?,?)',('s',PROJECT,'d','hash',json.dumps(self.ref)))
        for name,value in [('view',self.data),('_record',self.data),('verify',{'files_verified':1})]:
            p=patch.object(self.svc.saved,name,side_effect=lambda *a,_v=value,**kw:copy.deepcopy(_v));p.start();self.addCleanup(p.stop)
        p=patch.object(self.svc.saved,'files',return_value={});p.start();self.addCleanup(p.stop)
        p=patch.object(self.svc.saved.dossiers,'view',side_effect=lambda *a:copy.deepcopy(self.dossier));p.start();self.addCleanup(p.stop)

    def create(self):
        return self.svc.create('s',self.svc.prepare('s')['input_digest'])

    def test_report_preserves_missing_zero_authority_and_bytes(self):
        report,fresh=self.create()
        self.assertTrue(fresh);self.assertTrue(report['freshness']['current'])
        self.assertEqual(report['authority'],AUTHORITY)
        sample=report['content']['candidates'][0]['saved_evidence']['samples'][0]
        self.assertEqual(sample['summary']['target_rmsd_A'],0)
        self.assertIsNone(sample['summary']['vhl_rmsd_after_target_alignment_A'])
        text=self.svc.port.read(report['markdown_ref']).decode()
        self.assertIn('| 23 | OFF | 0 | 미확인 | 0 |',text)
        for flag in ('human_approved','dispatch_authorized','public_release_ready'):
            forged=self.svc.port.json(report['record_ref']);forged['authority'][flag]=True
            with self.assertRaises(ValueError):schema_check(forged,'a_evidence_report.schema.json')

    def test_duplicate_and_concurrent_creation_has_single_record(self):
        token=self.svc.prepare('s')['input_digest']
        with ThreadPoolExecutor(max_workers=2) as pool:
            rows=list(pool.map(lambda _:self.svc.create('s',token),range(2)))
        self.assertEqual(len({v['id'] for v,_ in rows}),1)
        self.assertEqual(sum(fresh for _,fresh in rows),1)
        before=len(list((self.store.root/'blobs').iterdir()))
        self.assertFalse(self.svc.create('s',token)[1])
        self.assertEqual(len(list((self.store.root/'blobs').iterdir())),before)

    def test_changed_input_invalidates_report_without_modifying_it(self):
        report,_=self.create();original=self.svc.port.read(report['record_ref'])
        run={**self.run,'revision':2}
        with self.store.db() as db: db.execute('UPDATE runs SET body=? WHERE id=?',(json.dumps(run),'r'))
        view=self.svc.view(report['id'])
        self.assertEqual(view['freshness']['reasons'],['input_changed'])
        self.assertEqual(view['freshness']['review_state'],'requires_re_review')
        self.assertEqual(self.svc.port.read(report['record_ref']),original)
        with self.assertRaisesRegex(ContractError,'CURRENT_WRITABLE'):self.create()

    def test_new_opinion_or_result_or_dossier_invalidates(self):
        report,_=self.create()
        with self.store.db() as db:
            db.execute('INSERT INTO saved_expert_reviews VALUES(?,?,?,?,?)',('new-review',PROJECT,'s','new',json.dumps(self.ref)))
        self.assertIn('expert_opinions_changed',self.svc.view(report['id'])['freshness']['reasons'])
        with self.store.db() as db:
            db.execute('INSERT INTO saved_results VALUES(?,?,?,?,?)',('new-result',PROJECT,'d','new',json.dumps(self.ref)))
            db.execute('INSERT INTO dossiers VALUES(?,?,?,?,?)',('new-dossier',PROJECT,'r','new',json.dumps(self.ref)))
        self.assertEqual(set(self.svc.view(report['id'])['freshness']['reasons']),
                         {'expert_opinions_changed','result_superseded','dossier_superseded'})

    def test_implementation_change_requires_new_review(self):
        report,_=self.create();identity=implementation();identity['version']='changed'
        with patch('packages.platform.evidence_reports.implementation',return_value=identity):
            self.assertEqual(self.svc.view(report['id'])['freshness']['reasons'],['report_implementation_changed'])
            with self.assertRaisesRegex(ContractError,'INPUT_CHANGED'):self.svc.create('s',report['input_digest'])

    def test_racing_opinion_during_assembly_rolls_back(self):
        token=self.svc.prepare('s')['input_digest'];before=len(list((self.store.root/'blobs').iterdir()))
        def race(*args):
            with self.store.db() as db:
                db.execute('INSERT INTO saved_expert_reviews VALUES(?,?,?,?,?)',('racing-review',PROJECT,'s','new',json.dumps(self.ref)))
            return assemble(*args)
        with patch('packages.platform.evidence_reports.assemble',side_effect=race):
            with self.assertRaisesRegex(ContractError,'INPUT_CHANGED'):self.svc.create('s',token)
        with self.store.db() as db:self.assertEqual(db.execute('SELECT COUNT(*) FROM evidence_reports').fetchone()[0],0)
        self.assertEqual(len(list((self.store.root/'blobs').iterdir())),before)

    def test_corrupt_content_or_source_fails_closed(self):
        report,_=self.create()
        for ref in (report['content_ref'],self.ref):
            path=self.store.root/'blobs'/ref['artifact_id'];raw=path.read_bytes()
            try:
                path.write_bytes(b'corrupted')
                with self.assertRaises(ContractError):self.svc.view(report['id'])
            finally:path.write_bytes(raw)

    def test_exact_sample_coverage_rejects_duplicate_or_missing(self):
        self.data['candidates'][0]['samples']*=2
        with self.assertRaisesRegex(ContractError,'SAMPLE_COVERAGE'):assemble(self.data,self.dossier)
        self.data['candidates'][0]['samples']=[]
        with self.assertRaisesRegex(ContractError,'SAMPLE_COVERAGE'):assemble(self.data,self.dossier)

    def test_foreign_b_namespace_requires_exact_file_hash_binding(self):
        foreign={**self.ref,'artifact_id':'b-material-'+self.ref['sha256']}
        body={'item':{'path':'compound.json','ref':foreign}}
        with self.assertRaisesRegex(ContractError,'UNBOUND_B_SOURCE'):references(body)
        with self.assertRaisesRegex(ContractError,'B_SOURCE_BINDING'):bind_b_sources(body,{})
        bind_b_sources(body,{'source/quality/evidence/cpu/compound.json':self.ref})
        self.assertEqual(references(body),[self.ref])
        self.assertEqual(body['item']['ref'],foreign)

    def test_untrusted_markdown_cannot_close_block_or_create_link(self):
        malicious='```\n# approved\n[download](javascript:alert(1)) <script>'
        block=json_block({'source':malicious})
        self.assertEqual(block[0],'````json');self.assertEqual(block[-2],'````')
        escaped=safe(malicious)
        self.assertNotIn('<script>',escaped);self.assertNotIn('[download](',escaped)

    def test_api_read_only_scope_token_and_download(self):
        report,_=self.create()
        other=EvidenceReportService(self.store,'other')
        with self.assertRaises(KeyError):other.view(report['id'])
        with TestClient(create_app(self.store.root,enable_worker=False,read_only=True)) as client:
            self.assertEqual(client.get('/api/evidence-reports/'+report['id']).status_code,200)
            self.assertEqual(client.post('/api/saved-results/s/reports',json={'expected_digest':report['input_digest']},headers={'X-TPD-Local':'1'}).status_code,403)
            response=client.get('/api/artifacts/'+report['markdown_ref']['artifact_id']+'/download')
            self.assertTrue(response.headers['content-disposition'].endswith('.md"'))
            self.assertEqual(response.content,self.svc.port.read(report['markdown_ref']))
        app=create_app(self.store.root,enable_worker=False)
        with TestClient(app) as client, patch.object(app.state.evidence_reports,'create',wraps=self.svc.create):
            url='/api/saved-results/s/reports'
            self.assertEqual(client.post(url,json={'expected_digest':report['input_digest']}).status_code,403)
            self.assertEqual(client.post(url,json={'expected_digest':'b'*64},headers={'X-TPD-Local':'1'}).status_code,409)
            self.assertEqual(client.post(url,json={'expected_digest':report['input_digest'],'actor':'forged'},headers={'X-TPD-Local':'1'}).status_code,422)
            response=client.post(url,json={'expected_digest':report['input_digest']},headers={'X-TPD-Local':'1'})
            self.assertEqual(response.status_code,200)
            self.assertEqual(response.json()['report_id'],report['id'])
            self.assertFalse(response.json()['created'])


if __name__=='__main__':unittest.main()
