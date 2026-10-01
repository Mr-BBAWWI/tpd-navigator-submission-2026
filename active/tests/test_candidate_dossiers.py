"""Offline consistency/admission tests. Synthetic source in temporary storage only.

The small committed B output package is an import fixture, not a new chemistry run.
No official raw-data cache, RDKit, GPU, API key or network is needed.
"""
import copy
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient
from apps.api.main import create_app,PROJECT
from packages.contracts import encoded,ContractError
from packages.platform.dossiers import DossierService,read_package,digest,ROOT
from packages.platform.store import Store
from packages.platform.orchestrator import Orchestrator
from packages.transport import Fetch
from test_i2 import FakeTransport,target_row

XML=b'<article><front><article-meta><title-group><article-title>Synthetic test only</article-title></title-group></article-meta></front><body><sec><p>PROTAC 1 ( 2 ); PDB ID: 6HAY; PROTAC 2 ( 3 ); PDB ID: 6HAX; piperazine ring. Synthetic test only, no new experiment.</p></sec></body></article>'


class DossierTransport(FakeTransport):
    def __init__(self):
        super().__init__(rows=[target_row('P51531')])

    def get(self,url,params=None,max_bytes=2000000):
        if 'uniprot' in url:
            return super().get(url,params,max_bytes)
        if url.endswith('/search'):
            data={'hitCount':1,'resultList':{'result':[{'source':'MED','id':'1','pmcid':'PMC6600871',
                'doi':'10.1038/s41589-019-0294-6','title':'Synthetic test only','isOpenAccess':'Y'}]}}
            return Fetch(encoded(data),url,200,{})
        return Fetch(XML,url,200,{})


class DossierTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.package=self.root/'package'
        shutil.copytree(ROOT/'outputs/smarca2_20260923',self.package)
        self.case=json.loads((ROOT/'cases/smarca2_farnaby2019.json').read_text('utf-8'))
        for s in self.case['sources']:
            if s['file']=='article.xml':s['sha256']=digest(XML)
        self.case_path=self.root/'case.json';self.case_path.write_bytes(encoded(self.case))
        (self.package/'case-input.json').write_bytes(encoded(self.case))
        sources=json.loads((self.package/'sources.json').read_text())
        for s in sources:
            if s['file']=='article.xml':s.update(sha256=digest(XML),bytes=len(XML))
        (self.package/'sources.json').write_bytes(encoded(sources))
        self.handoff=json.loads((self.package/'handoff.json').read_text())
        self.handoff['input_sha256']=digest(self.case_path.read_bytes());self.resign()
        links=json.loads((ROOT/'cases/smarca2_literature_links.json').read_text())
        for a in links['anchors']:a['xpath']='/article/body/sec/p'
        self.links_path=self.root/'links.json';self.links_path.write_bytes(encoded(links))
        self.store=Store(self.root/'store')
        self.run,_=self.store.create_run(PROJECT,'synthetic-dossier-test',{'query':'TESTGENE','query_kind':'gene_symbol','objective':'','cell_line':'','pdf_artifact_ids':[],'max_documents':2})
        Orchestrator(self.store,PROJECT,DossierTransport).execute(self.run['id'])
        self.service=DossierService(self.store,PROJECT)

    def resign(self):
        for a in self.handoff['artifacts']:
            raw=(self.package/a['path']).read_bytes();a.update(sha256=digest(raw),bytes=len(raw))
        (self.package/'handoff.json').write_bytes(encoded(self.handoff))

    def load(self):
        return self.service.import_case(self.run['id'],self.package,self.case_path,self.links_path)

    def test_import_history_idempotency_and_no_authority(self):
        a,fresh=self.load();b,again=self.load()
        self.assertTrue(fresh);self.assertFalse(again);self.assertEqual(a['id'],b['id'])
        self.assertTrue(a['current_input']);self.assertFalse(any(a['authority'].values()))
        self.assertEqual(a['candidates'][0]['attachment']['candidate_atom'],'N46')
        self.assertEqual(a['candidates'][1]['attachment']['candidate_atom'],'N41')
        self.assertEqual(a['planned_candidates'][0]['status'],'not_designed')
        self.assertEqual(len(self.service.list(self.run['id'])),1)

    def test_corrupted_artifact_rejected_before_registration(self):
        (self.package/'C01.smi').write_text('invalid')
        with self.assertRaises(ContractError):self.load()
        self.assertEqual(self.service.list(self.run['id']),[])

    def test_path_escape_and_duplicate_rejected(self):
        for name in ('../outside','C:\\outside','C01.smi:alternate'):
            h=copy.deepcopy(self.handoff);h['artifacts'][0]['path']=name
            (self.package/'handoff.json').write_bytes(encoded(h))
            with self.assertRaises(ContractError):self.load()
        h=copy.deepcopy(self.handoff);h['artifacts'].append(h['artifacts'][0])
        (self.package/'handoff.json').write_bytes(encoded(h))
        with self.assertRaises(ContractError):self.load()

    def test_paper_number_collision_rejected_even_with_valid_hashes(self):
        self.handoff['candidates'][0]['paper_compound_number']='1';self.resign()
        with self.assertRaisesRegex(ContractError,'IDENTITY'):self.load()

    def test_wrong_atom_pair_rejected_even_after_rehashing(self):
        p=self.package/'C02-atom-map.json';value=json.loads(p.read_text())
        for a in value['warhead_mapping']['atoms']:
            if a['source_ccd_atom']=='N18':a['candidate_ccd_atom']='N46'
        p.write_bytes(encoded(value));self.resign()
        with self.assertRaisesRegex(ContractError,'ATTACHMENT'):self.load()

    def test_forged_gpu_or_approval_rejected(self):
        self.handoff['candidates'][0]['prediction']['status']='completed';self.resign()
        with self.assertRaises(ContractError):self.load()

    def test_target_change_same_revision_invalidates_dossier(self):
        value,_=self.load()
        p=self.store.scope(PROJECT);r=self.store.get_run(PROJECT,self.run['id'])
        target=p.json(r['target_ref']);target['scope_reasons'].append('Synthetic changed context')
        ref=p.put_json(target,'TargetResolution','source')
        self.store.update(PROJECT,r['id'],{'target_ref':ref})
        self.assertFalse(self.service.view(value['id'])['current_input'])

    def test_changed_text_version_cannot_use_old_anchors(self):
        plan=json.loads(self.links_path.read_text());plan['anchors'][0]['quote']='PROTAC 1 ( 1 )'
        self.links_path.write_bytes(encoded(plan))
        with self.assertRaisesRegex(ContractError,'ANCHOR'):self.load()

    def test_scoped_read_and_tamper_detection(self):
        value,_=self.load()
        with self.assertRaises(KeyError):DossierService(self.store,'other').view(value['id'])
        ref=value['artifact_refs']['C01.sdf']
        (self.store.root/'blobs'/ref['artifact_id']).write_bytes(b'tampered')
        with self.assertRaises(ContractError):self.service.view(value['id'])

    def test_api_read_only_replay_and_no_import_endpoint(self):
        value,_=self.load()
        app=create_app(self.store.root,enable_worker=False)
        with TestClient(app) as client:
            self.assertEqual(client.get('/api/dossiers/'+value['id']).status_code,200)
            self.assertEqual(len(client.get('/api/runs/'+self.run['id']+'/dossiers').json()),1)
            self.assertEqual(client.post('/api/runs/'+self.run['id']+'/dossiers',json={},headers={'x-tpd-local':'1'}).status_code,405)
            self.assertFalse(client.get('/api/runs/'+self.run['id']+'/g1').json()['dispatch_authorized'])

    def test_fixture_run_cannot_be_imported_as_real_case(self):
        self.store.update(PROJECT,self.run['id'],{'data_mode':'test_fixture'})
        with self.assertRaisesRegex(ContractError,'REAL_COLLECTION'):self.load()
