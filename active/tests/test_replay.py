"""Closed-graph replay tests, synthetic text and committed CPU fixtures only."""
import io
import json
from pathlib import Path
import zipfile
import unittest

from fastapi.testclient import TestClient
from apps.api.main import create_app,PROJECT
from packages.contracts import ContractError
from packages.platform.replay import export_bytes,import_archive,inspect_archive
import test_candidate_dossiers as dossier_fixture
import test_literature_service as literature_fixture
from packages.platform.analysis import AnalysisService


class ReplayTests(unittest.TestCase):
    def setUp(self):
        self.fixture=dossier_fixture.DossierTests();self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.dossier,_=self.fixture.load()
        self.raw,self.manifest=export_bytes(self.fixture.store,PROJECT,self.fixture.run['id'])

    def test_fresh_store_read_only_roundtrip_and_graph(self):
        store,manifest=import_archive(self.raw,self.fixture.root/'replay',PROJECT)
        self.assertEqual(manifest['digest'],self.manifest['digest'])
        app=create_app(store.root,enable_worker=False)
        with TestClient(app) as client:
            response=client.get('/api/dossiers/'+self.dossier['id']);self.assertEqual(response.status_code,200)
            self.assertEqual(response.json()['candidates'],self.dossier['candidates'])
            self.assertTrue(client.get('/api/runs/'+self.fixture.run['id']).json()['replay_only'])
            self.assertFalse(client.get('/api/runs/'+self.fixture.run['id']).json()['llm_enabled'])

    def altered_zip(self, change):
        with zipfile.ZipFile(io.BytesIO(self.raw)) as archive: files={n:archive.read(n) for n in archive.namelist()}
        change(files);output=io.BytesIO()
        with zipfile.ZipFile(output,'w') as archive:
            for n,raw in files.items():archive.writestr(n,raw)
        return output.getvalue()

    def test_modified_and_missing_blob_rejected_before_destination(self):
        for missing in (False,True):
            def change(files):
                key=next(n for n in files if n.startswith('blobs/'))
                if missing:del files[key]
                else:files[key]=b'changed'
            dest=self.fixture.root/('bad-'+str(missing))
            with self.assertRaises((ContractError,KeyError)):import_archive(self.altered_zip(change),dest,PROJECT)
            self.assertFalse(dest.exists())

    def test_path_traversal_and_extra_files_rejected(self):
        for name in ('../escape','blobs/../../escape','unexpected.txt'):
            with self.assertRaises(ContractError):inspect_archive(self.altered_zip(lambda files:files.update({name:b'x'})))

    def test_existing_store_and_other_project_rejected(self):
        with self.assertRaises(ContractError):import_archive(self.raw,self.fixture.store.root,PROJECT)
        with self.assertRaises(ContractError):import_archive(self.raw,self.fixture.root/'other','other')

    def test_closed_graph_contains_all_candidate_files(self):
        manifest,blobs=inspect_archive(self.raw)
        for ref in self.dossier['artifact_refs'].values():
            self.assertIn(ref['artifact_id'],blobs)
        self.assertFalse(manifest['human_approval_transferred'])
        self.assertFalse(manifest['public_release_ready'])

    def test_completed_analysis_survives_fresh_import_without_app_initialization(self):
        fixture=literature_fixture.AnalysisTests();fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        completed=fixture.execute()
        self.assertEqual(completed['state'],'review_ready')
        # Deliberately sort the request key before the older one. Export must retain
        # insertion order rather than relying on the database's chosen index scan.
        later=fixture.start(key='aaa-newer-request');fixture.service.execute(later['id'])
        raw,_=export_bytes(fixture.store,PROJECT,fixture.run['id'])
        store,_=import_archive(raw,self.fixture.root/'analysis-replay',PROJECT)
        restored=AnalysisService(store,PROJECT).view(completed['id'])
        self.assertEqual(restored['assessment'],fixture.service.view(completed['id'])['assessment'])
        self.assertEqual(restored['review'],fixture.service.view(completed['id'])['review'])
        self.assertEqual([r['id'] for r in AnalysisService(store,PROJECT).list(fixture.run['id'])],
                         [r['id'] for r in fixture.service.list(fixture.run['id'])])
