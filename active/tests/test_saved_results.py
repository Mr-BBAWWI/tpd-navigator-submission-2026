"""Technical admission/rollback/opinion tests; synthetic data, no science claims."""
import copy
import io
import json
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

from fastapi.testclient import TestClient
from apps.api.main import create_app, PROJECT
from packages.contracts import ContractError, encoded
from packages.platform.store import Store
from packages.platform.saved_results import SavedResultsService, AUTHORITY
from packages.platform.expert_responses import extract_answers, performance_comparability
from packages.platform.dossiers import digest
from packages.science.evidence_common import schema_check


def document(text='synthetic opinion'):
    data=io.BytesIO()
    with zipfile.ZipFile(data,'w') as z:
        z.writestr('word/document.xml','<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body><w:p><w:r><w:t>Q01 답변</w:t></w:r></w:p><w:tbl><w:tr><w:tc><w:p><w:r><w:t>'+text+'</w:t></w:r></w:p></w:tc></w:tr></w:tbl></w:body></w:document>')
    return data.getvalue()


class SavedResultsTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.store=Store(self.tmp.name);self.svc=SavedResultsService(self.store,PROJECT)

    def test_failed_transaction_has_no_visible_rows_or_orphan_files(self):
        def build(put):
            put(b'first');put(b'second');raise RuntimeError('simulated disk failure')
        with self.assertRaisesRegex(RuntimeError,'simulated'):
            self.svc._atomic('saved_results','dossier_id','synthetic','hash',build,lambda db:None)
        with self.store.db() as db:
            self.assertEqual(db.execute('SELECT COUNT(*) FROM artifacts').fetchone()[0],0)
            self.assertEqual(db.execute('SELECT COUNT(*) FROM saved_results').fetchone()[0],0)
        self.assertEqual(list((self.store.root/'blobs').iterdir()),[])

    def test_duplicate_transaction_does_not_write_more_artifacts(self):
        build=lambda put:{'id':'synthetic','ref':put(b'test')}
        a=self.svc._atomic('saved_results','dossier_id','d','hash',build,lambda db:None)
        with patch.object(self.svc,'_write',side_effect=AssertionError('must not write')):
            b=self.svc._atomic('saved_results','dossier_id','d','hash',build,lambda db:None)
        self.assertEqual(a,('synthetic',True));self.assertEqual(b,('synthetic',False))

    def test_changed_snapshot_checked_before_duplicate_lookup(self):
        self.svc._atomic('saved_results','dossier_id','d','hash',lambda put:{'id':'synthetic'},lambda db:None)
        def changed(db):raise ContractError('stale input')
        with self.assertRaisesRegex(ContractError,'stale input'):
            self.svc._atomic('saved_results','dossier_id','d','hash',lambda put:{},changed)

    def test_read_only_api_and_cross_project_access(self):
        ref=self.store.scope('other').put_raw(b'private')
        doc=self.store.scope(PROJECT).put_raw(document(),'application/vnd.openxmlformats-officedocument.wordprocessingml.document','source')
        with TestClient(create_app(self.store.root,enable_worker=False,read_only=True)) as client:
            self.assertTrue(client.get('/api/health').json()['read_only'])
            self.assertEqual(client.post('/api/workflows',json={},headers={'X-TPD-Local':'1'}).status_code,403)
            self.assertEqual(client.get('/api/artifacts/'+ref['artifact_id']+'/download').status_code,404)
            download=client.get('/api/artifacts/'+doc['artifact_id']+'/download')
            self.assertTrue(download.headers['content-disposition'].endswith('.docx"'))
            self.assertEqual(digest(download.content),doc['sha256'])
            self.assertEqual(client.get('/api/saved-results/unknown').status_code,404)
            self.assertEqual(client.get('/api/saved-results/unknown/files').status_code,404)
            self.assertEqual(client.post('/api/saved-results',json={'path':'C:/private'},headers={'X-TPD-Local':'1'}).status_code,403)

    def test_document_tables_and_duplicate_question_detection(self):
        self.assertEqual(extract_answers(document()),{'Q01':'synthetic opinion'})
        with self.assertRaisesRegex(ContractError,'DUPLICATE_QUESTION'):
            extract_answers(document('Q01 답변'))

    def assessment(self):
        raw=document()
        return raw,{'version':'synthetic','source_sha256':digest(raw),'bundle_digest':'a'*64,
                    'reviewer_name':None,'review_date':None,'disposition_basis':'developer_mapping_of_free_text_not_a_signed_approval',
                    'questions':[{'question_id':'Q01','answer_original':'synthetic opinion','disposition':'opinion_received','applied':[],'remaining':[]}],
                    'literature_supplements':[],'authority':AUTHORITY.copy()}

    def test_opinion_cannot_claim_authority(self):
        raw,a=self.assessment();schema_check(a,'a_expert_response.schema.json')
        for k in ('human_approved','dispatch_authorized','public_release_ready'):
            forged=copy.deepcopy(a);forged['authority'][k]=True
            with self.assertRaises(ValueError):schema_check(forged,'a_expert_response.schema.json')

    def test_wrong_response_source_bundle_or_answer_rejected(self):
        raw,a=self.assessment();source=self.store.root/'answer.docx';source.write_bytes(raw)
        for field,value,error in [('source_sha256','b'*64,'SOURCE_HASH'),('bundle_digest','b'*64,'STALE_BUNDLE'),('answer_original','forged','ANSWER_SOURCE')]:
            bad=copy.deepcopy(a)
            if field=='answer_original':bad['questions'][0][field]=value
            else:bad[field]=value
            with patch.object(self.svc,'_record',return_value={'bundle_digest':'a'*64}):
                with self.assertRaisesRegex(ContractError,error):self.svc.import_review('id',source,bad)

    def test_missing_or_different_assay_conditions_cannot_rank(self):
        a={'target':'SMARCA2','cell_line':'test','exposure_time_h':24,'assay':'test','biological_replicates':2,'Dmax':65,'source_location':'fig1'}
        b=copy.deepcopy(a);b['exposure_time_h']=None
        self.assertEqual(performance_comparability(a,b)['status'],'conditions_missing')
        b['exposure_time_h']=8
        self.assertEqual(performance_comparability(a,b)['status'],'conditions_differ')
        r=performance_comparability(a,a)
        self.assertEqual(r['status'],'matched_conditions_pending_review');self.assertFalse(r['automatic_efficacy_ranking'])


if __name__=='__main__':unittest.main()
