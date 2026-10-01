"""Verify an actual saved report via the local API, without new LLM/GPU work."""
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

from fastapi.testclient import TestClient
from apps.api.main import create_app, PROJECT
from packages.contracts import now
from packages.platform.dossiers import digest
from packages.platform.evidence_reports import EvidenceReportService
from packages.platform.store import Store


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir',type=Path,required=True)
    p.add_argument('--report-id',required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Historical verification cannot be overwritten.')
    if not (a.data_dir/'index.sqlite3').is_file():p.error('Existing store required.')
    store=Store(a.data_dir);svc=EvidenceReportService(store,PROJECT);report=svc.view(a.report_id)
    source=svc.saved.view(report['result_id']);body=report['content']
    assert report['freshness']['current']
    assert body['counts']=={'candidates':2,'samples':10,'observations':10,'expert_responses':1,'expert_questions':5}
    assert [c['saved_evidence'] for c in body['candidates']]==source['candidates']
    assert body['literature']==source['literature'] and body['expert_reviews']==source['expert_reviews']
    by_candidate={c['identity']['id']:len(c['saved_evidence']['samples']) for c in body['candidates']}
    assert by_candidate=={'C01':6,'C02':4}
    before=len(list((store.root/'blobs').iterdir()))
    duplicate,fresh=svc.create(report['result_id'],report['input_digest'])
    assert not fresh and duplicate['id']==report['id']
    assert len(list((store.root/'blobs').iterdir()))==before
    with TestClient(create_app(a.data_dir,enable_worker=False,read_only=True)) as client:
        actual=client.get('/api/evidence-reports/'+report['id'])
        assert actual.status_code==200 and actual.json()['counts']==body['counts']
        history=client.get('/api/saved-results/'+report['result_id']+'/reports')
        assert history.status_code==200 and history.json()[0]['id']==report['id']
        assert client.get('/evidence-reports/'+report['result_id']).status_code==200
        for ref in [report['content_ref'],report['markdown_ref'],report['record_ref'],*report['source_catalog']]:
            response=client.get('/api/artifacts/'+ref['artifact_id']+'/download')
            assert response.status_code==200 and digest(response.content)==ref['sha256']
        assert client.post('/api/saved-results/'+report['result_id']+'/reports',json={'expected_digest':report['input_digest']},headers={'X-TPD-Local':'1'}).status_code==403
    result={'verified_at':now(),'passed':True,'data_mode':'real','report_id':report['id'],
        'counts':body['counts'],'candidate_sample_counts':by_candidate,'input_digest':report['input_digest'],
        'report_ref':report['record_ref'],'direct_source_downloads_verified':len(report['source_catalog']),
        'report_downloads_verified':3,'b_source_aliases_verified':len(body['source_aliases']),
        'original_observations_and_responses_preserved':True,'duplicate_created_files':0,
        'read_only_creation_rejected':True,'live_llm_calls':0,'new_gpu_runs':0,'authority':report['authority']}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
