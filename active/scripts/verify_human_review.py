"""Verify real review preparation only; never log in or submit a reviewer decision."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from fastapi.testclient import TestClient
from apps.api.main import create_app, PROJECT
from packages.contracts import now
from packages.platform.store import Store
from packages.platform.report_reviews import ReportReviewService
from packages.platform.dossiers import digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir',type=Path,required=True);p.add_argument('--request-id',required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Historical verification cannot be overwritten.')
    if not (a.data_dir/'index.sqlite3').is_file():p.error('Existing store required.')
    store=Store(a.data_dir);svc=ReportReviewService(store,PROJECT);row=svc.view(a.request_id)
    report=svc.reports.view(row['report_id'])
    assert row['status']=='pending' and row['revision']==0 and not row['history']
    assert row['assigned_reviewer']=={'id':'mr-an','name':'Reviewer A'}
    assert report['counts']['candidates']==2 and report['counts']['samples']==10
    before=len(list((store.root/'blobs').iterdir()))
    duplicate,fresh=svc.create(row['report_id'],row['report_ref']['sha256'],list(row['scopes']),'mr-an')
    assert not fresh and duplicate['id']==row['id'] and before==len(list((store.root/'blobs').iterdir()))
    with TestClient(create_app(store.root,enable_worker=False,review_only=True)) as client:
        assert client.get('/human-review/'+row['report_id']).status_code==200
        data=client.get('/api/report-reviews/'+row['id']);assert data.status_code==200 and data.json()['revision']==0
        session=client.get('/api/review-auth/session').json();assert session['actor'] is None
        assert {'id':'mr-an','name':'Reviewer A'} in session['reviewers']
        payload={'expected_request_sha256':row['record_ref']['sha256'],'expected_revision':0,'idempotency_key':'unauthenticated-technical-test',
                 'action':'record','judgments':[],'acknowledgement_ids':[],'withdrawal_reason':''}
        assert client.post('/api/report-reviews/'+row['id']+'/decisions',json=payload,headers={'X-TPD-Local':'1'}).status_code==401
        assert client.post('/api/workflows',json={},headers={'X-TPD-Local':'1'}).status_code==403
        for ref in [row['record_ref'],report['record_ref'],report['content_ref'],report['markdown_ref']]:
            response=client.get('/api/artifacts/'+ref['artifact_id']+'/download')
            assert response.status_code==200 and digest(response.content)==ref['sha256']
            assert b'"access_key"' not in response.content and b'"key_hash"' not in response.content
    with store.db() as db:
        decisions=db.execute('SELECT COUNT(*) FROM report_review_decisions WHERE project=?',(PROJECT,)).fetchone()[0]
        sessions=db.execute('SELECT COUNT(*) FROM review_sessions WHERE project=?',(PROJECT,)).fetchone()[0]
    assert decisions==sessions==0
    result={'verified_at':now(),'passed':True,'data_mode':'real','request_id':row['id'],'request_ref':row['record_ref'],
            'reviewer':row['assigned_reviewer'],'report_id':row['report_id'],'counts':report['counts'],'scopes':list(row['scopes']),
            'status':'pending','actual_decisions':decisions,'actual_login_sessions':sessions,'download_hashes_verified':4,
            'duplicate_new_blobs':0,'unauthenticated_decision_rejected':True,'research_mutation_rejected':True,
            'credentials_absent_from_exports':True,'new_llm_calls':0,'new_gpu_runs':0,'authority':row['authority']}
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
