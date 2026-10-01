"""Verify an isolated actual saved store. Mutated test copies never replace source data."""
import argparse
import copy
import json
import sys
import tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from fastapi.testclient import TestClient
from apps.api.main import create_app, PROJECT
from packages.contracts import ContractError, encoded, now
from packages.platform.saved_results import SavedResultsService
from packages.platform.store import Store
from packages.platform.dossiers import digest
from packages.science.evidence_common import seal


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir',type=Path,required=True);p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--result-id',required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();store=Store(a.data_dir);svc=SavedResultsService(store,PROJECT)
    result=svc.view(a.result_id);verified=svc.verify(a.result_id)
    before_blobs=set((store.root/'blobs').iterdir())
    duplicate,fresh=svc.import_archive(result['dossier_id'],a.archive)
    assert not fresh and duplicate['id']==result['id'] and set((store.root/'blobs').iterdir())==before_blobs
    assert result['summary']['candidate_count']==2 and result['summary']['sample_count']==10
    assert [len(c['samples']) for c in result['candidates']]==[6,4]
    assert len(result['literature']['records'])==10
    assert all(c['synthesis_assessment']['status']=='not_assessed' for c in result['candidates'])
    assert len(result['expert_reviews'][0]['assessment']['questions'])==5
    ref_downloads=0
    with TestClient(create_app(store.root,enable_worker=False,read_only=True)) as client:
        assert client.get('/saved-review/'+result['dossier_id']).status_code==200
        response=client.get('/api/saved-results/'+result['id']);assert response.status_code==200
        for c in response.json()['candidates']:
            for sample in c['samples']:
                for ref in sample['reports'].values():
                    raw=client.get('/api/artifacts/'+ref['artifact_id']+'/download')
                    assert raw.status_code==200 and digest(raw.content)==ref['sha256'];ref_downloads+=1
        ref=result['expert_reviews'][0]['source_ref']
        assert digest(client.get('/api/artifacts/'+ref['artifact_id']+'/download').content)==ref['sha256']
        assert client.get('/api/saved-results/'+result['id']+'/files').status_code==200
        assert client.post('/api/saved-results',json={'path':'../anywhere'},headers={'X-TPD-Local':'1'}).status_code==403
        # Tamper only a blob in the explicitly isolated verification store; restore in finally.
        ref=svc.files(result['id'])['README.md'];path=store.root/'blobs'/ref['artifact_id'];original=path.read_bytes()
        try:
            path.write_bytes(b'technical corruption test')
            assert client.get('/api/artifacts/'+ref['artifact_id']+'/download').status_code==409
            try:svc.verify(result['id']);raise AssertionError('accepted corrupted storage')
            except ContractError:pass
        finally:path.write_bytes(original)
    # Reject invalid top-level archives before registering anything. Small test copies
    # exercise the real reader's first admission checks; detailed binding covered by core tests.
    errors={}
    raw_manifest=(a.archive/'manifest.json').read_bytes();manifest=json.loads(raw_manifest)
    with tempfile.TemporaryDirectory(prefix='tpd-rejected-') as temp:
        folder=Path(temp)
        for label in ('missing','extra','unsafe_directory'):
            for f in folder.iterdir():f.unlink()
            m=copy.deepcopy(manifest)
            if label=='unsafe_directory':
                # Empty indexed set reaches path validation without copying the 300 MB source tree.
                m['files']={};m['preparation_directories']=['../outside'];m=seal({k:v for k,v in m.items() if k!='digest'})
            (folder/'manifest.json').write_bytes(encoded(m))
            if label=='extra':(folder/'unindexed.txt').write_bytes(b'test')
            try:svc.import_archive(result['dossier_id'],folder);raise AssertionError(label+' accepted')
            except ContractError as exc:errors[label]=str(exc)
    assert set((store.root/'blobs').iterdir())==before_blobs
    run=store.get_run(PROJECT,result['run_id'])
    try:
        store.update(PROJECT,run['id'],{'revision':run['revision']+1})
        assert not svc.view(result['id'])['current_input']
        try:svc.import_archive(result['dossier_id'],a.archive);raise AssertionError('stale import accepted')
        except ContractError as exc:errors['stale_input']=str(exc)
    finally:
        with store.db() as db:
            db.execute('UPDATE runs SET revision=?,state=?,body=? WHERE id=? AND project=?',
                       (run['revision'],run['state'],json.dumps(run),run['id'],PROJECT))
    svc.verify(result['id'])
    record={'verified_at':now(),'result_id':result['id'],'data_dir':str(store.root),'bundle_digest':result['bundle_digest'],
            'summary':result['summary'],**verified,'report_downloads_verified':ref_downloads,'expert_questions':5,
            'duplicate_import_no_new_blobs':True,'stored_tampering_rejected':True,'rejections':errors,
            'original_archive_modified':False,'new_llm_calls':0,'new_gpu_runs':0,'authority':result['authority']}
    a.output.write_bytes(encoded(record));print(json.dumps(record,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
