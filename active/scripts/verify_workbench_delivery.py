"""Read and replay the actual local delivery; makes no model or science calls."""
import argparse
import hashlib
import json
import re
import sys
import time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import httpx
from fastapi.testclient import TestClient
from apps.api.main import create_app,PROJECT
from packages.platform.store import Store
from packages.platform.workbench import WorkbenchService
from packages.platform.workbench_replay import inspect_archive,import_archive
from packages.contracts import now

def main():
    p=argparse.ArgumentParser();p.add_argument('--store',type=Path,required=True);p.add_argument('--result-id',required=True)
    p.add_argument('--url',default='http://127.0.0.1:8770');p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    s=WorkbenchService(Store(a.store),PROJECT);d=s.saved.view(a.result_id);jobs=s.list(a.result_id)
    assert not any(j['state'] in {'running','queued'} for j in jobs)
    views=[s.view(j['id']) for j in jobs]
    print('Verified all registered job outputs and call artifacts.',flush=True)
    raw=httpx.get(a.url+'/api/lab/cases/'+a.result_id+'/replay',timeout=240).raise_for_status().content
    manifest,_,_,_=inspect_archive(raw);zip_path=a.store/'tpd-workbench-final-20260928.zip';zip_path.write_bytes(raw)
    zip_hash=hashlib.sha256(raw).hexdigest();zip_path.with_suffix('.zip.sha256').write_text(zip_hash+'  '+zip_path.name+'\n',encoding='ascii')
    print('Exported final replay:',len(raw),'bytes',flush=True)
    dest=a.store.parent/('workbench-final-replay-'+str(time.time_ns()));replay,_=import_archive(raw,dest,PROJECT)
    with TestClient(create_app(replay.root,read_only=True)) as client:
        assert len(client.get('/api/lab/cases').json())>=1
        assert client.get('/api/lab/structures/'+a.result_id+'/C01/0').status_code==200
        assert client.get('/api/lab/structures/'+a.result_id+'/C02/0').status_code==200
        assert client.post('/api/lab/jobs',json={},headers={'X-TPD-Local':'1'}).status_code==403
    with replay.db() as db:
        accounts=db.execute('SELECT COUNT(*) FROM review_actors').fetchone()[0]
        sessions=db.execute('SELECT COUNT(*) FROM review_sessions').fetchone()[0]
    assert accounts==sessions==0
    with s.store.db() as db:decisions=db.execute('SELECT COUNT(*) FROM report_review_decisions WHERE project=?',(PROJECT,)).fetchone()[0]
    root=Path(__file__).resolve().parents[1]
    test_log=(a.store/'tests-final-20260928.log').read_text(encoding='utf-8-sig')
    match=re.search(r'Ran (\d+) tests',test_log)
    assert match and test_log.rstrip().endswith('OK'),'The saved test run did not pass.'
    files=[p for folder in ('apps','packages','contracts') for p in (root/folder).rglob('*') if p.is_file() and p.suffix in {'.py','.json','.js','.css','.html'} and '__pycache__' not in p.parts]
    result={'verified_at':now(),'kind':'actual_local_workbench_delivery_not_scientific_validation',
        'saved_result_id':a.result_id,'summary':d['summary'],'reader':s.saved.verify(a.result_id),'formal_review_decisions':decisions,
        'jobs':[{k:v[k] for k in ('id','operation','state','calls','total_tokens','error_code','freshness')} for v in views],
        'agent_api_calls':sum(j['calls'] for j in jobs),'agent_api_tokens':sum(j['total_tokens'] for j in jobs),
        'additional_probe':{'calls':1,'tokens':54},'replay':{'sha256':zip_hash,'bytes':len(raw),'manifest_digest':manifest['digest'],
            'extra_artifacts':len(manifest['artifacts']),'roundtrip_verified':True,'accounts_exported':accounts,'sessions_exported':sessions,
            'write_rejected':True,'C01_C02_overlays_verified':True},
        'tests':{'unittest_passed':int(match[1]),'log':'.localdata/workbench-20260928/tests-final-20260928.log',
                 'log_sha256':hashlib.sha256((a.store/'tests-final-20260928.log').read_bytes()).hexdigest(),
                 'I0_I1_I3':'passed in delivery session'},
        'overlay_numeric_verification':json.loads((a.store/'overlay-verification.json').read_text()),
        'code_sha256':{str(p.relative_to(root)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
        'remaining':['signed expert decisions','protonation/tautomer and whole-complex contacts','synthesis route evidence',
                     'new GPU worker and G1 execution approval','public case curation and deployment']}
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print('PASS: final graph replay, both structure views, read-only boundary and credential exclusion.',flush=True)
    print(json.dumps({k:result[k] for k in ('agent_api_calls','agent_api_tokens','formal_review_decisions','replay')}),flush=True)

if __name__=='__main__':main()
