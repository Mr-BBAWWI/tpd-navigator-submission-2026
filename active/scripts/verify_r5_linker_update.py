"""Record actual R5 CPU runs against an existing local case; never runs LLM/GPU."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from apps.api.main import PROJECT,create_app
from packages.platform.store import Store
from packages.platform.workbench import WorkbenchService
from packages.platform.saved_results import SavedResultsService
from packages.agents.saved_bundle import build_saved_bundle
from packages.science.linker_design import library
from fastapi.testclient import TestClient


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--store',type=Path,required=True)
    parser.add_argument('--response',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    store=Store(args.store);jobs=WorkbenchService(store,PROJECT);saved=SavedResultsService(store,PROJECT)
    with store.db() as db:
        rows=db.execute('SELECT id FROM saved_results WHERE project=? ORDER BY rowid DESC',(PROJECT,)).fetchall()
        assert len(rows)==1,'Choose a single-case store explicitly.'
        assert not db.execute("SELECT id FROM workbench_jobs WHERE project=? AND state IN ('queued','running')",(PROJECT,)).fetchone(),'A job is active.'
    result_id=rows[0]['id']
    previous={j['id']:hashlib.sha256(json.dumps(j,sort_keys=True).encode()).hexdigest() for j in jobs.list(result_id)}
    assessment=json.loads((ROOT/'cases/expert_response_20260930.json').read_text(encoding='utf-8'))
    review,fresh=saved.import_review(result_id,args.response,assessment)
    results=[]
    settings=[{'parent':'C01','linker_id':t['id'],'orientation':'forward'} for t in library()['templates']]
    settings.append({'parent':'C01','linker_id':'peg_alkyl','orientation':'reverse'})
    for parameters in settings:
        key='r5-20260930.1-'+parameters['parent']+'-'+parameters['linker_id']+'-'+parameters['orientation']
        job,created=jobs.create(result_id,'candidate',key,parameters)
        if created:jobs.execute(job['id'])
        value=jobs.view(job['id'])
        assert value['state']=='completed',(job['id'],value.get('error_code'))
        result=value['result']
        for output in value['outputs']:store.scope(PROJECT).read(output['ref'])
        assert not result['final_candidate'] and not result['authority']['human_approved']
        results.append({'job_id':job['id'],'candidate_id':result['candidate_id'],'parameters':parameters,
                        'identity':result['identity'],'comparison':result['linker_comparison'],
                        'protected_atom_count':len(result['transformation']['protected_atom_maps']),
                        'files':result['files']})
    report,created=jobs.create(result_id,'report','r5-20260930.1-report',{})
    if created:jobs.execute(report['id'])
    report=jobs.view(report['id']);assert report['state']=='completed',report
    supplement=next(x for x in report['outputs'] if x['name']=='CPU-실험-부록.md')
    text=store.scope(PROJECT).read(supplement['ref']).decode()
    assert '초기 linker perturbation 기능 검증용 prototype' in text
    assert 'Linker 설계 가설' in text
    view=saved.view(result_id)
    bundle=build_saved_bundle(view,saved.dossiers.view(view['dossier_id']),jobs.tool_snapshot(result_id))
    source=bundle['sources']['tools:cpu']['content']['prototypes']
    assert {p['design_class'] for p in source}>={'linker_perturbation_prototype','medchem_linker_hypothesis'}
    assert all(p['final_candidate'] is False for p in source)
    now={j['id']:hashlib.sha256(json.dumps(j,sort_keys=True).encode()).hexdigest() for j in jobs.list(result_id)}
    assert all(now[k]==v for k,v in previous.items()),'Existing job records changed.'
    with TestClient(create_app(args.store,enable_worker=False,read_only=True)) as client:
        catalog=client.get('/api/lab/linker-library');assert catalog.status_code==200
        assert len(catalog.json()['templates'])==6
        page=client.get('/');assert page.status_code==200 and 'design-linker' in page.text
        blocked=client.post('/api/lab/jobs',headers={'X-TPD-Local':'1'},json={
            'result_id':result_id,'operation':'candidate','request_key':'r5-readonly-rejection',
            'parameters':{'parent':'C01','linker_id':'peg3'}})
        assert blocked.status_code==403
    with store.db() as db:
        decisions=db.execute('SELECT count(*) FROM report_review_decisions WHERE project=?',(PROJECT,)).fetchone()[0]
    from packages.contracts import now as timestamp
    audit={'verified_at':timestamp(),'kind':'actual_CPU_linker_design_not_scientific_validation',
           'result_id':result_id,'opinion_id':review['id'],'opinion_created':fresh,
           'source_sha256':assessment['source_sha256'],'expert_opinion_count':len(view['expert_reviews']),
           'formal_decisions':decisions,'old_job_records_preserved':len(previous),
           'template_count':6,'generated_records':len(results),
           'unique_structures':len({r['candidate_id'] for r in results}),
           'candidate_results':results,'report_job_id':report['id'],
           'report_has_prototype_and_linker_hypothesis_labels':True,
           'AI_input_preserves_classification':True,'read_only_candidate_write_rejected':True,
           'baseline':view['summary'],'new_API_calls':0,'new_GPU_runs':0,
           'synthesis_route_search_performed':False,'full_complex_contacts_recomputed':False,
           'runtime':jobs.get(report['id'])['runtime']}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in audit.items() if k not in {'candidate_results','runtime','baseline'}},ensure_ascii=False))


if __name__=='__main__':main()
