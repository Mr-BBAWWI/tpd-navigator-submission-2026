"""Explicit live development probe, fixed prompts/labels, saved before the first call.

Does not mutate research state. No answer labels enter the agent input. No automatic
prompt evolution against this dataset; a new protocol/version is needed for changes.
"""
import argparse
import json
import random
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from packages.agents.provider import DaconProvider,load_key
from packages.agents.review_harness import ClaimReviewHarness,PROFILE_PATH,SCHEMA_PATH
from packages.platform.dossiers import DossierService,digest
from packages.platform.store import Store
from packages.contracts import encoded,parse_json,now


def evidence_cards(dossier, port):
    handoff=port.json(dossier['handoff_ref'])
    observations=port.json(dossier['artifact_refs']['literature-evidence.json'])['observations']
    start=handoff['starting_ligand']
    cards=[{'id':'START:identity','kind':'B_cpu_and_curated_paper_identity','value':{k:start[k] for k in ('paper_name','paper_compound_number','ccd','pdb','attachment_atom')},'source_sha256':dossier['handoff_ref']['sha256']}]
    for c in dossier['candidates']:
        cards.append({'id':c['id']+':identity','kind':'provisional_A_B_identity_link','value':{k:c[k] for k in ('id','paper_name','paper_compound_number','pdb','ccd','origin')},'source_sha256':dossier['handoff_ref']['sha256']})
        cards.append({'id':c['id']+':status','kind':'recorded_workflow_status','value':c['stages'],'source_sha256':dossier['dossier_ref']['sha256']})
        cards.append({'id':c['id']+':assay','kind':'manual_literature_transcription_pending_pharmacy_review','value':next(o for o in observations if o['candidate_id']==c['id']),'source_sha256':dossier['artifact_refs']['literature-evidence.json']['sha256']})
    cards.append({'id':'C03:status','kind':'recorded_workflow_status','value':dossier['planned_candidates'][0],'source_sha256':dossier['dossier_ref']['sha256']})
    return cards


def score(result, case):
    correct=result['status']=='reviewed_pending_human' and result['verdict']==case['expected']
    cited={e for r in result['reviews'] for e in r['evidence_ids']}
    return {'decision_correct':correct,'relevant_evidence_cited':bool(cited.intersection(case['relevant_evidence'])),
        'false_accept':case['expected']!='supported' and result['verdict']=='supported',
        'false_alarm':case['expected']=='supported' and result['verdict']!='supported'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir',type=Path,required=True)
    parser.add_argument('--dossier-id',required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--repeats',type=int,default=2,choices=range(1,6))
    parser.add_argument('--live',action='store_true',help='Required acknowledgement of actual API calls')
    args=parser.parse_args()
    if not args.live:parser.error('--live is required; offline tests do not call providers')
    key=load_key()
    store=Store(args.data_dir);port=store.scope('local-research')
    dossier=DossierService(store,port.project).view(args.dossier_id)
    if not dossier['current_input']:raise SystemExit('Stale dossier; evaluation was not started')
    dataset_path=ROOT/'cases/harness/claim_challenges.json'
    dataset=parse_json(dataset_path.read_bytes())
    facts=evidence_cards(dossier,port)
    profiles=list(parse_json(PROFILE_PATH.read_bytes())['profiles'])
    plan=[(p,c,r) for r in range(args.repeats) for c in dataset['cases'] for p in profiles]
    random.Random(20260923).shuffle(plan)
    args.output.mkdir(parents=True,exist_ok=False)
    def save(name,value):
        with (args.output/name).open('xb') as f:f.write(encoded(value))
    manifest={'created_at':now(),'mode':'live_api_synthetic_claims_over_real_saved_evidence','model':'gpt-5.6-sol',
        'source_dossier_ref':dossier['dossier_ref'],'dataset':dataset,'evidence':facts,
        'plan':[{'profile':p,'case_id':c['id'],'repeat':r} for p,c,r in plan],
        'frozen_files':{str(p.relative_to(ROOT)):digest(p.read_bytes()) for p in [dataset_path,PROFILE_PATH,SCHEMA_PATH,Path(__file__),ROOT/'packages/agents/review_harness.py']},
        'limits':{'max_calls_per_run':2,'max_tokens_per_run':40000,'experiment_observed_token_cap':600000},
        'limitations':['Development-set probe, not held-out or pharmacy-expert labelled.','Profiles share model/evidence/output cap and maximum run budget; actual call/token costs differ.',
            'Blind pair holds disagreements; no majority-vote truth or clinical efficacy claim.','Scorer checks verdict and cited IDs; rationale entailment still needs human evaluation.']}
    save('protocol.json',manifest)
    results=[];tokens=0;calls=0
    for profile,case,repeat in plan:
        if tokens+40000>600000:break
        provider=DaconProvider(key)
        try:
            folder=f'{profile}-{case["id"]}-{repeat}'
            runtime=ClaimReviewHarness(provider,args.output/folder)
            result=runtime.run(profile,{'proposition':case['proposition'],'context_note':case['context_note'],'evidence':facts})
        finally:provider.close()
        tokens+=result['total_tokens'];calls+=result['calls']
        row={'case_id':case['id'],'profile':profile,'repeat':repeat,'directory':folder,'status':result['status'],
            'verdict':result['verdict'],'expected':case['expected'],'calls':result['calls'],'tokens':result['total_tokens'],
            'seconds':result['elapsed_seconds'],'score':score(result,case),'error_code':result.get('error_code')}
        results.append(row);save(f'record-{len(results):03}.json',row)
        print(json.dumps({'completed':len(results),'planned':len(plan),**row},ensure_ascii=False),flush=True)
        # An uncertain provider charge or outage stops the experiment; no automatic retry storm.
        if result['usage_status']=='unreconciled':break
    aggregates=[]
    for p in profiles:
        rows=[r for r in results if r['profile']==p]
        aggregates.append({'profile':p,'attempted':len(rows),'decision_correct':sum(r['score']['decision_correct'] for r in rows),
            'relevant_evidence_cited':sum(r['score']['relevant_evidence_cited'] for r in rows),
            'false_accepts':sum(r['score']['false_accept'] for r in rows),'false_alarms':sum(r['score']['false_alarm'] for r in rows),
            'held':sum(r['status']=='held' for r in rows),'calls':sum(r['calls'] for r in rows),
            'tokens':sum(r['tokens'] for r in rows),'seconds':round(sum(r['seconds'] for r in rows),3)})
    summary={'status':'completed' if len(results)==len(plan) else 'incomplete','calls':calls,'tokens':tokens,
        'attempted':len(results),'planned':len(plan),'aggregates':aggregates,'results':results,'limitations':manifest['limitations']}
    save('summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k!='results'},ensure_ascii=False),flush=True)
    return 0 if summary['status']=='completed' else 1


if __name__=='__main__':raise SystemExit(main())
