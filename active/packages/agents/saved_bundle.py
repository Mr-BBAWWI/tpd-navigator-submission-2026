"""Bounded agent input from registered predictions, literature and expert opinions."""
from packages.agents.bundle import json_digest
from packages.platform.saved_results import AUTHORITY


VERSION = 'saved-agent-input/20260928.1'


def build_saved_bundle(data, dossier, tool_inputs=()):
    facts, sources = [], {}
    latest = data['expert_reviews'][-1]['assessment'] if data['expert_reviews'] else None
    opinions = {q['question_id']:q['answer_original'] for q in latest['questions']} if latest else {}
    sources['expert:literature'] = {'kind':'received_expert_opinion_not_signed_approval',
        'content':{'source_sha256':latest['source_sha256'] if latest else None,
                   'Q04':opinions.get('Q04'), 'authority':AUTHORITY}}
    sources['expert:structure'] = {'kind':'received_expert_opinion_not_signed_approval',
        'content':{'source_sha256':latest['source_sha256'] if latest else None,
                   **{q:opinions.get(q) for q in ('Q01','Q02','Q03','Q05')},'authority':AUTHORITY}}
    for candidate in data['candidates']:
        cid=candidate['compound_id']
        records=[r['record'] for r in data['literature']['records'] if r['record'].get('compound_id')==cid]
        # The exchange keeps the candidate name in query when no explicit compound_id exists.
        if not records:
            name=candidate['facts']['compound']['paper']['name']
            records=[r['record'] for r in data['literature']['records'] if r['record'].get('query',{}).get('paper_name')==name]
        compact=[]
        for r in records:
            observation=r['observation'];original=observation.get('A_original',{})
            compact.append({'record_id':r['record_id'],'query':r['query'],'source':r['source'],
                'observation':{k:v for k,v in observation.items() if k not in {'A_original','A_observation'}},
                'original_measurement':observation.get('A_observation'),
                'original_assay_context':original.get('assay_context'),
                'original_target':original.get('target'),'association':original.get('association')})
        records=compact
        sid='literature:'+cid
        sources[sid]={'kind':'registered_literature_observation_with_original_conditions','content':records}
        facts.append({'id':cid+':literature','candidate_id':cid,'domain':'literature','kind':'literature_annotation',
            'value':{'observation_count':len(records),'coverage':'Registered extracted passages only; not a whole-paper or SI absence claim.'},
            'source_ids':[sid,'expert:literature']})
        samples=[]
        for s in candidate['samples']:
            samples.append({k:s[k] for k in ('run_id','rank','condition_digest','summary','parts_rmsd_ranges_A',
                'status','bad_overlap_pairs') if k in s})
            samples[-1]['role_disagreement_atoms']=sorted({r['name'] for r in s.get('role_disagreements',[])})
            samples[-1]['settings']=s.get('reported_execution',{}).get('settings',{})
        sid='structure:'+cid
        sources[sid]={'kind':'registered_B_saved_prediction_and_quality','content':{'molecule_id':candidate['molecule_id'],
            'samples':samples,'synthesis_assessment':candidate.get('synthesis_assessment'),
            'interpretation':'Known reference reproduction only. OFF and ON are different conditions; contacts preliminary.'}}
        facts.append({'id':cid+':structure','candidate_id':cid,'domain':'molecule','kind':'computed_descriptive_geometry',
            'value':{'sample_count':len(samples),'potentials_off':sum(s['settings'].get('use_potentials') is False for s in samples),
                     'potentials_on':sum(s['settings'].get('use_potentials') is True for s in samples),'human_approved':False},
            'source_ids':[sid,'expert:structure']})
    sources['scope:C03']={'kind':'M2_scope_and_unfinished_work','content':{
        'planned_candidates':dossier.get('planned_candidates',[]),'authority':AUTHORITY,
        'notice':'C03 CPU prototypes, if any, are separately versioned tool results, not evidence in this frozen review.'}}
    facts.append({'id':'C03:scope','candidate_id':'C03','domain':'molecule','kind':'not_run',
        'value':{'new_candidate_in_this_bundle':False,'efficacy':'not_established','synthesis':'not_assessed'},
        'source_ids':['scope:C03','expert:structure']})
    if tool_inputs:
        prototypes=[{k:v for k,v in x['result'].items() if k not in {'fragments','parent_identity'}}
                    for x in tool_inputs if x['operation']=='candidate']
        from packages.science.linker_design import design_class
        for prototype in prototypes:
            prototype['design_class']=design_class(prototype)
            prototype['final_candidate']=False
        trials=[{k:v for k,v in x['result'].items() if k not in {'nitrogen_inventory','binding'}}
                for x in tool_inputs if x['operation']=='chemistry']
        # One structured source keeps the source-reading budget independent of experiment count.
        sources['tools:cpu']={'kind':'recorded_CPU_hypotheses_and_ligand_only_trials','content':{
            'prototypes':prototypes,'hydrogen_trials':trials,'authority':AUTHORITY}}
        if prototypes:
            facts[-1]={'id':'C03:scope','candidate_id':'C03','domain':'molecule','kind':'design_hypothesis',
                'value':{'prototype_count':sum(p['design_class']=='linker_perturbation_prototype' for p in prototypes),
                         'linker_hypothesis_count':sum(p['design_class']=='medchem_linker_hypothesis' for p in prototypes),
                         'new_candidate_in_this_bundle':True,
                         'efficacy':'not_established','synthesis':'not_assessed','bound_structure':'not_run'},
                'source_ids':['tools:cpu','expert:structure']}
        if trials:
            facts.append({'id':'CPU:hydrogen','candidate_id':'C01/C02','domain':'molecule','kind':'chemistry_trial',
                'value':{'trial_count':len(trials),'full_complex_contacts':'not_run'},'source_ids':['tools:cpu','expert:structure']})
        panels=[x['result'] for x in tool_inputs if x['operation']=='design_panel']
        if panels:
            sources['tools:design']={'kind':'actual_CPU_design_docking_dual_E3_hypotheses_not_approval','content':[
                {k:p[k] for k in ('candidate_id','summary','warhead','calibration','review_questions','limitations')} for p in panels]}
            facts.append({'id':'DESIGN:panel','candidate_id':'SMARCA2-dual-E3','domain':'molecule','kind':'design_hypothesis',
                'value':{'panel_count':len(panels),'summaries':[p['summary'] for p in panels],
                    'ternary_structure':'not_run','efficacy':'not_established','synthesis':'not_assessed'},
                'source_ids':['tools:design','expert:structure']})
    value={'version':VERSION,'saved_result_ref':data['record_ref'],'dossier_ref':data['dossier_ref'],
        'expert_review_refs':[r['record_ref'] for r in data['expert_reviews']], 'tool_inputs':[
            {k:v for k,v in x.items() if k!='result'} for x in tool_inputs], 'facts':facts,'sources':sources,
        'coverage':'Registered A observations, saved B predictions and quality summaries, and received expert opinions. '
                   'No new GPU, protonation prediction, protein hydrogen optimization or synthesis validation in this review.',
        'not_completed':['미확인 원문/보충자료 조건 보완','화학 상태·단백질 수소·접촉의 최종 검토',
                         '신규 후보의 과학 검토·실제 GPU 계산','정식 실행 승인과 효능·합성 검증']}
    value['input_digest']=json_digest(value)
    return value
