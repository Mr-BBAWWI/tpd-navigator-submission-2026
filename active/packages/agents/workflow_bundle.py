"""Read only the registered A/B evidence snapshot, without raw-cache dependencies."""
from packages.agents.bundle import json_digest
from packages.platform.analysis import AnalysisService
from packages.platform.dossiers import DossierService, require
from packages.platform.handoffs import HandoffService

VERSION = 'service-evidence-bundle-20260924.2'


def build_service_bundle(store, project, input_id):
    service = HandoffService(store, project)
    value = service.view(input_id)
    require(value['current_input'], 'WORKFLOW_INPUT_STALE')
    port = store.scope(project)
    dossier = DossierService(store, project).view(value['dossier_id'])
    material = port.json(value['files']['material.json'])
    analyses = AnalysisService(store, project)
    facts, sources = [], {}
    for candidate in dossier['candidates']:
        cid = candidate['id']
        require(candidate['literature_evidence'], 'CANDIDATE_LITERATURE_REQUIRED')
        records, passages, reviews = [], {}, {}
        for e in candidate['literature_evidence']:
            records.append({k:e[k] for k in ('evidence_id','statement','target','basis_kind','observations','assay_context','review')})
            a = analyses.view(e['analysis_id'])
            require(a['state']=='review_ready' and a['current_input'], 'WORKFLOW_ANALYSIS_NOT_CURRENT')
            reviews[e['analysis_id']] = a['review']['findings']
            for loc in e['source_locators']:
                segment = next(s for s in a['annotated_bundle']['segments'] if s['segment_id']==loc['block_id'])
                passages[segment['segment_id']] = {'locator':loc.get('xpath'), 'text':port.read(segment['text_artifact_ref']).decode('utf-8'),
                    'sha256':segment['text_artifact_ref']['sha256']}
        sid = 'A:'+cid
        sources[sid] = {'kind':'structured_A_extraction', 'content':{'paper_name':candidate['paper_name'],
            'paper_compound_number':candidate['paper_compound_number'],'records':records,
            'previous_literature_critic':reviews}}
        source_ids=[sid]
        for segment_id, passage in passages.items():
            paper_id='paper:'+segment_id
            sources[paper_id]={'kind':'original_article_passage',**passage}
            source_ids.append(paper_id)
        require(len(source_ids)<=4,'WORKFLOW_CITATION_SCOPE_TOO_LARGE')
        facts.append({'id':cid+':literature','candidate_id':cid,'domain':'literature','kind':'literature_annotation',
            'value':{'records':records,'human_review':'pending'},'source_ids':source_ids})
        compound = next(c for c in material['h1']['compounds'] if c['compound_id']==cid)
        sid = 'B:'+cid
        reconstruction = next(c for c in material['h2']['reference_reconstructions'] if c['compound_id']==cid)
        sources[sid] = {'kind':'B_registered_H1_H2', 'content':{'compound':compound,
            'mapping_summary':{k:reconstruction['warhead_mapping'][k] for k in ('equivalent_mapping_count','normalization_for_matching_only')},
            'boundary_bonds':reconstruction['boundary_bonds'],'meaning':reconstruction['meaning']}}
        facts.append({'id':cid+':reconstruction','candidate_id':cid,'domain':'molecule','kind':'reference_reconstruction',
            'value':{'molecule_id':compound['molecule_id'],'paper':compound['paper'],'structure':compound['structure'],
                'attachment_atom':compound['attachment_atom'],'prediction_status':'not_run','human_review':'pending'},'source_ids':[sid]})
    sources['B:START'] = {'kind':'B_registered_H2', 'content':{'geometry':material['h2']['descriptive_geometry'],
        'attachment_literature':material['h2']['attachment_literature'],'unresolved':material['h2']['unresolved']}}
    facts.append({'id':'START:geometry','candidate_id':'START','domain':'molecule','kind':'computed_descriptive_geometry',
        'value':material['h2']['descriptive_geometry'],'source_ids':['B:START']})
    sources['M2:status'] = {'kind':'M2_current_scope', 'content':{'pending_inputs':value['pending_inputs'],
        'authority':value['authority'],'planned_candidates':dossier['planned_candidates']}}
    facts.append({'id':'C03:not_designed','candidate_id':'C03','domain':'molecule','kind':'not_run',
        'value':{'structure':None,'status':'not_designed','gpu':'not_run','human_approved':False},'source_ids':['M2:status']})
    result = {'version':VERSION,'b_input_ref':value['input_ref'],'dossier_ref':dossier['dossier_ref'],
        'facts':facts,'sources':sources,
        'coverage':'First-case curated candidate links, supplied A extraction and original passages, B H1/H2 and CPU outputs. No whole-paper/image/supplement extraction or new chemistry/GPU execution.',
        'not_completed':['약학 검수와 실제 사람 승인','정식 Phase1/같은 후보 집합 CPU preview','계산 범위·비용·MSA/GPU 준비','실제 GPU 계산·구조 비교','C03 설계','최종 보고·공개 배포']}
    result['input_digest'] = json_digest(result)
    return result
