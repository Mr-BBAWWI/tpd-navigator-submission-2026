'use strict';
(function(){
  const api={};
  const n=(tag,cls,text)=>{const el=document.createElement(tag);if(cls)el.className=cls;if(text!==undefined)el.textContent=String(text);return el};
  const add=(parent,...children)=>{children.forEach(child=>child&&parent.appendChild(child));return parent};
  const val=(value,fallback='미보고')=>value===undefined||value===null||value===''?fallback:Array.isArray(value)?value.join(' · '):typeof value==='object'?JSON.stringify(value):String(value);
  const list=value=>Array.isArray(value)?value:value&&typeof value==='object'?(Array.isArray(value.records)?value.records:Array.isArray(value.cards)?value.cards:[]):[];
  const source=(evidence,name)=>evidence?.sources?.[name]||null;
  const safeUrl=value=>{if(typeof value!=='string'||!/^https?:\/\//i.test(value.trim()))return null;try{const url=new URL(value.trim());return url.protocol==='http:'||url.protocol==='https:'?url.href:null}catch{return null}};
  const title=(kicker,text,copy)=>{const head=n('div','de-heading');const left=n('div');add(left,n('p','kicker',kicker),n('h3','',text));if(copy)add(left,n('p','de-copy',copy));add(head,left);return head};
  const card=(heading,rows,cls='')=>{const el=n('article',`de-card ${cls}`.trim());add(el,n('h4','',heading));(rows||[]).forEach(([label,value])=>{const row=n('div','de-kv');add(row,n('span','',label),n('b','',val(value)));el.appendChild(row)});return el};
  const details=(summary,data,open=false)=>{const el=n('details','de-json');el.open=open;add(el,n('summary','',summary),n('pre','',JSON.stringify(data??null,null,2)));return el};
  const section=(root,kickerText,heading,copy,count)=>{const el=n('details','card de-section');const summary=n('summary','de-section-summary');add(summary,n('span','kicker',kickerText),n('strong','',heading));if(count!==undefined)add(summary,n('span','badge neutral',`${count}건`));el.appendChild(summary);if(copy)el.appendChild(n('p','de-copy',copy));root.appendChild(el);return el};
  const grid=parent=>add(parent,n('div','de-grid')).lastChild;
  const link=(label,url)=>{const href=safeUrl(url);if(!href)return null;const a=n('a','de-source-link',label);a.href=href;a.target='_blank';a.rel='noopener noreferrer';return a};
  const findValue=(obj,keys)=>{if(!obj||typeof obj!=='object')return undefined;for(const key of keys)if(obj[key]!==undefined)return obj[key];for(const value of Object.values(obj)){if(value&&typeof value==='object'){const found=findValue(value,keys);if(found!==undefined)return found}}};
  function renderScope(root,evidence,result){
    const war=source(evidence,'warhead_catalog.json');
    const data=war?.data;
    const cards=Array.isArray(data)?data:list(data);
    const scope=result?.parent_scope||{};
    const el=section(root,'Evidence scope','수집 카탈로그와 선택한 실제 설계 parent','카탈로그의 parent 기록은 각각 독립된 근거입니다. 이번 실행에는 선택한 SMARCA2 parent 하나의 구조와 그 parent에 연결된 SAR만 적용됩니다.');
    const g=grid(el);
    add(g,card('범위 구분',[
      ['수집 warhead 카드',cards.length||scope.collected_parent_catalog_count],
      ['실제 설계 parent',scope.actual_design_parent_count??(result?1:'실행 전')],
      ['설계 parent ID',scope.actual_design_parent_id||result?.warhead?.id],
      ['선택 입력 역할',scope.selected_parent_input_role],
      ['receptor frame',scope.receptor_frame],
      ['기술 ready = 과학 승인',scope.technical_ready_is_scientific_approval],
      ['source 상태',war?.status]
    ],'de-emphasis'));
    cards.forEach((item,index)=>{
      const versions=list(item.experimental_versions);
      const c=card(`근거 카드 ${index+1}`,[
        ['대표 record',item.representative_record_id],
        ['domain',item.domain_partition],
        ['실험 version',versions.length],
        ['co-crystal tier',item.evidence_tiers?.co_crystal],
        ['binding tier',item.evidence_tiers?.binding],
        ['PROTAC tier',item.evidence_tiers?.PROTAC]
      ]);
      const missing=item.missing_evidence||{};
      add(c,n('p','de-warning',`누락/미지: ${Object.entries(missing).map(([k,v])=>`${k} ${val(v)}`).join(' · ')||'별도 보고 없음'}`));
      if(versions.length)add(c,details(`experimental versions ${versions.length}건`,versions));
      g.appendChild(c);
    });
    if(!cards.length)add(el,details('warhead catalog 원문 또는 미구성 상태',data??war));
  }
  function renderMedchem(root,evidence,result){
    const src=source(evidence,'medchem_evidence.json');
    const data=src?.data;
    const records=Array.isArray(data)?data:list(data?.parent_catalog||data);
    const el=section(root,'Measured evidence','SMARCA2/4 parent 활성 근거','측정 IC50의 relation·값·원문 단위를 그대로 표시하며 Kd로 바꾸지 않습니다.',records.length);
    const g=grid(el);
    const measurement=(record,target)=>(record.measurements||[]).find(m=>m.target===target&&m.endpoint==='IC50');
    const formatted=m=>m?`${m.relation||''}${m.value} ${m.unit||''}`.trim():undefined;
    records.forEach((record,index)=>{
      const sm2=measurement(record,'SMARCA2'),sm4=measurement(record,'SMARCA4');
      const c=card(record.compound_id||record.parent_id||record.record_id||`parent record ${index+1}`,[
        ['SMARCA2 IC50',formatted(sm2)],['SMARCA4 IC50',formatted(sm4)],
        ['SMARCA2 조건',sm2?.conditions],['SMARCA4 조건',sm4?.conditions],
        ['근거 종류','IC50 · Kd 아님'],['source/SI',sm2?.assay_source_locator||sm4?.assay_source_locator]
      ]);
      const doi=record.source?.doi;
      const a=link('DOI 원문 열기',doi?`https://doi.org/${doi}`:null)||link('원문 출처 열기',record.source_url);if(a)c.appendChild(a);
      g.appendChild(c);
    });
    if(!records.length)add(el,card('측정 근거 상태',[['source',src?.status],['records',0],['해석','측정값이 구성되지 않았거나 알려지지 않은 schema입니다.']]));
    const joined=result?.selected_parent_measured_evidence;
    if(joined){
      const selected=section(root,'Selected parent measured evidence','선택 parent의 실제 측정 연결','정확한 heavy-atom graph와 formal charge가 모두 같은 경우만 연결합니다. atom map과 H만 제거하며 임의 중성화는 하지 않습니다.',joined.measurements?.length||0);
      const sg=grid(selected),record=joined.matched_record||{};
      add(sg,card('선택 parent join',[
        ['상태',joined.status],['compound',record.compound_id],['identity',record.identity_status],
        ['formal charge',record.formal_charge],['DOI',record.source?.doi],['source row',record.source?.record_number_1_based_including_header],
        ['신규 analog 측정 증명',joined.new_analog_measurement_status],['설계 policy 변경',joined.changes_design_policy]
      ],joined.status==='exact_measured_parent_join'?'':'de-pending'));
      (joined.measurements||[]).forEach(m=>sg.appendChild(card(`${m.target||'target'} ${m.endpoint||'measurement'}`,[
        ['reported',`${m.relation||''}${m.value} ${m.unit||''}`.trim()],['raw',m.raw_string],['assay',m.assay],
        ['conditions',m.conditions],['source locator',m.assay_source_locator],['source column',m.source_column_label]
      ])));
      (joined.source_sar||[]).forEach(pair=>sg.appendChild(card(`정보용 source pair ${pair.pair_id}`,[
        ['변화',pair.requested_change_label],['source site consensus',pair.source_site_consensus],
        ['selected parent maps',pair.selected_parent_affected_atom_maps],['full mapping unambiguous',pair.full_graph_mapping_unambiguous],
        ['해석',pair.interpretation]
      ])));
    }
  }
  function renderPipeline(root,result){
    const el=section(root,'Ordered pipeline','실행 funnel · 실제 count와 탈락 이유','각 단계의 configured 목표와 실제 결과를 분리합니다. 도킹을 끈 preview는 qualified가 0이며 조립 근거가 아닙니다.');
    if(!result){add(el,n('p','de-empty','결과를 불러오면 단계별 실제 count와 shortfall이 표시됩니다.'));return}
    const s=result.summary||{},counts=s.stage_counts||{};
    const stages=[
      ['generate','생성'],['cheap_filter','cheap filter'],['stereoisomer_expansion','입체 확장'],['eligibility','도킹 적격'],
      ['dock_all_eligible','모든 적격 analog 도킹'],['parent_redock_and_pose_filter','parent/pose filter'],['diverse_final_selection','다양성 최종 선택'],['assemble_qualified','CRBN·VHL 조립']
    ];
    const flow=n('ol','de-flow');
    stages.forEach(([key,label])=>{const li=n('li');add(li,n('b','',label),n('span','',counts[key]?Object.entries(counts[key]).map(([k,v])=>`${k} ${val(v)}`).join(' · '):'단계 count 미보고'));flow.appendChild(li)});
    el.appendChild(flow);
    const preview=s.selection_mode==='explicit_preview_not_qualified'||result.parameters?.dock===false;
    add(el,n('div',`de-callout ${s.shortfall?'warning':''}`,preview?`그래프 preview: qualified 0 · 선택 ${val(s.selected_analogs,0)} · shortfall ${val(s.shortfall,0)}`:`qualified ${val(s.qualified_analogs,0)} / 요청 ${val(s.panel_size_requested)} · shortfall ${val(s.shortfall,0)}`));
    const reasons={};(result.rejections||[]).forEach(row=>{const key=row.reason||(row.reasons||[]).join(', ')||row.stage||'미분류';reasons[key]=(reasons[key]||0)+1});
    const rejection=card('탈락 이유',Object.entries(reasons).map(([key,count])=>[key,count]));if(!Object.keys(reasons).length)add(rejection,n('p','de-empty','탈락 사유 없음 또는 미보고'));el.appendChild(rejection);
  }
  function renderScientific(root,evidence,result){
    const scientific=result?.chemical_states||{};
    const prep=source(evidence,'preparation.json');
    const prepData=prep?.data||{};
    const protein=prepData.protein_preparation||{};
    const pKas=protein.pKa_entries_from_actual_output||[];
    const counts=prepData.actual_interaction_counts||prepData.comparison?.actual_interaction_counts||{};
    const mapping=protein.heavy_atom_mapping||{};
    const shiftValues=[prepData.before_ligand_hydrogen_receipt?.final_heavy_max_displacement_A,prepData.after_ligand_hydrogen_receipt?.final_heavy_max_displacement_A,mapping.maximum_displacement_A];
    const maxShift=shiftValues.every(Number.isFinite)?Math.max(...shiftValues):null;
    const el=section(root,'Chemical state & preparation','미세상태·수소·접촉·단백질 준비','Ligand population/pKa는 예측하지 않았지만 protein pKa는 실제 PROPKA 출력값을 표시합니다.',pKas.length);
    const g=grid(el);
    const micro=scientific.microstates||{};
    const m=card('Parent microstates',[
      ['상태',micro.status],['열거 수',micro.count],['pH 맥락',findValue(micro,['pH_context','ph_context'])||'별도 pH population 미계산'],
      ['population / pKa','미보고 · 추정하지 않음'],['선택',micro.selection]
    ]);if(micro.records)add(m,details('열거된 microstate',micro.records));g.appendChild(m);
    const h=scientific.ligand_hydrogen_preparation||{};
    const hc=card('고정 heavy-atom H 준비',[['상태',h.status],['도킹 입력 사용',h.docking_input_used],['범위','parent reference · heavy atoms 고정']]);if(h.receipt)add(hc,details('수소 준비 receipt',h.receipt));g.appendChild(hc);
    const interactions=scientific.parent_interactions||{};
    const ic=card('방향성 접촉 검토',[['protein H 상태',interactions.protein_hydrogen_status],['before / after',`${val(counts.before)} / ${val(counts.after)}`],['gained / lost / retained',`${val(counts.gained)} / ${val(counts.lost)} / ${val(counts.retained)}`],['최대 heavy shift',maxShift===null?'미지':`${maxShift} Å`],['판정','사람 검토 필요']]);g.appendChild(ic);
    const pc=card('Protein preparation 근거',[['source',prep?.status],['상태',protein.status],['실제 protein pKa',`${pKas.length} entries · PROPKA`],['중원자 매핑',`${val(mapping.input_count)} → ${val(mapping.output_count)} · mapped ${val(mapping.mapped_count)}`],['좌표 최대 이동',`${val(mapping.maximum_displacement_A)} Å`],['추가 원자',`${val(mapping.extra_in_output)} · OXT retained, review 필요`]]);
    if(pKas.length)add(pc,details(`실제 PROPKA pKa ${pKas.length}건`,pKas.map(x=>({residue:x.residue,sequence_number:x.sequence_number,chain:x.chain,pKa:x.pKa}))));g.appendChild(pc);
  }
  function renderCalibration(root,evidence,result){
    const catalogCalibration=source(evidence,'crbn_calibration.json')?.data||source(evidence,'crbn_calibration')?.data||source(evidence,'crbn_calibration.json')||source(evidence,'crbn_calibration');
    const calibration=result?.calibration||catalogCalibration;
    const crbn=calibration?.CRBN||calibration||{};
    const seeds=Array.isArray(crbn.seed_receipts)?crbn.seed_receipts:Array.isArray(crbn.seeds)?crbn.seeds:[];
    const el=section(root,'CRBN calibration','Seed 재현성과 confidence 분리','세 seed의 실제 값을 모두 표시하며 임의 pass 판정을 만들지 않습니다. CRBN과 VHL raw score를 교차 순위화하지 않습니다.',seeds.length);
    if(!calibration){add(el,n('p','de-empty','catalog에 CRBN calibration source가 구성되지 않았습니다.'));return}
    const g=grid(el);
    add(g,card('Calibration 상태',[
      ['상태',crbn.status],['reproduction status',crbn.reported_reproduction_status||crbn.reproduction_status],
      ['confidence',crbn.confidence||findValue(crbn.source_record,['confidence'])||'별도 미보고'],['cross-E3 ranking',calibration.cross_e3_raw_score_ranking===false?'금지':'미보고'],
      ['new MSA baseline',calibration.new_MSA_baseline]
    ],'de-emphasis'));
    seeds.forEach((seed,index)=>{const metrics=seed.inspection?.comparison?.metrics||seed.comparison?.metrics||{};const confidence=seed.inspection?.model_confidence||seed.model_confidence||{};const c=card(`Seed ${val(seed.seed,index+1)}`,[
      ['execution status',seed.execution_status||seed.status||seed.inspection?.status||seed.reproduction_status],
      ['target Cα RMSD',metrics.target_CA_RMSD_A],['E3 Cα RMSD',metrics.e3_CA_RMSD_after_target_alignment_A],
      ['ligand RMSD',metrics.ligand_heavy_atom_RMSD_after_target_alignment_A],['contact Jaccard',metrics.contact_jaccard],
      ['model confidence',confidence.confidence_score]
    ]);g.appendChild(c)});
    if(!seeds.length)add(el,details('알려지지 않은 calibration schema',crbn.source_record??crbn));
    add(el,details('metrics · baseline · limitations',{metrics:crbn.metrics,new_MSA_baseline:calibration.new_MSA_baseline,limitations:crbn.limitations||findValue(crbn.source_record,['limitations'])}));
  }
  function renderRoutes(root,evidence,result){
    const route=source(evidence,'route_evidence.json');
    const routeData=route?.data;
    const records=Array.isArray(routeData)?routeData:Array.isArray(routeData?.evidence_records)?routeData.evidence_records:list(routeData);
    const el=section(root,'Linker & route review','Linker descriptor와 C01/C02 SI','C01/C02의 정확한 SI 수율·정제·특성화는 선례로만 표시하며 신규 후보의 합성 증명이 아닙니다. 요구 ternary geometry가 없으면 미지입니다.');
    const g=grid(el);
    records.filter(r=>['C01','C02'].includes(String(r.compound_id||r.id||'').toUpperCase())).forEach(record=>{
      const characterization=Array.isArray(record.characterization)?record.characterization:[];
      const locations=characterization.map(x=>`${x.type||'record'} ¶${x.paragraph_index_1_based||'?'}`).join(' · ');
      const c=card(record.compound_id||record.id,[['reported yield',record.yield?.reported_yield],['isolated mass',record.yield?.isolated_mass],['purification',record.purification?.summary],['characterization',`${characterization.length}건${locations?` · ${locations}`:''}`],['신규 합성 증명','아님']]);
      const a=link('SI 출처 열기',record.source_url);if(a)c.appendChild(a);g.appendChild(c);
    });
    const candidates=result?.protac_candidates||[];
    const unique=new Map();candidates.forEach(candidate=>{if(candidate.linker_id&&!unique.has(candidate.linker_id))unique.set(candidate.linker_id,candidate)});
    unique.forEach((candidate,id)=>{const a=candidate.linker_assessment||{};const c=card(`Linker ${id}`,[
      ['descriptor MW',a.descriptors?.molecular_weight_g_mol],['rotatable bonds',a.descriptors?.rotatable_bonds_rdkit_strict],
      ['sampled range Å',a.sampled_3d?.end_to_end_range_A?`${val(a.sampled_3d.end_to_end_range_A.min_observed)}–${val(a.sampled_3d.end_to_end_range_A.max_observed)}`:'미계산'],
      ['required geometry',a.ternary_geometry?.status||candidate.linker_geometry_status||'미지'],['합성 상태',candidate.synthetic_route_status||'미평가']
    ]);add(c,details('linker 평가 전체',a));g.appendChild(c)});
    if(!records.length&&!unique.size)add(el,details('route evidence 상태',route));
  }
  function renderLedger(root,evidence,result){
    const el=section(root,'Acceptance ledger','구성 · 실제 metric · 승인 대기','계산 완료는 과학적 승인이 아닙니다. Reviewer A의 명시적 검토 전에는 승인 상태로 바뀌지 않습니다.');
    const ledger=result?.acceptance_ledger||{};
    const s=result?.summary||{},g=grid(el);
    add(g,card('Configured',[
      ['evidence manifest',evidence?.manifest?.version],['evidence 상태',evidence?.status],['pipeline version',ledger.functional_implementation?.version],['요청 panel',s.panel_size_requested]
    ]),card('Actual metrics',[
      ['qualified',s.qualified_analogs],['selected',s.selected_analogs],['assembled',s.protac_count],['shortfall',s.shortfall]
    ]),card('Scientific approval',[
      ['상태',ledger.scientific_approval?.status||'not_approved'],['자동 승인','없음'],['검토자',ledger.scientific_approval?.reviewer||'Reviewer A 승인 대기'],['approved',ledger.scientific_approval?.approved??false]
    ],'de-pending'));
    if(result?.review_questions?.length){const q=n('div','de-questions');add(q,n('h4','','검토 질문'));const ul=n('ul');result.review_questions.forEach(item=>add(ul,n('li','',item)));add(q,ul);el.appendChild(q)}
  }
  api.render=function(catalog,result){
    const root=document.getElementById('designEvidenceContent');
    const panel=document.getElementById('designEvidencePanel');
    if(!root||!panel)return;
    while(root.firstChild)root.removeChild(root.firstChild);
    const evidence=result?.source_catalog||catalog?.acceptance_evidence||null;
    if(!evidence){panel.classList.remove('hidden');add(root,n('p','de-empty','Acceptance evidence가 아직 구성되거나 로드되지 않았습니다.'));return}
    panel.classList.remove('hidden');
    renderScope(root,evidence,result);
    renderMedchem(root,evidence,result);
    renderPipeline(root,result);
    renderScientific(root,evidence,result);
    renderCalibration(root,evidence,result);
    renderRoutes(root,evidence,result);
    renderLedger(root,evidence,result);
  };
  window.TPDDesignEvidence=api;
})();
