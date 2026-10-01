/* Saved observations and opinions only; no ranking or execution actions. */
function savedNumber(value){
 if(value===null||value===undefined)return '미확인';
 if(Array.isArray(value))return value.map(savedNumber).join('–');
 return typeof value==='number'?(Number.isInteger(value)?String(value):value.toPrecision(4)):String(value);
}
function savedDownload(label,ref){return link(label,'/api/artifacts/'+encodeURIComponent(ref.artifact_id)+'/download');}
function savedDetails(title,value){const d=node('details');d.append(node('summary',title),node('pre',JSON.stringify(value,null,2),'saved-json'));return d;}
function savedTable(headers,rows){
 const wrap=node('div',undefined,'saved-table-wrap'),table=node('table',undefined,'saved-table'),head=node('thead'),hr=node('tr');
 for(const h of headers){const th=node('th',h);th.scope='col';hr.append(th);}head.append(hr);table.append(head);
 const body=node('tbody');for(const values of rows){const tr=node('tr');for(const v of values){const td=node('td');td.append(v instanceof Node?v:document.createTextNode(savedNumber(v)));tr.append(td);}body.append(tr);}table.append(body);wrap.append(table);return wrap;
}
async function renderSavedResults(dossier,container){
 const section=node('section',undefined,'saved-results');container.append(section);
 section.append(node('h3','B 저장 계산 결과와 전문가 답변'));
 try{
  const history=await api('/api/dossiers/'+encodeURIComponent(dossier.id)+'/saved-results');
  if(!container.contains(section))return;
  if(!history.length){section.append(node('p','이 인계 자료에 연결된 저장 계산 결과가 없습니다.'));return;}
  const chooser=node('select');chooser.setAttribute('aria-label','B 저장 결과 버전');
  for(const h of history){const o=node('option',h.created_at+' · '+h.summary.sample_count+'개 · '+(h.current_input?'현재 입력':'과거 입력'));o.value=h.id;chooser.append(o);}section.append(chooser);
  const content=node('div');section.append(content);let revision=0;
  async function show(){const ticket=++revision;content.replaceChildren(node('p','저장 근거를 읽고 있습니다…'));
   try{const data=await api('/api/saved-results/'+encodeURIComponent(chooser.value));if(ticket!==revision)return;
    content.replaceChildren();if(!data.current_input)content.append(notice('입력이 변경된 과거 결과입니다. 현재 입력의 검토·계산 근거로 재사용하려면 다시 대조해야 합니다.'));
    content.append(node('p','후보 '+data.summary.candidate_count+'개 · 저장 샘플 '+data.summary.sample_count+'개 · '+data.file_count+'개 원파일'),
      notice('알려진 기준 후보의 구조 재현성 비교입니다. 접촉·수소결합·국소 충돌은 잠정 결과이며, 효능 예측이나 자동 합격·탈락에 사용하지 않습니다.'));
    content.append(savedDownload('저장 계산·문헌 요약 JSON',data.projection_ref),savedDetails('자료 버전과 권한',{bundle_digest:data.bundle_digest,authority:data.authority,source_authority:data.source_authority}));
    content.append(link('통합 보고서와 버전 확인','/evidence-reports/'+encodeURIComponent(data.id)));
    for(const candidate of data.candidates){
     const card=node('article',undefined,'document');card.append(node('h4',candidate.compound_id+' · '+candidate.facts.compound.paper.name),node('p','구조 재구성과 계산 물성 확보 · 별도 합성 평가 미실시'));
     const rows=candidate.samples.map((s,i)=>{const exec=s.reported_execution,settings=exec?.settings||{};const btn=node('button','샘플 '+(i+1)+' 상세');btn.type='button';btn.onclick=()=>{detail.replaceChildren();
      detail.append(node('h4',s.run_id+' / rank '+s.rank),node('p','B 측 과거 실행 기록 · 이 PC에서 새 계산한 결과가 아닙니다.'));
      detail.append(savedTable(['구간','RMSD 범위 (Å)'],[['Warhead',s.parts_rmsd_ranges_A.warhead],['Linker',s.parts_rmsd_ranges_A.linker],['E3 ligand',s.parts_rmsd_ranges_A.recruiter]]));
      detail.append(node('p','모든 대칭 원자 대응의 최솟값–최댓값입니다. 구간별 최솟값이 같은 원자 대응에서 나왔다고 가정하지 않습니다. 정렬은 표적 기준으로 동일하게 적용했습니다.','hint'));
      for(const [label,ref] of Object.entries(s.reports))if(ref)detail.append(savedDownload({comparison:'구조 비교 원본',prior_quality:'선행 품질 보고서',preparation:'준비 보고서',contacts:'접촉 보고서'}[label],ref));
      detail.append(savedDetails('원자별 역할 차이',s.role_disagreements),savedDetails('검토 한계와 사유',s.review_issues),savedDetails('설정·모델·MSA·원자 대응 버전',{condition_digest:s.condition_digest,binding:s.binding,execution:s.reported_execution,alignment:s.alignment_method}));};
      return [btn,settings.use_potentials===true?'ON':settings.use_potentials===false?'OFF':null,settings.seed,s.summary?.target_rmsd_A,s.summary?.vhl_rmsd_after_target_alignment_A,s.summary?.confidence_score,s.bad_overlap_pairs.xray,s.bad_overlap_pairs.nuclear];});
     card.append(savedTable(['샘플','Potentials','Seed','표적 Å','VHL Å','Confidence','접촉 X선','접촉 핵'],rows));
     card.append(node('p','접촉은 선택 범위의 bad-overlap pair 수입니다. 두 수소 길이 규약은 합산하지 않습니다. 전체 clashscore가 아닙니다.','hint'));
     const detail=node('div',undefined,'saved-sample-detail');card.append(detail);card.append(savedDetails('문헌·계산 물성 원기록',candidate.facts));content.append(card);
    }
    const lit=node('details');lit.append(node('summary','실제 A 문헌 '+data.literature.records.length+'개와 B 대조'));
    for(const row of data.literature.records){const r=row.record;lit.append(node('h4',r.query.paper_name+' · '+r.observation.endpoint+' = '+savedNumber(r.observation.value)+' '+(r.observation.unit||'')),node('p','원 A 기록 보존 · 실험 조건 대조 전 직접 성능 비교 보류'),savedDetails('원 관측·조건·출처와 대조 차이',row));}
    for(const missing of data.literature.missing)lit.append(notice(missing.compound_id+' DC50: 현재 검토한 구간에서 미확인. 논문 전체에 값이 없다는 뜻이 아닙니다.'));content.append(lit);
    content.append(node('h4','전문가 의견과 후속 확인'));
    if(!data.expert_reviews.length)content.append(node('p','아직 연결된 답변이 없습니다.'));
    for(const review of data.expert_reviews){content.append(savedDownload('받은 답변서 원본',review.source_ref),node('p','의견 수신 · 정식 승인 별도 · 검토자 이름/검토일 '+(review.assessment.reviewer_name||'미기재')+' / '+(review.assessment.review_date||'미기재'),'hint'));
     for(const q of review.assessment.questions){const d=node('details');d.append(node('summary',q.question_id+' · 의견 수신'),node('h4','개발 반영'));q.applied.forEach(t=>d.append(node('p',t)));d.append(node('h4','후속 확인'));q.remaining.forEach(t=>d.append(node('p',t)));d.append(node('h4','답변 원문'),node('pre',q.answer_original,'saved-json'));content.append(d);}
     for(const s of review.assessment.literature_supplements){content.append(node('p',s.compound_id+' 문헌 보완: biological replicates '+s.value+' · Figure 1 caption 출처 · A 원기록 보존'),savedDetails('보완 근거·적용 범위',s));}
    }
    const files=node('details');files.append(node('summary','전체 원파일 찾아서 내려받기'));let loaded=false;
    files.ontoggle=async()=>{if(!files.open||loaded)return;try{const refs=await api('/api/saved-results/'+encodeURIComponent(data.id)+'/files');loaded=true;const search=node('input');search.placeholder='파일명 또는 run ID';search.setAttribute('aria-label','원파일 검색');files.append(search);const list=node('div');files.append(list);function draw(){list.replaceChildren();const entries=Object.entries(refs).filter(([name])=>name.toLowerCase().includes(search.value.toLowerCase()));list.append(node('p',entries.length+'개 일치 · 처음 40개 표시'));for(const [name,ref] of entries.slice(0,40)){const row=node('p');row.append(savedDownload(name,ref));list.append(row);}}search.oninput=draw;draw();}catch(e){files.append(notice(e.message));}};content.append(files);
   }catch(e){if(ticket===revision)content.replaceChildren(notice('저장 결과를 확인하지 못했습니다: '+e.message));}
  }
  chooser.onchange=show;await show();
 }catch(e){section.append(notice(e.message));}
}
