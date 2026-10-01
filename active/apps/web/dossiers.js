let dossierRenderVersion=0;
const stageNames={literature:'문헌·물질 연결',structure:'구조 확인',attachment:'부착점',assembly:'후보 조립',gpu:'GPU 계산',pharmacy_review:'약학 검수',g1:'G1 승인',report:'보고서'};
const stageValues={linked_pending_human_review:'원문 위치 연결 · 사람 확인 대기',B_cpu_result_imported:'B의 CPU 결과 수신',B_mapping_and_paper_rationale_linked:'원자 대응과 문헌 근거 연결',known_reference_reconstructed_by_B:'알려진 후보 재구성 완료',not_run:'미실행',pending:'대기',not_approved:'미승인',internal_dossier_only:'내부 검토 자료까지'};
async function renderDossiers(run){
 const version=++dossierRenderVersion,panel=$('candidate-dossiers');panel.hidden=false;panel.replaceChildren(node('h2','후보 물질별 진행 상황과 저장 결과'));
 try{
  const items=await api('/api/runs/'+run.id+'/dossiers');if(version!==dossierRenderVersion)return;
  if(!items.length){panel.append(node('p','아직 이 탐색에 연결한 B 후보 자료가 없습니다.'));return;}
  const chooser=node('select');chooser.setAttribute('aria-label','후보 인계 이력');for(const d of items){const o=node('option',new Date(d.created_at).toLocaleString()+' · '+(d.current_input?'현재 입력':'이전 입력'));o.value=d.id;chooser.append(o);}panel.append(chooser);
  const body=node('div');panel.append(body);
  function show(){const d=items.find(x=>x.id===chooser.value);body.replaceChildren();
   if(!d.current_input)body.append(notice('입력 버전이 달라진 과거 인계 자료입니다. 현재 계산·승인에 사용할 수 없습니다.'));
   body.append(node('p','아래 인계 당시의 단계는 과거 기록입니다. 이후 계산 자료와 전문가 답변은 B 저장 결과에서 확인하세요.','hint'));
   body.append(node('p','C01·C02는 논문에 발표된 물질을 재구성한 기준 후보입니다. 새로운 후보 C03는 아직 설계하지 않았습니다.'));
   for(const c of d.candidates){const card=node('article',undefined,'document');card.append(node('h3',c.id+' · '+c.paper_name+' (논문 화합물 '+c.paper_compound_number+'번)'),node('p','PDB '+c.pdb+' · CCD '+c.ccd+' · '+c.identity.inchikey,'hint'));
    const list=node('dl',undefined,'candidate-stages');for(const k of ['literature','structure','attachment','assembly','pharmacy_review','g1','gpu','report']){const v=c.stages[k];list.append(node('dt',stageNames[k]),node('dd',stageValues[v]||v));}card.append(list);
    card.append(node('p','출발 부착 원자 '+c.attachment.starting_atom+' → '+c.attachment.candidate_atom+' · 대칭에 따른 대응 '+c.attachment.equivalent_mapping_count+'가지'));
    const details=node('details');details.append(node('summary','문헌 근거와 실제 구조 파일'));
    for(const key of c.anchors){const a=d.anchors[key];details.append(node('blockquote',a.quote),node('small',a.xpath),link('해당 원문 구간','/api/artifacts/'+encodeURIComponent(a.text_ref.artifact_id)+'/download'));}
    details.append(node('p','AI 분석에서 연결한 근거 기록 '+c.literature_evidence.length+'개. 값의 조건과 Critic 지적은 아래 저장된 분석에서 확인하세요.'));
    for(const name of [c.id+'.sdf',c.id+'.smi',c.id+'-atom-map.json',c.id+'-boltz-input.yaml']){details.append(link(name,'/api/artifacts/'+encodeURIComponent(c.files[name].artifact_id)+'/download'));}
    card.append(details,node('p','다음: '+c.missing.join(' → '),'hint'));body.append(card);
   }
   body.append(link('저장 결과 전용 화면 열기','/saved-review/'+encodeURIComponent(d.id)));
   renderSavedResults(d,body);
   d.limitations.forEach(x=>body.append(node('p',x,'hint')));body.append(link('현재 인계 자료 JSON 내려받기','/api/artifacts/'+encodeURIComponent(d.dossier_ref.artifact_id)+'/download'));
  }
  chooser.onchange=show;show();
 }catch(e){if(version===dossierRenderVersion)panel.append(notice('후보 인계 자료를 확인하지 못했습니다: '+e.message));}
}
