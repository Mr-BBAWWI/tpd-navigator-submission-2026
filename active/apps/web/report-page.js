const node=(tag,text,className)=>{const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(className)n.className=className;return n;};
const link=(text,url)=>{const a=node('a',text);a.href=url;a.target='_blank';a.rel='noopener noreferrer';return a;};
const notice=text=>node('div',text,'notice');
async function api(path,options){const response=await fetch(path,options);const value=await response.json();if(!response.ok)throw Error(value.detail||'자료 요청 실패');return value;}
const reportReasons={input_changed:'탐색 입력 변경',dossier_superseded:'새 A 문헌·후보 자료 등록',result_superseded:'새 B 결과 등록',expert_opinions_changed:'전문가 의견 추가·변경',report_implementation_changed:'보고서 작성 규칙 변경',storage_mode_changed:'저장소 모드 변경',dossier_changed:'A 자료 연결 변경',result_changed:'B 자료 연결 변경',replay_read_only:'재생 전용 사본'};
function reportStatus(freshness){return freshness.current?'자료 버전 일치 · 작성 당시 검토 초안':'재검토 필요 · '+freshness.reasons.map(r=>reportReasons[r]||r).join(', ');}
function renderReport(value,container){
 const body=value.content;container.replaceChildren(node('h2','저장된 보고서'),notice(reportStatus(value.freshness)),node('p','작성 '+value.created_at+' · 현재 버전 확인 '+value.freshness.checked_at));
 container.append(node('p',`후보 ${body.counts.candidates}개 · 샘플 ${body.counts.samples}개 · 문헌 관측 ${body.counts.observations}개 · 전문가 질문 ${body.counts.expert_questions}개`));
 container.append(savedDownload('보고서 내려받기 (.md)',value.markdown_ref),savedDownload('전체 근거 내려받기 (.json)',value.content_ref),savedDownload('버전·출처 목록 (.json)',value.record_ref));
 container.append(link('범위별 사람 검토와 판정 이력','/human-review/'+encodeURIComponent(value.id)));
 container.append(node('p','다운로드 파일은 작성 당시 기록입니다. 새 자료가 추가되면 이 화면에서 재검토 여부를 확인하세요.','hint'));
 const limits=node('details');limits.append(node('summary','해석 범위와 미완료 항목'));body.limitations.forEach(t=>limits.append(node('p',t)));container.append(limits);
 for(const c of body.candidates){const a=c.identity;container.append(node('h3',a.id+' · '+a.paper_name+' · 논문 화합물 '+a.paper_compound_number),node('p',a.pdb+' / '+a.ccd+' · 알려진 기준 후보'));
  const rows=c.saved_evidence.samples.map(s=>{const settings=s.reported_execution.settings,summary=s.summary;return [s.run_id+' / '+s.rank,settings.seed,settings.use_potentials===undefined?'미확인':settings.use_potentials?'ON':'OFF',summary.target_rmsd_A,summary.vhl_rmsd_after_target_alignment_A,summary.confidence_score,s.status.contacts];});
  container.append(savedTable(['실행 / rank','seed','potentials','표적 RMSD Å','표적 정렬 후 VHL RMSD Å','confidence','접촉 처리 상태'],rows));
  container.append(savedDetails(a.id+' 부위 비교·원 보고서 참조·미해결 항목',c));
 }
 container.append(node('h3','A 문헌 관측'),notice('표적·조건이 다른 관측을 직접 비교하지 않습니다. 미확인은 0이 아니며, 논문 전체에 값이 없다는 뜻도 아닙니다.'));
 container.append(savedTable(['논문 후보','관측 대상','관측','원값','단위','관계'],body.literature.records.map(row=>{const r=row.record,o=r.observation;return [r.query.paper_name,o.target,o.endpoint,o.value,o.unit,o.relation];})));
 container.append(savedDetails('원 관측·실험 조건·위치와 B 대조',body.literature),node('h3','전문가 의견과 후속 검토'));
 if(!body.expert_reviews.length)container.append(node('p','연결된 답변이 없습니다.'));
 for(const review of body.expert_reviews){container.append(savedDownload('받은 답변 원본',review.source_ref),node('p','검토자 '+(review.assessment.reviewer_name||'미기재')+' · 검토일 '+(review.assessment.review_date||'미기재')+' · 정식 승인 별도'));
  for(const q of review.assessment.questions){const d=node('details');d.append(node('summary',q.question_id+' · 의견과 반영 범위'),node('h4','개발 반영'));q.applied.forEach(t=>d.append(node('p',t)));d.append(node('h4','남은 검토'));q.remaining.forEach(t=>d.append(node('p',t)));d.append(node('h4','답변 원문'),node('pre',q.answer_original,'saved-json'));container.append(d);}
  container.append(savedDetails('출처가 확인된 문헌 보완',review.assessment.literature_supplements));
 }
 container.append(savedDetails('아직 설계하지 않은 후보',body.planned_candidates));
 const sources=node('details');sources.append(node('summary','근거 파일 '+value.source_catalog.length+'개와 버전'));for(const ref of value.source_catalog){const row=node('p');row.append(savedDownload(ref.artifact_id,ref),node('small','SHA256 '+ref.sha256));sources.append(row);}container.append(sources,savedDetails('보고서 입력 버전과 권한',{input_digest:value.input_digest,binding:value.binding,authority:value.authority}));
}
document.addEventListener('DOMContentLoaded',async()=>{
 const controls=document.getElementById('report-controls'),content=document.getElementById('report-content'),resultId=location.pathname.split('/').pop(),root='/api/saved-results/'+encodeURIComponent(resultId);
 let revision=0;
 async function load(selected){try{
  const [health,prep,history]=await Promise.all([api('/api/health'),api(root+'/report-preparation'),api(root+'/reports')]);controls.replaceChildren();
  controls.append(link('계산 결과 상세로 돌아가기','/saved-review/'+encodeURIComponent(prep.dossier_id)));
  if(prep.blocked_reasons.length)controls.append(notice('현재 자료로 새 보고서를 작성할 수 없습니다: '+prep.blocked_reasons.map(r=>reportReasons[r]||r).join(', ')));
  if(health.read_only||health.review_only)controls.append(node('p','보고서는 내부 작성 환경에서 생성합니다. 이 화면에서는 저장된 자료를 확인합니다.','hint'));
  else if(prep.can_create){const button=node('button','현재 자료로 보고서 저장');button.onclick=async()=>{button.disabled=true;try{const receipt=await api(root+'/reports',{method:'POST',headers:{'Content-Type':'application/json','X-TPD-Local':'1'},body:JSON.stringify({expected_digest:prep.input_digest})});await load(receipt.report_id);}catch(e){controls.append(notice('저장하지 못했습니다: '+e.message+' · 자료를 새로 확인한 뒤 다시 작성하세요.'));button.disabled=false;}};controls.append(button);}
  const refresh=node('button','자료 버전 새로 확인');refresh.onclick=()=>load(selected);controls.append(refresh);
  if(!history.length){content.replaceChildren(node('p','아직 저장된 보고서가 없습니다.'));return;}
  const chooser=node('select');chooser.setAttribute('aria-label','통합 보고서 버전');history.forEach(h=>{const option=node('option',h.created_at+' · '+reportStatus(h.freshness));option.value=h.id;chooser.append(option);});if(history.some(h=>h.id===selected))chooser.value=selected;controls.append(chooser);
  async function show(){const ticket=++revision;content.replaceChildren(node('p','보고서를 읽고 있습니다…'));try{const report=await api('/api/evidence-reports/'+encodeURIComponent(chooser.value));if(ticket===revision){renderReport(report,content);historyReplace(report.id);}}catch(e){if(ticket===revision)content.replaceChildren(notice(e.message));}}
  chooser.onchange=show;await show();
 }catch(e){controls.replaceChildren(notice(e.message));content.replaceChildren();}}
 function historyReplace(id){const url=new URL(location.href);url.searchParams.set('report',id);window.history.replaceState(null,'',url);}
 await load(new URL(location.href).searchParams.get('report'));
});
