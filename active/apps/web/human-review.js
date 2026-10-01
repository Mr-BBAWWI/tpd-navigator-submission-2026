const node=(tag,text,className)=>{const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(className)n.className=className;return n;};
const link=(text,url)=>{const a=node('a',text);a.href=url;a.target='_blank';a.rel='noopener noreferrer';return a;};
const notice=text=>node('div',text,'notice');
async function api(path,options){const r=await fetch(path,options),v=await r.json();if(!r.ok)throw Error(v.detail||'요청 실패');return v;}
const post=(path,body)=>api(path,{method:'POST',headers:{'Content-Type':'application/json','X-TPD-Local':'1'},body:JSON.stringify(body)});
const humanStatuses={pending:'판정 대기',unassigned:'검수자 배정 대기',accepted_for_scope:'요청 범위 내 수용',changes_requested:'수정 요청',deferred:'판단 보류',withdrawn:'판정 철회',requires_re_review:'재검토 필요'};
const verdictNames={accept_scope:'범위 내 수용',request_changes:'수정 요청',defer:'판단 보류'};
const reasonNames={input_changed:'입력 변경',dossier_superseded:'새 문헌·후보 자료 등록',result_superseded:'새 계산 결과 등록',expert_opinions_changed:'전문가 의견 변경',report_implementation_changed:'보고서 작성 규칙 변경',review_policy_changed:'검토 기록 규칙 변경',request_superseded:'다음 검토 요청으로 대체',reviewer_inactive:'검수자 계정 비활성',storage_mode_changed:'저장소 모드 변경',input_changed_during_read:'조회 중 입력 변경'};
const reviewDate=value=>new Date(value).toLocaleString('ko-KR');
document.addEventListener('DOMContentLoaded',()=>{
 const reportId=location.pathname.split('/').pop(),context=document.getElementById('review-context'),login=document.getElementById('review-login'),content=document.getElementById('human-review-content');
 let revision=0;
 async function load(){const ticket=++revision;try{
  const [report,requests,session,health,predecessors]=await Promise.all([api('/api/evidence-reports/'+encodeURIComponent(reportId)),api('/api/evidence-reports/'+encodeURIComponent(reportId)+'/review-requests'),api('/api/review-auth/session'),api('/api/health'),api('/api/evidence-reports/'+encodeURIComponent(reportId)+'/review-predecessors')]);
  if(ticket!==revision)return;
  context.replaceChildren(node('p','보고서 작성 '+reviewDate(report.created_at)+' · 후보 '+report.counts.candidates+'개 · 샘플 '+report.counts.samples+'개'),link('대상 보고서와 근거 확인','/evidence-reports/'+encodeURIComponent(report.result_id)+'?report='+encodeURIComponent(report.id)),notice('원문 의견 수신과 판정 기록은 별개입니다. 선택한 범위의 판단만 기록하며 전체 효능·합성 타당성을 자동 승인하지 않습니다.'));
  if(!report.freshness.current)context.append(notice('이 보고서의 자료가 변경되었습니다. 최신 보고서를 작성하고 검토를 다시 요청해 주세요.'));
  renderLogin(session,health);
  content.replaceChildren(node('h2','검토 요청과 판정 이력'));
  const refresh=node('button','최신 상태 확인');refresh.onclick=load;content.append(refresh);
  if(!requests.length)content.append(node('p','등록된 검토 요청이 없습니다.'));
  for(const request of requests)renderRequest(request,session,health);
  for(const prior of [...predecessors,...requests.filter(r=>r.recheck_reasons.includes('review_policy_changed')&&!r.superseded_by)]){
   const item=node('section',undefined,'document document-top');item.append(node('h3','이전 검토 요청에서 이어가기'),node('p',(prior.assigned_reviewer?.name||'미배정')+' · '+humanStatuses[prior.status]),link('이전 보고서의 검토 이력','/human-review/'+encodeURIComponent(prior.report_id)));
   if(session.actor&&prior.assigned_reviewer?.id===session.actor.id&&!health.read_only&&report.freshness.current){const replace=node('button','현재 보고서로 재검토 요청 만들기');replace.onclick=async()=>{replace.disabled=true;try{await post('/api/evidence-reports/'+reportId+'/review-requests',{expected_report_sha256:report.record_ref.sha256,scopes:Object.keys(prior.scopes),supersedes:prior.id});await load();}catch(error){item.append(notice(error.message));replace.disabled=false;}};item.append(replace);}content.append(item);
  }
  if(session.actor&&!health.read_only&&report.freshness.current){const create=node('button','이 보고서의 전체 범위 검토 요청 만들기');create.onclick=async()=>{create.disabled=true;try{await post('/api/evidence-reports/'+encodeURIComponent(report.id)+'/review-requests',{expected_report_sha256:report.record_ref.sha256,scopes:['chemical_state','structural_comparison','literature_conditions','synthesis_evidence'],supersedes:null});await load();}catch(e){content.append(notice(e.message));create.disabled=false;}};content.append(create);}
 }catch(e){if(ticket===revision)content.replaceChildren(notice(e.message));}}
 function renderLogin(session,health){login.replaceChildren(node('h2','검수자 확인'));
  if(health.read_only){login.append(node('p','열람 전용 화면입니다. 판정 기록은 검토 전용 실행 환경에서 가능합니다.'));return;}
  if(session.actor){login.append(node('p',session.actor.name+' 로그인 중 · 팀이 등록한 로컬 계정'));const logout=node('button','로그아웃');logout.onclick=async()=>{await post('/api/review-auth/logout',{});await load();};login.append(logout);return;}
  if(!session.reviewers.length){login.append(node('p','등록된 검수자가 없습니다. 담당자가 계정을 먼저 등록해야 합니다.'));return;}
  const form=node('form'),select=node('select');select.id='reviewer-account';select.name='reviewer';const label=node('label','검수자');label.htmlFor=select.id;
  for(const reviewer of session.reviewers){const option=node('option',reviewer.name);option.value=reviewer.id;select.append(option);}
  const key=node('input');key.id='reviewer-access-key';key.type='password';key.required=true;key.maxLength=128;key.autocomplete='current-password';const keyLabel=node('label','로그인 키');keyLabel.htmlFor=key.id;
  const submit=node('button','로그인');submit.type='submit';form.append(label,select,keyLabel,key,node('p','로그인 키는 별도로 전달받은 로컬 접속 정보입니다. 보고서나 검토 의견에 붙여 넣지 마세요.','hint'),submit);
  form.onsubmit=async e=>{e.preventDefault();submit.disabled=true;try{await post('/api/review-auth/login',{reviewer_id:select.value,access_key:key.value});key.value='';await load();}catch(error){key.value='';form.append(notice(error.message));submit.disabled=false;}};login.append(form);
 }
 function renderRequest(request,session,health){
  const section=node('section',undefined,'document document-top');content.append(section);
  section.append(node('h3',(request.assigned_reviewer?request.assigned_reviewer.name:'미배정')+' · '+humanStatuses[request.status]),node('p','요청 '+reviewDate(request.created_at)+' · 판정 이력 '+request.revision+'건'));
  for(const text of Object.values(request.scopes))section.append(node('p','검토 범위: '+text));
  if(request.recheck_reasons.length)section.append(notice(request.recheck_reasons.map(r=>reasonNames[r]||r).join(' · ')));
  if(request.superseded_by)section.append(node('p','후속 요청: '+request.superseded_by));
  section.append(savedDownload('고정된 검토 요청 내려받기',request.record_ref));
  if(!request.history.length)section.append(node('p','아직 제출된 판정이 없습니다.'));
  for(const decision of request.history){const details=node('details');details.append(node('summary','판정 '+decision.revision+' · '+decision.actor.name+' · '+reviewDate(decision.recorded_at)+(decision.action==='withdraw'?' · 철회':'')));
   if(decision.action==='withdraw')details.append(node('p',decision.withdrawal_reason));
   else for(const j of decision.judgments){details.append(node('h4',request.scopes[j.scope_id]+' · '+verdictNames[j.verdict]),node('p',j.reason));}
   details.append(savedDownload('판정 원기록',decision.record_ref));section.append(details);
  }
  if(health.read_only||!session.actor||!request.assigned_reviewer||session.actor.id!==request.assigned_reviewer.id){section.append(node('p','배정된 검수자가 로그인해야 판정을 기록할 수 있습니다.','hint'));return;}
  if(!request.recheck_reasons.length)renderDecisionForm(section,request);
  if(request.history.length&&request.history.at(-1).action==='record'){
   const form=node('form'),reason=node('textarea');reason.id=request.id+'-withdraw';reason.required=true;reason.maxLength=5000;const label=node('label','기존 판정 철회 사유');label.htmlFor=reason.id;const submit=node('button','판정 철회 기록');submit.type='submit';form.append(label,reason,submit);
   form.onsubmit=async e=>{e.preventDefault();submit.disabled=true;try{await post('/api/report-reviews/'+request.id+'/decisions',payload(request,'withdraw',[],[],reason.value));await load();}catch(error){form.append(notice(error.message+' · 최신 상태를 다시 확인해 주세요.'));submit.disabled=false;}};section.append(form);
  }
 }
 function payload(request,action,judgments,acknowledgements,reason=''){return {expected_request_sha256:request.record_ref.sha256,expected_revision:request.revision,idempotency_key:crypto.randomUUID(),action,judgments,acknowledgement_ids:acknowledgements,withdrawal_reason:reason};}
 function renderDecisionForm(section,request){
  const form=node('form'),fields=[],acks=[];form.append(node('h4',request.revision?'판정 수정 기록':'범위별 판정 입력'),node('p','각 범위의 판정과 이유를 직접 입력하세요. 미확인 사항은 판단 보류로 남길 수 있습니다.'));
  for(const [id,title] of Object.entries(request.scopes)){const group=node('fieldset'),legend=node('legend',title),select=node('select');select.required=true;select.id=request.id+'-'+id;select.setAttribute('aria-label',title+' 판정');const placeholder=node('option','판정을 선택하세요');placeholder.value='';select.append(placeholder);for(const [value,text] of Object.entries(verdictNames)){const option=node('option',text);option.value=value;select.append(option);}
   const reason=node('textarea');reason.required=true;reason.maxLength=5000;reason.id=select.id+'-reason';const label=node('label','판정 근거·조건·미확인 사항');label.htmlFor=reason.id;group.append(legend,select,label,reason);form.append(group);fields.push({id,select,reason});}
  for(const [id,text] of Object.entries(request.acknowledgements)){const label=node('label',undefined,'review-check'),checkbox=node('input');checkbox.type='checkbox';checkbox.value=id;label.append(checkbox,node('span',text));form.append(label);acks.push(checkbox);}
  const submit=node('button','판정과 근거 기록');submit.type='submit';form.append(submit);let lastBody=null;
  form.onsubmit=async e=>{e.preventDefault();const judgments=fields.map(f=>({scope_id:f.id,verdict:f.select.value,reason:f.reason.value})),ackIds=acks.filter(c=>c.checked).map(c=>c.value);if(judgments.some(j=>j.verdict==='accept_scope')&&ackIds.length!==acks.length){form.append(notice('범위 내 수용을 기록하려면 적용 범위 3개를 확인해 주세요.'));return;}
   const next=payload(request,'record',judgments,ackIds);if(lastBody&&JSON.stringify({...lastBody,idempotency_key:''})===JSON.stringify({...next,idempotency_key:''}))next.idempotency_key=lastBody.idempotency_key;lastBody=next;submit.disabled=true;
   try{await post('/api/report-reviews/'+request.id+'/decisions',next);await load();}catch(error){form.append(notice(error.message+' · 입력은 보존되어 있습니다. 최신 상태를 확인한 뒤 다시 제출해 주세요.'));submit.disabled=false;}};section.append(form);
 }
 load();
});
