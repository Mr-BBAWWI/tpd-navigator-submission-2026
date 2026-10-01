let workflowRenderVersion=0;
const workflowStates={queued:'검토 대기',running:'전문 Agent 검토 중',review_ready:'검토 초안 준비',held:'보완 필요',interrupted:'중단됨',cancelled:'사용자 중단',stale:'입력 변경으로 보류'};
const workflowErrors={AGENT_OUTPUT_SCHEMA:'Agent 답변 형식이 검토 계약과 맞지 않습니다.',STRUCTURED_CITATION_REQUIRES_NULL_QUOTE:'구조화 자료의 인용 방식을 보완해야 합니다.',CITATION_QUOTE_MISMATCH:'인용문이 저장된 원문과 일치하지 않습니다.',REQUIRED_SOURCE_NOT_CITED:'필수 근거의 인용이 누락됐습니다.',UNRESOLVED_CRITIC_FINDINGS:'한 차례 수정 후에도 Critic의 지적이 남았습니다.',WORKFLOW_INPUT_STALE:'입력이 바뀌어 이 검토를 중단했습니다.',WORKFLOW_RUNTIME_CHANGED:'검토 프로그램이 변경돼 이전 실행을 재개할 수 없습니다.',WORKFLOW_EXECUTION_FAILED:'실행을 마치지 못했습니다. 저장된 호출 기록을 확인해 주세요.'};
const pendingNames={phase1:'정식 초기 평가·변경 전후 비교',execution_scope:'계산 조건과 비용 범위',pharmacy_reviewer:'약학 검수 담당과 실제 검수',gpu_result:'실제 GPU 결과'};
async function renderWorkflows(run){
 const ticket=++workflowRenderVersion,panel=$('review-workflows');panel.hidden=false;panel.replaceChildren(node('h2','전문 검토와 연구자 검토 준비'));
 try{
  const dossiers=await api('/api/runs/'+run.id+'/dossiers');
  const all=await Promise.all(dossiers.map(d=>api('/api/dossiers/'+d.id+'/b-inputs')));
  if(ticket!==workflowRenderVersion)return;
  const inputs=all.flat();
  if(!inputs.length){panel.append(node('p','아직 이 탐색에 수용한 B H1/H2 자료가 없습니다.'));return;}
  if(run.replay_only)panel.append(notice('다른 환경에서 가져온 내부 재생 자료입니다. 저장 결과를 열람하며 새 AI 검토는 실행하지 않습니다.'));
  panel.append(node('p','계획 → 문헌·분자 전문 검토 → Critic 대조 → 필요 시 한 번 수정 → 검토 초안. 약학 검수와 계산 승인은 별도입니다.'));
  const chooser=node('select');chooser.setAttribute('aria-label','B 인계 자료 선택');
  for(const input of inputs){const option=node('option',new Date(input.created_at).toLocaleString()+' · '+(input.current_input?'현재 입력':'이전 입력'));option.value=input.id;chooser.append(option);}
  panel.append(chooser);const detail=node('div');panel.append(detail);
  async function showInput(){
   const input=inputs.find(i=>i.id===chooser.value);detail.replaceChildren();
   for(const missing of input.pending_inputs)detail.append(node('p',(pendingNames[missing.id]||missing.id)+' — 준비 대기 · 담당 '+missing.owner,'hint'));
   const button=node('button','이 자료로 전문 Agent 검토');button.disabled=!run.llm_enabled||!input.current_input||(run.workflows||[]).some(w=>['queued','running'].includes(w.state));
   button.onclick=async()=>{button.disabled=true;try{await api('/api/workflows',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({b_input_id:input.id,request_key:crypto.randomUUID()})});lastStamp=null;await refresh();}catch(e){detail.append(notice(e.message));button.disabled=false;}};
   detail.append(button,node('p','실행 시 저장된 원문·A/B 결과가 공모전 API로 전송됩니다. 최대 14회 호출, 180,000토큰, 600초의 단계 간 실행 한도입니다. 저장 결과 열람은 새로 호출하지 않습니다.','hint'));
   detail.append(link('A 문헌–B 관측 대조 JSON','/api/artifacts/'+encodeURIComponent(input.reconciled_ref.artifact_id)+'/download'));
   const jobs=(run.workflows||[]).filter(w=>w.b_input_id===input.id);
   for(const w of jobs){
    const card=node('article',undefined,'document');card.append(node('h3',workflowStates[w.state]||w.state),node('p','실제 요청 '+w.calls+'회 · 확인된 사용량 '+w.total_tokens.toLocaleString()+'토큰 · 실행 시도 '+w.attempts+'회'));
    if(w.usage_status==='unreconciled')card.append(notice('응답 또는 사용량이 확인되지 않은 요청이 있습니다. 자동 재호출하지 않습니다.'));
    if(!w.current_input)card.append(notice('이전 입력의 검토입니다. 현재 계산·승인에 사용할 수 없습니다.'));
    if(w.error_code)card.append(node('p','보류 사유: '+(workflowErrors[w.error_code]||'검토를 완료하지 못했습니다. 호출 기록에서 사유를 확인할 수 있습니다.')));
    if(w.result){const claims=new Map(w.result.claims.map(c=>[c.id,c]));for(const id of w.result.claim_order||[]){const c=claims.get(id);if(c)card.append(node('p',c.candidate_id+' · '+c.interpretation));}for(const finding of w.result.critique?.findings||[])card.append(notice(finding.problem));}
    const downloads=node('ul');for(const [field,label] of [['report_ref','Agent 검토 초안'],['packet_ref','연구자 검토 준비 자료'],['result_ref','검토 결과 JSON']])if(w[field]){const item=node('li');item.append(link(label,'/api/artifacts/'+encodeURIComponent(w[field].artifact_id)+'/download'));downloads.append(item);}if(downloads.childElementCount)card.append(downloads);
    const trace=node('details');trace.append(node('summary','호출 기록과 답변'));if(w.error_code)trace.append(node('p','기록 코드: '+w.error_code));for(const call of w.calls_detail||[]){trace.append(node('p',(call.ordinal+1)+'번째 요청 · '+call.usage_status));if(call.answer_ref)trace.append(link('저장된 답변','/api/artifacts/'+encodeURIComponent(call.answer_ref.artifact_id)+'/download'));}card.append(trace);
    for(const [action,label,enabled] of [['cancel','검토 중단',['queued','running'].includes(w.state)],['resume','저장된 답변부터 재개',w.can_resume&&run.llm_enabled]])if(enabled){const b=node('button',label);b.onclick=async()=>{b.disabled=true;try{await api('/api/workflows/'+w.id+'/'+action,{method:'POST'});lastStamp=null;await refresh();}catch(e){card.append(notice(e.message));}};card.append(b);}
    if(w===jobs[0])detail.append(card);else{const history=node('details');history.append(node('summary','이전 실행 · '+(workflowStates[w.state]||w.state)+' · '+w.calls+'회 호출'),card);detail.append(history);}
   }
  }
  chooser.onchange=showInput;await showInput();
  if(!(run.workflows||[]).some(w=>['queued','running'].includes(w.state))&&!(run.analyses||[]).some(a=>['queued','running'].includes(a.state)))panel.append(link('이 탐색의 내부 재생 패키지 받기','/api/runs/'+run.id+'/replay-package'));
  panel.append(node('p','재생 패키지에는 내부 원문과 답변이 포함됩니다. 공개 배포용 선별본은 아직 아닙니다.','hint'));
 }catch(e){if(ticket===workflowRenderVersion)panel.append(notice(e.message));}
}
