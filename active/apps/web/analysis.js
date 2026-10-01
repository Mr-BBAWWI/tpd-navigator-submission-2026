const selectedSegments = new Map();
const analysisLabels = {queued:'분석 대기',running:'AI 분석 중',review_ready:'연구자 검토 대기',held:'보완 필요',failed:'분석 보류',stale:'입력 변경으로 보류',interrupted:'중단됨'};
let analysisRenderVersion = 0;
function choices(run){if(!selectedSegments.has(run.id))selectedSegments.set(run.id,new Set());return selectedSegments.get(run.id);}
function updateSelectionCount(run){const count=$('selected-count');if(count)count.textContent=choices(run).size+'개 선택 / 최대 12개';const button=$('analyze-button');if(button)button.disabled=!run.llm_enabled||!choices(run).size||choices(run).size>12||(run.analyses||[]).some(a=>['queued','running'].includes(a.state));}
function addAnalysisSelection(part,segment,run){const label=node('label',undefined,'analysis-choice'),box=node('input');box.type='checkbox';box.checked=choices(run).has(segment.segment_id);box.disabled=!run.llm_enabled;box.onchange=()=>{box.checked?choices(run).add(segment.segment_id):choices(run).delete(segment.segment_id);updateSelectionCount(run);};label.append(box,document.createTextNode('이 구간을 AI 분석에 포함'));part.append(label);}
async function renderAnalysis(run){
 renderDossiers(run);
 renderWorkflows(run);
 const version=++analysisRenderVersion,panel=$('analysis-panel'),results=$('analysis-results');panel.replaceChildren();results.replaceChildren();panel.hidden=!run.bundle;results.hidden=true;if(!run.bundle)return;
 panel.append(node('h2','원문 근거 분석'),node('p','아래 문헌을 펼쳐 분석할 구간을 선택하세요. 문헌 분석가가 관측을 추출하고 별도 검토자가 원문과 대조합니다. 그림 속 구조와 표의 배치는 해석하지 않습니다.'));
 const count=node('p',undefined,'hint');count.id='selected-count';panel.append(count);
 const button=node('button','선택한 구간 분석');button.id='analyze-button';button.onclick=async()=>{button.disabled=true;try{await api('/api/runs/'+run.id+'/analyses',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({revision:run.revision,segment_ids:Array.from(choices(run)),request_key:crypto.randomUUID()})});lastStamp=null;await refresh();}catch(e){panel.append(notice(e.message));updateSelectionCount(run);}};panel.append(button);
 if(!run.llm_enabled)panel.append(notice('이 서버는 AI 분석이 꺼져 있습니다. 저장된 분석 결과는 계속 열람할 수 있습니다.'));
 panel.append(node('p','선택한 원문은 공모전 API로 전송됩니다. 실행당 최대 2회 호출하며, 결과를 다시 여는 동안에는 새로 호출하지 않습니다.','hint'));
 updateSelectionCount(run);
 const jobs=run.analyses||[];
 if(jobs.length){const latest=jobs[0];panel.append(node('p',analysisLabels[latest.state]||latest.state,'badge'));if(latest.stage==='Critic'&&latest.state==='running')panel.append(node('p','추출된 주장과 원문을 별도로 대조하고 있습니다.'));}
 const gate=node('div',undefined,'notice');gate.append(node('strong','G1 검토 준비'),node('p','B의 정식 후보·변경 전후 검토 결과, 실행 비용 범위와 연구자 검토 기록이 필요합니다. 현재 결과로 계산을 승인하거나 실행하지 않습니다.'));panel.append(gate);
 if(!jobs.length)return;
 results.hidden=false;
 const chooser=node('select');chooser.setAttribute('aria-label','저장된 분석 선택');for(const job of jobs){const option=node('option',(analysisLabels[job.state]||job.state)+' · '+new Date(job.created_at).toLocaleString());option.value=job.id;chooser.append(option);}results.append(node('h2','저장된 분석과 검토'),chooser);
 const body=node('div');results.append(body);
 let detailVersion=0;
 async function show(id){const ticket=++detailVersion;body.replaceChildren(node('p','저장된 결과를 읽고 있습니다.'));try{const result=await api('/api/analyses/'+id);if(version!==analysisRenderVersion||ticket!==detailVersion)return;body.replaceChildren(node('p',result.message));
  if(!result.current_input)body.append(notice('이 결과는 현재 입력과 다릅니다. 이전 결과로만 열람해 주세요.'));
  body.append(node('p','호출 '+result.calls+'회 · 확인된 사용량 '+result.total_tokens.toLocaleString()+' 토큰'+(result.usage_status==='unreconciled'?' · 일부 호출 사용량 미확인':''),'hint'));
  if(result.error_code)body.append(notice('보류 코드: '+result.error_code));
  if(result.assessment){const b=result.annotated_bundle;body.append(node('p','전체 '+b.segments.length+'개 중 AI 전달 '+b.coverage.llm_read_segment_count+'개 · 독립 실험 수 '+(b.coverage.independent_experiment_count.state==='known'?b.coverage.independent_experiment_count.value:'미확인')));
   const kinds={experiment:'문헌 실험',computation:'문헌 계산',author_interpretation:'저자 해석',researcher_hypothesis:'가설'};
   for(const claim of result.assessment.claims){const card=node('article',undefined,'document');card.append(node('small',kinds[claim.basis_kind]),node('h3',claim.statement));
    for(const citation of claim.citations){const segment=b.segments.find(s=>s.segment_id===citation.segment_id),doc=b.documents.find(d=>d.document_id===segment.locator.document_id);card.append(node('p',doc.title.value||'제목 미확인','hint'),node('blockquote',citation.quote),node('small',segment.locator.xpath||('페이지 '+segment.locator.page)));}
    for(const record of b.evidence_records.filter(r=>claim.evidence_ids.includes(r.evidence_id))){const names=b.compounds.filter(c=>record.compound_ids.includes(c.compound_id)).map(c=>c.source_label);if(names.length)card.append(node('p','원문 화합물: '+names.join(', ')+' · 구조 동일성은 B 검토 대기'));for(const o of record.observations)card.append(node('p',o.endpoint_label+': '+(o.measurement.raw_text||o.measurement.reason)));const ctx=record.assay_context;card.append(node('p','세포주: '+(ctx.cell_line.value||'미보고')+' · 노출 시간: '+(ctx.duration.raw_text||ctx.duration.reason),'hint'));}
    claim.limitations.forEach(x=>card.append(node('p',x,'hint')));body.append(card);
   }
   result.assessment.limitations.forEach(x=>body.append(notice(x)));
  }
  if(result.review){body.append(node('h3','별도 Critic 검토'));if(!result.review.findings.length)body.append(node('p','추가 지적이 없습니다. 약학적 정확성에 대한 사람 검수는 남아 있습니다.'));for(const f of result.review.findings)body.append(notice((f.severity==='blocking'?'보완 필요: ':'검토 참고: ')+f.message));result.review.limitations.forEach(x=>body.append(node('p',x,'hint')));}
  const records=node('details');records.append(node('summary','AI 답변 원문과 실행 기록'));for(const trace of result.trace){records.append(link(trace.role+' 답변 JSON','/api/artifacts/'+encodeURIComponent(trace.response_ref.artifact_id)+'/download'));}body.append(records);
 }catch(e){if(version===analysisRenderVersion)body.replaceChildren(notice(e.message));}}
 chooser.onchange=()=>show(chooser.value);await show(chooser.value);
}
