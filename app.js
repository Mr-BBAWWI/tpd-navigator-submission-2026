"use strict";

const CLOUD_ORIGIN = "https://tpd-navigator-live.wlgudyun.chatgpt.site";
const OPERATIONS = {
  candidate_compare: "CRBN·VHL 후보와 linker를 비교해 주세요.",
  failure_analysis: "도킹 실패와 pose 제외 원인을 분석해 주세요.",
  acceptance_plan: "14개 기준을 확인하고 후속 실행 계획을 제안해 주세요."
};
const EXPECTED_OPERATIONS = Object.keys(OPERATIONS);
const MAX_SELECTED = 3;
const MAX_VISIBLE = 50;

const state = {
  apiOrigin: null,
  ready: false,
  running: false,
  configured: false,
  candidates: [],
  selected: new Map(),
  controller: null,
  runToken: 0,
  latestBody: null,
  lastPayload: null
};

const refs = {
  form: document.querySelector("#run-form"),
  operationFieldset: document.querySelector("#operation-fieldset"),
  selectedQuery: document.querySelector("#selected-query"),
  connectionStatus: document.querySelector("#connection-status"),
  refreshButton: document.querySelector("#refresh-button"),
  parentSelect: document.querySelector("#parent-select"),
  e3Select: document.querySelector("#e3-select"),
  candidateDetails: document.querySelector("#candidate-details"),
  candidateSearch: document.querySelector("#candidate-search"),
  candidateList: document.querySelector("#candidate-list"),
  candidateCount: document.querySelector("#candidate-count"),
  catalogSummary: document.querySelector("#catalog-summary"),
  selectedCount: document.querySelector("#selected-count"),
  selectedPills: document.querySelector("#selected-pills"),
  runButton: document.querySelector("#run-button"),
  submitHelp: document.querySelector("#submit-help"),
  resultTitle: document.querySelector("#result-title"),
  resultContent: document.querySelector("#result-content"),
  resultAnnouncer: document.querySelector("#result-announcer"),
  cancelButton: document.querySelector("#cancel-button"),
  retryButton: document.querySelector("#retry-button"),
  downloadButton: document.querySelector("#download-button")
};

function make(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
}

function safeString(value) {
  if (value === null || value === undefined || value === "") return "미제공";
  if (typeof value === "string") return value;
  try {
    return JSON.stringify(value, null, 2);
  } catch (_error) {
    return String(value);
  }
}

function isLoopback(hostname) {
  return hostname === "localhost" || hostname === "127.0.0.1" || hostname === "[::1]";
}

function validateApiBase(raw) {
  let url;
  try {
    url = new URL(raw);
  } catch (_error) {
    throw new Error("public-config.json의 apiBase가 올바른 URL이 아닙니다.");
  }

  if (url.username || url.password || url.search || url.hash || url.pathname !== "/") {
    throw new Error("apiBase는 인증 정보, 경로, 쿼리 또는 해시가 없는 루트 origin이어야 합니다.");
  }

  const pageLoopback = isLoopback(window.location.hostname);
  const targetLoopback = isLoopback(url.hostname);
  const sameOrigin = window.location.origin !== "null" && url.origin === window.location.origin;
  const approvedCloud = url.origin === CLOUD_ORIGIN;
  const allowedLoopback = pageLoopback && targetLoopback && url.protocol === "http:";
  const secureAllowed = url.protocol === "https:" && (sameOrigin || approvedCloud);

  if (!secureAllowed && !allowedLoopback) {
    throw new Error("허용되지 않은 API origin입니다. 배포 origin 또는 로컬 개발 origin만 사용할 수 있습니다.");
  }
  return url.origin;
}

async function loadApiOrigin() {
  const controller = new AbortController();
  const timer = window.setTimeout(() => controller.abort(), 15000);
  try {
    const response = await fetch("./public-config.json", {
      cache: "no-store",
      credentials: "omit",
      signal: controller.signal
    });

  if (response.status === 404) {
    return validateApiBase(window.location.origin + "/");
  }
  if (!response.ok) {
    throw new Error(`public-config.json을 읽지 못했습니다. HTTP ${response.status}`);
  }

  let config;
  try {
    config = await response.json();
  } catch (error) {
    if (error && error.name === "AbortError") throw error;
    throw new Error("public-config.json이 올바른 JSON이 아닙니다.");
  }

  if (!config || typeof config !== "object" || Array.isArray(config)) {
    throw new Error("public-config.json은 객체여야 합니다.");
  }
  const keys = Object.keys(config);
  if (!Object.prototype.hasOwnProperty.call(config, "apiBase")) {
    if (keys.length === 0) return validateApiBase(window.location.origin + "/");
    throw new Error("public-config.json에 apiBase가 없습니다.");
  }
  if (typeof config.apiBase !== "string" || !config.apiBase.trim()) {
    throw new Error("public-config.json의 apiBase가 비어 있습니다.");
  }
  return validateApiBase(config.apiBase.trim());
  } catch (error) {
    if (error && error.name === "AbortError") {
      throw new Error("public-config.json 요청 시간이 15초를 초과했습니다.");
    }
    throw error;
  } finally {
    window.clearTimeout(timer);
  }
}

function endpoint(path) {
  return new URL(path, state.apiOrigin).toString();
}

async function fetchJsonWithTimeout(path, timeoutMs) {
  const controller = new AbortController();
  const timer = window.setTimeout(() => controller.abort(), timeoutMs);
  try {
    const response = await fetch(endpoint(path), {
      cache: "no-store",
      credentials: "omit",
      signal: controller.signal
    });
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    return await response.json();
  } catch (error) {
    if (error && error.name === "AbortError") throw new Error("연결 확인 시간이 초과되었습니다.");
    throw error;
  } finally {
    window.clearTimeout(timer);
  }
}

function validateHealth(data) {
  if (!data || data.status !== "ok" || data.service !== "tpd-live-demo") return false;
  if (!Array.isArray(data.operations) || typeof data.configured !== "boolean") return false;
  if (data.execution_mode !== "live_analysis_of_saved_evidence" || data.scientific_approval !== false) return false;
  return EXPECTED_OPERATIONS.every((operation) => data.operations.includes(operation));
}

function validateCatalog(data) {
  if (!data || data.status !== "ok" || data.scientific_approval !== false) return null;
  if (!Array.isArray(data.parents) || !Array.isArray(data.candidates)) return null;

  const parentIds = new Set();
  for (const parent of data.parents) {
    if (!parent || typeof parent.parent_id !== "string" || !parent.parent_id) return null;
    parentIds.add(parent.parent_id);
  }

  const records = [];
  const recordKeys = new Set();
  for (const candidate of data.candidates) {
    if (!candidate || typeof candidate.candidate_id !== "string" || !candidate.candidate_id) return null;
    if (typeof candidate.parent_id !== "string" || !parentIds.has(candidate.parent_id)) return null;
    if (candidate.e3_type !== "CRBN" && candidate.e3_type !== "VHL") return null;
    if (typeof candidate.linker_id !== "string" || !candidate.linker_id) return null;
    const key = JSON.stringify([candidate.parent_id, candidate.e3_type, candidate.linker_id, candidate.candidate_id]);
    if (recordKeys.has(key)) continue;
    recordKeys.add(key);
    records.push({
      key,
      candidate_id: candidate.candidate_id,
      parent_id: candidate.parent_id,
      e3_type: candidate.e3_type,
      linker_id: candidate.linker_id
    });
  }
  return { parentIds: Array.from(parentIds), records };
}

async function initialize() {
  if (state.running) return;
  state.ready = false;
  state.configured = false;
  updateAvailability();
  refs.refreshButton.disabled = true;
  refs.connectionStatus.className = "connection-status";
  refs.connectionStatus.textContent = "서버와 카탈로그를 확인하고 있습니다.";

  try {
    state.apiOrigin = await loadApiOrigin();
    const [health, catalogRaw] = await Promise.all([
      fetchJsonWithTimeout("/api/health", 15000),
      fetchJsonWithTimeout("/api/catalog", 15000)
    ]);
    const catalog = validateCatalog(catalogRaw);
    if (!validateHealth(health)) throw new Error("서버 상태 응답이 약속된 형식과 일치하지 않습니다.");
    if (!catalog) throw new Error("카탈로그 응답이 약속된 형식과 일치하지 않습니다.");

    state.configured = health.configured;
    state.candidates = catalog.records;
    state.selected.clear();
    populateParents(catalog.parentIds);
    renderCandidates();
    renderSelected();
    refs.catalogSummary.textContent = `Parent ${catalog.parentIds.length}개 · 후보 ${catalog.records.length}개`;

    if (health.configured) {
      state.ready = true;
      refs.connectionStatus.className = "connection-status is-good";
      refs.connectionStatus.textContent = "실행 서버와 저장 근거 카탈로그에 연결되었습니다.";
    } else {
      refs.connectionStatus.className = "connection-status is-waiting";
      refs.connectionStatus.textContent = "실행 서버 설정을 기다리고 있습니다.";
    }
  } catch (error) {
    state.candidates = [];
    state.selected.clear();
    renderCandidates();
    renderSelected();
    refs.catalogSummary.textContent = "카탈로그 연결 실패";
    refs.connectionStatus.className = "connection-status is-bad";
    refs.connectionStatus.textContent = `연결을 확인하지 못했습니다. ${error instanceof Error ? error.message : "알 수 없는 오류"}`;
  } finally {
    refs.refreshButton.disabled = state.running;
    updateAvailability();
  }
}

function populateParents(parentIds) {
  const current = refs.parentSelect.value;
  refs.parentSelect.replaceChildren();
  const all = make("option", "", "전체");
  all.value = "";
  refs.parentSelect.append(all);
  parentIds.sort((a, b) => a.localeCompare(b)).forEach((parentId) => {
    const option = make("option", "", parentId);
    option.value = parentId;
    refs.parentSelect.append(option);
  });
  if (parentIds.includes(current)) refs.parentSelect.value = current;
}

function candidateLabel(candidate) {
  return `${candidate.parent_id} · ${candidate.e3_type} · linker ${candidate.linker_id} · ${candidate.candidate_id}`;
}

function matchesFilters(candidate) {
  const parent = refs.parentSelect.value;
  const e3 = refs.e3Select.value;
  return (!parent || candidate.parent_id === parent) && (!e3 || candidate.e3_type === e3);
}

function renderCandidates() {
  refs.candidateList.replaceChildren();
  const query = refs.candidateSearch.value.trim().toLocaleLowerCase("ko");
  const filtered = state.candidates.filter((candidate) => {
    if (!matchesFilters(candidate)) return false;
    return !query || candidateLabel(candidate).toLocaleLowerCase("ko").includes(query);
  });
  const visible = filtered.slice(0, MAX_VISIBLE);
  refs.candidateCount.textContent = `조건과 검색에 맞는 후보 ${filtered.length}개 중 ${visible.length}개 표시 · 전체 ${state.candidates.length}개`;

  if (visible.length === 0) {
    refs.candidateList.append(make("p", "empty-list", "표시할 후보가 없습니다."));
    return;
  }

  const selectedIds = new Set(Array.from(state.selected.values(), (item) => item.candidate_id));
  visible.forEach((candidate) => {
    const label = make("label", "candidate-option");
    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.value = candidate.key;
    checkbox.checked = state.selected.has(candidate.key);
    const duplicateId = selectedIds.has(candidate.candidate_id) && !checkbox.checked;
    checkbox.disabled = state.running || duplicateId || (state.selected.size >= MAX_SELECTED && !checkbox.checked);
    checkbox.addEventListener("change", () => toggleCandidate(candidate, checkbox.checked));
    label.append(checkbox, make("span", "", candidateLabel(candidate)));
    refs.candidateList.append(label);
  });
}

function toggleCandidate(candidate, checked) {
  if (state.running) return;
  if (checked) {
    const duplicate = Array.from(state.selected.values()).some((item) => item.candidate_id === candidate.candidate_id);
    if (duplicate || state.selected.size >= MAX_SELECTED) {
      renderCandidates();
      return;
    }
    state.selected.set(candidate.key, candidate);
  } else {
    state.selected.delete(candidate.key);
  }
  renderSelected();
  renderCandidates();
}

function renderSelected() {
  refs.selectedPills.replaceChildren();
  refs.selectedCount.textContent = `${state.selected.size} / ${MAX_SELECTED}`;
  if (state.selected.size === 0) {
    refs.selectedPills.append(make("span", "empty-selection", "선택 사항입니다."));
    return;
  }

  state.selected.forEach((candidate) => {
    const pill = make("span", "candidate-pill");
    pill.append(make("span", "", candidateLabel(candidate)));
    const remove = make("button", "pill-remove", "제거");
    remove.type = "button";
    remove.disabled = state.running;
    remove.setAttribute("aria-label", `${candidateLabel(candidate)} 제거`);
    remove.addEventListener("click", () => {
      state.selected.delete(candidate.key);
      renderSelected();
      renderCandidates();
    });
    pill.append(remove);
    refs.selectedPills.append(pill);
  });
}

function reconcileSelection() {
  for (const [key, candidate] of state.selected) {
    if (!matchesFilters(candidate)) state.selected.delete(key);
  }
  renderSelected();
  renderCandidates();
}

function selectedOperation() {
  const checked = refs.form.querySelector('input[name="operation"]:checked');
  return checked ? checked.value : "candidate_compare";
}

function updateOperationCards() {
  const operation = selectedOperation();
  refs.selectedQuery.textContent = OPERATIONS[operation];
  document.querySelector("#query-scope").textContent = {
    candidate_compare: "부모·E3·후보 필터의 교집합에 해당하는 후보를 선택합니다.",
    failure_analysis: "필터로 부모 범위를 선택한 뒤 저장된 부모별 실패 카운터를 집계합니다. CRBN과 VHL의 실패 귀속이나 비교가 아닙니다.",
    acceptance_plan: "14개 기준 상태는 전역이며 필터로 바뀌지 않습니다. 필터는 선택 범위의 후보 수에만 영향을 줍니다.",
  }[operation];
  document.querySelectorAll(".operation-card").forEach((card) => {
    const radio = card.querySelector('input[type="radio"]');
    card.classList.toggle("is-selected", Boolean(radio && radio.checked));
  });
}

function buildPayload() {
  const filters = {};
  if (refs.parentSelect.value) filters.parent_id = refs.parentSelect.value;
  if (refs.e3Select.value) filters.e3_type = refs.e3Select.value;
  const candidateIds = Array.from(new Set(Array.from(state.selected.values(), (item) => item.candidate_id))).slice(0, MAX_SELECTED);
  if (candidateIds.length) filters.candidate_ids = candidateIds;
  return { operation: selectedOperation(), filters };
}

function setFormDisabled(disabled) {
  refs.form.querySelectorAll("input, select, button").forEach((control) => {
    control.disabled = disabled;
  });
  refs.candidateSearch.disabled = disabled;
  renderSelected();
  renderCandidates();
}

function updateAvailability() {
  if (state.running) {
    refs.runButton.disabled = true;
    refs.submitHelp.textContent = "현재 요청의 응답을 기다리고 있습니다.";
    return;
  }
  refs.runButton.disabled = !state.ready;
  refs.refreshButton.disabled = false;
  if (state.ready) refs.submitHelp.textContent = "선택한 작업과 필터만 서버로 전송합니다.";
  else if (!state.configured) refs.submitHelp.textContent = "실행 서버 설정을 기다리고 있습니다.";
  else refs.submitHelp.textContent = "서버 연결을 다시 확인해 주세요.";
}

function resetResultActions() {
  refs.cancelButton.hidden = true;
  refs.retryButton.hidden = true;
  refs.downloadButton.hidden = true;
}

function showWaiting(started) {
  refs.resultContent.replaceChildren();
  const box = make("div", "waiting-result");
  box.append(make("strong", "", "요청을 보냈습니다. 서버의 도구 분석과 모델 응답을 기다리고 있습니다."));
  const timer = make("p", "waiting-time", "브라우저 대기 시간: 0초");
  box.append(timer);
  refs.resultContent.append(box);
  refs.resultAnnouncer.textContent = "요청을 전송했습니다. 응답을 기다립니다.";
  refs.cancelButton.hidden = false;
  return window.setInterval(() => {
    timer.textContent = `브라우저 대기 시간: ${Math.floor((Date.now() - started) / 1000)}초`;
  }, 1000);
}

function appendDefinitionList(container, entries) {
  const list = make("dl", "metadata-grid");
  entries.forEach(([term, value]) => {
    list.append(make("dt", "", term), make("dd", "", safeString(value)));
  });
  container.append(list);
}

function elapsedFromServer(startedAt, completedAt) {
  const start = Date.parse(startedAt);
  const end = Date.parse(completedAt);
  if (!Number.isFinite(start) || !Number.isFinite(end) || end < start) return "미제공";
  return `${((end - start) / 1000).toFixed(2)}초`;
}

function addTrace(container, trace, open) {
  if (!Array.isArray(trace) || trace.length === 0) return;
  const details = make("details", "technical-details");
  details.open = Boolean(open);
  details.append(make("summary", "", `도구 추적 ${trace.length}건`));
  trace.forEach((entry, index) => {
    const card = make("article", "trace-card");
    card.append(make("h4", "", `${index + 1}. ${safeString(entry && entry.tool)}`));
    card.append(make("p", "trace-time", `계산 시각: ${safeString(entry && entry.computed_at)}`));
    card.append(make("h5", "", "입력"), make("pre", "trace-data", safeString(entry && entry.input)));
    card.append(make("h5", "", "출력"), make("pre", "trace-data", safeString(entry && entry.output)));
    details.append(card);
  });
  container.append(details);
}

function showSuccess(body, browserElapsed) {
  refs.resultContent.replaceChildren();
  const answerSection = make("article", "answer-card");
  answerSection.append(make("p", "result-kicker", "모델의 새 답변"));
  answerSection.append(make("div", "answer-text", body.answer.text));
  refs.resultContent.append(answerSection);

  if (Array.isArray(body.answer.evidence_refs) && body.answer.evidence_refs.length) {
    const evidence = make("section", "evidence-section");
    evidence.append(make("h3", "", "응답의 근거 참조"));
    const list = make("ul", "evidence-list");
    body.answer.evidence_refs.forEach((item) => list.append(make("li", "", safeString(item))));
    evidence.append(list);
    refs.resultContent.append(evidence);
  }

  if (Array.isArray(body.tool_trace) && body.tool_trace.length) {
    const summaries = make("section", "tool-summary");
    summaries.append(make("h3", "", "실제 도구 출력 요약"));
    const grid = make("div", "summary-grid");
    body.tool_trace.forEach((entry) => {
      const card = make("article", "summary-card");
      card.append(make("strong", "", safeString(entry && entry.tool)));
      card.append(make("pre", "", safeString(entry && entry.output)));
      grid.append(card);
    });
    summaries.append(grid);
    refs.resultContent.append(summaries);
  }

  const receipt = make("details", "technical-details receipt-details");
  receipt.append(make("summary", "", "실행 영수증 및 모델 메타데이터"));
  appendDefinitionList(receipt, [
    ["요청 ID", body.request_id],
    ["작업", body.operation],
    ["실행 모드", body.execution_mode],
    ["요청 모델", body.model && body.model.requested],
    ["반환 모델", body.model && body.model.returned],
    ["응답 ID", body.model && body.model.response_id],
    ["입력 토큰", body.usage && body.usage.input_tokens],
    ["출력 토큰", body.usage && body.usage.output_tokens],
    ["전체 토큰", body.usage && body.usage.total_tokens],
    ["근거 SHA-256", body.evidence_sha256],
    ["서버 시작", body.started_at],
    ["서버 완료", body.completed_at],
    ["서버 경과 시간", elapsedFromServer(body.started_at, body.completed_at)],
    ["브라우저 대기 시간", `${(browserElapsed / 1000).toFixed(2)}초`],
    ["과학적 승인", body.scientific_approval === false ? "수행하지 않음" : "미제공"]
  ]);
  refs.resultContent.append(receipt);
  addTrace(refs.resultContent, body.tool_trace, false);
  refs.downloadButton.hidden = false;
  refs.resultAnnouncer.textContent = "새 분석이 완료되었습니다.";
  refs.resultTitle.focus({ preventScroll: true });
}

function validCompleted(body, operation) {
  return body && body.status === "completed" && body.operation === operation &&
    body.execution_mode === "live_analysis_of_saved_evidence" &&
    body.live_model_called === true && body.scientific_approval === false &&
    body.answer && typeof body.answer.text === "string";
}

function showFailure(body, httpStatus, fallbackMessage) {
  refs.resultContent.replaceChildren();
  const card = make("div", "failure-result");
  card.append(make("strong", "", "요청을 완료하지 못했습니다."));

  const hasFailureBody = body && body.status === "failed" && body.error && typeof body.error === "object";
  const code = hasFailureBody ? body.error.code : (httpStatus ? `HTTP_${httpStatus}` : "CONNECTION_ERROR");
  const message = hasFailureBody && typeof body.error.message === "string" ? body.error.message : fallbackMessage;
  appendDefinitionList(card, [
    ["오류 코드", code],
    ["메시지", message || "서버 응답을 확인할 수 없습니다."],
    ["요청 ID", body && body.request_id],
    ["HTTP 상태", httpStatus || "미제공"]
  ]);

  if (httpStatus === 400) card.append(make("p", "error-guidance", "필터와 후보 선택을 수정한 뒤 실행하거나 같은 요청을 수동으로 다시 시도할 수 있습니다. AMBIGUOUS_CANDIDATE_ID라면 Parent 또는 E3 필터로 후보를 구분하세요."));
  if (httpStatus === 429) card.append(make("p", "error-guidance", "동시 실행 또는 요청 제한일 수 있습니다. 잠시 기다린 뒤 수동으로 다시 시도하세요."));
  refs.resultContent.append(card);

  if (body && body.usage) {
    const usage = make("details", "technical-details");
    usage.append(make("summary", "", "실패 응답 사용량"));
    appendDefinitionList(usage, [
      ["입력 토큰", body.usage.input_tokens],
      ["출력 토큰", body.usage.output_tokens],
      ["전체 토큰", body.usage.total_tokens]
    ]);
    refs.resultContent.append(usage);
  }
  addTrace(refs.resultContent, body && body.tool_trace, true);
  refs.retryButton.hidden = !state.lastPayload;
  refs.downloadButton.hidden = !body;
  refs.resultAnnouncer.textContent = "요청이 완료되지 않았습니다.";
  refs.resultTitle.focus({ preventScroll: true });
}

function showCanceled(timedOut) {
  refs.resultContent.replaceChildren();
  const card = make("div", timedOut ? "failure-result" : "canceled-result");
  if (timedOut) {
    card.append(make("strong", "", "브라우저의 180초 응답 대기 시간이 초과되었습니다."));
    card.append(make("p", "", "연결 대기는 중지했지만 서버 작업은 이미 진행 중일 수 있습니다. 자동으로 다시 요청하지 않습니다."));
  } else {
    card.append(make("strong", "", "응답 기다리기를 중지했습니다."));
    card.append(make("p", "", "서버 작업은 이미 진행 중일 수 있습니다."));
  }
  refs.resultContent.append(card);
  refs.retryButton.hidden = !state.lastPayload;
  refs.resultAnnouncer.textContent = timedOut ? "응답 대기 시간이 초과되었습니다." : "응답 기다리기를 중지했습니다.";
}

async function execute(payload) {
  if (state.running || !state.ready) return;
  state.running = true;
  state.lastPayload = JSON.parse(JSON.stringify(payload));
  state.latestBody = null;
  const token = ++state.runToken;
  const controller = new AbortController();
  state.controller = controller;
  let timedOut = false;
  const started = Date.now();
  const waitingInterval = showWaiting(started);
  resetResultActions();
  refs.cancelButton.hidden = false;
  setFormDisabled(true);
  refs.refreshButton.disabled = true;
  updateAvailability();

  const timeout = window.setTimeout(() => {
    timedOut = true;
    controller.abort();
  }, 180000);

  try {
    const response = await fetch(endpoint("/api/run"), {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
      cache: "no-store",
      credentials: "omit",
      signal: controller.signal
    });

    const raw = await response.text();
    let body = null;
    if (raw) {
      try {
        body = JSON.parse(raw);
      } catch (_error) {
        body = null;
      }
    }
    if (token !== state.runToken) return;
    state.latestBody = body;

    if (response.ok && validCompleted(body, payload.operation)) {
      showSuccess(body, Date.now() - started);
    } else if (body && body.status === "failed") {
      showFailure(body, response.status, `서버가 HTTP ${response.status}로 요청을 거절했습니다.`);
    } else if (!response.ok) {
      showFailure(body, response.status, `서버가 HTTP ${response.status} 응답을 반환했습니다. HTML 또는 일반 텍스트 오류 본문은 안전을 위해 표시하지 않습니다.`);
    } else {
      showFailure(body, response.status, "완료 응답이 약속된 형식 또는 실행 조건과 일치하지 않습니다.");
    }
  } catch (error) {
    if (token !== state.runToken) return;
    if (error && error.name === "AbortError") {
      showCanceled(timedOut);
    } else {
      showFailure(null, null, "브라우저가 실행 서버에 연결하지 못했습니다. 서버 실패로 단정할 수 없으며 자동 재시도하지 않습니다.");
    }
  } finally {
    window.clearTimeout(timeout);
    window.clearInterval(waitingInterval);
    if (token === state.runToken) {
      state.running = false;
      state.controller = null;
      refs.cancelButton.hidden = true;
      setFormDisabled(false);
      refs.refreshButton.disabled = false;
      updateAvailability();
    }
  }
}

refs.form.addEventListener("change", (event) => {
  if (event.target && event.target.name === "operation") updateOperationCards();
});
refs.form.addEventListener("submit", (event) => {
  event.preventDefault();
  execute(buildPayload());
});
refs.parentSelect.addEventListener("change", reconcileSelection);
refs.e3Select.addEventListener("change", reconcileSelection);
refs.candidateSearch.addEventListener("input", renderCandidates);
refs.refreshButton.addEventListener("click", initialize);
refs.cancelButton.addEventListener("click", () => {
  if (state.running && state.controller) state.controller.abort();
});
refs.retryButton.addEventListener("click", () => {
  if (!state.running && state.lastPayload) execute(JSON.parse(JSON.stringify(state.lastPayload)));
});
refs.downloadButton.addEventListener("click", () => {
  if (!state.latestBody) return;
  const blob = new Blob([JSON.stringify(state.latestBody, null, 2)], { type: "application/json;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = `tpd-run-${safeString(state.latestBody.request_id).replace(/[^a-zA-Z0-9._-]/g, "_")}.json`;
  document.body.append(anchor);
  anchor.click();
  anchor.remove();
  URL.revokeObjectURL(url);
});

updateOperationCards();
initialize();
