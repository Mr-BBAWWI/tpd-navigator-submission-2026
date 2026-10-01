import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { createWorker } from './worker.mjs';
import { computeEvidence, validateRequest } from './evidence-tools.mjs';

const receipt = JSON.parse(await readFile(new URL('./campaign-public.json', import.meta.url), 'utf8'));
const origin = 'https://demo.example';
const env = { ALLOWED_ORIGINS: origin, DACON_API_KEY: 'test-secret-value' };
const headers = { origin, 'content-type': 'application/json', 'cf-connecting-ip': '203.0.113.8' };
const run = (worker, body, extra = {}) => worker.fetch(new Request('https://worker.example/api/run', { method: 'POST', headers: { ...headers, ...(extra.headers || {}) }, body: typeof body === 'string' ? body : JSON.stringify(body) }), extra.env || env);
const answerFor = refs => JSON.stringify({ text: '저장된 근거만 해석했으며 과학적 승인을 의미하지 않습니다.', evidence_refs: [refs[0]] });
function successMock(calls, delay) {
  return async (url, init) => {
    calls.push({ url, init, payload: JSON.parse(init.body) });
    if (delay) await delay();
    const context = JSON.parse(JSON.parse(init.body).input);
    return new Response(JSON.stringify({ id: 'resp_test_1', model: 'gpt-5.6-sol', status: 'completed', usage: { input_tokens: 10, output_tokens: 5, total_tokens: 15 }, output: [{ type: 'message', role: 'assistant', content: [{ type: 'output_text', text: answerFor(context.allowed_evidence_refs) }] }] }), { status: 200 });
  };
}

test('full receipt numeric invariants are independently recomputed', () => {
  assert.equal(receipt.evidence.candidates.length, 324);
  assert.equal(receipt.evidence.candidates.filter(c => c.e3_type === 'CRBN').length, 162);
  assert.equal(receipt.evidence.candidates.filter(c => c.e3_type === 'VHL').length, 162);
  assert.equal(receipt.evidence.parents.length, 9);
  const sums = receipt.evidence.parents.reduce((a, p) => {
    const c = p.counts_recomputed_from_records;
    a.attempts += c.docking_attempts; a.completed += c.docking_completed; a.failed += c.failed_docking; a.pose += c.computed_but_pose_filter_rejected; return a;
  }, { attempts: 0, completed: 0, failed: 0, pose: 0 });
  assert.deepEqual(sums, { attempts: 234, completed: 231, failed: 3, pose: 91 });
  const modified = structuredClone(receipt);
  modified.evidence.campaign_counts_recomputed_from_records.docking_attempts = 999;
  const n = validateRequest({ operation: 'failure_analysis' }, modified);
  assert.equal(computeEvidence(modified, n).counts.docking_attempts, 234);
});

test('candidate compare is deterministic, lexical, limited to three and preserves complete identity', () => {
  const n = validateRequest({ operation: 'candidate_compare' }, receipt);
  const a = computeEvidence(receipt, n); const b = computeEvidence(receipt, n);
  assert.deepEqual(a, b); assert.equal(a.counts.matched_candidates, 324); assert.equal(a.displayed_candidates.length, 3);
  for (const c of a.displayed_candidates) {
    const source = receipt.evidence.candidates.find(x => x.candidate_key === c.candidate_key);
    assert.equal(c.canonical_smiles, source.canonical_smiles);
    assert.ok(['CRBN', 'VHL'].includes(c.e3_type));
  }
  assert.deepEqual(a.displayed_candidates.map(x => x.candidate_key), [...a.displayed_candidates.map(x => x.candidate_key)].sort());
});

test('validation rejects unknown, dangerous, duplicate, foreign and ambiguous IDs', () => {
  assert.throws(() => validateRequest({ operation: 'candidate_compare', url: 'https://evil' }, receipt));
  assert.throws(() => validateRequest({ operation: 'candidate_compare', filters: { command: 'run' } }, receipt));
  assert.throws(() => validateRequest({ operation: 'candidate_compare', filters: { candidate_ids: ['foreign'] } }, receipt));
  const id = receipt.evidence.candidates[0].candidate_id;
  assert.throws(() => validateRequest({ operation: 'candidate_compare', filters: { candidate_ids: [id, id] } }, receipt));
  const duplicate = structuredClone(receipt); duplicate.evidence.candidates.push({ ...duplicate.evidence.candidates[0], candidate_key: duplicate.evidence.candidates[0].candidate_key + '/duplicate' });
  assert.throws(() => validateRequest({ operation: 'candidate_compare', filters: { candidate_ids: [id] } }, duplicate), e => e.code === 'AMBIGUOUS_CANDIDATE_ID');
});

test('origin, methods, query, JSON and streamed body limits reject before fetch', async () => {
  let calls = 0; const worker = createWorker({ receipt, fetchImpl: async () => { calls++; throw new Error(); } });
  let r = await worker.fetch(new Request('https://worker.example/api/run', { method: 'POST', headers: { 'content-type': 'application/json' }, body: '{}' }), env); assert.equal(r.status, 403);
  r = await worker.fetch(new Request('https://worker.example/api/run', { method: 'GET', headers: { origin } }), env); assert.equal(r.status, 405);
  r = await worker.fetch(new Request('https://worker.example/api/run?x=1', { method: 'POST', headers, body: '{}' }), env); assert.equal(r.status, 400);
  r = await run(worker, '{bad'); assert.equal(r.status, 400);
  const stream = new ReadableStream({ start(c) { c.enqueue(new TextEncoder().encode('x'.repeat(4097))); c.close(); } });
  r = await worker.fetch(new Request('https://worker.example/api/run', { method: 'POST', headers, body: stream, duplex: 'half' }), env); assert.equal(r.status, 413);
  assert.equal(calls, 0);
});

test('success calls fixed provider with bounded verified context and returns trace, usage and false approval', async () => {
  const calls = []; let now = 1700000000000;
  const worker = createWorker({ receipt, fetchImpl: successMock(calls), now: () => ++now });
  const response = await run(worker, { operation: 'failure_analysis' });
  assert.equal(response.status, 200); const data = await response.json();
  assert.equal(calls.length, 1); assert.equal(calls[0].url, 'https://dacon-apim-hackathon-0903.azure-api.net/hackathon/openai/v1/responses');
  assert.equal(calls[0].payload.model, 'gpt-5.6-sol'); assert.equal(calls[0].payload.store, false); assert.equal(calls[0].init.redirect, 'manual');
  assert.equal(calls[0].init.headers['api-key'], env.DACON_API_KEY); assert.equal(calls[0].init.headers.Authorization, `Bearer ${env.DACON_API_KEY}`);
  assert.ok(new TextEncoder().encode(calls[0].payload.input).length <= 15000);
  const context = JSON.parse(calls[0].payload.input); assert.equal(context.verified_tool_output.counts.docking_attempts, 234);
  assert.equal(data.scientific_approval, false); assert.equal(data.live_model_called, true); assert.deepEqual(data.usage, { input_tokens: 10, output_tokens: 5, total_tokens: 15 });
  assert.equal(data.tool_trace.length, 1); assert.ok(data.answer.evidence_refs.every(x => context.allowed_evidence_refs.includes(x)));
});

test('same input has stable tool hash while different filters change it', async () => {
  const calls = []; let t = 1700000000000; const worker = createWorker({ receipt, fetchImpl: successMock(calls), now: () => t++ });
  const a = await (await run(worker, { operation: 'candidate_compare', filters: { e3_type: 'CRBN' } })).json();
  const b = await (await run(worker, { operation: 'candidate_compare', filters: { e3_type: 'CRBN' } })).json();
  const c = await (await run(worker, { operation: 'candidate_compare', filters: { e3_type: 'VHL' } }, { headers: { 'cf-connecting-ip': 'other' } })).json();
  assert.equal(calls.length, 3); assert.equal(a.evidence_sha256, b.evidence_sha256); assert.notEqual(a.evidence_sha256, c.evidence_sha256); assert.equal(a.snapshot_sha256, c.snapshot_sha256);
  assert.equal(a.evidence_hash_kind, 'canonical_tool_result'); assert.equal(a.snapshot_hash_kind, 'canonical_snapshot');
  assert.deepEqual(a.tool_trace[0].output, b.tool_trace[0].output); assert.notEqual(a.started_at, b.started_at);
});

test('provider errors are safe, incomplete preserves known usage, grounding and secret output are rejected', async () => {
  const rateWorker = createWorker({ receipt, fetchImpl: async () => new Response(env.DACON_API_KEY, { status: 429 }) });
  let d = await (await run(rateWorker, { operation: 'candidate_compare' })).json(); assert.equal(d.error.code, 'UPSTREAM_RATE_LIMIT'); assert.equal(JSON.stringify(d).includes(env.DACON_API_KEY), false); assert.equal(d.answer, undefined);
  const incomplete = createWorker({ receipt, fetchImpl: async () => new Response(JSON.stringify({ id: 'resp_x', model: 'gpt-5.6-sol', status: 'incomplete', usage: { input_tokens: 1, output_tokens: 2, total_tokens: 3 }, output: [] })) });
  d = await (await run(incomplete, { operation: 'candidate_compare' })).json(); assert.equal(d.error.code, 'UPSTREAM_INCOMPLETE'); assert.deepEqual(d.usage, { input_tokens: 1, output_tokens: 2, total_tokens: 3 });
  const invalid = createWorker({ receipt, fetchImpl: async () => new Response(JSON.stringify({ id: 'resp_x', model: 'gpt-5.6-sol', status: 'completed', usage: { input_tokens: 1, output_tokens: 1, total_tokens: 2 }, output: [{ type: 'message', role: 'assistant', content: [{ type: 'output_text', text: JSON.stringify({ text: '근거 해석입니다.', evidence_refs: ['tool:bad:E_X'] }) }] }] })) });
  d = await (await run(invalid, { operation: 'candidate_compare' })).json(); assert.equal(d.error.code, 'API_INVALID_GROUNDING');
  const secret = createWorker({ receipt, fetchImpl: async () => new Response(JSON.stringify({ id: 'resp_x', model: 'gpt-5.6-sol', status: 'completed', usage: { input_tokens: 1, output_tokens: 1, total_tokens: 2 }, output: [{ type: 'message', role: 'assistant', content: [{ type: 'output_text', text: JSON.stringify({ text: `근거 ${env.DACON_API_KEY}`, evidence_refs: ['tool:candidate_compare:E_COUNTS'] }) }] }] })) });
  d = await (await run(secret, { operation: 'candidate_compare' })).json(); assert.equal(d.error.code, 'API_INVALID_ANSWER'); assert.equal(JSON.stringify(d).includes(env.DACON_API_KEY), false);
});

test('timeout and deterministic concurrency release slots', async () => {
  let calls = 0;
  const never = (url, init) => {
    calls++;
    return new Promise((resolve, reject) => init.signal.addEventListener('abort', () => reject(new DOMException('Aborted', 'AbortError')), { once: true }));
  };
  const worker = createWorker({ receipt, fetchImpl: never, timeoutMs: 20 });
  const a = run(worker, { operation: 'candidate_compare' }, { headers: { 'cf-connecting-ip': 'a' } });
  const b = run(worker, { operation: 'candidate_compare' }, { headers: { 'cf-connecting-ip': 'b' } });
  while (calls < 2) await new Promise(resolve => setTimeout(resolve, 1));
  const third = await run(worker, { operation: 'candidate_compare' }, { headers: { 'cf-connecting-ip': 'c' } }); assert.equal(third.status, 429);
  assert.equal((await a).status, 504); assert.equal((await b).status, 504);
  const fourth = await run(worker, { operation: 'candidate_compare' }, { headers: { 'cf-connecting-ip': 'd' } }); assert.equal(fourth.status, 504);
});

test('per-IP four per minute, retry, and expiry', async () => {
  const calls = []; let clock = 1000; const worker = createWorker({ receipt, fetchImpl: successMock(calls), now: () => clock });
  for (let i = 0; i < 4; i++) assert.equal((await run(worker, { operation: 'candidate_compare' })).status, 200);
  let r = await run(worker, { operation: 'candidate_compare' }); assert.equal(r.status, 429); assert.ok(Number(r.headers.get('retry-after')) >= 1);
  clock += 60001; r = await run(worker, { operation: 'candidate_compare' }); assert.equal(r.status, 200);
});

test('health discloses no secret or endpoint and makes no upstream call', async () => {
  let calls = 0; const worker = createWorker({ receipt, fetchImpl: async () => { calls++; } });
  let r = await worker.fetch(new Request('https://worker.example/api/health'), { ALLOWED_ORIGINS: origin }); let text = await r.text(); assert.equal(r.status, 200); assert.equal(JSON.parse(text).configured, false);
  r = await worker.fetch(new Request('https://worker.example/api/health'), env); text = await r.text(); assert.equal(JSON.parse(text).configured, true); assert.equal(text.includes(env.DACON_API_KEY), false); assert.equal(text.includes('azure-api.net'), false); assert.equal(calls, 0);
  r = await worker.fetch(new Request('https://worker.example/api/health', { method: 'POST' }), env); assert.equal(r.status, 405);
  r = await worker.fetch(new Request('https://worker.example/api/health'), { ...env, DACON_API_KEY: ' bad key ' }); assert.equal((await r.json()).configured, false);
});

test('catalog exposes only sorted minimal receipt identity and no SMILES', async () => {
  const worker = createWorker({ receipt });
  const response = await worker.fetch(new Request('https://worker.example/api/catalog'), env);
  assert.equal(response.status, 200);
  const data = await response.json();
  assert.deepEqual(Object.keys(data).sort(), ['candidates', 'parents', 'scientific_approval', 'status']);
  assert.equal(data.status, 'ok'); assert.equal(data.scientific_approval, false);
  assert.equal(data.parents.length, 9); assert.equal(data.candidates.length, 324);
  for (const parent of data.parents) assert.deepEqual(Object.keys(parent), ['parent_id']);
  for (const candidate of data.candidates) assert.deepEqual(Object.keys(candidate).sort(), ['candidate_id', 'e3_type', 'linker_id', 'parent_id']);
  assert.equal(JSON.stringify(data).toLowerCase().includes('smiles'), false);
  assert.deepEqual(data.parents.map(x => x.parent_id), [...data.parents.map(x => x.parent_id)].sort());
  assert.equal((await worker.fetch(new Request('https://worker.example/api/catalog?x=1'), env)).status, 400);
  assert.equal((await worker.fetch(new Request('https://worker.example/api/catalog', { method: 'POST' }), env)).status, 405);
});

test('non-api GET and HEAD delegate to ASSETS while unknown API never does', async () => {
  const seen = [];
  const assets = { fetch: async request => { seen.push([request.method, new URL(request.url).pathname]); return new Response(request.method === 'HEAD' ? null : 'asset', { status: 200 }); } };
  const worker = createWorker({ receipt });
  assert.equal((await worker.fetch(new Request('https://worker.example/index.html'), { ...env, ASSETS: assets })).status, 200);
  assert.equal((await worker.fetch(new Request('https://worker.example/app.js', { method: 'HEAD' }), { ...env, ASSETS: assets })).status, 200);
  assert.equal((await worker.fetch(new Request('https://worker.example/api/unknown'), { ...env, ASSETS: assets })).status, 404);
  assert.deepEqual(seen, [['GET', '/index.html'], ['HEAD', '/app.js']]);
  assert.equal((await worker.fetch(new Request('https://worker.example/no-assets'), env)).status, 404);
});

test('incoming cancellation returns 499, aborts upstream, and releases concurrency', async () => {
  let calls = 0; let observedAbort = false;
  const mock = async (url, init) => {
    calls++;
    if (calls > 1) return successMock([])(url, init);
    return new Promise((resolve, reject) => init.signal.addEventListener('abort', () => { observedAbort = true; reject(new DOMException('Aborted', 'AbortError')); }, { once: true }));
  };
  const worker = createWorker({ receipt, fetchImpl: mock, timeoutMs: 1000 });
  const already = new AbortController(); already.abort();
  let response = await worker.fetch(new Request('https://worker.example/api/run', { method: 'POST', headers, body: JSON.stringify({ operation: 'candidate_compare' }), signal: already.signal }), env);
  let data = await response.json(); assert.equal(response.status, 499); assert.equal(data.error.code, 'REQUEST_ABORTED'); assert.equal(data.live_model_called, false); assert.equal(calls, 0);
  const controller = new AbortController();
  const pending = worker.fetch(new Request('https://worker.example/api/run', { method: 'POST', headers: { ...headers, 'cf-connecting-ip': 'abort-live' }, body: JSON.stringify({ operation: 'candidate_compare' }), signal: controller.signal }), env);
  while (calls < 1) await new Promise(resolve => setTimeout(resolve, 1));
  controller.abort();
  response = await pending; data = await response.json();
  assert.equal(response.status, 499); assert.equal(data.error.code, 'REQUEST_ABORTED'); assert.equal(data.live_model_called, true); assert.equal(data.usage, undefined); assert.equal(data.answer, undefined); assert.equal(observedAbort, true);
  response = await run(worker, { operation: 'candidate_compare' }, { headers: { 'cf-connecting-ip': 'after-abort' } }); assert.equal(response.status, 200);
});

test('strict Responses parsing concatenates text and rejects unknown output content', async () => {
  const split = createWorker({ receipt, fetchImpl: async (url, init) => {
    const context = JSON.parse(JSON.parse(init.body).input);
    const text = answerFor(context.allowed_evidence_refs);
    return new Response(JSON.stringify({ id: 'resp_split', model: 'gpt-5.6-sol', status: 'completed', usage: { input_tokens: 1, output_tokens: 1, total_tokens: 2 }, output: [{ type: 'reasoning' }, { type: 'message', role: 'assistant', content: [{ type: 'output_text', text: text.slice(0, 10) }, { type: 'output_text', text: text.slice(10) }] }] }));
  } });
  assert.equal((await run(split, { operation: 'candidate_compare' })).status, 200);
  const invalid = createWorker({ receipt, fetchImpl: async () => new Response(JSON.stringify({ id: 'resp_bad', model: 'gpt-5.6-sol', status: 'completed', usage: { input_tokens: 1, output_tokens: 1, total_tokens: 2 }, output: [{ type: 'tool_call' }] })) });
  const data = await (await run(invalid, { operation: 'candidate_compare' })).json();
  assert.equal(data.error.code, 'API_INVALID_ANSWER');
});

test('default fetch implementation preserves the global fetch receiver', { concurrency: false }, async () => {
  const originalFetch = globalThis.fetch;
  let capturedCalls = 0;
  try {
    globalThis.fetch = async function (url, init) {
      assert.equal(this, globalThis);
      capturedCalls++;
      const context = JSON.parse(JSON.parse(init.body).input);
      return new Response(JSON.stringify({ id: 'resp_receiver', model: 'gpt-5.6-sol', status: 'completed', usage: { input_tokens: 10, output_tokens: 5, total_tokens: 15 }, output: [{ type: 'message', role: 'assistant', content: [{ type: 'output_text', text: answerFor(context.allowed_evidence_refs) }] }] }), { status: 200 });
    };
    const worker = createWorker({ receipt });
    const response = await run(worker, { operation: 'candidate_compare' });
    const data = await response.json();
    assert.equal(response.status, 200);
    assert.equal(data.status, 'completed');
    assert.equal(data.live_model_called, true);
    assert.equal(capturedCalls, 1);
    assert.equal(data.usage.total_tokens, 15);
  } finally {
    globalThis.fetch = originalFetch;
  }
});

test('transport diagnostics log only fixed safe fetch metadata', { concurrency: false }, async () => {
  const originalConsoleError = console.error;
  const logs = [];
  const secret = 'transport-diagnostic-secret';
  try {
    console.error = entry => { logs.push(entry); };
    const worker = createWorker({
      receipt,
      fetchImpl: async () => { throw new TypeError('Illegal invocation ' + secret + ' https://private.example/path'); }
    });
    const response = await run(worker, { operation: 'candidate_compare' }, { headers: { 'cf-connecting-ip': 'transport-diagnostic' } });
    const data = await response.json();
    assert.equal(response.status, 502);
    assert.equal(data.error.code, 'UPSTREAM_TRANSPORT_ERROR');
    assert.deepEqual(logs, [{
      event: 'tpd_upstream_transport_error',
      request_id: data.request_id,
      phase: 'fetch',
      error_type: 'TypeError',
      category: 'ILLEGAL_INVOCATION'
    }]);
    const serialized = JSON.stringify(logs);
    assert.equal(serialized.includes(secret), false);
    assert.equal(serialized.includes('private.example'), false);
    assert.equal(serialized.includes('message'), false);
    assert.equal(serialized.includes('stack'), false);
  } finally {
    console.error = originalConsoleError;
  }
});

test('transport diagnostics classify redirect, egress, and hashed other errors safely', { concurrency: false }, async () => {
  const originalConsoleError = console.error;
  const logs = [];
  const secret = 'local-fake-transport-key';
  const cases = [
    {
      category: 'REDIRECT_UNSUPPORTED',
      message: 'redirect is not supported raw-redirect https://diagnostic.example/redirect'
    },
    {
      category: 'EGRESS_POLICY',
      message: 'outbound network access is blocked raw-egress https://diagnostic.example/egress'
    },
    {
      category: 'OTHER',
      message: `opaque raw-other ${secret} https://diagnostic.example/private stack-trace`
    }
  ];
  try {
    console.error = entry => { logs.push(entry); };
    for (const diagnostic of cases) {
      logs.length = 0;
      const worker = createWorker({
        receipt,
        fetchImpl: async () => { throw new TypeError(diagnostic.message); }
      });
      const response = await run(worker, { operation: 'candidate_compare' }, {
        headers: { 'cf-connecting-ip': `transport-${diagnostic.category.toLowerCase()}` },
        env: { ...env, DACON_API_KEY: secret }
      });
      const data = await response.json();
      assert.equal(response.status, 502);
      assert.equal(data.error.code, 'UPSTREAM_TRANSPORT_ERROR');
      assert.equal(data.answer, undefined);
      assert.equal(logs.length, 1);
      const entry = logs[0];
      assert.equal(entry.event, 'tpd_upstream_transport_error');
      assert.equal(entry.request_id, data.request_id);
      assert.equal(entry.phase, 'fetch');
      assert.equal(entry.error_type, 'TypeError');
      assert.equal(entry.category, diagnostic.category);
      const expectedKeys = diagnostic.category === 'OTHER'
        ? ['category', 'error_type', 'event', 'message_sha256', 'phase', 'request_id']
        : ['category', 'error_type', 'event', 'phase', 'request_id'];
      assert.deepEqual(Object.keys(entry).sort(), expectedKeys);
      if (diagnostic.category === 'OTHER') {
        const redacted = diagnostic.message.slice(0, 4096).split(secret).join('[REDACTED]');
        const digest = globalThis.crypto.subtle.digest.bind(globalThis.crypto.subtle);
        const expectedBuffer = await digest('SHA-256', new TextEncoder().encode(redacted));
        const expectedHash = Array.from(new Uint8Array(expectedBuffer), byte => byte.toString(16).padStart(2, '0')).join('');
        assert.match(entry.message_sha256, /^[0-9a-f]{64}$/);
        assert.equal(entry.message_sha256, expectedHash);
      }
      const serialized = JSON.stringify(logs);
      assert.equal(serialized.includes(secret), false);
      assert.equal(serialized.includes('diagnostic.example'), false);
      assert.equal(serialized.includes('raw-'), false);
      assert.equal(serialized.includes('stack-trace'), false);
      assert.equal(serialized.includes(diagnostic.message), false);
      assert.equal(serialized.includes('"stack"'), false);
    }
  } finally {
    console.error = originalConsoleError;
  }
});

test('upstream redirects are rejected without following or exposing redirect data', async () => {
  const calls = [];
  const rawLocation = 'https://untrusted.example/redirect';
  const localTestSecret = env.DACON_API_KEY;
  const worker = createWorker({ receipt, fetchImpl: async (url, init) => {
    calls.push({ url, init });
    return new Response(JSON.stringify({ diagnostic: `raw-${localTestSecret}`, location: rawLocation }), {
      status: 302,
      headers: { Location: rawLocation }
    });
  } });

  const response = await run(worker, { operation: 'candidate_compare' });
  const body = await response.json();
  const serialized = JSON.stringify(body);

  assert.equal(calls.length, 1);
  assert.equal(calls[0].url, 'https://dacon-apim-hackathon-0903.azure-api.net/hackathon/openai/v1/responses');
  assert.equal(calls[0].init.redirect, 'manual');
  assert.equal(calls.some(call => String(call.url).includes('untrusted.example')), false);
  assert.equal(response.status, 502);
  assert.equal(body.error.code, 'UPSTREAM_REDIRECT');
  assert.equal(body.live_model_called, true);
  assert.equal(body.answer, undefined);
  assert.equal(serialized.includes(rawLocation), false);
  assert.equal(serialized.includes('untrusted.example'), false);
  assert.equal(serialized.includes(localTestSecret), false);
});
