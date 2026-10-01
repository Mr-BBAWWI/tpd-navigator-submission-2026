import http from 'node:http';
import path from 'node:path';
import { lstat, readFile } from 'node:fs/promises';
import { fileURLToPath, pathToFileURL } from 'node:url';

const MAX_REQUEST_BODY = 4096;
const MODULE_DIRECTORY = path.dirname(fileURLToPath(import.meta.url));
const ALLOWED_ENV = ['DACON_API_KEY', 'ALLOWED_ORIGINS'];
const MIME_TYPES = new Map([
  ['.html', 'text/html; charset=utf-8'], ['.js', 'text/javascript; charset=utf-8'],
  ['.css', 'text/css; charset=utf-8'], ['.json', 'application/json; charset=utf-8'],
  ['.svg', 'image/svg+xml; charset=utf-8']
]);

function inside(root, candidate) {
  const relative = path.relative(root, candidate);
  return relative === '' || (!relative.startsWith(`..${path.sep}`) && relative !== '..' && !path.isAbsolute(relative));
}

function errorCode(error) {
  const text = String(error?.code ?? error?.name ?? 'WORKER_ERROR').toUpperCase();
  return /^[A-Z0-9_]{1,40}$/.test(text) ? text : 'WORKER_ERROR';
}

async function safeStaticPath(root, pathname) {
  let segments;
  try { segments = pathname.split('/').filter(Boolean).map(decodeURIComponent); }
  catch { return { error: 400 }; }
  if (segments.some((part) => part === '.' || part === '..' || part.includes('/') || part.includes('\\') || part.includes('\0'))) return { error: 400 };
  if (!segments.length) segments = ['index.html'];
  const candidate = path.resolve(root, ...segments);
  if (!inside(root, candidate)) return { error: 400 };
  let current = root;
  for (const segment of segments) {
    current = path.join(current, segment);
    const info = await lstat(current).catch(() => null);
    if (!info) return { error: 404 };
    if (info.isSymbolicLink()) return { error: 403 };
  }
  return (await lstat(candidate)).isFile() ? { file: candidate } : { error: 404 };
}

function createAssets(clientDirectory) {
  return { async fetch(request) {
    if (!['GET', 'HEAD'].includes(request.method)) return new Response('Method Not Allowed', { status: 405, headers: { Allow: 'GET, HEAD', 'X-Content-Type-Options': 'nosniff' } });
    const resolved = await safeStaticPath(clientDirectory, new URL(request.url).pathname);
    if (resolved.error) return new Response(resolved.error === 404 ? 'Not Found' : 'Bad Request', { status: resolved.error, headers: { 'X-Content-Type-Options': 'nosniff' } });
    const data = await readFile(resolved.file);
    const headers = { 'Content-Type': MIME_TYPES.get(path.extname(resolved.file).toLowerCase()) ?? 'application/octet-stream', 'X-Content-Type-Options': 'nosniff', 'Content-Length': String(data.length) };
    return new Response(request.method === 'HEAD' ? null : data, { headers });
  } };
}

async function readRequestBody(incoming, signal) {
  const declared = incoming.headers['content-length'];
  if (declared !== undefined) {
    const length = Number(declared);
    if (!Number.isSafeInteger(length) || length < 0) return { status: 400 };
    if (length > MAX_REQUEST_BODY) return { status: 413 };
  }
  const chunks = [];
  let total = 0;
  try {
    for await (const chunk of incoming) {
      if (signal.aborted) return { aborted: true };
      total += chunk.length;
      if (total > MAX_REQUEST_BODY) return { status: 413 };
      chunks.push(chunk);
    }
  } catch { return { aborted: true }; }
  return { body: Buffer.concat(chunks, total) };
}

function waitForDrain(outgoing, signal) {
  return new Promise((resolve, reject) => {
    const cleanup = () => {
      outgoing.off('drain', onDrain);
      outgoing.off('close', onClose);
      outgoing.off('error', onError);
      signal.removeEventListener('abort', onAbort);
    };
    const finish = (error) => { cleanup(); error ? reject(error) : resolve(); };
    const onDrain = () => finish();
    const onClose = () => finish(Object.assign(new Error('Response closed'), { code: 'CLIENT_CLOSED' }));
    const onError = (error) => finish(error);
    const onAbort = () => finish(Object.assign(new Error('Response aborted'), { code: 'ABORT_ERR' }));
    outgoing.once('drain', onDrain);
    outgoing.once('close', onClose);
    outgoing.once('error', onError);
    signal.addEventListener('abort', onAbort, { once: true });
    if (signal.aborted) onAbort();
    else if (outgoing.destroyed || outgoing.closed) onClose();
  });
}

async function writeWorkerResponse(outgoing, response, controller) {
  outgoing.statusCode = response.status;
  outgoing.statusMessage = response.statusText || outgoing.statusMessage;
  for (const [name, value] of response.headers) {
    try { outgoing.setHeader(name, value); } catch {}
  }
  if (!response.body) { outgoing.end(); return; }

  const reader = response.body.getReader();
  let complete = false;
  const cancel = () => reader.cancel(controller.signal.reason).catch(() => {});
  const onClose = () => {
    if (!outgoing.writableFinished) {
      controller.abort();
      void cancel();
    }
  };
  outgoing.once('close', onClose);
  try {
    while (!controller.signal.aborted) {
      const { done, value } = await reader.read();
      if (done) { complete = true; break; }
      if (!outgoing.write(Buffer.from(value))) await waitForDrain(outgoing, controller.signal);
    }
    if (!controller.signal.aborted && !outgoing.destroyed) outgoing.end();
  } finally {
    outgoing.off('close', onClose);
    if (!complete) await cancel();
    reader.releaseLock();
  }
}

function filteredEnvironment(source) {
  const result = {};
  for (const name of ALLOWED_ENV) if (typeof source?.[name] === 'string') result[name] = source[name];
  return result;
}

export async function createLocalServer({ distDirectory, env = process.env, port = 8788, quiet = false } = {}) {
  const dist = path.resolve(distDirectory ?? path.join(MODULE_DIRECTORY, 'dist'));
  const distInfo = await lstat(dist).catch(() => null);
  if (!distInfo?.isDirectory() || distInfo.isSymbolicLink()) throw new Error('Built dist directory is missing or unsafe; run npm run build first');
  const clientDirectory = path.join(dist, 'client');
  const workerFile = path.join(dist, 'server', 'index.js');
  const clientInfo = await lstat(clientDirectory).catch(() => null);
  const workerInfo = await lstat(workerFile).catch(() => null);
  if (!clientInfo?.isDirectory() || clientInfo.isSymbolicLink() || !workerInfo?.isFile() || workerInfo.isSymbolicLink()) throw new Error('Built client or Worker artifact is missing or unsafe');

  const worker = (await import(`${pathToFileURL(workerFile).href}?local=${Date.now()}-${Math.random()}`)).default;
  if (!worker || typeof worker.fetch !== 'function') throw new Error('Built Worker must default-export an object with fetch()');
  const runtimeEnv = { ...filteredEnvironment(env), ASSETS: createAssets(clientDirectory) };
  const numericPort = Number(port);
  if (!Number.isInteger(numericPort) || numericPort < 0 || numericPort > 65535) throw new Error('Port must be an integer from 0 through 65535');

  const server = http.createServer(async (incoming, outgoing) => {
    const controller = new AbortController();
    const onAborted = () => controller.abort();
    const onClose = () => { if (!outgoing.writableFinished) controller.abort(); };
    incoming.once('aborted', onAborted);
    outgoing.once('close', onClose);
    try {
      const address = server.address();
      const activePort = typeof address === 'object' && address ? address.port : numericPort;
      const expectedOrigin = `http://127.0.0.1:${activePort}`;
      let requestUrl;
      try {
        if (!incoming.url?.startsWith('/') || incoming.url.startsWith('//')) throw new Error('invalid target');
        requestUrl = new URL(incoming.url, expectedOrigin);
        if (requestUrl.origin !== expectedOrigin) throw new Error('origin changed');
      } catch {
        outgoing.writeHead(400, { 'Content-Type': 'text/plain; charset=utf-8', 'X-Content-Type-Options': 'nosniff' }).end('Bad Request');
        return;
      }

      const method = incoming.method ?? 'GET';
      const hasBody = method !== 'GET' && method !== 'HEAD';
      if (!hasBody && Number(incoming.headers['content-length'] ?? 0) > 0) {
        outgoing.writeHead(400, { 'Content-Type': 'text/plain; charset=utf-8' }).end('Bad Request');
        return;
      }
      let body;
      if (hasBody) {
        const read = await readRequestBody(incoming, controller.signal);
        if (read.status) {
          outgoing.writeHead(read.status, { 'Content-Type': 'text/plain; charset=utf-8', 'X-Content-Type-Options': 'nosniff' }).end(read.status === 413 ? 'Payload Too Large' : 'Bad Request');
          return;
        }
        if (read.aborted) return;
        body = read.body.length ? read.body : undefined;
      }

      const headers = new Headers();
      for (const [name, value] of Object.entries(incoming.headers)) {
        if (value === undefined || name === 'host') continue;
        if (Array.isArray(value)) value.forEach((item) => headers.append(name, item));
        else headers.set(name, value);
      }
      headers.set('host', `127.0.0.1:${activePort}`);
      const request = new Request(requestUrl, { method, headers, body, signal: controller.signal });
      const pending = [];
      const ctx = {
        waitUntil(promise) { pending.push(Promise.resolve(promise).catch((error) => { if (!quiet) console.error(`waitUntil error: ${errorCode(error)}`); })); },
        passThroughOnException() {}
      };
      const response = await worker.fetch(request, runtimeEnv, ctx);
      if (!(response instanceof Response)) throw Object.assign(new Error('Worker returned an invalid response'), { code: 'INVALID_RESPONSE' });
      await writeWorkerResponse(outgoing, response, controller);
      void Promise.allSettled(pending);
    } catch (error) {
      if (!quiet && !controller.signal.aborted) console.error(`Worker error: ${errorCode(error)}`);
      if (!outgoing.headersSent && !outgoing.destroyed) outgoing.writeHead(500, { 'Content-Type': 'text/plain; charset=utf-8', 'X-Content-Type-Options': 'nosniff' }).end('Internal Server Error');
      else if (!outgoing.writableEnded && !outgoing.destroyed) outgoing.destroy();
    } finally {
      incoming.off('aborted', onAborted);
      outgoing.off('close', onClose);
    }
  });

  await new Promise((resolve, reject) => {
    server.once('error', reject);
    server.listen(numericPort, '127.0.0.1', () => { server.off('error', reject); resolve(); });
  });
  if (!quiet) {
    console.log(`Local server: http://127.0.0.1:${server.address().port}`);
    console.log(`Runtime credentials configured: ${Boolean(runtimeEnv.DACON_API_KEY)}`);
  }
  return server;
}

function parseArguments(argv) {
  let port = 8788;
  for (let index = 0; index < argv.length; index += 1) {
    const argument = argv[index];
    if (argument === '--help' || argument === '-h') return { help: true };
    if (argument === '--port') port = Number(argv[++index]);
    else if (argument.startsWith('--port=')) port = Number(argument.slice(7));
    else throw new Error(`Unknown argument: ${argument}`);
  }
  if (!Number.isInteger(port) || port < 1 || port > 65535) throw new Error('Port must be an integer from 1 through 65535');
  return { port };
}

const invokedDirectly = process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url);
if (invokedDirectly) {
  try {
    const options = parseArguments(process.argv.slice(2));
    if (options.help) process.stdout.write('Usage: node local-server.mjs [--port PORT]\n');
    else await createLocalServer({ port: options.port });
  } catch (error) {
    process.stderr.write(`Local server failed: ${error.message}\n`);
    process.exitCode = 1;
  }
}
