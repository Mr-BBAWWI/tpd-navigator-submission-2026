import { createHash } from 'node:crypto';
import { createRequire } from 'node:module';
import { copyFile, lstat, mkdir, mkdtemp, readFile, readdir, rename, rm, stat, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const MAX_PUBLIC_BYTES = 40 * 1024 * 1024;
const MODULE_DIRECTORY = path.dirname(fileURLToPath(import.meta.url));
const require = createRequire(import.meta.url);

function inside(root, candidate) {
  const relative = path.relative(root, candidate);
  return relative === '' || (!relative.startsWith(`..${path.sep}`) && relative !== '..' && !path.isAbsolute(relative));
}

async function regularFile(file, root, label) {
  const absolute = path.resolve(file);
  if (!inside(root, absolute)) throw new Error(`${label} escapes the project directory`);
  let current = root;
  for (const part of path.relative(root, absolute).split(path.sep).filter(Boolean)) {
    current = path.join(current, part);
    const info = await lstat(current).catch(() => null);
    if (!info) throw new Error(`Missing required input: ${label}`);
    if (info.isSymbolicLink()) throw new Error(`Symbolic links are not allowed: ${label}`);
  }
  if (!(await lstat(absolute)).isFile()) throw new Error(`Expected a regular file: ${label}`);
}

async function safeDirectory(directory, root, label) {
  const absolute = path.resolve(directory);
  if (!inside(root, absolute)) throw new Error(`${label} escapes the project directory`);
  let current = root;
  for (const part of path.relative(root, absolute).split(path.sep).filter(Boolean)) {
    current = path.join(current, part);
    const info = await lstat(current).catch(() => null);
    if (!info) throw new Error(`Missing directory: ${label}`);
    if (info.isSymbolicLink()) throw new Error(`Symbolic links are not allowed: ${label}`);
    if (!info.isDirectory()) throw new Error(`Expected a directory: ${label}`);
  }
}

const sha256 = (data) => createHash('sha256').update(data).digest('hex');

function scan(data, label) {
  const configured = process.env.DACON_API_KEY;
  if (configured && data.includes(Buffer.from(configured))) throw new Error(`Runtime credential material was found in ${label}`);
  const sample = data.subarray(0, 8192);
  if (sample.includes(0) || [...sample].filter((byte) => byte < 9 || (byte > 13 && byte < 32)).length > sample.length * 0.02) return;
  const text = data.toString('utf8');
  const checks = [
    [/-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----/, 'private key material'],
    [/(?:^|[^A-Za-z0-9])(?:sk-[A-Za-z0-9_-]{20,}|gh[pousr]_[A-Za-z0-9]{30,})/, 'common access token'],
    [/\bBearer\s+[A-Za-z0-9._~-]{24,}/i, 'bearer token'],
    [/["']?DACON_API_KEY["']?\s*[:=]\s*["'][^"']+["']/, 'embedded DACON credential'],
    [/(?:^|[\s"'(])(?:\/Users\/[^/\s]+|\/home\/[^/\s]+|[A-Za-z]:\\Users\\[^\\\s]+)/, 'private local path']
  ];
  for (const [pattern, description] of checks) if (pattern.test(text)) throw new Error(`${description} was found in ${label}`);
}

async function loadEsbuild(explicitPath) {
  try {
    let value;
    if (!explicitPath) value = await import('esbuild');
    else {
      const absolute = path.resolve(explicitPath);
      const info = await stat(absolute);
      value = info.isDirectory() ? require(absolute) : await import(pathToFileURL(absolute).href);
    }
    const api = value.default?.build ? value.default : value;
    if (typeof api.build !== 'function') throw new Error('module does not export build()');
    if (api.version !== '0.25.11') throw new Error(`expected esbuild 0.25.11, found ${api.version ?? 'unknown'}`);
    return api;
  } catch (error) {
    throw new Error(`esbuild 0.25.11 is required; run npm install first (${error.message})`);
  }
}

async function collectPublic(sourceRoot, destinationRoot) {
  let total = 0;
  const records = [];
  async function visit(sourceDirectory, relativeDirectory = '') {
    await safeDirectory(sourceDirectory, sourceRoot, `public/${relativeDirectory}`);
    const entries = await readdir(sourceDirectory, { withFileTypes: true });
    entries.sort((a, b) => a.name.localeCompare(b.name));
    for (const entry of entries) {
      const relative = path.join(relativeDirectory, entry.name);
      const source = path.join(sourceDirectory, entry.name);
      const destination = path.join(destinationRoot, relative);
      const info = await lstat(source);
      if (!inside(destinationRoot, destination)) throw new Error(`Public path escapes destination: ${relative}`);
      if (info.isSymbolicLink()) throw new Error(`Symbolic links are not allowed in public: ${relative}`);
      if (info.isDirectory()) {
        await mkdir(destination);
        await visit(source, relative);
      } else if (info.isFile()) {
        total += info.size;
        if (total > MAX_PUBLIC_BYTES) throw new Error('Public files exceed the 40 MB build limit');
        const data = await readFile(source);
        const name = `public/${relative.split(path.sep).join('/')}`;
        scan(data, name);
        await copyFile(source, destination);
        records.push({ relative: name, data });
      } else throw new Error(`Only regular public files are allowed: ${relative}`);
    }
  }
  await visit(sourceRoot);
  return records;
}

async function removeOwnedTree(target, parent, prefixes) {
  const absolute = path.resolve(target);
  if (path.dirname(absolute) !== path.resolve(parent) || !prefixes.some((prefix) => path.basename(absolute).startsWith(prefix))) {
    throw new Error('Refusing to remove an unowned directory');
  }
  const rootInfo = await lstat(absolute).catch(() => null);
  if (!rootInfo) return;
  if (rootInfo.isSymbolicLink() || !rootInfo.isDirectory()) throw new Error('Refusing to remove an unsafe build path');
  async function inspect(directory) {
    for (const entry of await readdir(directory, { withFileTypes: true })) {
      const child = path.join(directory, entry.name);
      const info = await lstat(child);
      if (info.isSymbolicLink()) throw new Error('Refusing to remove a tree containing a symbolic link');
      if (info.isDirectory()) await inspect(child);
      else if (!info.isFile()) throw new Error('Refusing to remove a tree containing a special file');
    }
  }
  await inspect(absolute);
  await rm(absolute, { recursive: true });
}

export async function buildProject({ projectDir = MODULE_DIRECTORY, esbuildModulePath } = {}) {
  const project = path.resolve(projectDir);
  const projectInfo = await lstat(project).catch(() => null);
  if (!projectInfo?.isDirectory() || projectInfo.isSymbolicLink()) throw new Error('Project directory must be a regular directory');

  const worker = path.join(project, 'worker.mjs');
  const publicDirectory = path.join(project, 'public');
  const hostingSource = path.join(project, '.openai', 'hosting.json');
  await regularFile(worker, project, 'worker.mjs');
  await safeDirectory(publicDirectory, project, 'public');
  await regularFile(path.join(publicDirectory, 'index.html'), project, 'public/index.html');
  await regularFile(hostingSource, project, '.openai/hosting.json');

  const dist = path.join(project, 'dist');
  const existingDist = await lstat(dist).catch(() => null);
  if (existingDist && (existingDist.isSymbolicLink() || !existingDist.isDirectory())) throw new Error('Existing dist must be a regular directory');
  const esbuild = await loadEsbuild(esbuildModulePath ?? process.env.ESBUILD_MODULE_PATH);
  const staging = await mkdtemp(path.join(project, '.dist-build-'));
  const backup = path.join(project, `.dist-backup-${process.pid}-${Date.now()}`);
  let oldMoved = false;
  let installed = false;

  try {
    const client = path.join(staging, 'client');
    const server = path.join(staging, 'server');
    const openai = path.join(staging, '.openai');
    await Promise.all([mkdir(client), mkdir(server), mkdir(openai)]);
    const sources = await collectPublic(publicDirectory, client);

    const hostingData = await readFile(hostingSource);
    scan(hostingData, '.openai/hosting.json');
    let hosting;
    try { hosting = JSON.parse(hostingData); } catch { throw new Error('.openai/hosting.json must contain valid JSON'); }
    if (!hosting || typeof hosting !== 'object' || Array.isArray(hosting) || typeof hosting.project_id !== 'string' || !hosting.project_id) {
      throw new Error('.openai/hosting.json must contain a project_id');
    }
    const hostingOutput = Buffer.from(`${JSON.stringify({ d1: null, r2: null }, null, 2)}\n`);
    scan(hostingOutput, 'dist/.openai/hosting.json');
    await writeFile(path.join(openai, 'hosting.json'), hostingOutput, { flag: 'wx' });
    sources.push({ relative: '.openai/hosting.json', data: hostingData });

    let result;
    try {
      result = await esbuild.build({ absWorkingDir: project, entryPoints: ['worker.mjs'], bundle: true, platform: 'browser', format: 'esm', target: 'es2022', outfile: path.join(server, 'index.js'), metafile: true, write: true, sourcemap: false, loader: { '.json': 'json' }, logLevel: 'silent' });
    } catch (error) {
      throw new Error(`Worker bundle failed: ${error.errors?.[0]?.text ?? error.message}`);
    }

    for (const inputName of Object.keys(result.metafile.inputs)) {
      for (const part of inputName.split(/[\\/]/)) {
        const lower = part.toLowerCase();
        if (lower === '.env' || lower.startsWith('.env.') || lower.includes('credential')) throw new Error(`Unsafe build input: ${inputName}`);
      }
      const absolute = path.resolve(project, inputName);
      await regularFile(absolute, project, inputName);
      const data = await readFile(absolute);
      scan(data, inputName);
      const relative = inputName.split(path.sep).join('/');
      if (!sources.some((record) => record.relative === relative)) sources.push({ relative, data });
    }
    for (const output of Object.values(result.metafile.outputs)) {
      for (const imported of output.imports ?? []) if (imported.external) throw new Error(`External import is not allowed: ${imported.path}`);
    }

    const workerOutput = await readFile(path.join(server, 'index.js'));
    scan(workerOutput, 'dist/server/index.js');
    if (/\b(?:from\s*|import\s*\()\s*["'](?:node:|fs(?:\/|["'])|path["']|crypto["']|http["']|https["'])/.test(workerOutput.toString())) {
      throw new Error('Node built-in imports are not allowed in the Worker bundle');
    }

    const outputs = [];
    async function collect(directory, relative = '') {
      const entries = await readdir(directory, { withFileTypes: true });
      entries.sort((a, b) => a.name.localeCompare(b.name));
      for (const entry of entries) {
        const name = path.posix.join(relative, entry.name);
        const child = path.join(directory, entry.name);
        const info = await lstat(child);
        if (info.isSymbolicLink()) throw new Error(`Generated output contains a symbolic link: ${name}`);
        if (info.isDirectory()) await collect(child, name);
        else if (info.isFile()) {
          const data = await readFile(child);
          scan(data, `dist/${name}`);
          outputs.push({ path: name, sha256: sha256(data), bytes: data.length });
        } else throw new Error(`Generated output is not a regular file: ${name}`);
      }
    }
    await collect(staging);
    sources.sort((a, b) => a.relative.localeCompare(b.relative));
    outputs.sort((a, b) => a.path.localeCompare(b.path));
    const receipt = { version: 1, sources: sources.map(({ relative, data }) => ({ path: relative, sha256: sha256(data), bytes: data.length })), outputs };
    await writeFile(path.join(openai, 'build-receipt.json'), `${JSON.stringify(receipt, null, 2)}\n`, { flag: 'wx' });

    if (existingDist) { await rename(dist, backup); oldMoved = true; }
    try { await rename(staging, dist); installed = true; }
    catch (error) {
      if (oldMoved) {
        try { await rename(backup, dist); oldMoved = false; }
        catch { throw new Error(`Could not install the new build; the previous build remains at ${path.basename(backup)}`); }
      }
      throw error;
    }
    if (oldMoved) { await removeOwnedTree(backup, project, ['.dist-backup-']); oldMoved = false; }
    return { ok: true, outputs: outputs.length + 1 };
  } finally {
    if (!installed) await removeOwnedTree(staging, project, ['.dist-build-']).catch(() => {});
  }
}

async function writeFixture(project) {
  await mkdir(path.join(project, 'public'), { recursive: true });
  await mkdir(path.join(project, '.openai'), { recursive: true });
  await writeFile(path.join(project, 'package.json'), '{"private":true,"type":"module"}\n');
  await writeFile(path.join(project, 'campaign-public.json'), '{"name":"fixture-campaign"}\n');
  await writeFile(path.join(project, 'public', 'index.html'), '<!doctype html><title>Fixture</title><h1>Fixture site</h1>\n');
  await writeFile(path.join(project, 'public', 'app.js'), 'globalThis.fixture = true;\n');
  await writeFile(path.join(project, '.openai', 'hosting.json'), JSON.stringify({ project_id: 'appgprj_fixture_public_identifier', custom: 'preserved' }));
  await writeFile(path.join(project, 'worker.mjs'), `import campaign from './campaign-public.json';\nexport default { async fetch(request, env) { return new URL(request.url).pathname === '/api/demo' ? Response.json(campaign) : env.ASSETS.fetch(request); } };\n`);
}

async function expectReject(action, pattern, label) {
  try { await action(); }
  catch (error) {
    if (pattern.test(error.message)) return;
    throw new Error(`${label} rejected for the wrong reason: ${error.message}`);
  }
  throw new Error(`${label} was not rejected`);
}

async function runSelfTest() {
  const temporaryRoot = path.resolve(tmpdir());
  const parent = path.resolve(await mkdtemp(path.join(temporaryRoot, 'live-demo-build-test-')));
  if (path.dirname(parent) !== temporaryRoot || !path.basename(parent).startsWith('live-demo-build-test-')) throw new Error('Unsafe self-test temporary directory');
  const project = path.join(parent, 'fixture');
  await mkdir(project);
  try {
    await writeFixture(project);
    const missingIndex = path.join(parent, 'missing-index');
    await mkdir(missingIndex);
    await writeFixture(missingIndex);
    await rm(path.join(missingIndex, 'public', 'index.html'));
    await expectReject(() => buildProject({ projectDir: missingIndex }), /^Missing required input: public\/index\.html$/, 'Missing-index self-test');

    let symlink = 'unsupported';
    const link = path.join(project, 'public', 'linked.txt');
    try {
      const fs = await import('node:fs/promises');
      await fs.symlink(path.join(project, 'public', 'index.html'), link);
      await expectReject(() => buildProject({ projectDir: project }), /^Symbolic links are not allowed in public:/, 'Symlink self-test');
      symlink = 'passed';
      await rm(link, { force: true });
    } catch (error) {
      await rm(link, { force: true });
      if (!['EPERM', 'EACCES', 'ENOSYS'].includes(error.code)) throw error;
    }

    const workerFile = path.join(project, 'worker.mjs');
    const originalWorker = await readFile(workerFile);
    const previousSecret = process.env.DACON_API_KEY;
    const secret = 'x';
    try {
      process.env.DACON_API_KEY = secret;
      await writeFile(workerFile, `export default { fetch(){ return new Response('${secret}') } };\n`);
      await expectReject(() => buildProject({ projectDir: project }), /^Runtime credential material was found in /, 'Credential self-test');
    } finally {
      await writeFile(workerFile, originalWorker);
      if (previousSecret === undefined) delete process.env.DACON_API_KEY;
      else process.env.DACON_API_KEY = previousSecret;
    }

    const sourceHostingPath = path.join(project, '.openai', 'hosting.json');
    const sourceHostingBefore = await readFile(sourceHostingPath);
    await buildProject({ projectDir: project });
    const sourceHostingAfter = await readFile(sourceHostingPath);
    if (!sourceHostingAfter.equals(sourceHostingBefore)) throw new Error('Source hosting metadata self-test failed');
    const sourceHosting = JSON.parse(sourceHostingAfter);
    if (sourceHosting.project_id !== 'appgprj_fixture_public_identifier' || sourceHosting.custom !== 'preserved') throw new Error('Source hosting identity self-test failed');
    const hosting = JSON.parse(await readFile(path.join(project, 'dist', '.openai', 'hosting.json')));
    const hostingKeys = Object.keys(hosting).sort();
    if (hostingKeys.length !== 2 || hostingKeys[0] !== 'd1' || hostingKeys[1] !== 'r2' || hosting.d1 !== null || hosting.r2 !== null) throw new Error('Hosting metadata self-test failed');
    const module = await import(`${pathToFileURL(path.join(project, 'dist', 'server', 'index.js')).href}?selftest=${Date.now()}`);
    const response = await module.default.fetch(new Request('http://fixture.invalid/api/demo'), { ASSETS: { fetch: () => new Response('unused') } }, {});
    if ((await response.json()).name !== 'fixture-campaign') throw new Error('Bundled module self-test failed');

    const { createLocalServer } = await import('./local-server.mjs');
    const local = await createLocalServer({ distDirectory: path.join(project, 'dist'), env: {}, port: 0, quiet: true });
    try {
      const base = `http://127.0.0.1:${local.address().port}`;
      if (!(await (await fetch(`${base}/`)).text()).includes('Fixture site')) throw new Error('Local static self-test failed');
      if ((await (await fetch(`${base}/api/demo`)).json()).name !== 'fixture-campaign') throw new Error('Local API self-test failed');
      if ((await fetch(`${base}/unknown`)).status !== 404) throw new Error('Local 404 self-test failed');
      if ((await fetch(`${base}/api/demo`, { method: 'POST', body: 'x'.repeat(4097) })).status !== 413) throw new Error('Local size-limit self-test failed');
    } finally { await new Promise((resolve) => local.close(resolve)); }
    return { ok: true, result: 'PASS', symlink };
  } finally {
    await removeOwnedTree(parent, temporaryRoot, ['live-demo-build-test-']);
  }
}

const invokedDirectly = process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url);
if (invokedDirectly) {
  try {
    const result = process.argv.includes('--self-test') ? await runSelfTest() : await buildProject();
    process.stdout.write(`${JSON.stringify(result)}\n`);
  } catch (error) {
    process.stderr.write(`Build failed: ${error.message}\n`);
    process.exitCode = 1;
  }
}
