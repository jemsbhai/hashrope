import { mkdtempSync, mkdirSync, writeFileSync, readFileSync, existsSync, cpSync } from 'node:fs';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { createRequire } from 'node:module';
import { fileURLToPath, pathToFileURL } from 'node:url';
import vm from 'node:vm';
import { build } from 'esbuild';

const root = fileURLToPath(new URL('../', import.meta.url));
const require = createRequire(import.meta.url);
const archive = path.resolve(process.argv[2] ?? path.join(root, 'artifacts', `hashrope-${JSON.parse(readFileSync(path.join(root, 'package.json'), 'utf8')).version}.tgz`));
if (!existsSync(archive)) throw new Error(`Pack the release before consumer testing: ${archive}`);
const cleanRoot = mkdtempSync(path.join(tmpdir(), 'hashrope-consumers-'));
const npmCli = [process.env.npm_execpath,
  path.join(path.dirname(process.execPath), 'node_modules', 'npm', 'bin', 'npm-cli.js'),
  path.join(path.dirname(process.execPath), '..', 'lib', 'node_modules', 'npm', 'bin', 'npm-cli.js'),
].find(candidate => candidate && candidate.endsWith('npm-cli.js') && existsSync(candidate));
if (!npmCli) throw new Error('Cannot locate npm CLI; run this script through npm run');
const report = { startedAt: new Date().toISOString(), node: process.version, archive: path.basename(archive), archiveSha256: createHash('sha256').update(readFileSync(archive)).digest('hex'), cleanDirectory: cleanRoot, checks: [] };
function run(name, cwd, executable, args) {
  const r = spawnSync(executable, args, { cwd, encoding: 'utf8' });
  report.checks.push({ name, exitCode: r.status, output: (r.stdout ?? '') + (r.stderr ?? '') });
  if (r.error) throw r.error;
  if (r.status !== 0) throw new Error(`${name} failed (${r.status}): ${r.stdout ?? ''}${r.stderr ?? ''}`);
  process.stdout.write(`PASS ${name}\n`);
}
function consumer(name, type = 'module') {
  const dir = path.join(cleanRoot, name);
  mkdirSync(dir);
  writeFileSync(path.join(dir, 'package.json'), JSON.stringify({ name: `hashrope-clean-${name}`, private: true, type }));
  run(`${name}: clean offline tarball install`, dir, process.execPath, [npmCli, 'install', archive, '--ignore-scripts', '--offline', '--no-audit', '--no-fund', '--package-lock=false']);
  const manifest = JSON.parse(readFileSync(path.join(dir, 'node_modules/hashrope/package.json'), 'utf8'));
  if (Object.keys(manifest.dependencies ?? {}).length) throw new Error('Unexpected runtime dependencies');
  return dir;
}
function installedOracle(flavor, dir) {
  const target = path.join(dir, 'test');
  mkdirSync(target, { recursive: true });
  cpSync(path.join(root, 'test', 'fixtures'), path.join(target, 'fixtures'), { recursive: true });
  const source = readFileSync(path.join(root, 'test', 'oracle.test.mjs'), 'utf8');
  const original = "import * as api from '../dist/esm/index.js';";
  if (source.split(original).length !== 2) throw new Error('Oracle adapter import changed; review the one-line consumer adaptation');
  const replacement = flavor === 'esm' ? "import * as api from 'hashrope';" : "import { createRequire } from 'node:module';\nconst api = createRequire(import.meta.url)('hashrope');";
  const adapted = source.replace(original, replacement);
  writeFileSync(path.join(target, 'oracle.test.mjs'), adapted);
  const r = spawnSync(process.execPath, ['--test', '--test-reporter=tap', 'test/oracle.test.mjs'], { cwd: dir, encoding: 'utf8', maxBuffer: 64 * 1024 * 1024, timeout: 300000 });
  const output = (r.stdout ?? '') + (r.stderr ?? '') + (r.error ? String(r.error) : '');
  const log = `consumer-oracle-${process.version}-${flavor}.tap`;
  writeFileSync(path.join(root, 'evidence', log), output);
  const summary = Object.fromEntries([...output.matchAll(/^# (tests|suites|pass|fail|cancelled|skipped|todo) (\d+)$/gm)].map(m => [m[1], Number(m[2])]));
  report.checks.push({ name: `${flavor}: installed-tarball published-wheel oracle parity`, exitCode: r.status, summary, log, adapterSha256: createHash('sha256').update(adapted).digest('hex'), corpusSha256: createHash('sha256').update(readFileSync(path.join(target, 'fixtures', 'oracle.json'))).digest('hex') });
  if (r.error) throw r.error;
  if (r.status !== 0 || summary.tests !== 27565 || summary.pass !== 27565 || summary.fail !== 0 || summary.skipped !== 0) throw new Error(`${flavor} installed oracle failed: ${output.slice(-6000)}`);
  process.stdout.write(`PASS ${flavor}: all ${summary.pass} published-wheel oracle checks from installed tarball\n`);
}
function installedExact(flavor, dir) {
  const source = readFileSync(path.join(root, 'test', 'exact.test.mjs'), 'utf8');
  const original = "import * as api from '../dist/esm/index.js';";
  if (source.split(original).length !== 2) throw new Error('Exact-test adapter import changed');
  const replacement = flavor === 'esm' ? "import * as api from 'hashrope';" : "import { createRequire } from 'node:module';\nconst api = createRequire(import.meta.url)('hashrope');";
  const adapted = source.replace(original, replacement);
  writeFileSync(path.join(dir, 'test', 'exact.test.mjs'), adapted);
  run(`${flavor}: installed exact and profile regressions`, dir, process.execPath, ['--test', 'test/exact.test.mjs']);
  report.checks.at(-1).adapterSha256 = createHash('sha256').update(adapted).digest('hex');
}
function documentedExamples(dir) {
  const files = [];
  const sources = [];
  for (const name of ['README.md', 'docs/API.md']) {
    const markdown = readFileSync(path.join(dir, 'node_modules', 'hashrope', name), 'utf8');
    const snippets = [...markdown.matchAll(/```ts\s*\r?\n([\s\S]*?)```/g)];
    for (let index = 0; index < snippets.length; index++) {
      const file = `${name === 'README.md' ? 'readme' : 'api'}-${index + 1}.mts`;
      writeFileSync(path.join(dir, file), snippets[index][1]);
      files.push(file);
      sources.push({ source: name, snippet: index + 1, file, sourceSha256: createHash('sha256').update(snippets[index][1]).digest('hex') });
    }
  }
  if (files.length !== 6) throw new Error(`Expected all six packaged TS examples; found ${files.length}`);
  const c = spawnSync(process.execPath, [require.resolve('typescript/bin/tsc'), '--strict', '--target', 'ES2022', '--module', 'NodeNext', '--moduleResolution', 'NodeNext', '--lib', 'ES2022,DOM', '--outDir', 'examples-out', ...files], { cwd: dir, encoding: 'utf8', maxBuffer: 8 * 1024 * 1024 });
  const check = { name: 'All packaged README and API TypeScript examples', compiler: 'TypeScript 5.9.3', compileExitCode: c.status, compileOutput: (c.stdout ?? '') + (c.stderr ?? ''), examples: sources };
  report.checks.push(check);
  if (c.error) throw c.error;
  if (c.status !== 0) throw new Error(`Packaged TypeScript examples failed compilation: ${check.compileOutput}`);
  for (const example of sources) {
    const compiled = `./examples-out/${example.file.replace(/\.mts$/, '.mjs')}`;
    const program = `console.assert = (condition, ...details) => { if (!condition) throw new Error('Example assertion failed: ' + details.join(' ')); }; await import(${JSON.stringify(compiled)});`;
    const result = spawnSync(process.execPath, ['--input-type=module', '--eval', program], { cwd: dir, encoding: 'utf8' });
    example.exitCode = result.status;
    example.output = (result.stdout ?? '') + (result.stderr ?? '');
    if (result.error) throw result.error;
    if (result.status !== 0) throw new Error(`Packaged example ${example.file} failed: ${example.output}`);
  }
  check.exitCode = 0;
  process.stdout.write(`PASS all ${files.length} packaged TypeScript examples, compiled and independently executed\n`);
}
const names = 'PolynomialHash, MERSENNE_61, MERSENNE_127, Leaf, ropeConcat, ropeSplit, ropeRepeat, ropeSubstrHash, ropeLen, ropeHash, ropeToBytes, validateRope, SlidingWindow, HashProfile, ropeEqual, ropeLcp, ropeRangeEqual, ropeRangeLcp';
const smoke = `
const check = (x, message) => { if (!x) throw new Error(message); };
const same = (a, b) => a.length === b.length && a.every((v, i) => v === b[i]);
const h = new PolynomialHash();
const data = Uint8Array.of(97, 98);
const leaf = new Leaf(data, h);
data[0] = 0;
const exposed = leaf.data;
exposed[1] = 0;
check(same(ropeToBytes(leaf), Uint8Array.of(97, 98)), 'defensive ownership');
const combined = ropeConcat(leaf, leaf, h);
check(ropeHash(combined) === h.hash(Uint8Array.of(97, 98, 97, 98)), 'concat');
const [left, right] = ropeSplit(combined, 3n, h);
check(ropeLen(left) === 3n && ropeLen(right) === 1n, 'split');
validateRope(combined);
const q = 2n ** 54n + 3n;
const huge = ropeRepeat(leaf, q, h);
check(ropeLen(huge) === 2n * q, 'exact huge length');
check(ropeSubstrHash(huge, 2n * q - 3n, 3n, h) === h.hash(Uint8Array.of(98, 97, 98)), 'huge substring');
const [hl, hr] = ropeSplit(huge, 2n * q - 1n, h);
check(ropeLen(hl) === 2n * q - 1n && ropeLen(hr) === 1n, 'huge split');
for (const prime of [MERSENNE_61, MERSENNE_127]) {
  const profile = new PolynomialHash(prime, 131n);
  check(profile.hash(Uint8Array.of(0, 131)) === profile.hash(Uint8Array.of(1, 0)), 'explicit collision retained');
}
check(ropeEqual(huge, ropeRepeat(new Leaf(Uint8Array.of(97, 98), new PolynomialHash(MERSENNE_127)), q, new PolynomialHash(MERSENNE_127))), 'exact cross-profile compressed equality');
check(ropeLcp(huge, huge) === 2n * q, 'exact huge LCP');
check(ropeRangeEqual(huge, 1n, huge, 3n, 2n * q - 3n), 'exact shifted periodic range');
check(ropeRangeLcp(huge, 0, huge, 1, 20) === 0n, 'exact range mismatch');
check(!ropeEqual(new Leaf(Uint8Array.of(0, 131), h), new Leaf(Uint8Array.of(1, 0), h)), 'collision is not equality');
check(Object.isFrozen(new HashProfile()), 'immutable profile');
let rejected = false;
try { ropeConcat(null, leaf, new PolynomialHash(MERSENNE_127)); } catch (error) { rejected = error instanceof RangeError; }
check(rejected, 'mixed profile no-op rejected');
const sw = new SlidingWindow(8n, 4n);
sw.appendBytes(Uint8Array.of(97, 98));
sw.appendCopy(2n, 5n);
check(sw.currentHash() === h.hash(Uint8Array.of(97, 98, 97, 98, 97, 98, 97)), 'sliding overlap');
check(sw.pos === 7n, 'stream count');
globalThis.hashropeSmoke = { passed: true, hugeLength: ropeLen(huge).toString() };
`;
try {
  const esm = consumer('esm');
  writeFileSync(path.join(esm, 'smoke.mjs'), `import { ${names} } from 'hashrope';\n` + smoke);
  run('ESM import runtime', esm, process.execPath, ['smoke.mjs']);
  installedOracle('esm', esm);
  installedExact('esm', esm);
  documentedExamples(esm);
  const cjs = consumer('cjs', 'commonjs');
  writeFileSync(path.join(cjs, 'smoke.cjs'), `const { ${names} } = require('hashrope');\n` + smoke);
  run('CommonJS require runtime', cjs, process.execPath, ['smoke.cjs']);
  installedOracle('cjs', cjs);
  installedExact('cjs', cjs);
  const types = consumer('types');
  const typed = `import { PolynomialHash, Leaf, ropeLen, ropeSplit, ropeRepeat, ropeToBytes, SlidingWindow, HashProfile, ropeEqual, ropeLcp, ropeRangeEqual, ropeRangeLcp, type Node } from 'hashrope';
const h = new PolynomialHash();
const leaf: Node = new Leaf(new Uint8Array([1, 2]), h);
const big: Node = ropeRepeat(leaf, 2n ** 70n, h);
const length: bigint = ropeLen(big);
const split: readonly [Node, Node] = ropeSplit(big, length - 1n, h);
const bytes: Uint8Array = ropeToBytes(split[1]);
const stream = new SlidingWindow(8n, 4n);
stream.appendBytes(bytes);
stream.appendCopy(1n, 1n);
const hash: bigint = stream.currentHash();
const profile: HashProfile = h.profile;
const exact: boolean = ropeEqual(leaf, leaf) && ropeRangeEqual(big, 1n, big, 3n, 2n);
const prefix: bigint = ropeLcp(big, big) + ropeRangeLcp(big, 1n, big, 3n, 2n);
void profile; void exact; void prefix;
// @ts-expect-error integer positions cannot be strings
ropeSplit(big, '2', h);
// @ts-expect-error bytes must be explicitly encoded
new Leaf('bytes', h);
void hash;
`;
  writeFileSync(path.join(types, 'esm.mts'), typed);
  writeFileSync(path.join(types, 'cjs.cts'), typed);
  run('TypeScript NodeNext import and require declarations', types, process.execPath, [require.resolve('typescript/bin/tsc'), '--noEmit', '--strict', '--target', 'ES2022', '--module', 'NodeNext', '--moduleResolution', 'NodeNext', '--lib', 'ES2022', 'esm.mts', 'cjs.cts']);
  const browser = consumer('browser');
  const entry = path.join(browser, 'smoke.mjs');
  writeFileSync(entry, `import { ${names} } from 'hashrope';\n` + smoke);
  const bundle = await build({ entryPoints: [entry], bundle: true, write: false, platform: 'browser', format: 'iife', target: 'es2022', metafile: true });
  for (const meta of Object.values(bundle.metafile.outputs)) if (meta.imports.length !== 0) throw new Error('Browser bundle has unresolved runtime imports');
  const js = bundle.outputFiles[0].text;
  writeFileSync(path.join(browser, 'browser-bundle.js'), js);
  const sandbox = Object.create(null);
  vm.runInNewContext(js, sandbox, { timeout: 10000 });
  if (sandbox.hashropeSmoke?.passed !== true) throw new Error('Browser-style realm smoke did not pass');
  report.checks.push({ name: 'ES2022 browser bundle; VM realm without Node globals', exitCode: 0, output: JSON.stringify(sandbox.hashropeSmoke), browserBundleBytes: Buffer.byteLength(js), bundleSha256: createHash('sha256').update(js).digest('hex'), caveat: 'This is a browser-targeted bundle executed in a Node VM without Node globals, not a real browser-engine certification.' });
  process.stdout.write('PASS ES2022 browser bundle and isolated realm\n');
  const browserExe = ['C:/Program Files/Google/Chrome/Application/chrome.exe', 'C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe'].find(existsSync);
  if (browserExe) {
    const html = path.join(browser, 'index.html');
    writeFileSync(html, '<!doctype html><html><head><meta charset="utf-8"><title>Hashrope clean consumer</title></head><body><pre id="result">PENDING</pre><script src="browser-bundle.js"></script><script>document.getElementById("result").textContent = globalThis.hashropeSmoke?.passed ? "PASS " + JSON.stringify(globalThis.hashropeSmoke) : "FAIL";</script></body></html>');
    const profile = path.join(browser, 'fresh-browser-profile');
    const r = spawnSync(browserExe, ['--headless=new', '--no-first-run', '--disable-gpu', '--disable-extensions', '--disable-background-networking', `--user-data-dir=${profile}`, '--dump-dom', pathToFileURL(html).href], { encoding: 'utf8', timeout: 30000, windowsHide: true });
    const passed = r.status === 0 && /<pre id="result">PASS /.test(r.stdout ?? '');
    const version = spawnSync('powershell.exe', ['-NoProfile', '-NonInteractive', '-Command', `(Get-Item -LiteralPath '${browserExe.replaceAll("'", "''")}').VersionInfo.ProductVersion`], { encoding: 'utf8', windowsHide: true }).stdout?.trim();
    report.checks.push({ name: 'Real headless browser with fresh temporary profile', exitCode: r.status, executable: browserExe, version, passed, output: r.stdout ?? '', diagnostics: r.stderr ?? '' });
    if (!passed) throw new Error(`Headless browser failed: ${r.error ?? r.stderr ?? r.stdout}`);
    process.stdout.write(`PASS real headless browser ${version}\n`);
  } else {
    report.checks.push({ name: 'Real headless browser', status: 'not_run', reason: 'No supported local Chrome or Edge executable found; browser bundle and isolated realm were tested.' });
  }
  report.finishedAt = new Date().toISOString();
  report.passed = true;
} catch (e) {
  report.passed = false;
  report.error = String(e);
  throw e;
} finally {
  mkdirSync(path.join(root, 'evidence'), { recursive: true });
  for (const name of ['consumer-verification.json', `consumer-verification-${process.version}.json`]) writeFileSync(path.join(root, 'evidence', name), JSON.stringify(report, null, 2) + '\n');
}
