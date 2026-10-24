/** Read-only public-registry validation. Never publishes or overwrites a run. */
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { existsSync, mkdirSync, mkdtempSync, readFileSync, readdirSync, writeFileSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const root = fileURLToPath(new URL('../../', import.meta.url));
const version = '0.3.4';
const expectedSha256 = 'be10002475ccdb22738ff4430380130f30176a60d01a37eae71d4cae5b9ac580';
const registry = 'https://registry.npmjs.org/';
const metadataUrl = `${registry}hashrope/${version}`;
const tarballUrl = `${registry}hashrope/-/hashrope-${version}.tgz`;
const sha = (bytes, algorithm = 'sha256', encoding = 'hex') => createHash(algorithm).update(bytes).digest(encoding);
const output = process.argv[2];
assert.ok(output, 'Usage: node evidence/release-0.3.4/verify-public.mjs <new-report-directory>');
const reportDirectory = path.resolve(output);
assert.ok(!existsSync(reportDirectory), 'Use a new report directory; previous evidence is never overwritten');
mkdirSync(reportDirectory, { recursive: true });
mkdirSync(path.join(root, '.scratch'), { recursive: true });
const scratch = mkdtempSync(path.join(root, '.scratch', 'public-0.3.4-'));
const npmCli = [process.env.npm_execpath,
  path.join(path.dirname(process.execPath), 'node_modules/npm/bin/npm-cli.js'),
  path.join(path.dirname(process.execPath), '../lib/node_modules/npm/bin/npm-cli.js')
].find(candidate => candidate?.endsWith('npm-cli.js') && existsSync(candidate));
assert.ok(npmCli, 'Use a Node installation containing npm');
const inspection = JSON.parse(readFileSync(path.join(root, 'evidence/release-0.3.4/tarball-inspection.json'), 'utf8'));
assert.equal(inspection.sha256, expectedSha256);
assert.equal(sha(readFileSync(path.join(root, 'artifacts/0.3.4/hashrope-0.3.4.tgz'))), expectedSha256);
const report = { startedAt: new Date().toISOString(), version, node: process.version,
  verifierSha256: sha(readFileSync(fileURLToPath(import.meta.url))), expectedSha256,
  metadataUrl, tarballUrl, scratch, passed: false };
try {
  const metadataResponse = await fetch(metadataUrl, { signal: AbortSignal.timeout(30000) });
  report.metadataStatus = metadataResponse.status;
  assert.equal(metadataResponse.status, 200);
  assert.equal(metadataResponse.url, metadataUrl);
  const metadata = await metadataResponse.json();
  assert.equal(metadata.name, 'hashrope');
  assert.equal(metadata.version, version);
  assert.equal(metadata.dist.tarball, tarballUrl);
  assert.equal(metadata.dist.integrity, inspection.integrity);
  report.registryIntegrity = metadata.dist.integrity;
  // Canonical URL is required: a cache-bypass download alone cannot pass this gate.
  const response = await fetch(tarballUrl, { signal: AbortSignal.timeout(30000) });
  report.canonicalTarballStatus = response.status;
  report.responseHeaders = Object.fromEntries(['date', 'age', 'cf-cache-status', 'content-length'].map(key => [key, response.headers.get(key)]));
  assert.equal(response.status, 200);
  assert.equal(response.url, tarballUrl);
  const bytes = Buffer.from(await response.arrayBuffer());
  assert.equal(bytes.length, inspection.size);
  assert.equal(sha(bytes), expectedSha256);
  assert.equal(`sha512-${sha(bytes, 'sha512', 'base64')}`, metadata.dist.integrity);
  writeFileSync(path.join(scratch, `hashrope-${version}.tgz`), bytes, { flag: 'wx' });
  report.canonicalTarballSha256 = sha(bytes);

  const consumer = path.join(scratch, 'consumer');
  const cache = path.join(scratch, 'empty-cache');
  mkdirSync(consumer); mkdirSync(cache);
  const userconfig = path.join(scratch, 'empty-user.npmrc');
  const globalconfig = path.join(scratch, 'empty-global.npmrc');
  writeFileSync(userconfig, '', { flag: 'wx' }); writeFileSync(globalconfig, '', { flag: 'wx' });
  writeFileSync(path.join(consumer, 'package.json'), JSON.stringify({ private: true, type: 'module' }), { flag: 'wx' });
  const env = Object.fromEntries(Object.entries(process.env).filter(([key]) => !/^npm_config_/i.test(key) && !/^(npm_token|node_auth_token)$/i.test(key)));
  const args = [npmCli, 'install', `hashrope@${version}`, '--ignore-scripts', `--registry=${registry}`,
    '--cache', cache, '--userconfig', userconfig, '--globalconfig', globalconfig, '--no-audit', '--no-fund',
    '--package-lock=false', '--fetch-retries=0', '--update-notifier=false', '--loglevel=http'];
  report.installCommand = [process.execPath, ...args];
  report.emptyCacheBeforeInstall = readdirSync(cache).length === 0;
  assert.ok(report.emptyCacheBeforeInstall);
  const installed = spawnSync(process.execPath, args, { cwd: consumer, env, encoding: 'utf8', timeout: 120000 });
  report.install = { exitCode: installed.status, stdout: installed.stdout, stderr: installed.stderr };
  if (installed.error) throw installed.error;
  assert.equal(installed.status, 0, installed.stdout + installed.stderr);
  const hex = Buffer.from(metadata.dist.integrity.slice('sha512-'.length), 'base64').toString('hex');
  const cachedTarball = path.join(cache, '_cacache/content-v2/sha512', hex.slice(0, 2), hex.slice(2, 4), hex.slice(4));
  assert.ok(existsSync(cachedTarball), 'No downloaded artifact in the initially empty cache');
  report.installedDownloadSha256 = sha(readFileSync(cachedTarball));
  assert.equal(report.installedDownloadSha256, expectedSha256);
  for (const file of inspection.files) {
    assert.equal(sha(readFileSync(path.join(consumer, 'node_modules/hashrope', file.path))), file.sha256, file.path);
  }
  report.installedFilesVerified = inspection.files.length;
  const smoke = `
import assert from 'node:assert/strict';
import * as esm from 'hashrope';
import { createRequire } from 'node:module';
const cjs = createRequire(import.meta.url)('hashrope');
for (const api of [esm, cjs]) {
  assert.equal(api.VERSION, '0.3.4'); assert.equal(api.PYTHON_VERSION, '0.2.2');
  const h = new api.PolynomialHash(), other = new api.PolynomialHash(api.MERSENNE_127, 257n);
  const leaf = values => new api.Leaf(Uint8Array.from(values), h);
  const a = leaf([0,131,42]), b = leaf([1,0,42]);
  assert.equal(a.hashVal,b.hashVal); assert.equal(api.ropeEqual(a,b),false); assert.equal(api.ropeLcp(a,b),0n);
  assert.equal(api.ropeRangeEqual(a,2,b,2,1),true);
  assert.equal(api.ropeEqual(a,new api.Leaf(Uint8Array.of(0,131,42),other)),true);
  assert.ok(Object.isFrozen(h.profile));
  assert.throws(()=>api.ropeConcat(null,a,other),/Incompatible hash profile/);
  assert.throws(()=>api.ropeRangeLcp(a,4,a,4,0),RangeError);
  const q=1n<<100n,x=new api.RepeatNode(leaf([1,2]),q,h),y=new api.RepeatNode(leaf([1,2,1,2]),q/2n,h);
  assert.equal(api.ropeEqual(x,y),true); assert.equal(api.ropeLcp(x,y),x.len);
  assert.equal(api.ropeRangeLcp(x,1,y,3,x.len-3n),x.len-3n);
  let da=leaf([1,2]),db=leaf([1,2]);
  for(let i=0;i<100;i++){da=new api.Internal(da,da,h);db=new api.Internal(db,db,h);}
  assert.equal(api.ropeLcp(da,db),da.len);
  const custom=new api.PolynomialHash(17,2),ca=new api.Leaf(Uint8Array.of(0,13),custom);
  const cb=new api.Internal(new api.Leaf(Uint8Array.of(0),custom),new api.Leaf(Uint8Array.of(13),custom),custom);
  assert.notEqual(ca.hashVal,cb.hashVal); assert.equal(api.ropeEqual(ca,cb),true);
  const data=Uint8Array.from({length:65536},(_,i)=>i%251),changed=data.slice();changed[65535]^=1;
  assert.equal(api.ropeLcp(new api.Leaf(data,h),new api.Leaf(changed,h)),65535n);
}
console.log('PASS public registry ESM/CJS: version, collision, profile, bounds, compressed ranges/DAGs, legacy custom profile and late mismatch');
`;
  const smokePath = path.join(consumer, 'verify.mjs'); writeFileSync(smokePath, smoke, { flag: 'wx' });
  report.consumerSourceSha256 = sha(Buffer.from(smoke));
  const result = spawnSync(process.execPath, [smokePath], { cwd: consumer, env, encoding: 'utf8', timeout: 30000 });
  report.consumer = { exitCode: result.status, stdout: result.stdout, stderr: result.stderr };
  if (result.error) throw result.error;
  assert.equal(result.status, 0, result.stdout + result.stderr);
  report.passed = true;
} catch (error) { report.error = String(error); throw error; }
finally {
  report.finishedAt = new Date().toISOString();
  writeFileSync(path.join(reportDirectory, 'verification.json'), JSON.stringify(report, null, 2) + '\n', { flag: 'wx' });
  console.log(JSON.stringify({ passed: report.passed, reportDirectory, artifactSha256: report.canonicalTarballSha256 }));
}
