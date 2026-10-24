import { readFileSync, writeFileSync, mkdirSync, existsSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { gunzipSync } from 'node:zlib';
import { fileURLToPath } from 'node:url';
import path from 'node:path';
import { artifactDirectory, evidenceDirectory } from './release-paths.mjs';
const root = fileURLToPath(new URL('../', import.meta.url));
const npmCli = [process.env.npm_execpath,
  path.join(path.dirname(process.execPath), 'node_modules', 'npm', 'bin', 'npm-cli.js'),
  path.join(path.dirname(process.execPath), '..', 'lib', 'node_modules', 'npm', 'bin', 'npm-cli.js'),
].find(candidate => candidate && candidate.endsWith('npm-cli.js') && existsSync(candidate));
if (!npmCli) throw new Error('Cannot locate npm CLI; run this script through npm run');
function run(args) {
  const r = spawnSync(process.execPath, args, { cwd: root, encoding: 'utf8' });
  if (r.error) throw r.error;
  if (r.status !== 0) throw new Error((r.stdout ?? '') + (r.stderr ?? ''));
  return r.stdout;
}
run(['scripts/build.mjs']);
mkdirSync(artifactDirectory, { recursive: true });
const metadata = JSON.parse(run([npmCli, 'pack', '--ignore-scripts', '--json', '--pack-destination', artifactDirectory]))[0];
const filename = path.join(artifactDirectory, metadata.filename);
const archive = readFileSync(filename);
const tar = gunzipSync(archive);
const readString = (buf) => buf.toString('utf8').replace(/\0.*$/s, '');
const entries = [];
const contents = new Map();
for (let offset = 0; offset + 512 <= tar.length;) {
  const header = tar.subarray(offset, offset + 512);
  if (header.every(b => b === 0)) break;
  const shortName = readString(header.subarray(0, 100));
  const prefix = readString(header.subarray(345, 500));
  const name = prefix ? `${prefix}/${shortName}` : shortName;
  const size = Number.parseInt(readString(header.subarray(124, 136)).trim() || '0', 8);
  if (!Number.isSafeInteger(size) || size < 0 || offset + 512 + size > tar.length) throw new Error('Invalid tar entry size');
  const expectedChecksum = Number.parseInt(readString(header.subarray(148, 156)).trim(), 8);
  let checksum = 0;
  for (let i = 0; i < 512; i++) checksum += i >= 148 && i < 156 ? 32 : header[i];
  if (checksum !== expectedChecksum) throw new Error('Invalid tar header checksum');
  const type = String.fromCharCode(header[156]);
  if (!['0', '\0'].includes(type)) throw new Error(`Unexpected tar entry type: ${type} ${name}`);
  if (!name.startsWith('package/') || name.includes('..') || name.includes('\\')) throw new Error(`Invalid path: ${name}`);
  const rel = name.slice('package/'.length);
  const allow = ['package.json', 'README.md', 'LICENSE', 'docs/API.md', 'docs/COMPATIBILITY.md', 'docs/PARITY-MATRIX.md', 'docs/TDD.md', 'docs/VERIFICATION.md', 'docs/RELEASE.md', 'docs/ORACLE.md'].includes(rel) || /^src\/[a-zA-Z0-9_/-]+\.ts$/.test(rel) || /^dist\/(esm|cjs)\/[a-zA-Z0-9_./-]+\.(js|json|ts|map)$/.test(rel);
  if (!allow || /(^|\/)(\.env|\.npmrc|node_modules|test|tests|oracle|evidence|\.scratch)(\/|$)/i.test(rel)) throw new Error(`Disallowed package file: ${rel}`);
  if (contents.has(rel)) throw new Error(`Duplicate tar entry: ${rel}`);
  const bytes = tar.subarray(offset + 512, offset + 512 + size);
  // Every allowlisted payload is text; reject locale-encoded docs or source.
  new TextDecoder('utf-8', { fatal: true }).decode(bytes);
  contents.set(rel, bytes);
  entries.push({ path: rel, size, sha256: createHash('sha256').update(bytes).digest('hex') });
  offset += 512 + Math.ceil(size / 512) * 512;
}
const manifest = JSON.parse(contents.get('package.json').toString('utf8'));
if (Object.keys(manifest.dependencies ?? {}).length !== 0) throw new Error('Package has runtime dependencies');
for (const required of ['README.md','LICENSE','docs/API.md','docs/COMPATIBILITY.md','docs/PARITY-MATRIX.md','docs/TDD.md','docs/VERIFICATION.md','docs/RELEASE.md','docs/ORACLE.md','dist/esm/index.js','dist/esm/index.d.ts','dist/cjs/index.js','dist/cjs/index.d.ts','dist/cjs/package.json']) if (!contents.has(required)) throw new Error(`Missing release file: ${required}`);
if (JSON.parse(contents.get('dist/cjs/package.json').toString('utf8')).type !== 'commonjs') throw new Error('Missing CommonJS package scope');
const expected = metadata.files.map(f => f.path).sort();
if (JSON.stringify(expected) !== JSON.stringify(entries.map(f => f.path).sort())) throw new Error('Tar contents disagree with npm file manifest');
const integrity = 'sha512-' + createHash('sha512').update(archive).digest('base64');
if (metadata.integrity !== integrity) throw new Error('npm integrity mismatch');
const result = { checkedAt: new Date().toISOString(), filename: metadata.filename, sha256: createHash('sha256').update(archive).digest('hex'), integrity, size: archive.length, unpackedSize: entries.reduce((sum, e) => sum + e.size, 0), files: entries, runtimeDependencies: Object.keys(manifest.dependencies ?? {}), passed: true };
writeFileSync(path.join(evidenceDirectory, 'tarball-inspection.json'), JSON.stringify(result, null, 2) + '\n');
writeFileSync(path.join(artifactDirectory, 'SHA256SUMS'), `${result.sha256}  ${metadata.filename}\n`);
process.stdout.write(`Verified ${metadata.filename}: ${entries.length} files, ${archive.length} bytes, SHA256 ${result.sha256}\n`);
