import { mkdirSync, existsSync, readdirSync, readFileSync, copyFileSync, writeFileSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { fileURLToPath } from 'node:url';
import path from 'node:path';
const root = fileURLToPath(new URL('../', import.meta.url));
const [label, command, ...args] = process.argv.slice(2);
if (!label || !/^[a-z0-9][a-z0-9._-]*$/i.test(label) || !command) throw new Error('Usage: node scripts/tdd.mjs LABEL COMMAND [ARGS...]');
const base = path.join(root, 'evidence', 'tdd', label);
let dir = base;
if (existsSync(dir)) dir += '-' + Date.now();
mkdirSync(path.join(dir, 'snapshot'), { recursive: true });
const files = [];
function walk(relative) {
  const absolute = path.join(root, relative);
  if (!existsSync(absolute)) return;
  for (const entry of readdirSync(absolute, { withFileTypes: true })) {
    const rel = path.join(relative, entry.name);
    if (entry.isDirectory()) walk(rel);
    else if (entry.isFile()) files.push(rel);
  }
}
for (const dir of ['src', 'test', 'scripts']) walk(dir);
for (const f of ['package.json', 'package-lock.json', 'tsconfig.json', 'tsconfig.esm.json', 'tsconfig.cjs.json']) if (existsSync(path.join(root, f))) files.push(f);
const hashes = {};
for (const file of files) {
  const data = readFileSync(path.join(root, file));
  hashes[file.replaceAll('\\', '/')] = createHash('sha256').update(data).digest('hex');
  const dest = path.join(dir, 'snapshot', file);
  mkdirSync(path.dirname(dest), { recursive: true });
  copyFileSync(path.join(root, file), dest);
}
const startedAt = new Date().toISOString();
const r = spawnSync(command, args, { cwd: root, encoding: 'utf8', maxBuffer: 64 * 1024 * 1024, shell: process.platform === 'win32' && /\.(cmd|bat)$/i.test(command) });
const output = (r.stdout ?? '') + (r.stderr ?? '') + (r.error ? String(r.error) : '');
writeFileSync(path.join(dir, 'output.txt'), output);
writeFileSync(path.join(dir, 'run.json'), JSON.stringify({ label, command: [command, ...args], startedAt, finishedAt: new Date().toISOString(), exitCode: r.status, signal: r.signal, hashes }, null, 2) + '\n');
process.stdout.write(output);
process.stdout.write(`\nTDD evidence: ${path.relative(root, dir)}; exit=${r.status}\n`);
process.exit(r.status ?? 1);
