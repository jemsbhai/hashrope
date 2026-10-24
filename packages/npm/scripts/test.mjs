import { readdirSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
const root = fileURLToPath(new URL('../', import.meta.url));
for (const args of [['scripts/build.mjs'], ['--test', ...readdirSync(new URL('../test/', import.meta.url)).filter(f => f.endsWith('.test.mjs')).sort().map(f => `test/${f}`)]]) {
  const r = spawnSync(process.execPath, args, { cwd: root, stdio: 'inherit' });
  if (r.error) throw r.error;
  if (r.status !== 0) process.exit(r.status ?? 1);
}
