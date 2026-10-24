import { mkdirSync, writeFileSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';
const require = createRequire(import.meta.url);
const root = fileURLToPath(new URL('../', import.meta.url));
const tsc = require.resolve('typescript/bin/tsc');
for (const mode of ['esm', 'cjs']) {
  const result = spawnSync(process.execPath, [tsc, '-p', `tsconfig.${mode}.json`], { cwd: root, stdio: 'inherit' });
  if (result.error) throw result.error;
  if (result.status !== 0) process.exit(result.status ?? 1);
}
mkdirSync(new URL('../dist/cjs/', import.meta.url), { recursive: true });
writeFileSync(new URL('../dist/cjs/package.json', import.meta.url), '{"type":"commonjs"}\n');
