import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
const root = new URL('../', import.meta.url);
test('npm manifest is dual-mode, pinned, dependency-free, and allowlisted', () => {
  const p = JSON.parse(readFileSync(new URL('package.json', root), 'utf8'));
  assert.equal(p.name, 'hashrope');
  assert.equal(p.version, '0.3.3');
  assert.equal(p.type, 'module');
  assert.equal(p.license, 'MIT');
  assert.equal(p.engines.node, '>=20.11.1');
  assert.deepEqual(p.dependencies ?? {}, {});
  assert.equal(p.exports['.'].import.types, './dist/esm/index.d.ts');
  assert.equal(p.exports['.'].import.default, './dist/esm/index.js');
  assert.equal(p.exports['.'].require.types, './dist/cjs/index.d.ts');
  assert.equal(p.exports['.'].require.default, './dist/cjs/index.js');
  for (const version of Object.values(p.devDependencies)) assert.match(version, /^\d+\.\d+\.\d+$/);
  assert.deepEqual(p.files, ['dist/', 'src/', 'README.md', 'LICENSE', 'docs/API.md', 'docs/COMPATIBILITY.md', 'docs/PARITY-MATRIX.md', 'docs/TDD.md', 'docs/VERIFICATION.md', 'docs/RELEASE.md', 'docs/ORACLE.md']);
  assert.equal(p.sideEffects, false);
});
test('compiler configurations require exact BigInt-capable output', () => {
  const base = JSON.parse(readFileSync(new URL('tsconfig.json', root), 'utf8'));
  assert.equal(base.compilerOptions.target, 'ES2022');
  assert.equal(base.compilerOptions.strict, true);
  assert.equal(base.compilerOptions.declaration, true);
  assert.equal(base.compilerOptions.noEmitOnError, true);
});
