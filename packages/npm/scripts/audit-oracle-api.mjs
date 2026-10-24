/** Audit every published Python export and public class member against the JS mapping. */
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync } from 'node:fs';
import * as api from '../dist/esm/index.js';
import path from 'node:path';
import { evidenceDirectory } from './release-paths.mjs';
const manifest = JSON.parse(readFileSync(new URL('../test/fixtures/manifest.json', import.meta.url)));
const camel = name => name.replace(/_([a-z])/g, (_, c) => c.toUpperCase());
const mappings = {};
for (const name of manifest.pythonExports) {
  if (name === 'Node') {
    const declarations = readFileSync(new URL('../dist/esm/index.d.ts', import.meta.url), 'utf8');
    assert.match(declarations, /export type \{ Node, NonEmptyNode \}/);
    mappings.Node = { js: 'Node', kind: 'TypeScript type-only union; null is empty' };
  } else {
    const js = /^[A-Z]/.test(name) ? name : camel(name);
    assert.ok(js in api, `Missing JS mapping for Python ${name}: ${js}`);
    assert.equal(typeof api[js], name.startsWith('MERSENNE_') ? 'bigint' : 'function');
    mappings[name] = { js, kind: typeof api[js] };
  }
}
assert.equal(api.VERSION, JSON.parse(readFileSync(new URL('../package.json', import.meta.url))).version);
assert.equal(api.PYTHON_VERSION, manifest.pythonVersion);
assert.equal(api.ALPHA_NUM, BigInt(manifest.balanceConstants.ALPHA_NUM));
assert.equal(api.ALPHA_DEN, BigInt(manifest.balanceConstants.ALPHA_DEN));
const h = new api.PolynomialHash(), leaf = new api.Leaf(new Uint8Array([97]), h);
const instances = { PolynomialHash: h, Leaf: leaf, Internal: new api.Internal(leaf, leaf, h), RepeatNode: new api.RepeatNode(leaf, 2n, h), SlidingWindow: new api.SlidingWindow() };
const publicMembers = {};
for (const [name, members] of Object.entries(manifest.classPublicMembers)) {
  publicMembers[name] = {};
  for (const member of members) {
    const js = camel(member);
    assert.ok(js in instances[name], `Missing public ${name}.${member} -> ${js}`);
    publicMembers[name][member] = js;
  }
}
assert.equal(h.prime, api.MERSENNE_61); assert.equal(h.base, 131n);
assert.equal(instances.SlidingWindow.h.prime, h.prime);
assert.equal(instances.SlidingWindow.h.base, h.base);
const result = { pythonVersion: manifest.pythonVersion, exportedNames: Object.keys(mappings).length, mappings, publicMembers, balanceConstantsChecked: true, versionsChecked: true, result: 'all mapped exports and public members present' };
writeFileSync(path.join(evidenceDirectory, 'oracle-api-audit.json'), JSON.stringify(result, null, 2) + '\n');
console.log(JSON.stringify({ exportedNames: result.exportedNames, classes: Object.keys(publicMembers).length, result: result.result }));
