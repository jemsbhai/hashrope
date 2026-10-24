/** Published-wheel differential tests. No expected value is calculated by the JS implementation. */
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import * as api from '../dist/esm/index.js';
const bytes = readFileSync(new URL('./fixtures/oracle.json', import.meta.url));
const corpus = JSON.parse(bytes);
const manifest = JSON.parse(readFileSync(new URL('./fixtures/manifest.json', import.meta.url)));
const sha256 = value => createHash('sha256').update(value).digest('hex');
const PIN = '710112e3aa14a991918a15331c1f418ce9d1c9d93c34a9da09e22be33c68a2d2';
const fromHex = hex => Uint8Array.from(Buffer.from(hex, 'hex'));
const toHex = bytes => Buffer.from(bytes).toString('hex');
const hashers = new Map();
function hasher(profile) {
  if (!hashers.has(profile)) {
    const p = corpus.profiles[profile];
    hashers.set(profile, new api.PolynomialHash(BigInt(p.prime), BigInt(p.base)));
  }
  return hashers.get(profile);
}
function make(spec, h) {
  if (spec === null) return null;
  switch (spec[0]) {
    case 'leaf': return api.ropeFromBytes(fromHex(spec[1]), h);
    case 'directLeaf': return new api.Leaf(fromHex(spec[1]), h);
    case 'concat': return api.ropeConcat(make(spec[1], h), make(spec[2], h), h);
    case 'internal': return new api.Internal(make(spec[1], h), make(spec[2], h), h);
    case 'repeat': return api.ropeRepeat(make(spec[1], h), BigInt(spec[2]), h);
    case 'directRepeat': return new api.RepeatNode(make(spec[1], h), BigInt(spec[2]), h);
    default: throw new Error(`Unknown fixture node: ${spec[0]}`);
  }
}
function shape(n) {
  if (n === null) return null;
  if (n instanceof api.Leaf) return ['L', toHex(n.data)];
  if (n instanceof api.Internal) return ['I', shape(n.left), shape(n.right)];
  assert.ok(n instanceof api.RepeatNode);
  return ['R', shape(n.child), n.reps.toString()];
}
function state(n) {
  let valid = true;
  try { api.validateRope(n); } catch (error) {
    assert.equal(error.name, 'AssertionError', 'Only structural assertion failures represent invalid oracle nodes');
    valid = false;
  }
  const len = api.ropeLen(n);
  assert.equal(typeof len, 'bigint');
  assert.equal(typeof api.ropeHash(n), 'bigint');
  if (n !== null) { assert.equal(typeof n.weight, 'bigint'); assert.equal(typeof n.hashVal, 'bigint'); }
  return {
    len: len.toString(), hash: api.ropeHash(n).toString(), weight: (n?.weight ?? 0n).toString(),
    height: api.ropeHeight(n).toString(), bytes: len <= 8192n ? toHex(api.ropeToBytes(n)) : null,
    valid, shapeSha256: sha256(JSON.stringify(shape(n))),
  };
}
function streamState(sw) {
  return {
    hPrefix: sw.hPrefix.toString(), lPrefix: sw.lPrefix.toString(), pos: sw.pos.toString(),
    windowLen: sw.windowLen.toString(), W: sw.W.toString(), dMax: sw.dMax.toString(), mMax: sw.mMax.toString(),
    currentHash: sw.currentHash().toString(), finalHash: sw.finalHash().toString(), window: state(sw.rWindow),
  };
}
function checkError(error, expected) {
  assert.ok(error, `Expected Python ${expected.type}: ${expected.message}`);
  const mapped = { ValueError: 'RangeError', ZeroDivisionError: 'RangeError', AttributeError: 'TypeError', TypeError: 'TypeError', AssertionError: 'AssertionError' }[expected.type];
  assert.ok(mapped, `Unmapped Python exception ${expected.type}`);
  assert.equal(error.name, mapped);
  // VM diagnostics differ; explicit upstream library messages and inversion errors remain exact.
  if (expected.type === 'ValueError' || expected.type === 'AssertionError') assert.equal(error.message, expected.message);
}
const traces = new Map();
function trace(name, profile) {
  const key = `${profile}/${name}`;
  if (!traces.has(key)) {
    const desc = corpus.streams[name]; const h = hasher(profile);
    const sw = new api.SlidingWindow(BigInt(desc.dMax), BigInt(desc.mMax), h.prime, h.base);
    const snapshots = [];
    for (const op of desc.steps) {
      let error = null;
      try {
        if (op[0] === 'bytes') sw.appendBytes(fromHex(op[1]));
        else sw.appendCopy(BigInt(op[1]), BigInt(op[2]));
      } catch (caught) { error = caught; }
      snapshots.push({ error, state: streamState(sw) });
    }
    traces.set(key, snapshots);
  }
  return traces.get(key);
}
function run(row) {
  const h = hasher(row.profile), a = row.args;
  switch (row.kind) {
    case 'hash': return h.hash(fromHex(a[0])).toString();
    case 'power': return h.power(BigInt(a[0])).toString();
    case 'phi': return api.phi(...a.map(BigInt)).toString();
    case 'mod': return api.mersenneMod(...a.map(BigInt)).toString();
    case 'mul': return api.mersenneMul(...a.map(BigInt)).toString();
    case 'concatHash': return h.hashConcat(...a.map(BigInt)).toString();
    case 'repeatHash': return h.hashRepeat(...a.map(BigInt)).toString();
    case 'overlapHash': return h.hashOverlap(...a.map(BigInt)).toString();
    case 'rope': case 'nodeError': return state(make(corpus.ropes[a[0]], h));
    case 'split': return api.ropeSplit(make(corpus.ropes[a[0]], h), BigInt(a[1]), h).map(state);
    case 'substr': return api.ropeSubstrHash(make(corpus.ropes[a[0]], h), BigInt(a[1]), BigInt(a[2]), h).toString();
    case 'validate': api.validateRope(make(corpus.ropes[a[0]], h)); return null;
    case 'constructorError': return new api.PolynomialHash(...a.map(BigInt));
    case 'windowConstructorError': return new api.SlidingWindow(...a.map(BigInt));
    case 'sliding': {
      const snapshot = trace(a[0], row.profile)[a[1]];
      if (snapshot.error) throw snapshot.error;
      return snapshot.state;
    }
    case 'slidingEdge': {
      const snapshot = trace(a[0], row.profile)[a[1]];
      if (row.expected.operationError) checkError(snapshot.error, row.expected.operationError);
      else assert.equal(snapshot.error, null);
      return { operationError: row.expected.operationError, state: snapshot.state };
    }
    default: throw new Error(`Unexercised fixture kind: ${row.kind}`);
  }
}
test('oracle manifest pins the published wheel, corpus bytes, exports and complete record inventory', () => {
  assert.equal(manifest.wheel.sha256, PIN);
  assert.equal(manifest.pythonVersion, '0.2.2');
  assert.equal(manifest.corpus.sha256, sha256(bytes));
  assert.equal(manifest.corpus.records, corpus.rows.length);
  assert.equal(manifest.pythonExports.length, 21);
  assert.equal(new Set(corpus.rows.map(row => row.id)).size, corpus.rows.length);
  assert.ok(Object.values(manifest.importOrigins).every(origin => origin.includes('.whl/')));
  const kinds = {}, countErrors = corpus.rows.filter(row => row.error).length;
  for (const row of corpus.rows) kinds[row.kind] = (kinds[row.kind] ?? 0) + 1;
  assert.deepEqual(kinds, manifest.corpus.kinds);
  assert.equal(countErrors, manifest.corpus.errors);
});
for (const row of corpus.rows) {
  test(`Python 0.2.2 ${row.id} [${row.profile}]`, () => {
    if (row.error) {
      let caught;
      try { run(row); } catch (error) { caught = error; }
      checkError(caught, row.error);
    } else assert.deepEqual(run(row), row.expected);
  });
}
