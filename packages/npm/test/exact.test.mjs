import test from 'node:test';
import assert from 'node:assert/strict';
import * as api from '../dist/esm/index.js';

const { PolynomialHash, HashProfile, MERSENNE_61, MERSENNE_127, Leaf, Internal, RepeatNode,
  ropeEqual, ropeLcp, ropeRangeEqual, ropeRangeLcp, contentEquals, ropeConcat, ropeSplit,
  ropeRepeat, ropeSubstrHash, ropeInsert, ropeDelete, ropeToBytes, validateRope, SlidingWindow } = api;
const h = new PolynomialHash();
const leaf = (bytes, hash = h) => new Leaf(Uint8Array.from(bytes), hash);
const lcp = (a, b) => { let i = 0; while (i < Math.min(a.length, b.length) && a[i] === b[i]) i++; return BigInt(i); };

test('exact: default collision is rejected in both Mersenne profiles and tree shapes', () => {
  for (const prime of [MERSENNE_61, MERSENNE_127]) {
    const hash = new PolynomialHash(prime);
    const a = leaf([0, 131, 42], hash), b = leaf([1, 0, 42], hash);
    assert.equal(a.hashVal, b.hashVal);
    for (const left of [a, new Internal(leaf([0], hash), leaf([131, 42], hash), hash)]) {
      assert.equal(ropeEqual(left, b), false);
      assert.equal(contentEquals(left, b), false);
      assert.equal(ropeLcp(left, b), 0n);
      assert.equal(ropeRangeEqual(left, 2, b, 2, 1), true);
    }
  }
});

test('exact: legacy custom profiles may give unequal hashes for equal bytes', () => {
  const hash = new PolynomialHash(17, 2);
  const a = leaf([0, 13], hash), b = new Internal(leaf([0], hash), leaf([13], hash), hash);
  assert.notEqual(a.hashVal, b.hashVal);
  assert.equal(ropeEqual(a, b), true);
  assert.equal(ropeLcp(a, b), 2n);
  assert.equal(ropeEqual(a, leaf([0, 13])), true);
});

test('exact: comparison never invokes a hasher or a defensive leaf-data getter', () => {
  const hash = new PolynomialHash();
  class NoGetterLeaf extends Leaf { get data() { throw new Error('No copied buffers'); } }
  const a = new NoGetterLeaf(Uint8Array.of(1, 2, 3), hash);
  const b = new Internal(leaf([1], hash), leaf([2, 3], hash), hash);
  for (const name of ['hash', 'power', 'hashConcat', 'hashRepeat', 'hashOverlap']) hash[name] = () => { throw new Error('No hashing'); };
  assert.equal(ropeEqual(a, b), true);
  assert.equal(ropeLcp(a, b), 3n);
  assert.equal(ropeRangeEqual(a, 1, b, 1, 2), true);
});

test('exact: differently sized periods and shifted huge repeats use finite witnesses', () => {
  const q = 1n << 110n;
  const a = new RepeatNode(leaf([1, 2]), 2n * q, h);
  const b = new RepeatNode(leaf([1, 2, 1, 2]), q, h);
  assert.equal(ropeEqual(a, b), true);
  assert.equal(ropeLcp(a, b), 4n * q);
  assert.equal(ropeRangeLcp(a, 1n, b, 3n, 4n * q - 3n), 4n * q - 3n);
  const shifted = new Internal(leaf([2]), new RepeatNode(leaf([1, 2]), 2n * q, h), h);
  assert.equal(ropeRangeLcp(a, 1n, shifted, 0n, a.len - 1n), a.len - 1n);
  const late = new Internal(b, leaf([9]), h);
  const good = new Internal(a, leaf([1]), h);
  assert.equal(ropeLcp(good, late), a.len);
  assert.equal(ropeEqual(good, late), false);
});

test('exact: independent shared DAGs are memoized without expanding logical bytes', () => {
  let a = leaf([7, 8]), b = leaf([7, 8]);
  for (let i = 0; i < 100; i++) { a = new Internal(a, a, h); b = new Internal(b, b, h); }
  assert.equal(ropeLcp(a, b), 2n << 100n);
  assert.equal(contentEquals(a, b), true);
  assert.equal(ropeLcp(new Internal(a, leaf([4]), h), new Internal(b, leaf([5]), h)), a.len);
  assert.equal(ropeRangeLcp(a, 17n, b, 19n, a.len - 19n), a.len - 19n);
});

test('exact: deterministic random trees and ranges agree with materialized byte oracles', () => {
  let state = 91731;
  const random = n => { state = (Math.imul(state, 1664525) + 1013904223) >>> 0; return state % n; };
  const alternate = new PolynomialHash(MERSENNE_127, 257);
  const shape = (data, hash) => {
    if (data.length === 0) return null;
    if (data.length < 2 || random(3) === 0) return leaf(data, hash);
    const cut = 1 + random(data.length - 1);
    return new Internal(shape(data.slice(0, cut), hash), shape(data.slice(cut), hash), hash);
  };
  for (let trial = 0; trial < 250; trial++) {
    const data = Uint8Array.from({length: 1 + random(65)}, () => random(256));
    const reps = 2 + random(4);
    const a = new RepeatNode(shape(data, h), reps, h);
    const bytes = ropeToBytes(a), changed = bytes.slice();
    if (trial % 3) changed[random(changed.length)] ^= 1;
    const b = shape(changed, alternate);
    assert.equal(ropeLcp(a, b), lcp(bytes, changed));
    assert.equal(ropeEqual(a, b), Buffer.from(bytes).equals(changed));
    for (let range = 0; range < 12; range++) {
      const sa = random(bytes.length + 1), sb = random(changed.length + 1);
      const size = random(Math.min(bytes.length - sa, changed.length - sb) + 1);
      const expected = lcp(bytes.subarray(sa, sa + size), changed.subarray(sb, sb + size));
      assert.equal(ropeRangeLcp(a, sa, b, sb, size), expected);
      assert.equal(ropeRangeEqual(a, sa, b, sb, size), expected === BigInt(size));
    }
  }
});

test('exact: strict range validation precedes identical and empty shortcuts', () => {
  const a = leaf([1, 2]);
  assert.equal(ropeEqual(null, null), true);
  assert.equal(ropeEqual(a, null), false);
  assert.equal(ropeLcp(a, null), 0n);
  assert.equal(ropeRangeLcp(null, 0, null, 0, 0), 0n);
  assert.equal(ropeRangeLcp(a, 2, a, 2, 0), 0n);
  for (const args of [[a, -1, a, -1, 0], [a, 3, a, 3, 0], [a, 0, a, 0, 3], [a, 0, a, 0, -1], [null, 1, null, 1, 0], [a, Number.MAX_SAFE_INTEGER + 1, a, 0, 0]]) {
    assert.throws(() => ropeRangeLcp(...args), RangeError);
    assert.throws(() => ropeRangeEqual(...args), RangeError);
  }
  assert.throws(() => ropeRangeLcp(a, '0', a, 0, 0), TypeError);
});

test('exact: word scans agree with byte oracles at every mismatch and unaligned tail', () => {
  for (const size of [0, 1, 7, 8, 15, 16, 17, 31, 32, 33, 63, 64, 65, 71, 72, 79, 127, 128, 129]) {
    for (let offset = 0; offset < 8; offset++) {
      const storage = Uint8Array.from({length:size + 32}, (_, i) => (i * 31) & 255);
      const data = storage.subarray(offset, offset + size + 8);
      const a = leaf(data);
      for (let mismatch = 0; mismatch <= size; mismatch++) {
        const other = new Uint8Array(size + 16);
        other.set(data.subarray(3, size + 3), 5);
        if (mismatch < size) other[5 + mismatch] ^= 1;
        const b = leaf(other);
        assert.equal(ropeRangeLcp(a, 3, b, 5, size), BigInt(mismatch));
        assert.equal(ropeRangeEqual(a, 3, b, 5, size), mismatch === size);
      }
    }
  }
});

test('profile: immutable metadata, equivalent hashers, and cross-profile exact equality', () => {
  const a = leaf([1, 2]);
  const same = new PolynomialHash(), alternate = new PolynomialHash(MERSENNE_127, 257);
  assert.equal(a.profile, h.profile);
  assert.deepEqual(h.profile, new HashProfile());
  assert.ok(Object.isFrozen(a.profile));
  assert.throws(() => { a.profile.base = 257n; }, TypeError);
  assert.throws(() => { a.profile = alternate.profile; }, TypeError);
  assert.throws(() => { h.profile = alternate.profile; }, TypeError);
  assert.equal(ropeConcat(null, a, same), a);
  assert.equal(ropeEqual(a, leaf([1, 2], alternate)), true);
  assert.equal(ropeLcp(a, leaf([1, 2], alternate)), 2n);
});

test('profile: mismatch is rejected before no-op and cached-fingerprint shortcuts', () => {
  const a = leaf([1, 2]), other = new PolynomialHash(MERSENNE_127, 257);
  const b = leaf([3], other);
  const calls = [() => new Internal(a, b, h), () => new RepeatNode(a, 2, other),
    () => ropeConcat(null, a, other), () => ropeConcat(a, null, other),
    () => ropeSplit(a, 0, other), () => ropeSplit(a, a.len, other),
    () => ropeRepeat(a, 0, other), () => ropeRepeat(a, 1, other),
    () => ropeSubstrHash(a, 0, 0, other), () => ropeSubstrHash(a, 0, a.len, other),
    () => ropeInsert(a, 0, b, h), () => ropeDelete(a, 0, 0, other), () => validateRope(a, other)];
  for (const call of calls) assert.throws(call, /Incompatible hash profile/);
});

test('profile: stream cannot combine its retained/prefix state under a replacement profile', () => {
  const stream = new SlidingWindow(2, 1);
  stream.appendBytes(Uint8Array.of(1, 2, 3, 4));
  const pos = stream.pos;
  stream.h = new PolynomialHash(MERSENNE_127, 257);
  for (const operation of [() => stream.currentHash(), () => stream.finalHash(), () => stream.appendBytes(new Uint8Array()), () => stream.appendCopy(1, 0)]) {
    assert.throws(operation, /Incompatible hash profile/);
    assert.equal(stream.pos, pos);
  }
  stream.h = new PolynomialHash();
  assert.equal(stream.currentHash(), h.hash(Uint8Array.of(1, 2, 3, 4)));
  stream.rWindow = leaf([1], new PolynomialHash(MERSENNE_127, 257));
  assert.throws(() => stream.currentHash(), /Incompatible hash profile/);
});
