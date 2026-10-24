import test from 'node:test';
import assert from 'node:assert/strict';
import { PolynomialHash, MERSENNE_127 } from '../dist/esm/polynomial.js';
import { Leaf, Internal, RepeatNode, ropeToBytes, ropeLen, validateRope } from '../dist/esm/nodes.js';
import { contentEquals, ropeSplit, ropeSubstrHash } from '../dist/esm/rope.js';
const h = new PolynomialHash();

test('review: typed-array hashing reads stored bytes even with a custom iterator', () => {
  const data = Uint8Array.of(1, 2);
  data[Symbol.iterator] = function* () { yield 255; yield 255; };
  assert.equal(h.hash(data), 265n);
  assert.equal(new Leaf(data, h).hashVal, h.hash(data));
  class CustomIterableBytes extends Uint8Array {
    *[Symbol.iterator]() { yield 10000; }
  }
  const subclass = new CustomIterableBytes([1, 2]);
  assert.equal(h.hash(subclass), 265n);
});

test('review: exact content comparison ignores tree shape and different coherent hash profiles', () => {
  const alternate = new PolynomialHash(MERSENNE_127, 257n);
  let state = 42023;
  const rnd = () => { state = (Math.imul(state, 1664525) + 1013904223) >>> 0; return state; };
  const shape = (data, profile = (rnd() & 1) ? h : alternate) => {
    if (data.length < 2 || rnd() % 3 === 0) return new Leaf(data, profile);
    const at = 1 + rnd() % (data.length - 1);
    return new Internal(shape(data.slice(0, at), profile), shape(data.slice(at), profile), profile);
  };
  for (let i = 0; i < 100; i++) {
    const data = Uint8Array.from({length:1 + rnd() % 80}, () => rnd() & 255);
    const left = shape(data), right = shape(data);
    assert.equal(contentEquals(left, right), true);
    const changed = data.slice(); changed[rnd() % data.length] ^= 1;
    assert.equal(contentEquals(left, shape(changed)), false);
    assert.deepEqual(ropeToBytes(left), data);
  }
});

test('review: 150-bit repetition offsets split and query without materializing the prefix', () => {
  const count = (1n << 150n) + 1n;
  const pattern = new Leaf(Uint8Array.of(17,45,255), h);
  const repeated = new RepeatNode(pattern, count, h);
  const position = repeated.len - 4n;
  const [left, right] = ropeSplit(repeated, position, h);
  assert.equal(ropeLen(left), position);
  assert.equal(ropeLen(right), 4n);
  const tail = Uint8Array.of(255,17,45,255);
  assert.deepEqual(ropeToBytes(right), tail);
  assert.equal(ropeSubstrHash(repeated, position, 4n, h), h.hash(tail));
  validateRope(left); validateRope(right);
});

test('review: detached returned arrays and retained callback buffers do not change leaves', () => {
  const backing = Uint8Array.of(99,1,2,99);
  let observed;
  class CapturingHash extends PolynomialHash {
    hash(data) { observed = data; return super.hash(data); }
  }
  const leaf = new Leaf(backing.subarray(1,3), new CapturingHash());
  backing.fill(0); observed.fill(0);
  const returned = leaf.data;
  structuredClone(returned.buffer, {transfer:[returned.buffer]});
  assert.equal(returned.byteLength, 0);
  assert.deepEqual(leaf.data, Uint8Array.of(1,2));
  assert.equal(leaf.hashVal, 265n);
  assert.deepEqual(ropeToBytes(new RepeatNode(leaf,3n,h)), Uint8Array.of(1,2,1,2,1,2));
});

test('review: malformed leaf slicing normalizes huge BigInts before narrowing', () => {
  const leaf = new Leaf(Uint8Array.of(10,20,30,40),h);
  const huge = 1n << 200n;
  assert.equal(ropeSubstrHash(leaf, -huge, huge + 2n, h), h.hash(Uint8Array.of(10,20)));
  assert.equal(ropeSubstrHash(leaf, huge, -huge + 2n, h), 0n);
  assert.equal(ropeSubstrHash(leaf, -1n, 1n, h), 0n); // Python slice [-1:0] is empty.
});
