import test from 'node:test';
import assert from 'node:assert/strict';
import { SlidingWindow } from '../dist/esm/sliding.js';
import { PolynomialHash, MERSENNE_61, MERSENNE_127 } from '../dist/esm/polynomial.js';
import { ropeToBytes, validateRope, ropeHeight } from '../dist/esm/nodes.js';

const bytes = (text) => new TextEncoder().encode(text);

test('defaults expose exact empty mutable streaming state', () => {
  const s = new SlidingWindow();
  assert.equal(s.dMax, 32768n);
  assert.equal(s.mMax, 258n);
  assert.equal(s.W, 33026n);
  assert.equal(s.h.prime, MERSENNE_61);
  assert.equal(s.h.base, 131n);
  assert.equal(s.hPrefix, 0n);
  assert.equal(s.lPrefix, 0n);
  assert.equal(s.rWindow, null);
  assert.equal(s.pos, 0n);
  assert.equal(s.windowLen, 0n);
  assert.equal(s.currentHash(), 0n);
  assert.equal(s.finalHash(), 0n);
});

test('constructor validates capacities in Python order with explicit messages', () => {
  assert.throws(() => new SlidingWindow(0, 0), {name:'RangeError', message:'d_max must be >= 1, got 0'});
  assert.throws(() => new SlidingWindow(1, -2), {name:'RangeError', message:'m_max must be >= 1, got -2'});
  assert.throws(() => new SlidingWindow(1, 1, 31, 1), {name:'RangeError', message:'Base must be in [2, p-1], got 1'});
});

test('custom profiles and logical capacities retain BigInt precision', () => {
  const huge = (1n << 105n) + 9n;
  const s = new SlidingWindow(huge, 7, MERSENNE_127, 257);
  assert.equal(s.dMax, huge);
  assert.equal(s.mMax, 7n);
  assert.equal(s.W, huge + 7n);
  assert.equal(s.h.prime, MERSENNE_127);
  assert.equal(s.h.base, 257n);
  assert.throws(() => new SlidingWindow(Number.MAX_SAFE_INTEGER + 1, 1), RangeError);
  assert.throws(() => new SlidingWindow(1.5, 1), RangeError);
  assert.throws(() => new SlidingWindow('1', 1), TypeError);
});

test('current hash composes assignable state and finalHash dispatches to override', () => {
  const s = new SlidingWindow(4, 2);
  s.hPrefix = 123n;
  s.lPrefix = 99n;
  s.pos = 888n;
  assert.equal(s.currentHash(), 123n);
  class Override extends SlidingWindow { currentHash() { return 987n; } }
  assert.equal(new Override().finalHash(), 987n);
  class CountingHash extends PolynomialHash {
    calls = [];
    hashConcat(a, n, b) { this.calls.push([a,n,b]); return super.hashConcat(a,n,b); }
  }
  const h = new CountingHash();
  s.h = h;
  assert.equal(s.currentHash(), 123n);
  assert.deepEqual(h.calls, [[123n,0n,0n]]);
});

function equalWindow(s, text) {
  assert.deepEqual(ropeToBytes(s.rWindow), bytes(text));
  assert.equal(s.windowLen, BigInt(bytes(text).length));
  validateRope(s.rWindow);
}

test('literal append owns bytes and persistent snapshots survive later appends', () => {
  const s = new SlidingWindow(20, 10);
  const input = bytes('hello ');
  s.appendBytes(input);
  const snapshot = s.rWindow;
  input.fill(0);
  s.appendBytes(bytes('world'));
  equalWindow(s, 'hello world');
  assert.deepEqual(ropeToBytes(snapshot), bytes('hello '));
  assert.equal(s.pos, 11n);
  assert.equal(s.currentHash(), s.h.hash(bytes('hello world')));
  assert.equal(s.finalHash(), s.currentHash());
  s.appendBytes(bytes('!'));
  equalWindow(s, 'hello world!');
});

test('literal eviction happens strictly above capacity and composes exact prefix', () => {
  const s = new SlidingWindow(4, 2);
  s.appendBytes(bytes('abcdef'));
  equalWindow(s, 'abcdef');
  assert.equal(s.lPrefix, 0n);
  s.appendBytes(bytes('g'));
  equalWindow(s, 'defg');
  assert.equal(s.lPrefix, 3n);
  assert.equal(s.hPrefix, s.h.hash(bytes('abc')));
  assert.equal(s.currentHash(), s.h.hash(bytes('abcdefg')));
  s.appendBytes(bytes('hijk'));
  equalWindow(s, 'hijk');
  assert.equal(s.lPrefix, 7n);
  assert.equal(s.hPrefix, s.h.hash(bytes('abcdefg')));
  assert.equal(s.currentHash(), s.h.hash(bytes('abcdefghijk')));
});

test('a single oversized literal evicts without limiting its input size', () => {
  const s = new SlidingWindow(4, 2);
  s.appendBytes(bytes('0123456789'));
  equalWindow(s, '6789');
  assert.equal(s.pos, 10n);
  assert.equal(s.lPrefix, 6n);
  assert.equal(s.hPrefix, s.h.hash(bytes('012345')));
  assert.equal(s.currentHash(), s.h.hash(bytes('0123456789')));
});

test('public capacity assignments do not silently recalculate other fields', () => {
  const s = new SlidingWindow(2, 2);
  s.appendBytes(bytes('abcd'));
  s.dMax = 1n;
  s.mMax = 100n;
  assert.equal(s.W, 4n);
  s.W = 2n;
  const snapshot = s.rWindow;
  s.appendBytes(new Uint8Array());
  assert.equal(s.rWindow, snapshot);
  assert.equal(s.windowLen, 4n);
  s.appendBytes(bytes('e'));
  equalWindow(s, 'e');
  assert.equal(s.pos, 5n);
  assert.equal(s.lPrefix, 4n);
  assert.equal(s.currentHash(), s.h.hash(bytes('abcde')));
});

test('literal API rejects foreign types and keeps empty Uint8Array a no-op', () => {
  const s = new SlidingWindow();
  for (const invalid of [null, undefined, '', [], new Int8Array(), new DataView(new ArrayBuffer(1))]) {
    assert.throws(() => s.appendBytes(invalid), TypeError);
  }
  s.appendBytes(new Uint8Array());
  assert.equal(s.rWindow, null);
  assert.equal(s.pos, 0n);
});


function unchangedState(s, action, expectedError) {
  const before = [s.rWindow, s.pos, s.hPrefix, s.lPrefix];
  assert.throws(action, expectedError);
  assert.deepEqual([s.rWindow, s.pos, s.hPrefix, s.lPrefix], before);
}

test('copy supports nonoverlap, overlap, partial final periods and old snapshots', () => {
  const s = new SlidingWindow(100, 20);
  s.appendBytes(bytes('abcde'));
  const old = s.rWindow;
  s.appendCopy(3, 2);
  equalWindow(s, 'abcdecd');
  s.appendCopy(2, 7);
  equalWindow(s, 'abcdecdcdcdcdc');
  assert.equal(s.currentHash(), s.h.hash(bytes('abcdecdcdcdcdc')));
  assert.deepEqual(ropeToBytes(old), bytes('abcde'));
});

test('copy bounds check decoded length before retained length even for zero output', () => {
  const s = new SlidingWindow(4, 2);
  s.appendBytes(bytes('0123456789'));
  unchangedState(s, () => s.appendCopy(11, 0), {
    name: 'RangeError', message: 'Invalid copy: offset 11 exceeds decoded length 10',
  });
  unchangedState(s, () => s.appendCopy(5, 0), {
    name: 'RangeError', message: 'Copy offset 5 exceeds window size 4. Increase d_max (currently 4) to at least 5.',
  });
});

test('copy distance and length are not capped separately by dMax and mMax', () => {
  const s = new SlidingWindow(4, 2);
  s.appendBytes(bytes('abcdef'));
  s.appendCopy(6, 1);
  equalWindow(s, 'defa');
  s.appendCopy(1, 10);
  equalWindow(s, 'aaaa');
  assert.equal(s.pos, 17n);
  assert.equal(s.currentHash(), s.h.hash(bytes('abcdef' + 'a'.repeat(11))));
});

test('zero distance follows Python division failure without mutating stream state', () => {
  const empty = new SlidingWindow();
  empty.appendCopy(0, 0);
  assert.equal(empty.pos, 0n);
  unchangedState(empty, () => empty.appendCopy(0, 1), {
    name: 'RangeError', message: 'integer division or modulo by zero',
  });
  const s = new SlidingWindow();
  s.appendBytes(bytes('abc'));
  const before = s.rWindow;
  s.appendCopy(0, 0);
  assert.equal(s.rWindow, before);
  unchangedState(s, () => s.appendCopy(0, 5), RangeError);
});

test('signed-copy quirks preserve published Python counters and bytes', () => {
  for (const [offset, length, pos] of [[1,-1,2],[-1,2,5],[-1,-2,1],[-2,-1,2]]) {
    const s = new SlidingWindow(4,2);
    s.appendBytes(bytes('abc'));
    s.appendCopy(offset, length);
    equalWindow(s, 'abc');
    assert.equal(s.pos, BigInt(pos));
    assert.equal(s.hPrefix, 0n);
    assert.equal(s.lPrefix, 0n);
    assert.equal(s.currentHash(), 1694847n);
  }
});

test('copy integer conversions reject lossy Number inputs before mutation', () => {
  const s = new SlidingWindow();
  s.appendBytes(bytes('abc'));
  for (const invalid of [1.5, NaN, Infinity, Number.MAX_SAFE_INTEGER + 1]) {
    unchangedState(s, () => s.appendCopy(invalid, 1), RangeError);
    unchangedState(s, () => s.appendCopy(1, invalid), RangeError);
  }
  unchangedState(s, () => s.appendCopy('1', 1), TypeError);
  unchangedState(s, () => s.appendCopy(1, null), TypeError);
});

test('compressed copy beyond 2^100 agrees with independent published-wheel constants', () => {
  const s = new SlidingWindow(4,2);
  s.appendBytes(bytes('ab'));
  s.appendCopy(2n, (1n << 100n) + 3n);
  assert.equal(s.pos, 1267650600228229401496703205381n);
  assert.equal(s.lPrefix, 1267650600228229401496703205377n);
  assert.equal(s.hPrefix, 885451250975102620n);
  equalWindow(s, 'baba');
  assert.equal(s.rWindow.hashVal, 224255854n);
  assert.equal(s.currentHash(), 189876086087698542n);
  assert.equal(ropeHeight(s.rWindow), 2);
});

test('copy indexes UTF-8 bytes and can copy within a multibyte character', () => {
  const s = new SlidingWindow();
  const encoded = bytes('A\u00e9');
  s.appendBytes(encoded);
  s.appendCopy(1, 2);
  assert.deepEqual(ropeToBytes(s.rWindow), Uint8Array.of(65,195,169,169,169));
  assert.equal(s.pos, 5n);
});

function referenceHash(data, prime, base) {
  let result = 0n;
  for (const byte of data) result = (result * base + BigInt(byte) + 1n) % prime;
  return result;
}

test('seeded streaming agrees with independent expanded-byte and modular models', () => {
  for (const [prime, base] of [[MERSENNE_61,131n],[MERSENNE_127,257n],[31n,3n],[127n,17n]]) {
    const s = new SlidingWindow(17,9,prime,base);
    let seed = 0x534c4944;
    const random = (limit) => {
      seed ^= seed << 13; seed ^= seed >>> 17; seed ^= seed << 5;
      return (seed >>> 0) % limit;
    };
    const full = [];
    let prefixLength = 0;
    for (let step = 0; step < 240; step++) {
      const snapshot = s.rWindow;
      const oldWindow = ropeToBytes(snapshot);
      if (s.windowLen === 0n || random(3) === 0) {
        const value = Uint8Array.from({length: 1 + random(10)}, () => random(256));
        s.appendBytes(value);
        full.push(...value);
      } else {
        const offset = 1 + random(Number(s.windowLen));
        const length = random(20);
        s.appendCopy(offset, length);
        for (let i=0; i<length; i++) full.push(full[full.length-offset]);
      }
      if (full.length - prefixLength > 26) prefixLength = full.length - 17;
      assert.equal(s.pos, BigInt(full.length));
      assert.equal(s.lPrefix, BigInt(prefixLength));
      assert.equal(s.windowLen, BigInt(full.length-prefixLength));
      assert.deepEqual(ropeToBytes(s.rWindow), Uint8Array.from(full.slice(prefixLength)));
      assert.deepEqual(ropeToBytes(snapshot), oldWindow);
      assert.equal(s.hPrefix, referenceHash(full.slice(0,prefixLength),prime,base));
      assert.equal(s.currentHash(), referenceHash(full,prime,base));
      validateRope(s.rWindow);
    }
  }
});

