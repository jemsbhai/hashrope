import test from 'node:test';
import assert from 'node:assert/strict';
import * as p from '../dist/esm/polynomial.js';

test('integer adapter accepts exact integers and rejects lossy Numbers and duck types', () => {
  assert.equal(p.toBigInt(123), 123n);
  assert.equal(p.toBigInt(-(2n ** 200n)), -(2n ** 200n));
  for (const bad of [1.5, NaN, Infinity, Number.MAX_SAFE_INTEGER + 1]) assert.throws(() => p.toBigInt(bad), RangeError);
  for (const bad of ['1', true, null, undefined]) assert.throws(() => p.toBigInt(bad), TypeError);
});

test('constants and Mersenne arithmetic match independent mathematical remainders', () => {
  assert.equal(p.MERSENNE_61, (1n << 61n) - 1n);
  assert.equal(p.MERSENNE_127, (1n << 127n) - 1n);
  for (const prime of [31n, p.MERSENNE_61, p.MERSENNE_127]) {
    for (const a of [0n, 1n, prime, prime + 1n, (1n << 401n) + 251n, 10n ** 200n]) {
      assert.equal(p.mersenneMod(a, prime), a % prime);
      assert.equal(p.mersenneMul(a, a + 23n, prime), a * (a + 23n) % prime);
    }
  }
});

test('reference signed folding and non-Mersenne masks remain observable', () => {
  // Values generated from the SHA256-pinned published Python 0.2.2 wheel.
  for (const [a, prime, expected] of [
    [-100n, 31n, 24n], [-1n, 31n, 30n],
    [-(10n ** 50n), 31n, -97656250000000000000000000000000000000000000000n],
    [1n, 30n, 0n], [-1n, 30n, 28n], [10n ** 50n, 30n, 14n],
  ]) assert.equal(p.mersenneMod(a, prime), expected);
  for (const prime of [0n, -1n]) {
    assert.throws(() => p.mersenneMod(4n, prime), RangeError);
    assert.throws(() => p.mersenneMul(2n, 3n, prime), RangeError);
  }
});

test('division uses Python floor quotient and divisor-signed remainder', () => {
  assert.deepEqual(p.floorDivmod(8n, 3n), [2n, 2n]);
  assert.deepEqual(p.floorDivmod(-8n, 3n), [-3n, 1n]);
  assert.deepEqual(p.floorDivmod(8n, -3n), [-3n, -1n]);
  assert.deepEqual(p.floorDivmod(-8n, -3n), [2n, -2n]);
  assert.throws(() => p.floorDivmod(1n, 0n), RangeError);
});


test('PolynomialHash exposes immutable profile values and validates base interval', () => {
  const h = new p.PolynomialHash();
  assert.equal(h.prime, p.MERSENNE_61);
  assert.equal(h.base, 131n);
  assert.throws(() => { h.prime = 31n; }, TypeError);
  assert.throws(() => { h.base = 3n; }, TypeError);
  for (const base of [1n, -1n, 31n]) assert.throws(() => new p.PolynomialHash(31n, base), {name:'RangeError', message:`Base must be in [2, p-1], got ${base}`});
  assert.equal(new p.PolynomialHash(30, 3).prime, 30n);
  assert.throws(() => new p.PolynomialHash(2n ** 61n, Number.MAX_SAFE_INTEGER + 1), RangeError);
});

test('byte hashing includes all unsigned bytes and collision regression under both moduli', () => {
  const h = new p.PolynomialHash();
  assert.equal(h.hash(new Uint8Array()), 0n);
  assert.equal(h.hash(Uint8Array.of(0)), 1n);
  assert.equal(h.hash(Uint8Array.of(255)), 256n);
  assert.equal(h.hash(Uint8Array.from({length:256}, (_,i) => i)), 773597766054883837n);
  const a = Uint8Array.of(0, 131), b = Uint8Array.of(1, 0);
  assert.notDeepEqual(a, b);
  for (const prime of [p.MERSENNE_61, p.MERSENNE_127]) {
    const sameBase = new p.PolynomialHash(prime, 131n);
    assert.equal(sameBase.hash(a), 263n);
    assert.equal(sameBase.hash(b), 263n);
  }
  for (const bad of [[1,2], 'abc', null, new Int8Array(2), new DataView(new ArrayBuffer(2))]) assert.throws(() => h.hash(bad), TypeError);
});

test('power supports caching, negative modular inverses and huge exact exponents', () => {
  const h = new p.PolynomialHash();
  for (const [n, expected] of [[-2n,1348221476280531735n],[-1n,1372944692508917009n],[0n,1n],[1n,131n],[17n,1688172759669340850n],[(1n<<100n)+3n,1347858949099060711n]]) {
    assert.equal(h.power(n), expected);
    assert.equal(h.power(n), expected);
  }
  assert.throws(() => new p.PolynomialHash(15, 3).power(-1), {name:'RangeError', message:'base is not invertible for the given modulus'});
  assert.throws(() => h.power(Number.MAX_SAFE_INTEGER + 1), RangeError);
});


test('phi agrees with an independently summed geometric series and pinned signed cases', () => {
  for (const prime of [31n, p.MERSENNE_61, p.MERSENNE_127]) {
    for (const alpha of [0n, 1n, 3n, prime - 1n]) {
      let sum = 0n, power = 1n;
      for (let q = 0n; q <= 128n; q++) {
        assert.equal(p.phi(q, alpha, prime), sum);
        sum = (sum + power) % prime;
        power = power * alpha % prime;
      }
    }
  }
  for (const [q, expected] of [[-9n,30n],[-4n,9n],[-3n,13n],[-2n,4n],[-1n,1n],[(1n<<100n)+3n,21n]]) assert.equal(p.phi(q, 3n, 31n), expected);
  const huge = (1n << 200n) + 937n;
  assert.equal(p.phi(huge, 1n), huge % p.MERSENNE_61);
  assert.equal(p.phi(0n, 3n, 0n), 0n);
  assert.equal(p.phi(1n, 3n, 0n), 1n);
  assert.throws(() => p.phi(2n, 3n, 0n), RangeError);
  assert.throws(() => p.phi(Number.MAX_SAFE_INTEGER + 1, 3), RangeError);
});

test('composition and repetition match independently expanded byte sequences', () => {
  for (const prime of [31n, p.MERSENNE_61, p.MERSENNE_127]) {
    const h = new p.PolynomialHash(prime, 3n);
    const a = Uint8Array.of(0, 255, 7), b = Uint8Array.of(42, 1);
    assert.equal(h.hashConcat(h.hash(a), b.length, h.hash(b)), h.hash(Uint8Array.of(...a,...b)));
    for (let q = 0; q < 40; q++) {
      assert.equal(h.hashRepeat(h.hash(a), a.length, q), h.hash(Uint8Array.from(Array(q).fill([...a]).flat())));
    }
    assert.throws(() => h.hashRepeat(1n, 2n, Number.MAX_SAFE_INTEGER + 1), RangeError);
  }
});

test('overlap uses floor division, exact prefix and pinned signed results', () => {
  const h = new p.PolynomialHash();
  const pattern = Uint8Array.of(97, 98, 99);
  for (let length = 0; length < 40; length++) {
    const prefix = pattern.slice(0, length % pattern.length);
    const expected = Uint8Array.from({length}, (_, i) => pattern[i % 3]);
    assert.equal(h.hashOverlap(h.hash(pattern), pattern.length, length, h.hash(prefix)), h.hash(expected));
  }
  for (const [d,l,expected] of [[3n,8n,1620339286111n],[3n,-8n,27806642488836253n],[-3n,8n,1368391352313939690n],[-3n,-8n,258791227259697484n]]) assert.equal(h.hashOverlap(42n,d,l,7n),expected);
  assert.throws(() => h.hashOverlap(0n,0n,0n,0n), RangeError);
});

test('composition dispatches power exactly like the Python public methods', () => {
  class Observed extends p.PolynomialHash {
    calls = [];
    power(n) { this.calls.push(n); return super.power(n); }
  }
  const h = new Observed();
  h.hashConcat(0n, 7n, 0n);
  h.hashRepeat(0n, 3n, 0n);
  h.hashOverlap(42n, 3n, 8n, 7n);
  assert.deepEqual(h.calls, [7n, 3n, 3n, 2n]);
  assert.throws(() => new p.PolynomialHash(15n,3n).hashRepeat(0n,-1n,0n), RangeError);
});
