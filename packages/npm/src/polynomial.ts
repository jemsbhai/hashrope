/** Exact integers; Number inputs must be safe integers. */
export type Integer = bigint | number;

export const MERSENNE_61 = (1n << 61n) - 1n;
export const MERSENNE_127 = (1n << 127n) - 1n;

/** Convert only demonstrably lossless Number inputs. */
export function toBigInt(value: Integer, name = 'value'): bigint {
  if (typeof value === 'bigint') return value;
  if (typeof value !== 'number') throw new TypeError(`${name} must be a bigint or safe integer number`);
  if (!Number.isSafeInteger(value)) throw new RangeError(`${name} must be a safe integer number; use bigint for larger values`);
  return BigInt(value);
}

/** Python's floor division and signed remainder, without Number conversion. */
export function floorDivmod(a: bigint, b: bigint): readonly [bigint, bigint] {
  if (b === 0n) throw new RangeError('integer division or modulo by zero');
  let q = a / b;
  let r = a % b;
  if (r !== 0n && (r < 0n) !== (b < 0n)) { q -= 1n; r += b; }
  return [q, r];
}

// The string already exists in addressable memory; its length is an exact Number.
function bitLength(value: bigint): number {
  return value === 0n ? 0 : (value < 0n ? -value : value).toString(2).length;
}

/** Python 0.2.2 bit-folding behavior, including positive non-Mersenne masks. */
export function mersenneMod(value: Integer, modulus: Integer = MERSENNE_61): bigint {
  let a = toBigInt(value, 'a');
  const p = toBigInt(modulus, 'p');
  if (p <= 0n) throw new RangeError('p must be positive');
  const k = p === MERSENNE_61 ? 61n : p === MERSENNE_127 ? 127n : BigInt(bitLength(p));
  const threshold = 1n << (2n * k);
  while (a >= threshold) a = (a >> k) + (a & p);
  a = (a >> k) + (a & p);
  a = (a >> k) + (a & p);
  return a >= p ? a - p : a;
}

/** Exact arbitrary-precision multiplication followed by reference folding. */
export function mersenneMul(a: Integer, b: Integer, p: Integer = MERSENNE_61): bigint {
  return mersenneMod(toBigInt(a, 'a') * toBigInt(b, 'b'), p);
}


/** Mathematical positive-modulus remainder (distinct from reference folding). */
function modulo(a: bigint, p: bigint): bigint {
  const r = a % p;
  return r < 0n ? r + p : r;
}

function inverse(value: bigint, p: bigint): bigint {
  let a = modulo(value, p), b = p;
  let x = 1n, y = 0n;
  while (b !== 0n) {
    const q = a / b;
    [a, b] = [b, a - q * b];
    [x, y] = [y, x - q * y];
  }
  if (a !== 1n) throw new RangeError('base is not invertible for the given modulus');
  return modulo(x, p);
}

/** Immutable identity of the arithmetic used to construct stored fingerprints. */
export class HashProfile {
  readonly prime: bigint;
  readonly base: bigint;
  constructor(prime: Integer = MERSENNE_61, base: Integer = 131n) {
    this.prime = toBigInt(prime, 'prime');
    this.base = toBigInt(base, 'base');
    if (this.base < 2n || this.base >= this.prime) {
      throw new RangeError(`Base must be in [2, p-1], got ${this.base}`);
    }
    Object.freeze(this);
  }
}
export function sameProfile(a: HashProfile, b: HashProfile): boolean {
  return a.prime === b.prime && a.base === b.base;
}

/** Polynomial fingerprint, not a cryptographic hash or content-equality proof. */
export class PolynomialHash {
  readonly #profile: HashProfile;
  readonly #prime: bigint;
  readonly #base: bigint;
  readonly #powerCache = new Map<bigint, bigint>();

  constructor(prime: Integer = MERSENNE_61, base: Integer = 131n) {
    this.#prime = toBigInt(prime, 'prime');
    this.#base = toBigInt(base, 'base');
    if (this.#base < 2n || this.#base >= this.#prime) {
      throw new RangeError(`Base must be in [2, p-1], got ${this.#base}`);
    }
    this.#profile = new HashProfile(this.#prime, this.#base);
    this.#powerCache.set(0n, 1n);
    this.#powerCache.set(1n, this.#base);
  }

  get profile(): HashProfile { return this.#profile; }
  get prime(): bigint { return this.#prime; }
  get base(): bigint { return this.#base; }

  /** Cached modular power; negative exponents require an invertible base. */
  power(exponent: Integer): bigint {
    const n = toBigInt(exponent, 'n');
    const cached = this.#powerCache.get(n);
    if (cached !== undefined) return cached;
    let remaining = n < 0n ? -n : n;
    let factor = n < 0n ? inverse(this.#base, this.#prime) : this.#base;
    let result = 1n;
    while (remaining !== 0n) {
      if ((remaining & 1n) !== 0n) result = result * factor % this.#prime;
      remaining >>= 1n;
      factor = factor * factor % this.#prime;
    }
    this.#powerCache.set(n, result);
    return result;
  }

  /** Hash bytes with unsigned byte + 1 coefficients. No implicit encoding. */
  hash(data: Uint8Array): bigint {
    if (!(data instanceof Uint8Array)) throw new TypeError('data must be a Uint8Array');
    let result = 0n;
    // Read typed-array storage, not a caller-defined iterable representation.
    for (const byte of Uint8Array.prototype.values.call(data)) {
      result = mersenneMod(result * this.#base + BigInt(byte) + 1n, this.#prime);
    }
    return result;
  }
  /** H(A || B) from fingerprints and the exact byte length of B. */
  hashConcat(hA: Integer, lenB: Integer, hB: Integer): bigint {
    const power = this.power(lenB);
    return mersenneMod(toBigInt(hA, 'hA') * power + toBigInt(hB, 'hB'), this.#prime);
  }

  /** H(S repeated q times), without expanding the repeated bytes. */
  hashRepeat(hS: Integer, d: Integer, q: Integer): bigint {
    const power = this.power(d);
    return mersenneMul(hS, phi(q, power, this.#prime), this.#prime);
  }

  /** Compose repeated periods and a caller-supplied partial-prefix hash. */
  hashOverlap(hPattern: Integer, d: Integer, length: Integer, hPrefix: Integer): bigint {
    const period = toBigInt(d, 'd');
    const [q, r] = floorDivmod(toBigInt(length, 'length'), period);
    const periodPower = this.power(period);
    const remainderPower = this.power(r);
    const repeated = mersenneMul(hPattern, phi(q, periodPower, this.#prime), this.#prime);
    return mersenneMod(repeated * remainderPower + toBigInt(hPrefix, 'hPrefix'), this.#prime);
  }
}

/** Inverse-free geometric accumulator; signed counts retain Python's bit scan. */
export function phi(count: Integer, alphaValue: Integer, modulus: Integer = MERSENNE_61): bigint {
  const q = toBigInt(count, 'q');
  const alpha = toBigInt(alphaValue, 'alpha');
  const p = toBigInt(modulus, 'p');
  if (q === 0n) return 0n;
  if (q === 1n) return 1n;
  if (p <= 0n) throw new RangeError('p must be positive');
  const bits = bitLength(q);
  let result = 1n;
  let alphaPower = alpha;
  for (let i = bits - 2; i >= 0; i--) {
    result = mersenneMul(result, modulo(1n + alphaPower, p), p);
    alphaPower = mersenneMul(alphaPower, alphaPower, p);
    if (((q >> BigInt(i)) & 1n) !== 0n) {
      result = modulo(mersenneMul(result, alpha, p) + 1n, p);
      alphaPower = mersenneMul(alphaPower, alpha, p);
    }
  }
  return result;
}
