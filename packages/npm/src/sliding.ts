import { PolynomialHash, MERSENNE_61, toBigInt, floorDivmod, type Integer, type HashProfile, sameProfile } from './polynomial.js';
import { Leaf, ropeLen, ropeHash, requireProfile, type Node } from './nodes.js';
import { ropeConcat, ropeSplit, ropeRepeat } from './rope.js';

/**
 * Mutable streaming state; offsets and lengths count bytes, never characters.
 * Public fields map Python slots. Assignments do not repair other fields.
 */
export class SlidingWindow {
  readonly #profile: HashProfile;
  h: PolynomialHash;
  dMax: bigint;
  mMax: bigint;
  W: bigint;
  hPrefix = 0n;
  lPrefix = 0n;
  rWindow: Node = null;
  pos = 0n;

  constructor(dMax: Integer = 32768n, mMax: Integer = 258n, prime: Integer = MERSENNE_61, base: Integer = 131n) {
    this.dMax = toBigInt(dMax, 'dMax');
    if (this.dMax < 1n) throw new RangeError(`d_max must be >= 1, got ${this.dMax}`);
    this.mMax = toBigInt(mMax, 'mMax');
    if (this.mMax < 1n) throw new RangeError(`m_max must be >= 1, got ${this.mMax}`);
    this.h = new PolynomialHash(prime, base);
    this.#profile = this.h.profile;
    this.W = this.dMax + this.mMax;
  }

  #checkProfile(): void {
    if (!sameProfile(this.#profile, this.h.profile)) throw new RangeError('Incompatible hash profile: create a new SlidingWindow');
    requireProfile(this.rWindow, this.h);
  }

  get windowLen(): bigint { return ropeLen(this.rWindow); }

  /** Compose the evicted prefix with retained bytes. Does not finalize state. */
  currentHash(): bigint {
    this.#checkProfile();
    const length = this.windowLen;
    const hash = ropeHash(this.rWindow);
    return this.h.hashConcat(this.hPrefix, length, hash);
  }

  /** Alias that deliberately dispatches to overridden currentHash methods. */
  finalHash(): bigint { return this.currentHash(); }

  /** Append owned literal bytes. Empty arrays are a no-op. */
  appendBytes(data: Uint8Array): void {
    this.#checkProfile();
    if (!(data instanceof Uint8Array)) throw new TypeError('data must be a Uint8Array');
    if (data.length === 0) return;
    const leaf = new Leaf(data, this.h);
    this.rWindow = ropeConcat(this.rWindow, leaf, this.h);
    this.pos += BigInt(data.length);
    this.#evict();
  }

  /**
   * Append a back-reference, retaining compressed repetitions when it overlaps.
   * Signed-input quirks are preserved from Python 0.2.2; use nonnegative lengths
   * and positive offsets for valid decoded-stream operations.
   */
  appendCopy(offsetValue: Integer, lengthValue: Integer): void {
    this.#checkProfile();
    const offset = toBigInt(offsetValue, 'offset');
    const length = toBigInt(lengthValue, 'length');
    const windowLength = this.windowLen;
    if (offset > this.pos) {
      throw new RangeError(`Invalid copy: offset ${offset} exceeds decoded length ${this.pos}`);
    }
    if (offset > windowLength) {
      throw new RangeError(`Copy offset ${offset} exceeds window size ${windowLength}. Increase d_max (currently ${this.dMax}) to at least ${offset}.`);
    }
    const start = windowLength - offset;
    const [, tail] = ropeSplit(this.rWindow, start, this.h);
    if (offset >= length) {
      const [source] = ropeSplit(tail, length, this.h);
      this.rWindow = ropeConcat(this.rWindow, source, this.h);
    } else {
      const [pattern] = ropeSplit(tail, offset, this.h);
      const [repetitions, remainder] = floorDivmod(length, offset);
      let repeated: Node = repetitions >= 1n ? ropeRepeat(pattern, repetitions, this.h) : null;
      if (remainder > 0n) {
        const [partial] = ropeSplit(pattern, remainder, this.h);
        repeated = ropeConcat(repeated, partial, this.h);
      }
      this.rWindow = ropeConcat(this.rWindow, repeated, this.h);
    }
    this.pos += length;
    this.#evict();
  }

  #evict(): void {
    const windowLength = ropeLen(this.rWindow);
    // Match the reference threshold: exactly W bytes are retained unchanged.
    if (windowLength <= this.W) return;
    const excess = windowLength - this.dMax;
    const [old, keep] = ropeSplit(this.rWindow, excess, this.h);
    const oldLength = ropeLen(old);
    const oldHash = ropeHash(old);
    this.hPrefix = this.h.hashConcat(this.hPrefix, oldLength, oldHash);
    this.lPrefix += excess;
    this.rWindow = keep;
  }
}
