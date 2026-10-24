import { PolynomialHash, phi, mersenneMul, toBigInt, type Integer, type HashProfile, sameProfile } from './polynomial.js';

const ownedLeafBytes = new WeakMap<Leaf, Uint8Array>();

function leafBytes(leaf: Leaf): Uint8Array {
  const data = ownedLeafBytes.get(leaf);
  if (data === undefined) throw new TypeError('Expected an initialized Leaf');
  return data;
}

function assertBytes(data: unknown): asserts data is Uint8Array {
  if (!(data instanceof Uint8Array)) throw new TypeError('data must be a Uint8Array');
}

/** A persistent, immutable rope leaf. Byte arrays are defensively owned. */
export class Leaf {
  readonly len: bigint;
  readonly weight: bigint = 1n;
  readonly hashVal: bigint;
  readonly profile: HashProfile;

  constructor(data: Uint8Array, h: PolynomialHash) {
    assertBytes(data);
    if (data.length === 0) throw new RangeError('Leaf cannot be empty');
    const owned = new Uint8Array(data);
    ownedLeafBytes.set(this, owned);
    this.profile = h.profile;
    this.len = BigInt(owned.length);
    this.hashVal = h.hash(new Uint8Array(owned));
    Object.freeze(this);
  }

  get data(): Uint8Array { return new Uint8Array(leafBytes(this)); }
}

/** Direct immutable join. Use ropeConcat when rebalancing is required. */
export class Internal {
  readonly left: NonEmptyNode;
  readonly right: NonEmptyNode;
  readonly len: bigint;
  readonly weight: bigint;
  readonly hashVal: bigint;
  readonly profile: HashProfile;
  constructor(left: NonEmptyNode, right: NonEmptyNode, h: PolynomialHash) {
    assertNonEmptyNode(left);
    assertNonEmptyNode(right);
    requireProfile(left, h);
    requireProfile(right, h);
    this.profile = h.profile;
    this.left = left;
    this.right = right;
    this.len = left.len + right.len;
    this.weight = left.weight + right.weight;
    this.hashVal = h.hashConcat(left.hashVal, right.len, right.hashVal);
    Object.freeze(this);
  }
}

/** Shares one child while retaining an exact, unexpanded repetition count. */
export class RepeatNode {
  readonly child: NonEmptyNode;
  readonly reps: bigint;
  readonly len: bigint;
  readonly weight: bigint;
  readonly hashVal: bigint;
  readonly profile: HashProfile;
  constructor(child: NonEmptyNode, reps: Integer, h: PolynomialHash) {
    const count = toBigInt(reps, 'reps');
    if (count < 2n) throw new RangeError(`RepeatNode reps must be >= 2, got ${count}`);
    assertNonEmptyNode(child);
    requireProfile(child, h);
    this.profile = h.profile;
    this.child = child;
    this.reps = count;
    this.len = child.len * count;
    this.weight = child.weight * count;
    const power = h.power(child.len);
    this.hashVal = mersenneMul(child.hashVal, phi(count, power, h.prime), h.prime);
    Object.freeze(this);
  }
}

function assertNonEmptyNode(node: unknown): asserts node is NonEmptyNode {
  if (!(node instanceof Leaf || node instanceof Internal || node instanceof RepeatNode)) {
    throw new TypeError('Expected a Leaf, Internal, or RepeatNode');
  }
}

export type NonEmptyNode = Leaf | Internal | RepeatNode;
export type Node = NonEmptyNode | null;

/** Internal guard shared by operations that consume stored fingerprints. */
export function requireProfile(node: Node, h: PolynomialHash): void {
  if (node === null) return;
  assertNonEmptyNode(node);
  if (!sameProfile(node.profile, h.profile)) throw new RangeError('Incompatible hash profile');
}

/** Internal byte verifier; never exposes the owned backing arrays. */
export function leafRangeLcp(a: Leaf, aStart: number, b: Leaf, bStart: number, length: number): number {
  const x = leafBytes(a), y = leafBytes(b);
  // A short byte scout keeps early mismatches cheap, with no view allocation.
  let i = 0;
  const scout = length < 64 ? length : 16;
  for (; i < scout; i++) if (x[aStart + i] !== y[bStart + i]) return i;
  if (i === length) return length;
  // DataView supports unaligned ranges and browser runtimes. These are views
  // of private owned buffers, never copies or public mutable buffer exposure.
  const xv = new DataView(x.buffer, x.byteOffset, x.byteLength);
  const yv = new DataView(y.buffer, y.byteOffset, y.byteLength);
  for (; i + 8 <= length; i += 8) {
    if (xv.getUint32(aStart + i, true) !== yv.getUint32(bStart + i, true)
      || xv.getUint32(aStart + i + 4, true) !== yv.getUint32(bStart + i + 4, true)) {
      // Word equality is an exact byte test; locate the first differing byte
      // only inside the unequal word pair, without a hash or encoding step.
      for (let j = 0; j < 8; j++) if (x[aStart + i + j] !== y[bStart + i + j]) return i + j;
    }
  }
  for (; i < length; i++) if (x[aStart + i] !== y[bStart + i]) return i;
  return length;
}

export function ropeLen(node: Node): bigint { return node === null ? 0n : node.len; }
export function ropeHash(node: Node): bigint { return node === null ? 0n : node.hashVal; }
export function ropeHeight(node: Node): number {
  if (node === null || node instanceof Leaf) return 0;
  if (node instanceof RepeatNode) return 1 + ropeHeight(node.child);
  assertNonEmptyNode(node);
  return 1 + Math.max(ropeHeight(node.left), ropeHeight(node.right));
}
export function ropeFromBytes(data: Uint8Array, h: PolynomialHash): Node {
  assertBytes(data);
  return data.length === 0 ? null : new Leaf(data, h);
}

/** A failed structural invariant; independent of runtime assertion flags. */
export class RopeValidationError extends Error {
  constructor(message: string) { super(message); this.name = 'AssertionError'; }
}
/** Materialize fresh bytes. The default allocation limit is 64 MiB. */
export function ropeToBytes(node: Node, maxBytes: Integer = 64 * 1024 * 1024): Uint8Array {
  const limit = toBigInt(maxBytes, 'maxBytes');
  if (limit < 0n || limit > BigInt(Number.MAX_SAFE_INTEGER)) {
    throw new RangeError('maxBytes must be between zero and Number.MAX_SAFE_INTEGER');
  }
  if (node !== null) assertNonEmptyNode(node);
  const length = ropeLen(node);
  if (length > limit) throw new RangeError(`Rope length ${length} exceeds materialization limit ${limit}`);
  // The bound above proves this conversion is exact; native allocation may still fail.
  const output = new Uint8Array(Number(length));
  if (node !== null) writeBytes(node, output, 0);
  return output;
}

function writeBytes(node: NonEmptyNode, output: Uint8Array, offset: number): number {
  if (node instanceof Leaf) {
    const data = leafBytes(node);
    output.set(data, offset);
    return offset + data.length;
  }
  if (node instanceof Internal) {
    return writeBytes(node.right, output, writeBytes(node.left, output, offset));
  }
  // Reuse the first materialized period; work and storage still scale with output bytes.
  const end = offset + Number(node.len);
  let written = writeBytes(node.child, output, offset) - offset;
  while (offset + written < end) {
    const take = Math.min(written, end - offset - written);
    output.copyWithin(offset + written, offset, offset + take);
    written += take;
  }
  return end;
}
/** Checks structural and profile consistency; does not rehash stored fingerprints. */
export function validateRope(node: Node, h?: PolynomialHash): void {
  if (node === null) return;
  assertNonEmptyNode(node);
  if (h !== undefined) requireProfile(node, h);
  if (node instanceof Leaf) {
    if (node.weight !== 1n) throw new RopeValidationError(`Leaf weight must be 1, got ${node.weight}`);
    if (node.len !== BigInt(leafBytes(node).length)) throw new RopeValidationError('Leaf len mismatch');
    return;
  }
  if (node instanceof Internal) {
    if (!sameProfile(node.profile, node.left.profile) || !sameProfile(node.profile, node.right.profile)) throw new RangeError('Incompatible hash profile');
    if (node.len !== node.left.len + node.right.len) throw new RopeValidationError('Internal len mismatch');
    if (node.weight !== node.left.weight + node.right.weight) throw new RopeValidationError('Internal weight mismatch');
    if (node.weight > 2n) {
      if (2n * node.weight > 7n * node.left.weight) {
        throw new RopeValidationError(`Left child too light: ${node.left.weight}/${node.weight}`);
      }
      if (7n * node.left.weight > 5n * node.weight) {
        throw new RopeValidationError(`Left child too heavy: ${node.left.weight}/${node.weight}`);
      }
    }
    validateRope(node.left);
    validateRope(node.right);
    return;
  }
  if (!sameProfile(node.profile, node.child.profile)) throw new RangeError('Incompatible hash profile');
  if (node.len !== node.child.len * node.reps) throw new RopeValidationError('RepeatNode len mismatch');
  if (node.weight !== node.child.weight * node.reps) throw new RopeValidationError('RepeatNode weight mismatch');
  if (node.reps < 2n) throw new RopeValidationError(`RepeatNode reps must be >= 2, got ${node.reps}`);
  validateRope(node.child);
}
/** Python-style structural value equality. Equal fingerprints alone never suffice. */
export function structuralEquals(a: Node, b: Node): boolean {
  const pending: Array<readonly [Node, Node]> = [[a, b]];
  while (pending.length !== 0) {
    const [left, right] = pending.pop()!;
    if (left !== null) assertNonEmptyNode(left);
    if (right !== null) assertNonEmptyNode(right);
    if (left === right) continue;
    if (left === null || right === null) return false;
    if (Object.getPrototypeOf(left) !== Object.getPrototypeOf(right)) return false;
    if (left.len !== right.len || left.weight !== right.weight || left.hashVal !== right.hashVal) return false;
    if (left instanceof Leaf && right instanceof Leaf) {
      const x = leafBytes(left);
      const y = leafBytes(right);
      if (x.length !== y.length) return false;
      for (let i = 0; i < x.length; i += 1) if (x[i] !== y[i]) return false;
    } else if (left instanceof Internal && right instanceof Internal) {
      pending.push([left.right, right.right], [left.left, right.left]);
    } else if (left instanceof RepeatNode && right instanceof RepeatNode) {
      if (left.reps !== right.reps) return false;
      pending.push([left.child, right.child]);
    } else {
      return false;
    }
  }
  return true;
}
