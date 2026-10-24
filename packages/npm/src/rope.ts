import { ropeEqual } from './comparison.js';
import { PolynomialHash, type Integer, toBigInt, floorDivmod, phi, mersenneMod, mersenneMul } from './polynomial.js';
import { Leaf, Internal, RepeatNode, type Node, type NonEmptyNode, ropeLen, requireProfile } from './nodes.js';

export const ALPHA_NUM = 2n;
export const ALPHA_DEN = 7n;

function balanced(wl: bigint, wr: bigint): boolean {
  const total = wl + wr;
  return total <= 2n || (ALPHA_NUM * total <= ALPHA_DEN * wl && ALPHA_DEN * wl <= (ALPHA_DEN - ALPHA_NUM) * total);
}
function makeRepeat(child: NonEmptyNode, q: bigint, h: PolynomialHash): Node {
  if (q === 0n) return null;
  if (q === 1n) return child;
  return new RepeatNode(child, q, h);
}
function decompose(node: NonEmptyNode, h: PolynomialHash): readonly [NonEmptyNode, NonEmptyNode] {
  if (node instanceof Internal) return [node.left, node.right];
  if (node instanceof RepeatNode) {
    const half = node.reps / 2n;
    return [makeRepeat(node.child, half, h)!, makeRepeat(node.child, node.reps - half, h)!];
  }
  throw new TypeError('Cannot decompose Leaf in balance');
}
function balance(left: NonEmptyNode, right: NonEmptyNode, h: PolynomialHash): NonEmptyNode {
  const wl = left.weight, wr = right.weight;
  if (balanced(wl, wr)) return new Internal(left, right, h);
  if (ALPHA_NUM * (wl + wr) > ALPHA_DEN * wl) {
    const [rl, rr] = decompose(right, h);
    if (balanced(wl, rl.weight) && balanced(wl + rl.weight, rr.weight)) {
      return new Internal(new Internal(left, rl, h), rr, h);
    }
    const [rll, rlr] = decompose(rl, h);
    return balance(balance(left, rll, h), balance(rlr, rr, h), h);
  }
  const [ll, lr] = decompose(left, h);
  if (balanced(lr.weight, wr) && balanced(ll.weight, lr.weight + wr)) {
    return new Internal(ll, new Internal(lr, right, h), h);
  }
  const [lrl, lrr] = decompose(lr, h);
  return balance(balance(ll, lrl, h), balance(lrr, right, h), h);
}
function join(left: NonEmptyNode, right: NonEmptyNode, h: PolynomialHash): NonEmptyNode {
  if (balanced(left.weight, right.weight)) return new Internal(left, right, h);
  if (left.weight > right.weight) {
    if (left instanceof Internal) return balance(left.left, join(left.right, right, h), h);
    if (left instanceof RepeatNode) {
      const [ll, lr] = decompose(left, h);
      return join(ll, join(lr, right, h), h);
    }
  } else {
    if (right instanceof Internal) return balance(join(left, right.left, h), right.right, h);
    if (right instanceof RepeatNode) {
      const [rl, rr] = decompose(right, h);
      return join(join(left, rl, h), rr, h);
    }
  }
  return new Internal(left, right, h);
}
/** Persistent BB[2/7] concatenation. Empty sides preserve root identity. */
export function ropeConcat(left: Node, right: Node, h: PolynomialHash): Node {
  requireProfile(left, h);
  requireProfile(right, h);
  if (left === null) return right;
  if (right === null) return left;
  return join(left, right, h);
}
/** Split at an exact byte offset; out-of-range endpoints are clamped. */
export function ropeSplit(node: Node, position: Integer, h: PolynomialHash): readonly [Node, Node] {
  requireProfile(node, h);
  const pos = toBigInt(position, 'pos');
  if (node === null) return [null, null];
  if (pos <= 0n) return [null, node];
  if (pos >= node.len) return [node, null];
  return split(node, pos, h);
}
function split(node: NonEmptyNode, pos: bigint, h: PolynomialHash): readonly [Node, Node] {
  if (node instanceof Leaf) {
    // Only an interior index of an already allocated Uint8Array reaches this conversion.
    const data = node.data, index = Number(pos);
    return [new Leaf(data.subarray(0, index), h), new Leaf(data.subarray(index), h)];
  }
  if (node instanceof Internal) {
    const ll = node.left.len;
    if (pos === ll) return [node.left, node.right];
    if (pos < ll) {
      const [l1, l2] = split(node.left, pos, h);
      return [l1, ropeConcat(l2, node.right, h)];
    }
    const [r1, r2] = split(node.right, pos - ll, h);
    return [ropeConcat(node.left, r1, h), r2];
  }
  if (node instanceof RepeatNode) {
    const d = node.child.len, m = pos / d, r = pos % d;
    if (r === 0n) return [makeRepeat(node.child, m, h), makeRepeat(node.child, node.reps - m, h)];
    const [cl, cr] = split(node.child, r, h);
    return [ropeConcat(makeRepeat(node.child, m, h), cl, h), ropeConcat(cr, makeRepeat(node.child, node.reps - m - 1n, h), h)];
  }
  throw new TypeError('Unknown node type');
}
/** Compact repetition; zero is empty and one preserves the root reference. */
export function ropeRepeat(node: Node, repetitions: Integer, h: PolynomialHash): Node {
  requireProfile(node, h);
  const q = toBigInt(repetitions, 'q');
  if (q === 0n) return null;
  if (q === 1n) return node;
  return new RepeatNode(node as NonEmptyNode, q, h);
}
/** Hash a byte interval without constructing rope nodes. Use valid ranges for ordinary queries. */
export function ropeSubstrHash(node: Node, start: Integer, length: Integer, h: PolynomialHash): bigint {
  requireProfile(node, h);
  const s = toBigInt(start, 'start'), len = toBigInt(length, 'length');
  if (node === null || len === 0n) return 0n;
  return hashRange(node, s, len, h);
}
function sliceIndex(index: bigint, length: bigint): number {
  let i = index < 0n ? index + length : index;
  if (i < 0n) i = 0n;
  if (i > length) i = length;
  // length belongs to a real typed array; normalization precedes conversion.
  return Number(i);
}
function hashRange(node: NonEmptyNode, start: bigint, length: bigint, h: PolynomialHash): bigint {
  if (start === 0n && length === node.len) return node.hashVal;
  if (node instanceof Leaf) {
    const data = node.data;
    return h.hash(data.subarray(sliceIndex(start, node.len), sliceIndex(start + length, node.len)));
  }
  if (node instanceof Internal) {
    const ll = node.left.len;
    if (start + length <= ll) return hashRange(node.left, start, length, h);
    if (start >= ll) return hashRange(node.right, start - ll, length, h);
    const lLen = ll - start, rLen = length - lLen;
    return h.hashConcat(hashRange(node.left, start, lLen, h), rLen, hashRange(node.right, 0n, rLen, h));
  }
  if (node instanceof RepeatNode) {
    const d = node.child.len;
    const [first, startInCopy] = floorDivmod(start, d);
    const [last, endInCopy] = floorDivmod(start + length - 1n, d);
    if (first === last) return hashRange(node.child, startInCopy, length, h);
    const tailLen = d - startInCopy;
    const hTail = hashRange(node.child, startInCopy, tailLen, h);
    const fullCopies = last - first - 1n;
    let hFull = 0n;
    if (fullCopies > 0n) hFull = mersenneMul(node.child.hashVal, phi(fullCopies, h.power(d), h.prime), h.prime);
    const fullLen = fullCopies * d, headLen = endInCopy + 1n;
    const hHead = hashRange(node.child, 0n, headLen, h);
    return mersenneMod(h.power(fullLen + headLen) * hTail + h.power(headLen) * hFull + hHead, h.prime);
  }
  throw new TypeError('Unknown node type');
}
/** JavaScript convenience: insert a rope at a clamped byte offset. */
export function ropeInsert(node: Node, pos: Integer, value: Node, h: PolynomialHash): Node {
  requireProfile(node, h);
  requireProfile(value, h);
  const [left, right] = ropeSplit(node, pos, h);
  return ropeConcat(ropeConcat(left, value, h), right, h);
}
/** JavaScript convenience: remove up to length bytes from a clamped start. */
export function ropeDelete(node: Node, start: Integer, length: Integer, h: PolynomialHash): Node {
  requireProfile(node, h);
  const len = toBigInt(length, 'length');
  if (len < 0n) throw new RangeError('length must be nonnegative');
  const [left, rest] = ropeSplit(node, start, h);
  const [, right] = ropeSplit(rest, len, h);
  return ropeConcat(left, right, h);
}
/** Exact byte equality. Compatibility alias for ropeEqual; no fingerprints are used. */
export function contentEquals(a: Node, b: Node): boolean { return ropeEqual(a, b); }
