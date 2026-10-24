import { type Integer, toBigInt } from './polynomial.js';
import { Leaf, Internal, RepeatNode, type Node, type NonEmptyNode, ropeLen, leafRangeLcp } from './nodes.js';

function checkNode(node: Node): void {
  if (node !== null && !(node instanceof Leaf || node instanceof Internal || node instanceof RepeatNode)) {
    throw new TypeError('Expected a Leaf, Internal, RepeatNode, or null');
  }
}

/** Exact content equality, independent of tree shape, profile, and fingerprints. */
export function ropeEqual(left: Node, right: Node): boolean {
  checkNode(left); checkNode(right);
  if (left === right) return true;
  return ropeLen(left) === ropeLen(right) && ropeLcp(left, right) === ropeLen(left);
}

/** Exact byte LCP. Shared ranges, verified subproblems, and periodic proofs can
 * skip compressed content. Unshared bytes may require a full scan; this is not
 * a universal logarithmic bound. Query work includes traversal, memoization and
 * arbitrary-precision integer arithmetic. Neither rope nor hash cache changes.
 */
export function ropeLcp(left: Node, right: Node): bigint {
  checkNode(left); checkNode(right);
  const length = ropeLen(left) < ropeLen(right) ? ropeLen(left) : ropeLen(right);
  if (left === right || length === 0n) return length;
  return rangeLcp(left!, 0n, right!, 0n, length);
}

function checkRange(node: Node, start: bigint, length: bigint): void {
  checkNode(node);
  if (start < 0n || length < 0n || start > ropeLen(node) || length > ropeLen(node) - start) {
    throw new RangeError('Comparison range is out of bounds');
  }
}

/** Exact comparison of equally sized ranges; bounds are checked even for
 * empty and identical ranges. Hash profiles may differ across the two roots.
 */
export function ropeRangeEqual(left: Node, leftStart: Integer, right: Node, rightStart: Integer, length: Integer): boolean {
  const size = toBigInt(length, 'length');
  return ropeRangeLcp(left, leftStart, right, rightStart, size) === size;
}

/** Exact range LCP relative to the starts, without splitting or rehashing. */
export function ropeRangeLcp(left: Node, leftStart: Integer, right: Node, rightStart: Integer, length: Integer): bigint {
  const a = toBigInt(leftStart, 'leftStart'), b = toBigInt(rightStart, 'rightStart'), size = toBigInt(length, 'length');
  checkRange(left, a, size); checkRange(right, b, size);
  if (size === 0n || (left === right && a === b)) return size;
  return rangeLcp(left!, a, right!, b, size);
}

function normalize(node: NonEmptyNode, start: bigint, size: bigint): readonly [NonEmptyNode, bigint] {
  for (;;) {
    if (node instanceof Internal) {
      if (start >= node.left.len) { start -= node.left.len; node = node.right; }
      else if (start + size <= node.left.len) node = node.left;
      else return [node, start];
    } else if (node instanceof RepeatNode) {
      start %= node.child.len;
      if (start + size <= node.child.len) node = node.child;
      else return [node, start];
    } else return [node, start];
  }
}

function gcd(a: bigint, b: bigint): bigint {
  while (b !== 0n) [a, b] = [b, a % b];
  return a;
}

type Range = { kind: 'range'; a: NonEmptyNode; aStart: bigint; b: NonEmptyNode; bStart: bigint; size: bigint };
type Task = Range | { kind: 'skip'; size: bigint } | { kind: 'commit'; key: string };

function rangeLcp(left: NonEmptyNode, leftStart: bigint, right: NonEmptyNode, rightStart: bigint, length: bigint): bigint {
  const pending: Task[] = [{ kind: 'range', a: left, aStart: leftStart, b: right, bStart: rightStart, size: length }];
  const verified = new Set<string>();
  const ids = new WeakMap<NonEmptyNode, number>();
  let nextId = 0;
  const id = (node: NonEmptyNode): number => {
    let value = ids.get(node);
    if (value === undefined) { value = nextId++; ids.set(node, value); }
    return value;
  };
  let matched = 0n;
  while (pending.length !== 0) {
    const task = pending.pop()!;
    if (task.kind === 'skip') { matched += task.size; continue; }
    if (task.kind === 'commit') { verified.add(task.key); continue; }
    const size = task.size;
    const [a, aStart] = normalize(task.a, task.aStart, size);
    const [b, bStart] = normalize(task.b, task.bStart, size);
    if (a === b && aStart === bStart) { matched += size; continue; }
    if (a instanceof Leaf && b instanceof Leaf) {
      // Both normalized ranges fit real Uint8Arrays; narrowing is now exact.
      const common = BigInt(leafRangeLcp(a, Number(aStart), b, Number(bStart), Number(size)));
      matched += common;
      if (common !== size) return matched;
      continue;
    }
    const key = `${id(a)}:${aStart}:${id(b)}:${bStart}:${size}`;
    if (verified.has(key)) { matched += size; continue; }
    // Commit only after every subrange has been verified. A mismatch returns
    // before this marker, so hashes or unproved descendants cannot enter memo.
    pending.push({ kind: 'commit', key });
    if (a instanceof RepeatNode && b instanceof RepeatNode) {
      const p = a.child.len, q = b.child.len;
      const witness = p + q - gcd(p, q);
      if (witness < size) {
        // Fine-Wilf applies to the two rotated periodic streams. Equality of
        // this entire witness proves their remaining common-length suffix.
        pending.push({ kind: 'skip', size: size - witness });
        pending.push({ kind: 'range', a, aStart, b, bStart, size: witness });
        continue;
      }
    }
    let cut = size;
    // Internal boundaries take precedence: splitting at every period first
    // would expand a huge repeated prefix next to one edited suffix.
    if (a instanceof Internal) cut = a.left.len - aStart < cut ? a.left.len - aStart : cut;
    if (b instanceof Internal) cut = b.left.len - bStart < cut ? b.left.len - bStart : cut;
    if (cut === size) {
      if (a instanceof RepeatNode) cut = a.child.len - aStart < cut ? a.child.len - aStart : cut;
      if (b instanceof RepeatNode) cut = b.child.len - bStart < cut ? b.child.len - bStart : cut;
    }
    if (cut <= 0n || cut >= size) throw new Error('Exact comparison failed to make progress');
    pending.push({ kind: 'range', a, aStart: aStart + cut, b, bStart: bStart + cut, size: size - cut });
    pending.push({ kind: 'range', a, aStart, b, bStart, size: cut });
  }
  return matched;
}
