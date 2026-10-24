# API examples and Python name map

Import from `hashrope` (ES modules) or `require('hashrope')` (CommonJS). Choose one module form per node graph; the two builds have separate constructor identities. Type declarations ship for both. Integer arguments use `Integer = bigint | number`; number values must be safe integers. Return values representing hashes, lengths, weights, counts and positions are `bigint`. Structural height is a `number` bounded in practice by the runtime stack.

## Polynomial arithmetic

```ts
import { PolynomialHash, MERSENNE_61, MERSENNE_127,
  mersenneMod, mersenneMul, phi } from 'hashrope';
const h = new PolynomialHash(MERSENNE_61, 131n);
const a = Uint8Array.of(0, 1, 255);
const ha = h.hash(a);
h.power(100000000000000000000n);
h.hashConcat(ha, 3n, ha);             // hash(a || a)
h.hashRepeat(ha, 3n, 20n);            // hash(a repeated twenty times)
h.hashOverlap(ha, 3n, 8n, h.hash(a.subarray(0, 2)));
mersenneMod(123456n); mersenneMul(123n, 456n);
phi(20n, h.power(3n));
const wider = new PolynomialHash(MERSENNE_127, 131n);
```

`prime`, `base` and `profile` are read-only. `profile` is a frozen `HashProfile` with prime/base fields; nodes share their originating hasher's profile. Fingerprint operations reject incompatible profiles even at no-op boundaries; exact comparison allows different profiles across its roots. The constructor checks `2 <= base < prime`; it does not prove primality or the Mersenne shape. `power` supports negative exponents when an inverse exists. Helpers preserve Python's bit-folding behavior, which is not the ordinary mathematical remainder for arbitrary non-Mersenne masks or negative operands. Nonpositive helper moduli are rejected to prevent reference nontermination (with `phi`'s zero/one short circuits preserved). See compatibility notes before relying on malformed arithmetic inputs.

## Nodes and operations

```ts
import { Leaf, Internal, RepeatNode, ropeFromBytes, ropeConcat, ropeSplit,
  ropeRepeat, ropeSubstrHash, ropeLen, ropeHash, ropeHeight, ropeToBytes,
  validateRope, PolynomialHash } from 'hashrope';
import type { Node, NonEmptyNode } from 'hashrope';
const h = new PolynomialHash();
const a = new Leaf(Uint8Array.of(97, 98), h);
const b = new Leaf(Uint8Array.of(99), h);
const direct = new Internal(a, b, h); // direct metadata join; does not rebalance
const repeated = new RepeatNode(a, 100n, h);
const root: Node = ropeConcat(repeated, b, h); // balanced persistent join
const [left, right] = ropeSplit(root, 37n, h);
const repeatedAgain = ropeRepeat(a, 1n << 100n, h);
const fingerprint = ropeSubstrHash(repeatedAgain, (1n << 99n) + 1n, 17n, h);
console.log(ropeLen(root), ropeHash(root), ropeHeight(root), fingerprint);
validateRope(root); // structural/profile invariants; does not rehash
const result = ropeToBytes(root); // fresh Uint8Array, default 64 MiB budget
const resultWithBudget = ropeToBytes(root, 1024n);
const empty = ropeFromBytes(new Uint8Array(), h); // null
```

All nonempty nodes expose read-only `len`, `weight`, `hashVal`, `profile`. `Leaf.data` returns a fresh byte copy; `Internal.left/right` and `RepeatNode.child/reps` expose immutable structure. Leaf weight is one independent of byte length; repeat weight counts expanded child leaves. Nodes are frozen at runtime. `Node` is a TypeScript union with `null` (not a runtime constructor); `NonEmptyNode` excludes null.

`ropeSplit` clamps endpoint positions. `ropeRepeat(node,0,h)` returns null, count one returns the same root, and counts at least two create a repeat node. Repeating null at least twice throws. Invalid negative counts throw. Valid substring ranges satisfy `0 <= start`, `0 <= length`, `start + length <= ropeLen(node)`. Full coverage returns stored metadata. The inherited malformed-range behavior is retained for parity but is shape-dependent and is not a clipping contract.

`validateRope(node, optionalHasher)` throws `RopeValidationError` (whose `.name` is `AssertionError`) on structural violations. The optional hasher must match the root profile. Child profiles must match their parent; mismatch throws `RangeError`. Validation does not recalculate fingerprints.

## Persistent edits and equality (JS additions)

```ts
import { ropeInsert, ropeDelete, contentEquals, structuralEquals, ropeEqual, ropeLcp, ropeRangeEqual, ropeRangeLcp,
  Leaf, PolynomialHash } from 'hashrope';
const h = new PolynomialHash();
const a = new Leaf(Uint8Array.of(1, 2, 3), h);
const inserted = ropeInsert(a, 1n, new Leaf(Uint8Array.of(9), h), h);
const deleted = ropeDelete(inserted, 1n, 1n, h);
console.assert(contentEquals(a, deleted));
console.log(structuralEquals(a, deleted)); // structure may differ after edits
const x = new Leaf(Uint8Array.of(0, 131), h);
const y = new Leaf(Uint8Array.of(1, 0), h);
console.assert(x.hashVal === y.hashVal && !contentEquals(x, y));
console.assert(!ropeEqual(x, y) && ropeLcp(x, y) === 0n);
console.assert(ropeRangeEqual(a, 1n, deleted, 1n, 2n));
console.assert(ropeRangeLcp(a, 0n, deleted, 0n, 3n) === 3n);
```

`ropeInsert(node, pos, value, h)` clamps its insertion offset using `ropeSplit`. `ropeDelete(node, start, length, h)` clamps start to the rope, rejects negative length, and removes up to the available length. Both preserve prior roots. Python has no corresponding exported functions: these are compositions of the ported split and concat APIs.

`structuralEquals` compares concrete node type, structure, bytes and stored metadata. `contentEquals` compares byte sequences regardless of shape or originating hash profile; it never treats matching fingerprints as proof or mismatching fingerprints as a safe rejection across profiles. It is an alias for `ropeEqual`. `ropeLcp` returns an exact BigInt byte count. `ropeRangeEqual(a,aStart,b,bStart,length)` and `ropeRangeLcp` validate both ranges strictly, even for empty or identical ranges, and do not split or rehash. Shared ranges, verified-range memoization and finite periodic witnesses avoid expanding common compressed cases; arbitrary unshared content can require a full byte scan. Query storage and BigInt arithmetic are not free.

## SlidingWindow

```ts
import { SlidingWindow } from 'hashrope';
const sw = new SlidingWindow(8n, 4n); // dMax, mMax, optional prime, optional base
sw.appendBytes(new TextEncoder().encode('abcdefgh'));
sw.appendCopy(3n, 10n); // repeating 'fgh' with overlap
console.log(sw.pos, sw.windowLen, sw.currentHash(), sw.finalHash());
const immutableSnapshot = sw.rWindow;
```

Mutable Python slots map to `h`, `dMax`, `mMax`, `W`, `hPrefix`, `lPrefix`, `rWindow`, `pos`. `windowLen` is a read-only getter. Direct state assignment is possible and does not repair invariants. Valid input operations keep `pos-lPrefix === windowLen`; inherited negative-copy quirks can break that relation and should not be used for ordinary decoded streams. Null/non-byte append inputs are rejected; an empty Uint8Array is a no-op. `currentHash`/`finalHash` cover all ingested bytes, including evictions.

## Constants and representation

`VERSION` is the npm implementation version; `PYTHON_VERSION` names the pinned compatibility target (`0.3.3` and `0.2.2`, respectively). `ALPHA_NUM=2n`, `ALPHA_DEN=7n` describe the balance condition. Python snake_case is mapped to camelCase; no accidental snake_case aliases are claimed. JSON cannot encode bigint directly: use an explicit decimal-string schema when serializing metadata.
