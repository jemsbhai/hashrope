# Changelog

## 0.2.2 — 2026-06-10

**Performance fix: `rope_substr_hash` now realizes the Theorem 9 `O(k · log w)` bound.**

`_hash_range` was descending into subtrees fully covered by the queried range and re-hashing every byte through `PolynomialHash.hash`, instead of reading the hash metadata stored at every node. The Invariant I3 values were written by all constructors but never consulted by range queries, making `rope_substr_hash` effectively O(length).

### What was wrong

The canonical-decomposition recursion had no early exit for full coverage. A range query decomposes into O(log w) fully covered subtrees plus at most two partial leaves at the ends — but each fully covered subtree was recursed to its leaves and re-hashed byte-by-byte. On a 500 KB rope (4 KB leaves), a half-range query cost ~120 ms and an LCP binary search over prefix hashes cost ~8 s.

### What changed

- **`_hash_range` fast path** — a node fully covered by the query (`start == 0 and length == node.len`) returns its stored `hash_val` in O(1). Eight lines added (guard + comment); nothing else touched.

Measured on CPython 3.12 (Linux, x86-64): half-range query ~120 ms → ~1 ms; prefix-hash LCP ~8 s → ~28 ms, size-independent from 500 KB to 16 MB.

### Correctness

For a fully covered node, the stored `hash_val` equals the recomputed hash by Invariants I1/I3/I6, maintained by every constructor (structural induction). Verified against a byte-level oracle on 425+ randomized ranges, including spans across RepeatNode tail/full-copies/head decompositions.

### Impact on existing users

- **No API changes.** All public signatures are identical.
- **Hash values are unchanged** for every input — bit-for-bit cross-language consistency with the Rust crate is preserved.
- **Return values are unchanged** for all inputs, including out-of-contract ones (the guard requires exact full coverage).
- The only observable difference is speed.

### Testing

146 tests pass (141 pre-existing with no regressions + 5 new in `tests/test_substr_hash_theorem9.py`):
- Byte-oracle property tests over random ranges and boundary cases (empty, single byte, full range, leaf-boundary straddles)
- RepeatNode span coverage (tail / full copies / head)
- Deterministic `CountingHash` complexity guard asserting a query byte-rehashes at most the two partial end-leaves — fails on 0.2.1, passes on 0.2.2
- Exact-LCP regression via prefix-hash binary search

## 0.2.1 — 2026-04-09

**Critical bugfix: BB[2/7] balance violations in `_join` and `_rebalance`.**

Two bugs in the balance-maintenance code caused `_join` (and therefore `rope_concat`) to produce trees that violate the BB[2/7] weight-balance invariant. Trees constructed with `rope_concat`, `rope_split`+`rope_concat`, or operations involving `RepeatNode` could contain Internal nodes with child-weight ratios outside the required [2/7, 5/7] range.

This is a port of the same fix applied to the Rust crate in v0.2.2.

### What was wrong

**Bug A — `_rebalance` could not handle `RepeatNode`.** When `_rebalance` needed to rotate a subtree rooted at a `RepeatNode`, the `isinstance(x, Internal)` check failed silently, falling through to return an unbalanced node. This affected any `rope_concat` or `rope_split`+`rope_concat` sequence involving `RepeatNode` (including all LZ77 back-reference patterns).

**Bug B — Rotations created unbalanced children.** `_rotate_left` and `_rotate_right` called `Internal(a, b, h)` on the inner child pair without checking whether the result satisfied BB[2/7]. In the `_join` context — where subtrees of arbitrary weight ratio are combined — the classical single-insertion rotation analysis does not apply, and the unchecked inner nodes could violate the invariant. This bug existed even in trees containing only Internal and Leaf nodes (no `RepeatNode` required).

### What changed

- **`_decompose(node, h)`** — new function that splits any non-Leaf node into two children. Internal nodes yield their existing children; RepeatNodes are split by halving the reps count. This eliminates the `isinstance` failure of Bug A.
- **`_balance(left, right, h)`** — replaces `_rebalance`. Adams-style checked rotations: single rotation is used only when both the inner and outer pairings satisfy BB[2/7]; otherwise double rotation decomposes further with recursive balancing. This eliminates the unchecked inner nodes of Bug B.
- **`_join` RepeatNode case** — now splits by reps count directly instead of calling `rope_split` at the byte midpoint, avoiding unnecessary intermediate node allocation.
- **Removed** `_rotate_left`, `_rotate_right`, and `_rebalance`.

### Termination and correctness

The `_balance` function terminates for all valid inputs. The ranking function M = max(weight(left), weight(right)) strictly decreases in every recursive call, with worst-case contraction ratio 223/245 ≈ 0.918, giving O(log W) recursion depth. Verified computationally (see Rust crate's `BALANCE_TERMINATION_PROOF.md`) and by an arithmetic sweep test in the Python test suite.

### Impact on existing users

- **No API changes.** All public function signatures are identical.
- **Hash values are unchanged.** The same input produces the same polynomial hash.
- **Tree shapes may differ.** The new balancing produces different (correctly balanced) tree structures, so `rope_height()` and internal structure may change for the same input sequence.
- **`validate_rope()` will no longer raise** on trees that previously triggered "Left child too heavy/light" errors.

### Testing

141 tests pass, including:
- All 50 pre-existing rope tests (no regressions)
- 44 polynomial hash and sliding window tests (no regressions)
- 22 bug-reproduction tests covering both bugs across all node-type combinations
- 25 comprehensive balance tests: direct `_balance` probes, exhaustive small-weight validation, termination verification, cross-language hash consistency, property-based weight sweeps, recursion depth safety

## 0.2.0

Initial release. Rope with polynomial hash metadata, RepeatNode, SlidingWindow.
