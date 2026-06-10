# Changelog

## 0.3.0 — 2026-04-10

*(Entry backfilled 2026-06-10; this release shipped without changelog notes.)*

**Feature release: lazy hash materialization (Definitions 19–20).**

Hash values can now be deferred at construction and computed on demand. A lazy `Arena` skips all polynomial-hash work while building ropes — leaves, internal nodes, and repeat nodes are written with `LAZY_SENTINEL` in place of a hash — making construction O(t·log t) instead of O(k·c_F·t) (Theorem 44). Hashes are materialized transparently the first time a query needs them and cached in place; subsequent reads are O(1). Consumers that never read hashes (CD-Radix sort, Corollary 8) pay zero hashing cost end-to-end.

### New API (additive only)

- **`Arena::new_lazy()`** — construct a lazy arena (Definition 19).
- **`Arena::is_lazy()`** — query the arena's mode.
- **`Arena::ensure_hash(id)`** — explicitly materialize a node's hash, children first, caching results (Definition 20). Idempotent; a no-op on eager arenas and on already-materialized nodes.
- **`NodeInner::hash_val()`** — accessor for the stored hash of any node variant.
- **`LAZY_SENTINEL: u64`** (= `u64::MAX`, re-exported from the crate root) — sentinel marking unmaterialized hashes; always greater than p = 2^61 − 1, so it can never collide with a valid hash.

### Semantics and caveats

- `substr_hash` on a lazy rope materializes on demand: only the fully covered canonical nodes a query touches are computed, then cached. Lazy and eager arenas produce **identical hash values** for the same content (verified by cross-mode equivalence tests).
- **Caveat:** `Arena::hash()` / `rope_hash()` take `&self` and cannot materialize. On a lazy rope whose root has not yet been materialized, they return `LAZY_SENTINEL`, not a hash. Call `ensure_hash(root)` first, or use `substr_hash(node, 0, len)`, which materializes as needed.
- Eager arenas (the default) are completely unaffected: every hash is computed at construction exactly as in 0.2.2.

### Impact on existing users

- **No breaking changes.** All existing signatures are identical; the additions above are new surface only.
- **Eager behavior is byte-identical to 0.2.2** — same hashes, same tree shapes, same costs.
- Minor version bump (0.2.2 → 0.3.0) signals the new construction mode rather than any incompatibility.

### Testing

107 tests in the published crate (74 in 0.2.2). The 33 additions cover lazy mode: sentinel state after construction, on-demand materialization via `substr_hash` (full, partial, and multi-range sweeps), RepeatNode materialization, split/rejoin on lazy ropes, concat of independently lazy ropes, idempotent `ensure_hash`, and lazy-vs-eager hash equivalence across all node-type combinations.

## 0.2.2 — 2026-04-09

**Critical bugfix: BB[2/7] balance violations in `join` and `rebalance`.**

Two bugs in the balance-maintenance code caused `join` (and therefore `concat`) to produce trees that violate the BB[2/7] weight-balance invariant. Trees constructed with `concat`, `split`+`concat`, or operations involving `RepeatNode` could contain Internal nodes with child-weight ratios outside the required [2/7, 5/7] range.

### What was wrong

**Bug A — `rebalance` could not handle `RepeatNode`.** When `rebalance` needed to rotate a subtree rooted at a `RepeatNode`, the `if let Internal { .. }` pattern match failed silently, falling through to an unchecked `make_internal` that created an unbalanced node. This affected any `concat` or `split`+`concat` sequence involving `RepeatNode` (including all LZ77 back-reference patterns in cdh-sort).

**Bug B — Rotations created unbalanced children.** `rotate_left` and `rotate_right` called `make_internal` on the inner child pair without checking whether the result satisfied BB[2/7]. In the `join` context — where subtrees of arbitrary weight ratio are combined — the classical single-insertion rotation analysis does not apply, and the unchecked inner nodes could violate the invariant. This bug existed even in trees containing only Internal and Leaf nodes (no `RepeatNode` required).

### What changed

The rotation-based `rebalance` function has been replaced by a new `balance` function using Adams-style checked rotations with a `decompose` helper that transparently handles all node types:

- **`decompose(node)`** splits any non-Leaf node into two children. Internal nodes yield their existing children; RepeatNodes are split by halving the reps count. This eliminates the pattern-match failure of Bug A.
- **`balance(L, R)`** performs single rotation only when both the inner and outer pairings satisfy BB[2/7]; otherwise it decomposes further via double rotation with recursive balancing. This eliminates the unchecked inner nodes of Bug B.
- **`join`'s RepeatNode case** now splits by reps count directly instead of calling `split_inner` at the byte midpoint, avoiding unnecessary intermediate node allocation.

The old `rotate_left`, `rotate_right`, and rotation-based `rebalance` functions have been removed.

### Termination and correctness

The new `balance` function terminates for all valid inputs. The ranking function M = max(weight(L), weight(R)) strictly decreases in every recursive call, with a worst-case contraction ratio of 223/245 ≈ 0.918, giving O(log W) recursion depth. This has been verified computationally over 9.2 million (weight, decomposition) combinations and proven analytically. A full proof is included in `BALANCE_TERMINATION_PROOF.md`.

### Impact on existing users

- **No API changes.** All public function signatures are identical.
- **Hash values are unchanged.** The same input produces the same polynomial hash.
- **Tree shapes may differ.** The new balancing produces different (correctly balanced) tree structures, so `node_count()`, `height()`, and internal node IDs may change for the same input sequence.
- **Benchmarks from v0.2.1 are invalid.** Unbalanced trees had pessimistic heights, making `substr_hash` and `byte_at` traversals slower than they should be. All benchmarks should be re-collected on v0.2.2.
- **`validate()` will no longer panic** on trees that previously triggered "Left child too heavy/light" errors.

### Testing

63 tests pass, including:
- All 42 pre-existing tests (no regressions)
- 9 bug-reproduction tests covering both bugs across all node-type combinations
- 8 diagnostic tests directly probing `rebalance` and `balance` with exhaustive weight sweeps
- 1 termination verification test (9.2M combinations, contraction ratio < 0.92)
- cdh-sort LZ77 back-reference pattern simulation
- distance-1 run-length chain (the exact production trigger)

## 0.2.1

Initial arena-based implementation.

## 0.2.0

Migration from `Arc`-based to arena-based node storage.
