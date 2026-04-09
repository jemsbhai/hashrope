# Changelog

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
