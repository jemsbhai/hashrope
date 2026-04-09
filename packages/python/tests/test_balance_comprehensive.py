"""
Comprehensive balance tests beyond bug reproduction (Step 4).

Covers:
- Direct _balance probes across all node-type combinations
- Exhaustive small-weight validation
- Termination verification (ported from Rust)
- Cross-language hash consistency
- Property-based weight sweeps
- Recursion depth safety
"""

import math
import pytest
from hashrope.polynomial_hash import PolynomialHash, MERSENNE_61
from hashrope.rope import (
    Leaf, Internal, RepeatNode,
    rope_concat, rope_split, rope_repeat,
    rope_len, rope_hash, rope_height,
    rope_to_bytes, rope_from_bytes,
    validate_rope,
    _balance, _decompose, _is_balanced_wt, _weight, _make_repeat,
    ALPHA_NUM, ALPHA_DEN,
)


@pytest.fixture
def h():
    return PolynomialHash(prime=MERSENNE_61, base=131)


# =========================================================================
# Direct _balance probes: every node-type combination
# =========================================================================

class TestBalanceDirectProbes:
    """Call _balance(L, R, h) directly with all node-type combos."""

    def _make_weighted_leaf_tree(self, weight, h):
        """Build a balanced tree of given weight using single-byte leaves."""
        if weight <= 0:
            return None
        node = Leaf(b"\x01", h)
        for i in range(1, weight):
            node = rope_concat(node, Leaf(bytes([i % 256]), h), h)
        return node

    def test_leaf_leaf(self, h):
        for wl in range(1, 6):
            for wr in range(1, 6):
                left = self._make_weighted_leaf_tree(wl, h)
                right = self._make_weighted_leaf_tree(wr, h)
                result = _balance(left, right, h)
                validate_rope(result, h)

    def test_internal_internal(self, h):
        for wl in range(2, 10):
            for wr in range(2, 10):
                left = self._make_weighted_leaf_tree(wl, h)
                right = self._make_weighted_leaf_tree(wr, h)
                result = _balance(left, right, h)
                validate_rope(result, h)

    def test_repeat_leaf(self, h):
        child = Leaf(b"x", h)
        leaf = Leaf(b"y", h)
        for q in range(2, 30):
            rep = RepeatNode(child, q, h)
            result = _balance(rep, leaf, h)
            validate_rope(result, h)
            result2 = _balance(leaf, rep, h)
            validate_rope(result2, h)

    def test_repeat_internal(self, h):
        child = Leaf(b"x", h)
        for q in range(2, 20):
            rep = RepeatNode(child, q, h)
            for wr in range(2, 10):
                right = self._make_weighted_leaf_tree(wr, h)
                result = _balance(rep, right, h)
                validate_rope(result, h)
                result2 = _balance(right, rep, h)
                validate_rope(result2, h)

    def test_repeat_repeat(self, h):
        child = Leaf(b"a", h)
        for q1 in range(2, 15):
            for q2 in range(2, 15):
                r1 = RepeatNode(child, q1, h)
                r2 = RepeatNode(child, q2, h)
                result = _balance(r1, r2, h)
                validate_rope(result, h)

    def test_nested_repeat(self, h):
        """RepeatNode whose child is also a RepeatNode."""
        inner = RepeatNode(Leaf(b"z", h), 3, h)
        outer = RepeatNode(inner, 4, h)  # weight = 12
        leaf = Leaf(b"w", h)
        result = _balance(outer, leaf, h)
        validate_rope(result, h)


# =========================================================================
# Exhaustive small-weight validation
# =========================================================================

class TestExhaustiveSmallWeight:
    """Build every reachable tree shape up to moderate weight via operations."""

    def test_all_concat_pairs_up_to_weight_20(self, h):
        """For every (wl, wr) with wl+wr <= 20, build and validate."""
        trees = {}
        for w in range(1, 21):
            trees[w] = Leaf(bytes([w % 256]), h)
            if w >= 2:
                # Also build via concat
                node = None
                for i in range(w):
                    node = rope_concat(node, Leaf(bytes([i % 256]), h), h)
                trees[w] = node

        for wl in range(1, 11):
            for wr in range(1, 11):
                result = rope_concat(trees[wl], trees[wr], h)
                validate_rope(result, h)

    def test_split_rejoin_all_positions_weight_15(self, h):
        """Weight-15 tree: split at every position, rejoin, validate."""
        node = None
        for i in range(15):
            node = rope_concat(node, Leaf(bytes([65 + i]), h), h)
        orig_hash = rope_hash(node)
        for pos in range(rope_len(node) + 1):
            left, right = rope_split(node, pos, h)
            rejoined = rope_concat(left, right, h)
            assert rope_hash(rejoined) == orig_hash
            validate_rope(rejoined, h)

    def test_repeat_then_split_all_positions(self, h):
        """RepeatNode: split at every byte position, validate both halves."""
        child = Leaf(b"abc", h)
        for q in range(2, 10):
            rep = RepeatNode(child, q, h)
            for pos in range(1, rep.len):
                left, right = rope_split(rep, pos, h)
                if left is not None:
                    validate_rope(left, h)
                if right is not None:
                    validate_rope(right, h)
                rejoined = rope_concat(left, right, h)
                assert rope_hash(rejoined) == rep.hash_val


# =========================================================================
# Termination verification (ported from Rust)
# =========================================================================

class TestBalanceTermination:
    """Verify that _balance terminates: M = max(wl, wr) strictly decreases."""

    def test_is_balanced_wt_consistency(self, h):
        """_is_balanced_wt must agree with the BB[2/7] definition."""
        for total in range(1, 50):
            for wl in range(0, total + 1):
                wr = total - wl
                result = _is_balanced_wt(wl, wr)
                if total <= 2:
                    assert result is True
                else:
                    expected = (ALPHA_NUM * total <= ALPHA_DEN * wl <=
                                (ALPHA_DEN - ALPHA_NUM) * total)
                    assert result == expected, f"wl={wl}, wr={wr}"

    def test_decompose_weight_conservation(self, h):
        """_decompose must preserve total weight."""
        child = Leaf(b"x", h)
        # Internal
        left = Leaf(b"a", h)
        right = Leaf(b"b", h)
        internal = Internal(left, right, h)
        dl, dr = _decompose(internal, h)
        assert _weight(dl) + _weight(dr) == internal.weight

        # RepeatNode
        for q in range(2, 20):
            rep = RepeatNode(child, q, h)
            dl, dr = _decompose(rep, h)
            assert _weight(dl) + _weight(dr) == rep.weight, f"q={q}"

    def test_decompose_halving_for_repeat(self, h):
        """RepeatNode decompose splits reps as floor(q/2), ceil(q/2)."""
        child = Leaf(b"x", h)
        for q in range(2, 50):
            rep = RepeatNode(child, q, h)
            dl, dr = _decompose(rep, h)
            assert _weight(dl) == q // 2
            assert _weight(dr) == q - q // 2

    def test_max_weight_decreases_arithmetic(self):
        """
        Pure arithmetic verification: for every unbalanced (wl, wr) pair
        with total <= 200, verify that all recursive _balance calls have
        strictly smaller max(wl, wr).

        This is the Python port of the Rust termination test.
        """
        violations = 0
        for total in range(3, 201):
            for wl in range(0, total + 1):
                wr = total - wl
                if _is_balanced_wt(wl, wr):
                    continue
                M = max(wl, wr)

                if ALPHA_NUM * total > ALPHA_DEN * wl:
                    # Left too light, decompose right (wr)
                    # Worst case: wr splits as (ceil(wr/2), floor(wr/2))
                    wrl = (wr + 1) // 2
                    wrr = wr // 2
                    # Single rotation case
                    if _is_balanced_wt(wl, wrl) and _is_balanced_wt(wl + wrl, wrr):
                        # No recursion, just Internal construction
                        pass
                    else:
                        # Double: decompose wrl -> (wrl//2, wrl - wrl//2)
                        wrll = wrl // 2
                        wrlr = wrl - wrl // 2
                        # Recursive calls: balance(wl, wrll), balance(wrlr, wrr), balance(result1, result2)
                        M1 = max(wl, wrll)
                        M2 = max(wrlr, wrr)
                        M3 = max(wl + wrll, wrlr + wrr)
                        if M1 >= M or M2 >= M or M3 >= M:
                            violations += 1
                else:
                    # Left too heavy, decompose left (wl)
                    wll = (wl + 1) // 2
                    wlr = wl // 2
                    if _is_balanced_wt(wlr, wr) and _is_balanced_wt(wll, wlr + wr):
                        pass
                    else:
                        wlrl = wlr // 2
                        wlrr = wlr - wlr // 2
                        M1 = max(wll, wlrl)
                        M2 = max(wlrr, wr)
                        M3 = max(wll + wlrl, wlrr + wr)
                        if M1 >= M or M2 >= M or M3 >= M:
                            violations += 1

        assert violations == 0, f"{violations} termination violations found"


# =========================================================================
# Cross-language hash consistency
# =========================================================================

class TestCrossLanguageConsistency:
    """Verify Python hashes match the Rust reference values."""

    def test_hello_world_hash(self, h):
        """hash('hello world') must match the Rust reference value."""
        expected = 430229793999670395
        assert h.hash(b"hello world") == expected

    def test_hello_world_rope_hash(self, h):
        """Rope-based hash of 'hello world' must match."""
        node = Leaf(b"hello world", h)
        assert rope_hash(node) == 430229793999670395

    def test_concat_hash_matches_direct(self, h):
        """concat(Leaf('hello '), Leaf('world')).hash == h.hash('hello world')"""
        a = Leaf(b"hello ", h)
        b = Leaf(b"world", h)
        result = rope_concat(a, b, h)
        assert rope_hash(result) == h.hash(b"hello world")

    def test_repeat_hash_matches_direct(self, h):
        """RepeatNode hash must match naive repetition hash."""
        for q in [2, 5, 10, 100, 1000]:
            child = Leaf(b"abc", h)
            rep = RepeatNode(child, q, h)
            assert rope_hash(rep) == h.hash(b"abc" * q), f"q={q}"

    def test_split_repeat_rejoin_hash(self, h):
        """Split a repeat, rejoin — hash must match original."""
        child = Leaf(b"xy", h)
        rep = RepeatNode(child, 50, h)
        for pos in [1, 10, 25, 49, 50, 75, 99]:
            left, right = rope_split(rep, pos, h)
            rejoined = rope_concat(left, right, h)
            assert rope_hash(rejoined) == rep.hash_val, f"pos={pos}"


# =========================================================================
# Property-based weight sweeps
# =========================================================================

class TestPropertyBasedBalance:
    """For random (wl, wr) pairs, verify _balance output satisfies BB[2/7]."""

    def test_weight_sweep_balance_output_valid(self, h):
        """Sweep wl, wr in [1..30]: _balance output must pass validate_rope."""
        child = Leaf(b"x", h)
        for wl in range(1, 31):
            for wr in range(1, 31):
                left = _build_weighted(wl, h)
                right = _build_weighted(wr, h)
                result = _balance(left, right, h)
                validate_rope(result, h)
                assert result.weight == wl + wr

    def test_balance_preserves_content(self, h):
        """_balance must not lose or reorder bytes."""
        for wl in range(1, 15):
            for wr in range(1, 15):
                left_data = bytes(range(wl))
                right_data = bytes(range(100, 100 + wr))
                left = _build_from_bytes_seq(left_data, h)
                right = _build_from_bytes_seq(right_data, h)
                result = _balance(left, right, h)
                assert rope_to_bytes(result) == left_data + right_data

    def test_balance_preserves_hash(self, h):
        """_balance must produce correct hash."""
        for wl in range(1, 15):
            for wr in range(1, 15):
                left_data = bytes(range(wl))
                right_data = bytes(range(100, 100 + wr))
                left = _build_from_bytes_seq(left_data, h)
                right = _build_from_bytes_seq(right_data, h)
                result = _balance(left, right, h)
                assert rope_hash(result) == h.hash(left_data + right_data)


# =========================================================================
# Recursion depth safety
# =========================================================================

class TestRecursionDepthSafety:
    """Verify operations don't blow Python's recursion limit."""

    def test_large_repeat_concat(self, h):
        """RepeatNode with large reps + Leaf must not stack overflow."""
        child = Leaf(b"x", h)
        rep = RepeatNode(child, 10000, h)
        result = rope_concat(rep, Leaf(b"y", h), h)
        validate_rope(result, h)
        assert rope_len(result) == 10001

    def test_large_repeat_split(self, h):
        """Splitting a large RepeatNode must not stack overflow."""
        child = Leaf(b"ab", h)
        rep = RepeatNode(child, 5000, h)
        left, right = rope_split(rep, 5000, h)
        validate_rope(left, h)
        validate_rope(right, h)

    def test_deep_sequential_concat(self, h):
        """2000 sequential concats must stay within recursion limit."""
        node = None
        for i in range(2000):
            node = rope_concat(node, Leaf(bytes([i % 256]), h), h)
        validate_rope(node, h)
        height = rope_height(node)
        bound = math.ceil(math.log(2000) / math.log(7 / 5))
        assert height <= bound

    def test_height_bound_universal(self, h):
        """Height bound O(log w / log(7/5)) holds for mixed operations."""
        node = None
        for i in range(1, 50):
            if i % 5 == 0:
                child = Leaf(bytes([i]), h)
                rep = RepeatNode(child, 3, h)
                node = rope_concat(node, rep, h)
            else:
                node = rope_concat(node, Leaf(bytes([i]), h), h)
            w = node.weight
            height = rope_height(node)
            bound = math.ceil(math.log(max(w, 2)) / math.log(7 / 5))
            assert height <= bound, f"step {i}: h={height}, bound={bound}, w={w}"


# =========================================================================
# Helpers
# =========================================================================

def _build_weighted(w, h):
    """Build a valid tree of exact weight w."""
    if w == 1:
        return Leaf(b"\x01", h)
    node = None
    for i in range(w):
        node = rope_concat(node, Leaf(bytes([i % 256]), h), h)
    return node


def _build_from_bytes_seq(data, h):
    """Build a tree with one Leaf per byte (weight = len(data))."""
    node = None
    for b in data:
        node = rope_concat(node, Leaf(bytes([b]), h), h)
    return node
