"""
Tests that reproduce BB[2/7] balance bugs in rope.py (Bug A and Bug B).

These tests are written BEFORE the fix. They MUST FAIL on the current
v0.2.0 code and PASS after the fix is applied. This is strict TDD.

Bug A: _rebalance cannot handle RepeatNode — isinstance(x, Internal)
       check fails, returning unbalanced trees silently.
Bug B: _rotate_left/_rotate_right create unbalanced inner children
       because they don't check balance of the newly created Internal.
"""

import math
import pytest
from hashrope import (
    PolynomialHash,
    Leaf,
    Internal,
    RepeatNode,
    rope_concat,
    rope_split,
    rope_repeat,
    rope_len,
    rope_hash,
    rope_height,
    rope_to_bytes,
    validate_rope,
    MERSENNE_61,
)


@pytest.fixture
def h():
    return PolynomialHash(prime=MERSENNE_61, base=131)


# =========================================================================
# Bug A: _rebalance fails on RepeatNode
# =========================================================================


class TestBugA_RepeatNodeRebalance:
    """Bug A: concat(RepeatNode, Leaf) produces unbalanced trees."""

    def test_repeat3_concat_leaf(self, h):
        """RepeatNode(child, 3) + Leaf -> weight 4, must be balanced."""
        rep = RepeatNode(Leaf(b"x", h), 3, h)
        leaf = Leaf(b"y", h)
        result = rope_concat(rep, leaf, h)
        validate_rope(result, h)

    def test_repeat5_concat_leaf(self, h):
        """RepeatNode(child, 5) + Leaf -> weight 6, must be balanced."""
        rep = RepeatNode(Leaf(b"x", h), 5, h)
        leaf = Leaf(b"y", h)
        result = rope_concat(rep, leaf, h)
        validate_rope(result, h)

    def test_repeat9_concat_leaf(self, h):
        """RepeatNode(child, 9) + Leaf -> weight 10, must be balanced."""
        rep = RepeatNode(Leaf(b"a", h), 9, h)
        leaf = Leaf(b"b", h)
        result = rope_concat(rep, leaf, h)
        validate_rope(result, h)

    def test_repeat20_concat_leaf(self, h):
        """RepeatNode(child, 20) + Leaf -> weight 21."""
        rep = RepeatNode(Leaf(b"z", h), 20, h)
        leaf = Leaf(b"w", h)
        result = rope_concat(rep, leaf, h)
        validate_rope(result, h)

    def test_sweep_repeat_concat_leaf(self, h):
        """Sweep q in 2..100: RepeatNode(child, q) + Leaf must always validate."""
        child = Leaf(b"x", h)
        leaf = Leaf(b"y", h)
        for q in range(2, 101):
            rep = RepeatNode(child, q, h)
            result = rope_concat(rep, leaf, h)
            validate_rope(result, h), f"Failed at q={q}"

    def test_leaf_concat_repeat_sweep(self, h):
        """Sweep: Leaf + RepeatNode(child, q) must always validate."""
        child = Leaf(b"x", h)
        leaf = Leaf(b"y", h)
        for q in range(2, 101):
            rep = RepeatNode(child, q, h)
            result = rope_concat(leaf, rep, h)
            validate_rope(result, h), f"Failed at q={q}"

    def test_repeat_concat_repeat(self, h):
        """RepeatNode + RepeatNode with different reps must validate."""
        child = Leaf(b"a", h)
        for q1, q2 in [(2, 5), (3, 10), (7, 3), (15, 2), (4, 20)]:
            rep1 = RepeatNode(child, q1, h)
            rep2 = RepeatNode(child, q2, h)
            result = rope_concat(rep1, rep2, h)
            validate_rope(result, h), f"Failed at q1={q1}, q2={q2}"

    def test_accumulate_small_repeats(self, h):
        """Accumulate many small RepeatNodes via concat, validate at every step."""
        node = None
        for i in range(1, 50):
            child = Leaf(bytes([i % 256]), h)
            reps = 2 + (i % 5)
            rep = RepeatNode(child, reps, h)
            node = rope_concat(node, rep, h)
            validate_rope(node, h), f"Failed at step {i}"

    def test_repeat_concat_preserves_content(self, h):
        """Content correctness after concat involving RepeatNode."""
        rep = RepeatNode(Leaf(b"ab", h), 4, h)  # "abababab"
        leaf = Leaf(b"cd", h)
        result = rope_concat(rep, leaf, h)
        assert rope_to_bytes(result) == b"ababababcd"
        validate_rope(result, h)

    def test_repeat_concat_preserves_hash(self, h):
        """Hash correctness after concat involving RepeatNode."""
        rep = RepeatNode(Leaf(b"ab", h), 4, h)
        leaf = Leaf(b"cd", h)
        result = rope_concat(rep, leaf, h)
        expected = h.hash(b"ababababcd")
        assert rope_hash(result) == expected
        validate_rope(result, h)


# =========================================================================
# Bug B: Rotations create unbalanced inner children
# =========================================================================


class TestBugB_RotationInnerBalance:
    """Bug B: split-and-reverse exposes unchecked inner nodes after rotation."""

    def test_4_leaves_split_at_1_reverse(self, h):
        """Build 4 leaves, split@1, concat reversed -> must be balanced."""
        node = None
        for c in b"ABCD":
            node = rope_concat(node, Leaf(bytes([c]), h), h)
        left, right = rope_split(node, 1, h)
        result = rope_concat(right, left, h)
        validate_rope(result, h)

    def test_4_leaves_split_at_3_reverse(self, h):
        """Build 4 leaves, split@3, concat reversed -> must be balanced."""
        node = None
        for c in b"ABCD":
            node = rope_concat(node, Leaf(bytes([c]), h), h)
        left, right = rope_split(node, 3, h)
        result = rope_concat(right, left, h)
        validate_rope(result, h)

    def test_8_leaves_split_at_2_reverse(self, h):
        """Build 8 leaves, split@2, concat reversed."""
        node = None
        for i in range(8):
            node = rope_concat(node, Leaf(bytes([65 + i]), h), h)
        left, right = rope_split(node, 2, h)
        result = rope_concat(right, left, h)
        validate_rope(result, h)

    def test_split_reverse_sweep(self, h):
        """For n in [4..20], split at every position and reverse-concat."""
        for n in range(4, 21):
            node = None
            for i in range(n):
                node = rope_concat(node, Leaf(bytes([i % 256]), h), h)
            for pos in range(1, n):
                left, right = rope_split(node, pos, h)
                result = rope_concat(right, left, h)
                validate_rope(result, h), f"Failed at n={n}, pos={pos}"

    def test_split_reverse_preserves_content(self, h):
        """Content must be correct after split-reverse."""
        data = b"ABCD"
        node = None
        for c in data:
            node = rope_concat(node, Leaf(bytes([c]), h), h)
        left, right = rope_split(node, 1, h)
        result = rope_concat(right, left, h)
        assert rope_to_bytes(result) == b"BCDA"
        validate_rope(result, h)


# =========================================================================
# Combined: LZ77-style workloads (triggers both bugs)
# =========================================================================


class TestLZ77WorkloadBugs:
    """LZ77-style split+repeat+concat chains — the real-world trigger."""

    def test_lz77_chain_validate_every_step(self, h):
        """Split, repeat, concat — validate at EVERY step."""
        node = Leaf(b"abcdef", h)
        for i in range(20):
            nlen = rope_len(node)
            split_pos = max(1, nlen // 3)
            left, right = rope_split(node, split_pos, h)
            if right is None:
                break
            rep = rope_repeat(right, 3, h) if right else right
            node = rope_concat(left, rep, h)
            validate_rope(node, h), f"Failed at LZ77 step {i}"

    def test_lz77_content_integrity(self, h):
        """LZ77 chain must produce correct byte content."""
        node = Leaf(b"abc", h)
        # Split at 1: left=a, right=bc
        left, right = rope_split(node, 1, h)
        # Repeat right 3x: bcbcbc
        rep = rope_repeat(right, 3, h)
        # Concat: a + bcbcbc = abcbcbc
        node = rope_concat(left, rep, h)
        assert rope_to_bytes(node) == b"abcbcbc"
        assert rope_hash(node) == h.hash(b"abcbcbc")
        validate_rope(node, h)

    def test_lz77_hash_integrity(self, h):
        """Hash must match naive computation through LZ77 chain."""
        node = Leaf(b"hello world", h)
        for i in range(10):
            nlen = rope_len(node)
            split_pos = max(1, nlen // 4)
            left, right = rope_split(node, split_pos, h)
            if right is None:
                break
            rep = rope_repeat(right, 2, h) if right else right
            node = rope_concat(left, rep, h)
            # Hash must match reconstructed bytes
            actual_bytes = rope_to_bytes(node)
            assert rope_hash(node) == h.hash(actual_bytes), f"Hash mismatch at step {i}"
            validate_rope(node, h), f"Validation failed at step {i}"


# =========================================================================
# Stress tests (will fail due to balance violations)
# =========================================================================


class TestStressBugs:
    """Stress patterns that expose the bugs under volume."""

    def test_1000_sequential_concats_all_valid(self, h):
        """1000 sequential Leaf concats — validate final tree."""
        node = None
        for i in range(1000):
            node = rope_concat(node, Leaf(bytes([i % 256]), h), h)
        validate_rope(node, h)

    def test_100_repeat_accumulations_all_valid(self, h):
        """100 RepeatNode accumulations — validate at every step."""
        node = None
        for i in range(1, 101):
            child = Leaf(bytes([i % 256]), h)
            reps = 2 + (i % 7)
            rep = RepeatNode(child, reps, h)
            node = rope_concat(node, rep, h)
            validate_rope(node, h), f"Failed at step {i}"

    def test_split_at_every_position_and_rejoin(self, h):
        """Split at every byte position and rejoin — hash must be preserved."""
        node = None
        for i in range(20):
            node = rope_concat(node, Leaf(bytes([65 + i]), h), h)
        original_hash = rope_hash(node)
        original_bytes = rope_to_bytes(node)
        for pos in range(rope_len(node) + 1):
            left, right = rope_split(node, pos, h)
            rejoined = rope_concat(left, right, h)
            assert rope_hash(rejoined) == original_hash, f"Hash mismatch at pos={pos}"
            assert rope_to_bytes(rejoined) == original_bytes, f"Content mismatch at pos={pos}"
            validate_rope(rejoined, h), f"Validation failed at pos={pos}"

    def test_height_bound_after_repeat_concats(self, h):
        """Height must satisfy O(log w) bound even with RepeatNodes."""
        node = None
        for i in range(1, 30):
            child = Leaf(bytes([i]), h)
            reps = 2 + (i % 10)
            rep = RepeatNode(child, reps, h)
            node = rope_concat(node, rep, h)
        w = node.weight
        height = rope_height(node)
        bound = math.ceil(math.log(w) / math.log(7 / 5))
        validate_rope(node, h)
        assert height <= bound, f"height {height} exceeds bound {bound} for w={w}"
