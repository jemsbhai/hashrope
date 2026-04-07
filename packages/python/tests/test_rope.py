"""Tests for hashrope.rope module."""

import pytest
from hashrope import (
    PolynomialHash,
    Leaf,
    Internal,
    RepeatNode,
    rope_concat,
    rope_split,
    rope_repeat,
    rope_substr_hash,
    rope_len,
    rope_hash,
    rope_from_bytes,
    rope_to_bytes,
    validate_rope,
)


class TestLeaf:
    def test_basic(self, h):
        leaf = Leaf(b"hello", h)
        assert leaf.len == 5
        assert leaf.weight == 1
        assert leaf.hash_val == h.hash(b"hello")
        assert leaf.data == b"hello"

    def test_empty_raises(self, h):
        with pytest.raises(ValueError):
            Leaf(b"", h)

    def test_immutable(self, h):
        leaf = Leaf(b"test", h)
        with pytest.raises(AttributeError):
            leaf.len = 10


class TestInternal:
    def test_basic(self, h):
        left = Leaf(b"hello ", h)
        right = Leaf(b"world", h)
        node = Internal(left, right, h)
        assert node.len == 11
        assert node.weight == 2
        assert node.hash_val == h.hash(b"hello world")

    def test_hash_invariant(self, h):
        """I3: hash = left.hash * x^(right.len) + right.hash."""
        left = Leaf(b"abc", h)
        right = Leaf(b"def", h)
        node = Internal(left, right, h)
        expected = h.hash_concat(left.hash_val, right.len, right.hash_val)
        assert node.hash_val == expected


class TestRepeatNode:
    def test_basic(self, h):
        child = Leaf(b"ab", h)
        node = RepeatNode(child, 5, h)
        assert node.len == 10
        assert node.weight == 5
        assert node.hash_val == h.hash(b"ab" * 5)

    def test_reps_lt_2_raises(self, h):
        child = Leaf(b"x", h)
        with pytest.raises(ValueError):
            RepeatNode(child, 1, h)

    def test_large_reps(self, h):
        child = Leaf(b"xyz", h)
        node = RepeatNode(child, 10000, h)
        assert node.len == 30000
        assert node.hash_val == h.hash(b"xyz" * 10000)


class TestRopeConcat:
    def test_none_left(self, h):
        leaf = Leaf(b"test", h)
        assert rope_concat(None, leaf, h) is leaf

    def test_none_right(self, h):
        leaf = Leaf(b"test", h)
        assert rope_concat(leaf, None, h) is leaf

    def test_both_none(self, h):
        assert rope_concat(None, None, h) is None

    def test_two_leaves(self, h):
        a = Leaf(b"hello ", h)
        b = Leaf(b"world", h)
        ab = rope_concat(a, b, h)
        assert rope_len(ab) == 11
        assert rope_hash(ab) == h.hash(b"hello world")

    def test_associativity(self, h):
        """H((A || B) || C) == H(A || (B || C))."""
        a = Leaf(b"aaa", h)
        b = Leaf(b"bbb", h)
        c = Leaf(b"ccc", h)
        ab_c = rope_concat(rope_concat(a, b, h), c, h)
        a_bc = rope_concat(a, rope_concat(b, c, h), h)
        assert rope_hash(ab_c) == rope_hash(a_bc)
        assert rope_hash(ab_c) == h.hash(b"aaabbbccc")

    def test_many_concats_balanced(self, h):
        """Concatenating many leaves should produce a balanced tree."""
        node = None
        for i in range(50):
            leaf = Leaf(bytes([i % 256]), h)
            node = rope_concat(node, leaf, h)
        validate_rope(node)
        assert rope_len(node) == 50


class TestRopeSplit:
    def test_split_at_zero(self, h):
        leaf = Leaf(b"hello", h)
        left, right = rope_split(leaf, 0, h)
        assert left is None
        assert right is leaf

    def test_split_at_end(self, h):
        leaf = Leaf(b"hello", h)
        left, right = rope_split(leaf, 5, h)
        assert left is leaf
        assert right is None

    def test_split_leaf(self, h):
        leaf = Leaf(b"hello", h)
        left, right = rope_split(leaf, 2, h)
        assert rope_to_bytes(left) == b"he"
        assert rope_to_bytes(right) == b"llo"

    def test_split_preserves_hash(self, h):
        """Split and re-concat should preserve the hash."""
        data = b"the quick brown fox"
        node = Leaf(data, h)
        for pos in range(len(data) + 1):
            left, right = rope_split(node, pos, h)
            rejoined = rope_concat(left, right, h)
            assert rope_hash(rejoined) == h.hash(data), f"Failed at pos={pos}"

    def test_split_internal(self, h):
        a = Leaf(b"hello ", h)
        b = Leaf(b"world", h)
        ab = rope_concat(a, b, h)
        left, right = rope_split(ab, 6, h)
        assert rope_to_bytes(left) == b"hello "
        assert rope_to_bytes(right) == b"world"

    def test_split_none(self, h):
        left, right = rope_split(None, 0, h)
        assert left is None
        assert right is None

    def test_split_repeat_on_boundary(self, h):
        child = Leaf(b"abc", h)
        rep = RepeatNode(child, 4, h)  # "abcabcabcabc"
        left, right = rope_split(rep, 6, h)
        assert rope_to_bytes(left) == b"abcabc"
        assert rope_to_bytes(right) == b"abcabc"

    def test_split_repeat_within(self, h):
        child = Leaf(b"abcd", h)
        rep = RepeatNode(child, 3, h)  # "abcdabcdabcd"
        left, right = rope_split(rep, 5, h)
        assert rope_to_bytes(left) == b"abcda"
        assert rope_to_bytes(right) == b"bcdabcd"
        # Hashes must match
        assert rope_hash(rope_concat(left, right, h)) == rep.hash_val


class TestRopeRepeat:
    def test_repeat_zero(self, h):
        leaf = Leaf(b"x", h)
        assert rope_repeat(leaf, 0, h) is None

    def test_repeat_one(self, h):
        leaf = Leaf(b"x", h)
        assert rope_repeat(leaf, 1, h) is leaf

    def test_repeat_many(self, h):
        leaf = Leaf(b"ab", h)
        rep = rope_repeat(leaf, 100, h)
        assert rope_len(rep) == 200
        assert rope_hash(rep) == h.hash(b"ab" * 100)


class TestRopeSubstrHash:
    def test_full_string(self, h):
        leaf = Leaf(b"hello", h)
        assert rope_substr_hash(leaf, 0, 5, h) == h.hash(b"hello")

    def test_prefix(self, h):
        leaf = Leaf(b"hello", h)
        assert rope_substr_hash(leaf, 0, 3, h) == h.hash(b"hel")

    def test_suffix(self, h):
        leaf = Leaf(b"hello", h)
        assert rope_substr_hash(leaf, 2, 3, h) == h.hash(b"llo")

    def test_middle(self, h):
        leaf = Leaf(b"hello world", h)
        assert rope_substr_hash(leaf, 3, 5, h) == h.hash(b"lo wo")

    def test_spanning_internal(self, h):
        a = Leaf(b"hello ", h)
        b = Leaf(b"world", h)
        ab = rope_concat(a, b, h)
        assert rope_substr_hash(ab, 4, 4, h) == h.hash(b"o wo")

    def test_within_repeat(self, h):
        child = Leaf(b"abc", h)
        rep = RepeatNode(child, 5, h)  # "abcabcabcabcabc"
        # Entirely within one copy
        assert rope_substr_hash(rep, 1, 2, h) == h.hash(b"bc")

    def test_spanning_repeat(self, h):
        child = Leaf(b"abcd", h)
        rep = RepeatNode(child, 4, h)  # "abcdabcdabcdabcd"
        # Spanning multiple copies
        assert rope_substr_hash(rep, 2, 8, h) == h.hash(b"cdabcdab")

    def test_empty_range(self, h):
        leaf = Leaf(b"test", h)
        assert rope_substr_hash(leaf, 0, 0, h) == 0

    def test_none_node(self, h):
        assert rope_substr_hash(None, 0, 0, h) == 0


class TestRopeHelpers:
    def test_from_bytes_empty(self, h):
        assert rope_from_bytes(b"", h) is None

    def test_from_bytes_nonempty(self, h):
        node = rope_from_bytes(b"test", h)
        assert isinstance(node, Leaf)
        assert rope_to_bytes(node) == b"test"

    def test_to_bytes_none(self):
        assert rope_to_bytes(None) == b""

    def test_roundtrip(self, h):
        data = b"round trip test data"
        node = rope_from_bytes(data, h)
        assert rope_to_bytes(node) == data

    def test_to_bytes_complex_tree(self, h):
        """Build a complex tree and verify byte reconstruction."""
        a = Leaf(b"hello ", h)
        b = Leaf(b"world", h)
        ab = rope_concat(a, b, h)
        rep = rope_repeat(Leaf(b"! ", h), 3, h)
        full = rope_concat(ab, rep, h)
        assert rope_to_bytes(full) == b"hello world! ! ! "


class TestValidateRope:
    def test_valid_leaf(self, h):
        validate_rope(Leaf(b"test", h))

    def test_valid_internal(self, h):
        a = Leaf(b"a", h)
        b = Leaf(b"b", h)
        validate_rope(Internal(a, b, h))

    def test_valid_repeat(self, h):
        child = Leaf(b"abc", h)
        validate_rope(RepeatNode(child, 5, h))

    def test_valid_none(self, h):
        validate_rope(None)

    def test_complex_tree(self, h):
        """Build a complex tree via operations and validate."""
        node = None
        for i in range(30):
            leaf = Leaf(bytes([i % 256] * (i % 5 + 1)), h)
            node = rope_concat(node, leaf, h)
        validate_rope(node)
