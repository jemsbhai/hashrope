"""Exact ranges must agree with byte slices without allocating rope nodes."""

import random

import pytest

from hashrope import (
    Leaf, PolynomialHash, rope_concat, rope_from_bytes, rope_range_equal,
    rope_range_lcp, rope_repeat,
)


def test_collision_inside_ranges_and_cross_profile_equality():
    h = PolynomialHash()
    left = Leaf(b"a\x00\x83*z", h)
    right = Leaf(b"bb\x01\x00*z", h)
    assert not rope_range_equal(left, 1, right, 2, 3)
    assert rope_range_lcp(left, 1, right, 2, 3) == 0
    other = Leaf(b"xx\x00\x83*", PolynomialHash(base=257))
    assert rope_range_equal(left, 1, other, 2, 3)


def test_same_node_different_offsets_and_shifted_huge_repeat():
    h = PolynomialHash()
    left = rope_repeat(Leaf(b"ab", h), 10**12, h)
    right = rope_repeat(Leaf(b"baba", h), 10**12, h)
    assert rope_range_lcp(left, 0, left, 1, 100) == 0
    assert rope_range_equal(left, 1, right, 0, 10**12)
    assert rope_range_lcp(left, 1, right, 0, 10**12) == 10**12


@pytest.mark.parametrize("start,length", [(-1, 0), (0, -1), (4, 0), (2, 2)])
def test_invalid_ranges_rejected_on_either_side(start, length):
    node = Leaf(b"abc", PolynomialHash())
    for compare in (rope_range_equal, rope_range_lcp):
        with pytest.raises(ValueError):
            compare(node, start, node, 0, length)
        with pytest.raises(ValueError):
            compare(node, 0, node, start, length)
        with pytest.raises(ValueError):
            compare(node, start, node, start, length)


def test_empty_ranges_still_validate_bounds():
    node = Leaf(b"abc", PolynomialHash())
    assert rope_range_equal(None, 0, node, 3, 0)
    assert rope_range_lcp(node, 3, None, 0, 0) == 0
    with pytest.raises(ValueError):
        rope_range_equal(None, 1, None, 1, 0)


def test_random_ranges_across_shapes_match_byte_oracle():
    rng = random.Random(20261024)
    h = PolynomialHash()
    for _ in range(100):
        data = bytes(rng.randrange(4) for _ in range(100))
        left = rope_from_bytes(data, h)
        right = None
        for offset in range(0, len(data), 7):
            right = rope_concat(right, Leaf(data[offset:offset + 7], h), h)
        a, b = rng.randrange(101), rng.randrange(101)
        length = rng.randrange(min(100 - a, 100 - b) + 1)
        x, y = data[a:a + length], data[b:b + length]
        expected = next((i for i, (u, v) in enumerate(zip(x, y)) if u != v), length)
        assert rope_range_lcp(left, a, right, b, length) == expected
        assert rope_range_equal(left, a, right, b, length) == (x == y)
