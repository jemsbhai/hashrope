"""Independent exact-comparison checks against bytes, shapes, and repetition."""

import random
import sys
from contextlib import contextmanager

import pytest

import hashrope.rope as rope_module
from hashrope import (
    Internal,
    Leaf,
    MERSENNE_127,
    PolynomialHash,
    RepeatNode,
    rope_concat,
    rope_equal,
    rope_from_bytes,
    rope_lcp,
    rope_repeat,
    rope_split,
    rope_to_bytes,
)


def _oracle_lcp(a, b):
    return next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))


def _chunked(data, width, h):
    node = None
    for start in range(0, len(data), width):
        node = rope_concat(node, rope_from_bytes(data[start:start + width], h), h)
    return node


@contextmanager
def _bounded_comparison_work(limit=50_000):
    """Fail promptly if a compressed comparison enumerates enormous repeats.

    A deterministic source-line budget avoids platform-dependent timing limits.
    Construction is outside this guard; only the requested comparisons count.
    """
    previous_trace = sys.gettrace()
    remaining = limit

    def trace(frame, event, arg):
        nonlocal remaining
        if event == "line" and frame.f_code.co_filename == rope_module.__file__:
            remaining -= 1
            assert remaining >= 0, "comparison expanded a compressed repetition"
        return trace

    sys.settrace(trace)
    try:
        yield
    finally:
        sys.settrace(previous_trace)


@pytest.mark.parametrize("prefix", [b"", b"same prefix", b"x" * 257])
def test_collision_after_shared_prefix_is_not_equality(prefix):
    h = PolynomialHash()
    shared = rope_from_bytes(prefix, h)
    a = rope_concat(shared, Leaf(b"\x00\x83*", h), h)
    b = rope_concat(shared, Leaf(b"\x01\x00*", h), h)
    assert a.hash_val == b.hash_val
    assert not rope_equal(a, b)
    assert rope_lcp(a, b) == len(prefix)
    assert rope_lcp(b, a) == len(prefix)


@pytest.mark.parametrize("data", [b"", b"a", bytes(range(256)) * 3])
def test_equal_content_across_profiles_and_chunk_boundaries(data):
    h = PolynomialHash()
    other_h = PolynomialHash(prime=MERSENNE_127, base=257)
    a = _chunked(data, 7, h)
    b = _chunked(data, 19, other_h)
    assert rope_equal(a, b)
    assert rope_lcp(a, b) == len(data)
    longer = rope_concat(b, Leaf(b"!", other_h), other_h)
    assert not rope_equal(a, longer)
    assert rope_lcp(a, longer) == len(data)
    assert rope_lcp(longer, a) == len(data)


def test_independent_equal_repeats_with_different_period_lengths():
    h = PolynomialHash()
    repetitions = 10**12
    a = rope_repeat(Leaf(b"ab", h), repetitions, h)
    b = rope_repeat(Leaf(b"abab", h), repetitions // 2, h)
    with _bounded_comparison_work():
        assert rope_equal(a, b)
        assert rope_lcp(a, b) == 2 * repetitions


def test_periodic_witness_checks_its_last_byte():
    h = PolynomialHash()
    # Periods 3 and 4 agree for five bytes, then differ. The exact periodic
    # equality witness is 3 + 4 - gcd(3, 4) = 6, including that sixth byte.
    a = rope_repeat(Leaf(b"aab", h), 10**12, h)
    b = rope_repeat(Leaf(b"aaba", h), 10**12, h)
    with _bounded_comparison_work():
        assert rope_lcp(a, b) == 5
        assert rope_lcp(b, a) == 5


def test_shifted_repeat_phase_and_mismatch_after_huge_prefix():
    h = PolynomialHash()
    repetitions = 10**12
    a = rope_repeat(Leaf(b"ab", h), repetitions + 1, h)
    middle = rope_repeat(Leaf(b"ba", h), repetitions, h)
    b = rope_concat(rope_concat(Leaf(b"a", h), middle, h), Leaf(b"b", h), h)
    different_tail = rope_concat(rope_concat(Leaf(b"a", h), middle, h), Leaf(b"c", h), h)
    with _bounded_comparison_work():
        assert rope_equal(a, b)
        assert rope_lcp(a, b) == 2 * (repetitions + 1)
        assert not rope_equal(a, different_tail)
        assert rope_lcp(a, different_tail) == 2 * repetitions + 1


def test_nested_repeats_and_shorter_prefix_without_materializing():
    h = PolynomialHash()
    other_h = PolynomialHash(base=257)
    repetitions = 1 << 70
    # Direct RepeatNode construction preserves the nested representation.
    a = RepeatNode(RepeatNode(Leaf(b"abc", h), 7, h), repetitions, h)
    b = rope_repeat(Leaf(b"abcabcabc", other_h), (7 * repetitions) // 3, other_h)
    assert a.len > b.len
    with _bounded_comparison_work():
        assert not rope_equal(a, b)
        assert rope_lcp(a, b) == b.len
        assert rope_lcp(b, a) == b.len


def test_independent_equal_self_concat_dags_do_not_expand():
    h = PolynomialHash()
    a = Leaf(b"ab", h)
    b = Leaf(b"ab", h)
    for _ in range(40):
        a = Internal(a, a, h)
        b = Internal(b, b, h)
    with _bounded_comparison_work():
        assert rope_equal(a, b)
        assert rope_lcp(a, b) == 2 << 40


@pytest.mark.parametrize("position", [0, 1, 2, 3, 4, 5])
def test_mismatch_in_repeated_period(position):
    h = PolynomialHash()
    unit = b"abcabc"
    changed = unit[:position] + b"!" + unit[position + 1:]
    a = rope_repeat(Leaf(unit, h), 10**12, h)
    b = rope_repeat(_chunked(changed, 2, h), 10**12, h)
    with _bounded_comparison_work():
        assert not rope_equal(a, b)
        assert rope_lcp(a, b) == position
        assert rope_lcp(b, a) == position


def test_randomized_nested_shapes_against_byte_oracle():
    rng = random.Random(20261023)
    h = PolynomialHash()
    alternate_h = PolynomialHash(base=257)

    def tree(depth):
        if depth == 0 or rng.randrange(4) == 0:
            data = rng.randbytes(rng.randrange(1, 10))
            return Leaf(data, h)
        if rng.randrange(2):
            return rope_repeat(tree(depth - 1), rng.randrange(2, 5), h)
        return rope_concat(tree(depth - 1), tree(depth - 1), h)

    for _ in range(100):
        a = tree(3)
        data = rope_to_bytes(a)
        b = _chunked(data, rng.randrange(1, 13), alternate_h)
        assert rope_equal(a, b)
        assert rope_lcp(a, b) == len(data)

        position = rng.randrange(len(data) + 1)
        left, right = rope_split(a, position, h)
        edited = rope_concat(rope_concat(left, Leaf(b"\xff", h), h), right, h)
        expected = data[:position] + b"\xff" + data[position:]
        assert rope_lcp(a, edited) == _oracle_lcp(data, expected)
        assert rope_lcp(edited, b) == _oracle_lcp(expected, data)
        assert not rope_equal(a, edited)

        unrelated = tree(3)
        unrelated_bytes = rope_to_bytes(unrelated)
        assert rope_lcp(a, unrelated) == _oracle_lcp(data, unrelated_bytes)
        assert rope_equal(a, unrelated) == (data == unrelated_bytes)
