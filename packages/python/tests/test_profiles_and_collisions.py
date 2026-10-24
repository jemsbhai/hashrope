"""Regression tests for exact content decisions and hash-profile provenance."""

import struct

import pytest

from hashrope import (
    HashProfile, Internal, Leaf, MERSENNE_61, MERSENNE_127, PolynomialHash,
    RepeatNode, SlidingWindow, rope_concat, rope_equal, rope_hash, rope_lcp,
    rope_repeat, rope_split, rope_substr_hash, rope_to_bytes, validate_rope,
)


@pytest.mark.parametrize("prime", [MERSENNE_61, MERSENNE_127])
def test_legacy_collision_fingerprints_remain_compatible_but_comparison_is_exact(prime):
    h = PolynomialHash(prime=prime)
    a, b = Leaf(b"\x00\x83", h), Leaf(b"\x01\x00", h)
    assert rope_hash(a) == rope_hash(b) == 263
    assert rope_substr_hash(a, 0, 2, h) == rope_substr_hash(b, 0, 2, h) == 263
    assert a != b
    assert not rope_equal(a, b)
    assert rope_lcp(a, b) == 0

    # The extra suffix changes the old hash-only binary-search path.
    a, b = Leaf(a.data + b"*", h), Leaf(b.data + b"*", h)
    assert rope_hash(a) == rope_hash(b)
    assert not rope_equal(a, b)
    assert rope_lcp(a, b) == 0

    a, b = Leaf(struct.pack("<I", 33536), h), Leaf(struct.pack("<I", 1), h)
    assert rope_hash(a) == rope_hash(b)
    assert not rope_equal(a, b)
    assert rope_lcp(a, b) // 4 == 0


@pytest.mark.parametrize("other", [PolynomialHash(base=257), PolynomialHash(prime=MERSENNE_127)])
def test_mixed_profiles_rejected_at_construction_and_query_boundaries(other):
    h = PolynomialHash()
    a, b = Leaf(b"AB", h), Leaf(b"CD", other)
    operations = [
        lambda: Internal(a, b, h), lambda: Internal(a, b, other),
        lambda: RepeatNode(b, 2, h), lambda: rope_concat(a, b, h),
        lambda: rope_concat(None, b, h), lambda: rope_concat(b, None, h),
        lambda: rope_split(b, 0, h), lambda: rope_split(b, b.len, h),
        lambda: rope_split(b, 1, h), lambda: rope_repeat(b, 0, h),
        lambda: rope_repeat(b, 1, h), lambda: rope_repeat(b, 2, h),
        lambda: rope_substr_hash(b, 0, 0, h),
        lambda: rope_substr_hash(b, 0, b.len, h),
        lambda: rope_substr_hash(b, 0, 1, h), lambda: validate_rope(b, h),
    ]
    for operation in operations:
        with pytest.raises(ValueError, match="Incompatible hash profile"):
            operation()
    assert rope_to_bytes(a) == b"AB"
    assert rope_to_bytes(b) == b"CD"


def test_equal_profiles_from_distinct_hashers_are_compatible():
    h, other = PolynomialHash(), PolynomialHash()
    assert h.profile == other.profile == HashProfile(MERSENNE_61, 131)
    assert h.profile is not other.profile
    a, b = Leaf(b"AB", h), Leaf(b"CD", other)
    joined = rope_concat(a, b, other)
    assert rope_hash(joined) == h.hash(b"ABCD")
    left, right = rope_split(joined, 1, h)
    restored = rope_concat(left, right, other)
    assert rope_equal(restored, joined)
    assert rope_substr_hash(restored, 1, 2, h) == h.hash(b"BC")
    validate_rope(restored, h)


def test_hash_parameters_are_immutable_but_power_cache_still_works():
    h = PolynomialHash()
    for attribute, value in [("_p", MERSENNE_127), ("_x", 257), ("_profile", HashProfile(MERSENNE_61, 257))]:
        with pytest.raises(AttributeError, match="immutable"):
            setattr(h, attribute, value)
        with pytest.raises(AttributeError, match="immutable"):
            delattr(h, attribute)
    with pytest.raises(AttributeError):
        h.profile.base = 257
    with pytest.raises(AttributeError):
        h.__init__(base=257)
    assert h.power(19) == pow(131, 19, MERSENNE_61)
    assert h.power(19) == pow(131, 19, MERSENNE_61)


def test_stream_profile_replacement_rejected_before_mutating_state():
    stream = SlidingWindow(d_max=2, m_max=1)
    stream.append_bytes(b"abcdef")  # Exercise an already-evicted prefix.
    state = (stream.pos, stream.l_prefix, stream.h_prefix, stream.r_window)
    original = stream.h
    expected = stream.current_hash()
    stream.h = PolynomialHash(base=257)
    for operation in [stream.current_hash, stream.final_hash,
                      lambda: stream.append_bytes(b"!"), lambda: stream.append_bytes(b""),
                      lambda: stream.append_copy(1, 1)]:
        with pytest.raises(ValueError, match="Incompatible hash profile"):
            operation()
    assert (stream.pos, stream.l_prefix, stream.h_prefix, stream.r_window) == state
    stream.h = PolynomialHash()
    assert stream.current_hash() == expected == original.hash(b"abcdef")
    stream.append_bytes(b"!")
    assert stream.current_hash() == original.hash(b"abcdef!")


def test_exact_apis_accept_empty_and_prefix_inputs():
    h = PolynomialHash()
    a, b = Leaf(b"abc", h), Leaf(b"abcd", h)
    assert rope_equal(None, None)
    assert not rope_equal(None, a)
    assert not rope_equal(a, b)
    assert rope_lcp(None, a) == rope_lcp(a, None) == 0
    assert rope_lcp(a, a) == rope_lcp(a, b) == rope_lcp(b, a) == 3


def test_exact_equality_does_not_depend_on_nonstandard_hash_arithmetic():
    # The legacy constructor accepts non-Mersenne moduli, whose folding does
    # not obey the polynomial composition law. Preserve that behavior while
    # keeping exact comparison independent of cached fingerprints altogether.
    h = PolynomialHash(prime=17, base=2)
    a = Leaf(b"\x00\x0d", h)
    b = Internal(Leaf(b"\x00", h), Leaf(b"\x0d", h), h)
    assert a.hash_val != b.hash_val
    assert rope_equal(a, b)
    assert rope_lcp(a, b) == 2
