"""Exact mismatch localization across native-comparison chunk boundaries."""
import pytest

from hashrope import Leaf, PolynomialHash, rope_equal, rope_lcp, rope_range_lcp


@pytest.mark.parametrize("length", [1, 2, 3, 7, 16, 63, 64, 65, 127])
def test_every_mismatch_position_in_short_leaves(length):
    data = bytes((i * 17) % 256 for i in range(length))
    h = PolynomialHash()
    left = Leaf(data, h)
    for position in range(length):
        changed = bytearray(data)
        changed[position] ^= 1
        right = Leaf(bytes(changed), h)
        assert rope_lcp(left, right) == position
        assert not rope_equal(left, right)


@pytest.mark.parametrize("position", [0, 1, 2, 32767, 32768, 65534, 65535,
                                      65536, 65537, 98303, 131071, 131072])
def test_long_leaf_mismatch_and_offset_ranges(position):
    data = bytes(range(256)) * 512 + b"z"
    changed = bytearray(data)
    changed[position] ^= 1
    h = PolynomialHash()
    left = Leaf(b"left" + data + b"tail", h)
    right = Leaf(b"right-prefix" + bytes(changed) + b"other", h)
    assert rope_range_lcp(left, 4, right, 12, len(data)) == position
    if position:
        assert rope_range_lcp(left, 4, right, 12, position) == position


def test_later_mismatches_do_not_hide_the_first():
    data = b"x" * 131073
    changed = bytearray(data)
    for position in (129, 32768, 65535, 65536, 131072):
        changed[position] = ord("y")
    h = PolynomialHash()
    assert rope_lcp(Leaf(data, h), Leaf(bytes(changed), h)) == 129
