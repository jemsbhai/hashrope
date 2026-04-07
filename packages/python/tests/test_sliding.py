"""Tests for hashrope.sliding module."""

import pytest
from hashrope import PolynomialHash, SlidingWindow


class TestSlidingWindowBasic:
    def test_empty(self):
        sw = SlidingWindow()
        assert sw.current_hash() == 0
        assert sw.pos == 0
        assert sw.window_len == 0

    def test_append_bytes(self):
        h = PolynomialHash()
        sw = SlidingWindow()
        sw.append_bytes(b"hello")
        assert sw.current_hash() == h.hash(b"hello")
        assert sw.pos == 5

    def test_append_bytes_incremental(self):
        h = PolynomialHash()
        sw = SlidingWindow()
        sw.append_bytes(b"hello ")
        sw.append_bytes(b"world")
        assert sw.current_hash() == h.hash(b"hello world")

    def test_append_empty(self):
        sw = SlidingWindow()
        sw.append_bytes(b"")
        assert sw.pos == 0
        assert sw.current_hash() == 0

    def test_final_hash_equals_current(self):
        sw = SlidingWindow()
        sw.append_bytes(b"test")
        assert sw.final_hash() == sw.current_hash()


class TestSlidingWindowCopy:
    def test_non_overlapping_copy(self):
        h = PolynomialHash()
        sw = SlidingWindow()
        sw.append_bytes(b"hello world")
        sw.append_copy(offset=5, length=5)  # copy "world"
        assert sw.current_hash() == h.hash(b"hello worldworld")

    def test_overlapping_copy(self):
        h = PolynomialHash()
        sw = SlidingWindow()
        sw.append_bytes(b"ab")
        sw.append_copy(offset=2, length=6)  # "ab" repeated -> "abababab"
        assert sw.current_hash() == h.hash(b"abababab")

    def test_overlapping_copy_partial(self):
        h = PolynomialHash()
        sw = SlidingWindow()
        sw.append_bytes(b"abc")
        sw.append_copy(offset=3, length=7)  # "abcabca" -> "abcabcabca"
        assert sw.current_hash() == h.hash(b"abcabcabca")

    def test_copy_single_byte(self):
        h = PolynomialHash()
        sw = SlidingWindow()
        sw.append_bytes(b"x")
        sw.append_copy(offset=1, length=4)  # "xxxxx"
        assert sw.current_hash() == h.hash(b"xxxxx")

    def test_copy_offset_exceeds_pos_raises(self):
        sw = SlidingWindow()
        sw.append_bytes(b"hi")
        with pytest.raises(ValueError, match="exceeds decoded length"):
            sw.append_copy(offset=10, length=1)

    def test_copy_offset_exceeds_window_raises(self):
        sw = SlidingWindow(d_max=4, m_max=4)
        # Fill enough to trigger eviction, then try a large offset
        sw.append_bytes(b"a" * 20)
        with pytest.raises(ValueError, match="exceeds window size"):
            sw.append_copy(offset=15, length=1)


class TestSlidingWindowEviction:
    def test_eviction_preserves_hash(self):
        """Hash must be correct even after eviction."""
        h = PolynomialHash()
        sw = SlidingWindow(d_max=8, m_max=4)
        data = b"the quick brown fox jumps over the lazy dog"
        for byte in data:
            sw.append_bytes(bytes([byte]))
        assert sw.current_hash() == h.hash(data)

    def test_window_bounded(self):
        sw = SlidingWindow(d_max=16, m_max=8)
        sw.append_bytes(b"a" * 100)
        # Window should be bounded: at most d_max after eviction
        assert sw.window_len <= sw.d_max

    def test_eviction_with_copies(self):
        """Mix of bytes and copies with a small window."""
        h = PolynomialHash()
        sw = SlidingWindow(d_max=16, m_max=8)
        sw.append_bytes(b"abcdef")
        sw.append_copy(offset=3, length=6)  # "defdef"
        sw.append_bytes(b"ghi")
        # Total: "abcdefdefdefghi"
        expected = b"abcdef" + b"def" * 2 + b"ghi"
        assert sw.current_hash() == h.hash(expected)

    def test_large_stream(self):
        """Stream a lot of data through a small window."""
        h = PolynomialHash()
        sw = SlidingWindow(d_max=32, m_max=16)
        chunks = [bytes([i % 256] * 10) for i in range(50)]
        full_data = b"".join(chunks)
        for chunk in chunks:
            sw.append_bytes(chunk)
        assert sw.current_hash() == h.hash(full_data)


class TestSlidingWindowInit:
    def test_invalid_d_max(self):
        with pytest.raises(ValueError):
            SlidingWindow(d_max=0)

    def test_invalid_m_max(self):
        with pytest.raises(ValueError):
            SlidingWindow(m_max=0)

    def test_custom_params(self):
        sw = SlidingWindow(d_max=1024, m_max=512, base=257)
        assert sw.d_max == 1024
        assert sw.m_max == 512
        assert sw.W == 1536
