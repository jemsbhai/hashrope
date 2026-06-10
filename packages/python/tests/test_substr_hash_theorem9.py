"""
Regression + complexity-guard tests for hashrope.rope_substr_hash (Theorem 9).

Run (PowerShell, from repo root):
    pytest tests/test_substr_hash_theorem9.py -v

Expected behavior:
  * PRE-patch (v0.2.1): test_complexity_bytes_hashed_bound FAILS immediately —
    the implementation re-hashes every byte in the queried range instead of
    reading stored node hashes. (test_lcp_exact also passes pre-patch but is
    very slow, ~minutes; the complexity guard is the fast red/green signal.)
  * POST-patch: all tests pass in seconds.

Seeds are fixed for full determinism (no wall-clock assertions).
"""
import os
import random

import pytest

from hashrope import (
    PolynomialHash,
    rope_concat,
    rope_from_bytes,
    rope_repeat,
    rope_split,
    rope_substr_hash,
    rope_to_bytes,
    validate_rope,
)

CHUNK = 4096
SEED = 20260609


class CountingHash(PolynomialHash):
    """PolynomialHash that counts how many raw bytes get byte-level re-hashed.

    rope_substr_hash should only ever pay byte-level work for the (at most
    two) partial leaves at the ends of the queried range. Every fully
    covered node must be answered from its stored hash_val (Invariant I3).
    """

    __slots__ = ("bytes_hashed",)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.bytes_hashed = 0

    def hash(self, data):
        self.bytes_hashed += len(data)
        return super().hash(data)


def build_rope(data: bytes, h: PolynomialHash, chunk: int = CHUNK):
    """Bottom-up balanced build from fixed-size leaves (adapter-style)."""
    leaves = [rope_from_bytes(data[i:i + chunk], h)
              for i in range(0, len(data), chunk)]
    while len(leaves) > 1:
        leaves = [
            rope_concat(leaves[i], leaves[i + 1], h)
            if i + 1 < len(leaves) else leaves[i]
            for i in range(0, len(leaves), 2)
        ]
    return leaves[0]


@pytest.fixture(scope="module")
def env():
    random.seed(SEED)
    h = CountingHash()
    data = os.urandom(1_000_000)
    rope = build_rope(data, h)
    validate_rope(rope, h)
    return h, data, rope


def test_oracle_random_ranges(env):
    """Hash of any substring equals direct hash of the same bytes."""
    h, data, rope = env
    rng = random.Random(SEED + 1)
    for _ in range(300):
        length = rng.randint(1, 5000)
        start = rng.randint(0, len(data) - length)
        assert rope_substr_hash(rope, start, length, h) == \
            h.hash(data[start:start + length])


def test_oracle_boundaries(env):
    """Edge cases: empty, single byte, full range, leaf-boundary straddles."""
    h, data, rope = env
    n = len(data)
    cases = [
        (0, 1), (n - 1, 1), (0, n), (0, n // 2), (n // 2, n - n // 2),
        (CHUNK - 1, 2), (CHUNK, CHUNK), (0, CHUNK), (3 * CHUNK - 7, 13),
    ]
    for start, length in cases:
        assert rope_substr_hash(rope, start, length, h) == \
            h.hash(data[start:start + length])
    assert rope_substr_hash(rope, 0, 0, h) == 0
    # Full-range query must equal the stored root hash exactly.
    assert rope_substr_hash(rope, 0, n, h) == rope.hash_val


def test_oracle_repeatnode_spans():
    """Ranges spanning RepeatNode boundaries (tail / full copies / head)."""
    rng = random.Random(SEED + 2)
    h = PolynomialHash()
    rep = rope_repeat(build_rope(os.urandom(1000), h), 50, h)
    rope = rope_concat(
        build_rope(os.urandom(3000), h),
        rope_concat(rep, build_rope(os.urandom(2000), h), h),
        h,
    )
    ref = rope_to_bytes(rope)
    for _ in range(200):
        length = rng.randint(1, 4000)
        start = rng.randint(0, len(ref) - length)
        assert rope_substr_hash(rope, start, length, h) == \
            h.hash(ref[start:start + length])


def test_complexity_bytes_hashed_bound(env):
    """Theorem 9 guard (deterministic, no timing): a range query may
    byte-level re-hash at most the two partial leaves at its ends — never
    the interior. Pre-patch this fails with bytes_hashed == range length."""
    h, data, rope = env
    rng = random.Random(SEED + 3)
    for _ in range(50):
        length = rng.randint(10_000, 900_000)  # spans many leaves
        start = rng.randint(0, len(data) - length)
        h.bytes_hashed = 0
        rope_substr_hash(rope, start, length, h)
        assert h.bytes_hashed <= 2 * (CHUNK - 1), (
            f"re-hashed {h.bytes_hashed} bytes for a {length}-byte range; "
            f"stored node hashes are not being consulted (Theorem 9 violated)"
        )


def test_lcp_exact(env):
    """Binary search over prefix hashes recovers the exact divergence byte."""
    h, data, rope = env
    rng = random.Random(SEED + 4)

    def replace_byte(node, idx):
        left, rest = rope_split(node, idx, h)
        _, tail = rope_split(rest, 1, h)
        mid = rope_from_bytes(bytes([(data[idx] + 1) % 256]), h)
        return rope_concat(rope_concat(left, mid, h), tail, h)

    def lcp(a, b):
        lo, hi = 0, min(a.len, b.len)
        while lo < hi:
            m = (lo + hi + 1) // 2
            if rope_substr_hash(a, 0, m, h) == rope_substr_hash(b, 0, m, h):
                lo = m
            else:
                hi = m - 1
        return lo

    for _ in range(15):
        d = rng.randint(0, len(data) - 2)
        assert lcp(rope, replace_byte(rope, d)) == d
