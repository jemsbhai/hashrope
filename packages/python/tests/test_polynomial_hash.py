"""Tests for hashrope.polynomial_hash module."""

import pytest
from hashrope import (
    PolynomialHash,
    phi,
    mersenne_mod,
    mersenne_mul,
    MERSENNE_61,
    MERSENNE_127,
)


class TestMersenneArithmetic:
    """Tests for Mersenne prime modular arithmetic (Lemma 14)."""

    def test_mersenne_mod_small(self):
        assert mersenne_mod(0) == 0
        assert mersenne_mod(1) == 1
        assert mersenne_mod(MERSENNE_61) == 0
        assert mersenne_mod(MERSENNE_61 + 1) == 1

    def test_mersenne_mod_large(self):
        # Very large value should reduce correctly
        val = MERSENNE_61 * 1000 + 42
        assert mersenne_mod(val) == 42

    def test_mersenne_mul_commutative(self):
        a, b = 12345, 67890
        assert mersenne_mul(a, b) == mersenne_mul(b, a)

    def test_mersenne_mul_identity(self):
        assert mersenne_mul(42, 1) == 42

    def test_mersenne_mul_zero(self):
        assert mersenne_mul(42, 0) == 0

    def test_mersenne_mod_127(self):
        assert mersenne_mod(MERSENNE_127, MERSENNE_127) == 0
        assert mersenne_mod(MERSENNE_127 + 5, MERSENNE_127) == 5


class TestPhi:
    """Tests for geometric accumulator Phi (Definition 3, Theorem 3)."""

    def test_phi_base_cases(self):
        assert phi(0, 131) == 0
        assert phi(1, 131) == 1

    def test_phi_two(self):
        # Phi(2, alpha) = 1 + alpha
        alpha = 131
        assert phi(2, alpha) == (1 + alpha) % MERSENNE_61

    def test_phi_three(self):
        # Phi(3, alpha) = 1 + alpha + alpha^2
        alpha = 131
        expected = (1 + alpha + alpha * alpha) % MERSENNE_61
        assert phi(3, alpha) == expected

    def test_phi_degenerate_alpha_one(self):
        # Phi(q, 1) = q (mod p) -- Lemma 3
        assert phi(100, 1) == 100
        assert phi(1000, 1) == 1000

    def test_phi_large_q(self):
        # Phi(q, alpha) = sum_{i=0}^{q-1} alpha^i
        alpha = 5
        q = 20
        expected = sum(pow(alpha, i, MERSENNE_61) for i in range(q)) % MERSENNE_61
        assert phi(q, alpha) == expected

    def test_phi_power_of_two(self):
        alpha = 131
        q = 64
        expected = sum(pow(alpha, i, MERSENNE_61) for i in range(q)) % MERSENNE_61
        assert phi(q, alpha) == expected

    def test_phi_large_alpha(self):
        alpha = MERSENNE_61 - 1  # largest valid
        q = 10
        expected = sum(pow(alpha, i, MERSENNE_61) for i in range(q)) % MERSENNE_61
        assert phi(q, alpha) == expected


class TestPolynomialHash:
    """Tests for PolynomialHash class (Definition 2, Theorems 1-2, 5)."""

    def test_empty_string(self, h):
        assert h.hash(b"") == 0

    def test_single_byte(self, h):
        # H([c]) = c + 1
        assert h.hash(b"\x00") == 1
        assert h.hash(b"\x01") == 2

    def test_deterministic(self, h):
        data = b"hello world"
        assert h.hash(data) == h.hash(data)

    def test_different_data_different_hash(self, h):
        assert h.hash(b"hello") != h.hash(b"world")

    def test_base_validation(self):
        with pytest.raises(ValueError):
            PolynomialHash(base=0)
        with pytest.raises(ValueError):
            PolynomialHash(base=1)
        with pytest.raises(ValueError):
            PolynomialHash(base=MERSENNE_61)

    def test_hash_concat_theorem1(self, h):
        """Theorem 1: H(A || B) = H(A) * x^|B| + H(B)."""
        a = b"hello "
        b_ = b"world"
        h_a = h.hash(a)
        h_b = h.hash(b_)
        h_ab = h.hash(a + b_)
        assert h.hash_concat(h_a, len(b_), h_b) == h_ab

    def test_hash_concat_empty_right(self, h):
        data = b"test"
        assert h.hash_concat(h.hash(data), 0, 0) == h.hash(data)

    def test_hash_concat_empty_left(self, h):
        data = b"test"
        assert h.hash_concat(0, len(data), h.hash(data)) == h.hash(data)

    def test_hash_repeat_theorem2(self, h):
        """Theorem 2: H(S^q) = H(S) * Phi(q, x^d)."""
        s = b"abc"
        for q in [1, 2, 5, 10, 100]:
            expected = h.hash(s * q)
            computed = h.hash_repeat(h.hash(s), len(s), q)
            assert computed == expected, f"Failed for q={q}"

    def test_hash_overlap_theorem5(self, h):
        """Theorem 5: overlapping back-reference hash."""
        # Pattern "ab" (d=2), total length 7 -> "abababa"
        pattern = b"ab"
        d = len(pattern)
        l = 7
        q, r = divmod(l, d)  # q=3, r=1
        expected = h.hash(b"abababa")
        h_prefix = h.hash(pattern[:r])  # H("a")
        computed = h.hash_overlap(h.hash(pattern), d, l, h_prefix)
        assert computed == expected

    def test_power_caching(self, h):
        # First call computes, second uses cache
        p1 = h.power(100)
        p2 = h.power(100)
        assert p1 == p2
        assert p1 == pow(131, 100, MERSENNE_61)

    def test_different_bases_different_hashes(self):
        h1 = PolynomialHash(base=131)
        h2 = PolynomialHash(base=257)
        data = b"test data"
        assert h1.hash(data) != h2.hash(data)

    def test_mersenne_127(self):
        h = PolynomialHash(prime=MERSENNE_127, base=131)
        assert h.hash(b"test") > 0
        # Concat still works
        a, b_ = b"hello ", b"world"
        assert h.hash_concat(h.hash(a), len(b_), h.hash(b_)) == h.hash(a + b_)
