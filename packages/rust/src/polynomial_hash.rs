//! Polynomial hash over Z/pZ with Mersenne prime arithmetic.
//!
//! Implements Definitions 1-3, Lemmas 1-2, 13-14, Theorems 1-3, 5
//! from the compressed-domain hashing framework.
//!
//! All arithmetic is performed modulo a Mersenne prime p = 2^61 - 1
//! using bit-shift reduction (Lemma 14) to avoid expensive division.
//!
//! Power cache uses a fixed inline array for exponents 0..255 (O(1) lookup)
//! with BTreeMap fallback for larger exponents. This is `no_std + alloc`
//! compatible.


use alloc::collections::BTreeMap;

/// 2^61 - 1
pub const MERSENNE_61: u64 = (1u64 << 61) - 1;

/// Number of precomputed powers in the inline array.
const INLINE_CACHE_SIZE: usize = 256;

// ---------------------------------------------------------------------------
// Mersenne prime arithmetic (Lemma 14)
// ---------------------------------------------------------------------------

/// Reduce a non-negative integer `a` modulo Mersenne prime `p = 2^k - 1`.
///
/// Uses bit-shift reduction: `a mod (2^k - 1) = (a >> k) + (a & p)`
/// with at most one conditional subtraction.
#[inline]
pub fn mersenne_mod(a: u128, p: u64) -> u64 {
    let k = 64 - p.leading_zeros(); // p.bit_length()
    let mask = p as u128;

    // Fold
    let mut v = (a >> k) + (a & mask);
    // One more fold (handles overflow from first fold)
    v = (v >> k) + (v & mask);
    // Conditional subtraction
    let r = v as u64;
    if r >= p {
        r - p
    } else {
        r
    }
}

/// Compute `(a * b) mod p` for Mersenne prime `p`.
#[inline]
pub fn mersenne_mul(a: u64, b: u64, p: u64) -> u64 {
    mersenne_mod(a as u128 * b as u128, p)
}

// ---------------------------------------------------------------------------
// Geometric accumulator Phi (Definition 3, Theorem 3)
// ---------------------------------------------------------------------------

/// Compute `Phi(q, alpha) = sum_{i=0}^{q-1} alpha^i (mod p)` via
/// inverse-free repeated doubling.
///
/// Recurrence (Theorem 3):
///   - Phi(0, alpha) = 0
///   - Phi(1, alpha) = 1
///   - Phi(2k, alpha) = Phi(k, alpha) * (1 + alpha^k)
///   - Phi(2k+1, alpha) = Phi(2k, alpha) * alpha + 1
///
/// Time: O(log q) multiplications in Z/pZ.
pub fn phi(q: u64, alpha: u64, p: u64) -> u64 {
    if q == 0 {
        return 0;
    }
    if q == 1 {
        return 1;
    }

    let bits = 64 - q.leading_zeros(); // q.bit_length()

    let mut phi_val: u64 = 1;
    let mut alpha_power: u64 = alpha;

    // Scan from second-most-significant bit down to bit 0
    for i in (0..bits - 1).rev() {
        // Even step: Phi(2k) = Phi(k) * (1 + alpha^k)
        let sum = mersenne_mod(1u128 + alpha_power as u128, p);
        phi_val = mersenne_mul(phi_val, sum, p);
        // alpha^k -> alpha^(2k)
        alpha_power = mersenne_mul(alpha_power, alpha_power, p);

        // Odd step if bit is set: Phi(2k+1) = Phi(2k) * alpha + 1
        if (q >> i) & 1 == 1 {
            phi_val = mersenne_mod(
                mersenne_mul(phi_val, alpha, p) as u128 + 1,
                p,
            );
            alpha_power = mersenne_mul(alpha_power, alpha, p);
        }
    }

    phi_val
}

// ---------------------------------------------------------------------------
// Polynomial Hash (Definition 2, Theorems 1-2, 5)
// ---------------------------------------------------------------------------

/// Polynomial hash function over Z/pZ.
///
/// `H(s) = sum_{i=0}^{n-1} (s_i + 1) * x^(n-1-i) (mod p)`
///
/// The +1 offset (Definition 2) ensures no byte maps to zero in Z/pZ,
/// preventing the zero-padding collision (Lemma 1).
///
/// Power cache: inline `[u64; 256]` for exponents 0..255 (O(1) lookup),
/// `BTreeMap` fallback for larger exponents.
pub struct PolynomialHash {
    p: u64,
    x: u64,
    inline_cache: [u64; INLINE_CACHE_SIZE],
    overflow_cache: BTreeMap<u64, u64>,
}

impl PolynomialHash {
    /// Create a new polynomial hash function.
    ///
    /// Precomputes `x^0` through `x^255` in the inline cache.
    ///
    /// # Panics
    /// Panics if `base < 2` or `base >= prime`.
    pub fn new(prime: u64, base: u64) -> Self {
        assert!(base >= 2 && base < prime, "Base must be in [2, p-1]");
        let mut inline_cache = [0u64; INLINE_CACHE_SIZE];
        inline_cache[0] = 1;
        for i in 1..INLINE_CACHE_SIZE {
            inline_cache[i] = mersenne_mul(inline_cache[i - 1], base, prime);
        }
        Self {
            p: prime,
            x: base,
            inline_cache,
            overflow_cache: BTreeMap::new(),
        }
    }

    /// Create with default parameters (p = 2^61 - 1, base = 131).
    pub fn default_hash() -> Self {
        Self::new(MERSENNE_61, 131)
    }

    /// The prime modulus.
    #[inline]
    pub fn prime(&self) -> u64 {
        self.p
    }

    /// The hash base.
    #[inline]
    pub fn base(&self) -> u64 {
        self.x
    }

    /// Compute `x^n mod p` with caching.
    ///
    /// O(1) for n < 256 (inline array), O(log n) for cache misses on larger n.
    #[inline]
    pub fn power(&mut self, n: u64) -> u64 {
        if (n as usize) < INLINE_CACHE_SIZE {
            return self.inline_cache[n as usize];
        }
        if let Some(&v) = self.overflow_cache.get(&n) {
            return v;
        }
        // Binary exponentiation
        let mut result = 1u64;
        let mut base = self.x;
        let mut exp = n;
        let p = self.p;
        while exp > 0 {
            if exp & 1 == 1 {
                result = mersenne_mul(result, base, p);
            }
            base = mersenne_mul(base, base, p);
            exp >>= 1;
        }
        self.overflow_cache.insert(n, result);
        result
    }

    /// Compute `H(data)`.
    ///
    /// For the empty slice, returns 0 (Definition 2).
    pub fn hash(&self, data: &[u8]) -> u64 {
        if data.is_empty() {
            return 0;
        }
        let mut h: u64 = 0;
        for &byte in data {
            // h = h * x + (byte + 1)
            h = mersenne_mod(h as u128 * self.x as u128 + byte as u128 + 1, self.p);
        }
        h
    }

    /// Compute `H(A || B)` from `H(A)`, `|B|`, `H(B)` via Theorem 1.
    ///
    /// `H(A || B) = H(A) * x^|B| + H(B) (mod p)`
    pub fn hash_concat(&mut self, h_a: u64, len_b: u64, h_b: u64) -> u64 {
        let x_pow = self.power(len_b);
        mersenne_mod(h_a as u128 * x_pow as u128 + h_b as u128, self.p)
    }

    /// Compute `H(S^q)` from `H(S)` and `|S| = d` via Theorem 2.
    ///
    /// `H(S^q) = H(S) * Phi(q, x^d) (mod p)`
    pub fn hash_repeat(&mut self, h_s: u64, d: u64, q: u64) -> u64 {
        let x_d = self.power(d);
        let phi_val = phi(q, x_d, self.p);
        mersenne_mul(h_s, phi_val, self.p)
    }

    /// Compute `H(W)` for overlapping back-reference via Theorem 5.
    ///
    /// `W = P^q || P[0..r-1]` where `q = l / d`, `r = l % d`.
    pub fn hash_overlap(&mut self, h_p: u64, d: u64, l: u64, h_prefix: u64) -> u64 {
        let q = l / d;
        let r = l % d;
        let x_d = self.power(d);
        let x_r = self.power(r);

        let phi_val = phi(q, x_d, self.p);
        let h_repeated = mersenne_mul(h_p, phi_val, self.p);

        mersenne_mod(h_repeated as u128 * x_r as u128 + h_prefix as u128, self.p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec::Vec;

    #[test]
    fn test_mersenne_mod_basic() {
        assert_eq!(mersenne_mod(0, MERSENNE_61), 0);
        assert_eq!(mersenne_mod(1, MERSENNE_61), 1);
        assert_eq!(mersenne_mod(MERSENNE_61 as u128, MERSENNE_61), 0);
        assert_eq!(mersenne_mod(MERSENNE_61 as u128 + 1, MERSENNE_61), 1);
    }

    #[test]
    fn test_mersenne_mod_large() {
        let val = MERSENNE_61 as u128 * 1000 + 42;
        assert_eq!(mersenne_mod(val, MERSENNE_61), 42);
    }

    #[test]
    fn test_mersenne_mul_commutative() {
        let a = 12345u64;
        let b = 67890u64;
        assert_eq!(mersenne_mul(a, b, MERSENNE_61), mersenne_mul(b, a, MERSENNE_61));
    }

    #[test]
    fn test_phi_base_cases() {
        assert_eq!(phi(0, 131, MERSENNE_61), 0);
        assert_eq!(phi(1, 131, MERSENNE_61), 1);
    }

    #[test]
    fn test_phi_two() {
        let alpha = 131u64;
        let expected = mersenne_mod(1 + alpha as u128, MERSENNE_61);
        assert_eq!(phi(2, alpha, MERSENNE_61), expected);
    }

    #[test]
    fn test_phi_degenerate_alpha_one() {
        assert_eq!(phi(100, 1, MERSENNE_61), 100);
        assert_eq!(phi(1000, 1, MERSENNE_61), 1000);
    }

    #[test]
    fn test_phi_large_q() {
        let alpha = 5u64;
        let q = 20u64;
        let mut expected = 0u64;
        let mut power = 1u64;
        for _ in 0..q {
            expected = mersenne_mod(expected as u128 + power as u128, MERSENNE_61);
            power = mersenne_mul(power, alpha, MERSENNE_61);
        }
        assert_eq!(phi(q, alpha, MERSENNE_61), expected);
    }

    #[test]
    fn test_hash_empty() {
        let h = PolynomialHash::default_hash();
        assert_eq!(h.hash(b""), 0);
    }

    #[test]
    fn test_hash_deterministic() {
        let h = PolynomialHash::default_hash();
        assert_eq!(h.hash(b"hello world"), h.hash(b"hello world"));
    }

    #[test]
    fn test_hash_different_data() {
        let h = PolynomialHash::default_hash();
        assert_ne!(h.hash(b"hello"), h.hash(b"world"));
    }

    #[test]
    fn test_hash_concat_theorem1() {
        let mut h = PolynomialHash::default_hash();
        let a = b"hello ";
        let b_ = b"world";
        let h_a = h.hash(a);
        let h_b = h.hash(b_);
        let h_ab = h.hash(b"hello world");
        assert_eq!(h.hash_concat(h_a, b_.len() as u64, h_b), h_ab);
    }

    #[test]
    fn test_hash_repeat_theorem2() {
        let mut h = PolynomialHash::default_hash();
        let s = b"abc";
        for q in [1, 2, 5, 10, 100] {
            let repeated: Vec<u8> = s.iter().copied().cycle().take(s.len() * q).collect();
            let expected = h.hash(&repeated);
            let computed = h.hash_repeat(h.hash(s), s.len() as u64, q as u64);
            assert_eq!(computed, expected, "Failed for q={}", q);
        }
    }

    #[test]
    fn test_hash_overlap_theorem5() {
        let mut h = PolynomialHash::default_hash();
        let pattern = b"ab";
        let d = pattern.len() as u64;
        let l = 7u64;
        let r = l % d;
        let expected = h.hash(b"abababa");
        let h_prefix = h.hash(&pattern[..r as usize]);
        let computed = h.hash_overlap(h.hash(pattern), d, l, h_prefix);
        assert_eq!(computed, expected);
    }

    #[test]
    fn test_python_cross_validation() {
        // Exact values from Python: PolynomialHash(base=131, p=2^61-1)
        let h = PolynomialHash::default_hash();
        assert_eq!(h.hash(b"hello world"), 430229793999670395);
        let ab100: Vec<u8> = b"ab".iter().copied().cycle().take(200).collect();
        assert_eq!(h.hash(&ab100), 463276214322025987);
        assert_eq!(h.hash(b"\x00"), 1);
    }

    #[test]
    fn test_inline_cache_correctness() {
        // Verify inline cache matches binary exponentiation
        let mut h = PolynomialHash::default_hash();
        for n in 0..256u64 {
            let cached = h.power(n);
            let computed = {
                let mut r = 1u64;
                let mut b = h.base();
                let mut e = n;
                let p = h.prime();
                while e > 0 {
                    if e & 1 == 1 { r = mersenne_mul(r, b, p); }
                    b = mersenne_mul(b, b, p);
                    e >>= 1;
                }
                r
            };
            assert_eq!(cached, computed, "Mismatch at n={}", n);
        }
    }

    #[test]
    fn test_overflow_cache() {
        // Verify powers beyond 255 are computed and cached correctly
        let mut h = PolynomialHash::default_hash();
        let p1 = h.power(1000);
        let p2 = h.power(1000);
        assert_eq!(p1, p2);
        assert_ne!(p1, 0);
    }
}

