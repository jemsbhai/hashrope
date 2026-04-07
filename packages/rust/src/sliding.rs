//! Sliding window hash rope for bounded-memory streaming hash computation.
//!
//! Provides a generic sliding window over a hash rope without any
//! compression-format dependencies.

use crate::polynomial_hash::PolynomialHash;
use crate::rope::*;

/// Bounded-memory sliding window over a hash rope.
///
/// Maintains a window of at most `W = d_max + m_max` bytes as a rope,
/// plus a prefix hash covering all evicted bytes.
pub struct SlidingWindow {
    h: PolynomialHash,
    d_max: u64,
    pub m_max: u64,
    w: u64,
    h_prefix: u64,
    l_prefix: u64,
    r_window: Node,
    pos: u64,
}

impl SlidingWindow {
    /// Create a new sliding window.
    ///
    /// # Panics
    /// Panics if `d_max < 1` or `m_max < 1`.
    pub fn new(d_max: u64, m_max: u64, prime: u64, base: u64) -> Self {
        assert!(d_max >= 1, "d_max must be >= 1, got {}", d_max);
        assert!(m_max >= 1, "m_max must be >= 1, got {}", m_max);
        Self {
            h: PolynomialHash::new(prime, base),
            d_max,
            m_max,
            w: d_max + m_max,
            h_prefix: 0,
            l_prefix: 0,
            r_window: None,
            pos: 0,
        }
    }

    /// Create with default parameters.
    pub fn default_window() -> Self {
        Self::new(32768, 258, crate::polynomial_hash::MERSENNE_61, 131)
    }

    /// Current window size in bytes.
    #[inline]
    pub fn window_len(&self) -> u64 {
        rope_len(&self.r_window)
    }

    /// Total bytes ingested.
    #[inline]
    pub fn pos(&self) -> u64 {
        self.pos
    }

    /// Compute `H(T[0..pos-1])` from sliding state.
    pub fn current_hash(&mut self) -> u64 {
        let wl = self.window_len();
        let wh = rope_hash(&self.r_window);
        self.h.hash_concat(self.h_prefix, wl, wh)
    }

    /// Alias for `current_hash()`.
    pub fn final_hash(&mut self) -> u64 {
        self.current_hash()
    }

    /// Append literal bytes to the stream.
    pub fn append_bytes(&mut self, data: &[u8]) {
        if data.is_empty() {
            return;
        }
        let leaf = rope_from_bytes(data, &self.h);
        self.r_window = rope_concat(&self.r_window, &leaf, &mut self.h);
        self.pos += data.len() as u64;
        self.evict();
    }

    /// Append a copy-from-window reference.
    ///
    /// Copies `length` bytes starting `offset` bytes back from the
    /// current end of the window.
    ///
    /// # Panics
    /// Panics if offset exceeds decoded length or window size.
    pub fn append_copy(&mut self, offset: u64, length: u64) {
        let win_len = self.window_len();

        assert!(
            offset <= self.pos,
            "Invalid copy: offset {} exceeds decoded length {}",
            offset, self.pos
        );
        assert!(
            offset <= win_len,
            "Copy offset {} exceeds window size {}. Increase d_max (currently {}) to at least {}.",
            offset, win_len, self.d_max, offset
        );

        let start = win_len - offset;

        if offset >= length {
            // Non-overlapping
            let (_, tmp) = rope_split(&self.r_window, start, &mut self.h);
            let (source, _) = rope_split(&tmp, length, &mut self.h);
            self.r_window = rope_concat(&self.r_window, &source, &mut self.h);
        } else {
            // Overlapping: extract pattern of length `offset`, repeat
            let (_, tmp) = rope_split(&self.r_window, start, &mut self.h);
            let (pattern, _) = rope_split(&tmp, offset, &mut self.h);

            let q = length / offset;
            let r = length % offset;

            let mut rep = if q >= 1 {
                rope_repeat(&pattern, q, &mut self.h)
            } else {
                None
            };

            if r > 0 {
                let (partial, _) = rope_split(&pattern, r, &mut self.h);
                rep = rope_concat(&rep, &partial, &mut self.h);
            }

            self.r_window = rope_concat(&self.r_window, &rep, &mut self.h);
        }

        self.pos += length;
        self.evict();
    }

    fn evict(&mut self) {
        let win_len = rope_len(&self.r_window);
        if win_len <= self.w {
            return;
        }

        let excess = win_len - self.d_max;
        let (r_old, r_keep) = rope_split(&self.r_window, excess, &mut self.h);

        let old_len = rope_len(&r_old);
        let old_hash = rope_hash(&r_old);
        self.h_prefix = self.h.hash_concat(self.h_prefix, old_len, old_hash);

        self.l_prefix += excess;
        self.r_window = r_keep;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;
    use alloc::vec::Vec;
    use crate::polynomial_hash::MERSENNE_61;

    fn ph() -> PolynomialHash {
        PolynomialHash::default_hash()
    }

    #[test]
    fn test_empty() {
        let mut sw = SlidingWindow::default_window();
        assert_eq!(sw.current_hash(), 0);
        assert_eq!(sw.pos(), 0);
    }

    #[test]
    fn test_append_bytes() {
        let ph = ph();
        let mut sw = SlidingWindow::default_window();
        sw.append_bytes(b"hello");
        assert_eq!(sw.current_hash(), ph.hash(b"hello"));
    }

    #[test]
    fn test_append_incremental() {
        let ph = ph();
        let mut sw = SlidingWindow::default_window();
        sw.append_bytes(b"hello ");
        sw.append_bytes(b"world");
        assert_eq!(sw.current_hash(), ph.hash(b"hello world"));
    }

    #[test]
    fn test_non_overlapping_copy() {
        let ph = ph();
        let mut sw = SlidingWindow::default_window();
        sw.append_bytes(b"hello world");
        sw.append_copy(5, 5);
        assert_eq!(sw.current_hash(), ph.hash(b"hello worldworld"));
    }

    #[test]
    fn test_overlapping_copy() {
        let ph = ph();
        let mut sw = SlidingWindow::default_window();
        sw.append_bytes(b"ab");
        sw.append_copy(2, 6);
        assert_eq!(sw.current_hash(), ph.hash(b"abababab"));
    }

    #[test]
    fn test_eviction_preserves_hash() {
        let ph = ph();
        let mut sw = SlidingWindow::new(8, 4, MERSENNE_61, 131);
        let data = b"the quick brown fox jumps over the lazy dog";
        for &byte in data.iter() {
            sw.append_bytes(&[byte]);
        }
        assert_eq!(sw.current_hash(), ph.hash(data));
    }

    #[test]
    fn test_large_stream() {
        let ph = ph();
        let mut sw = SlidingWindow::new(32, 16, MERSENNE_61, 131);
        let mut full_data = Vec::new();
        for i in 0..50u8 {
            let chunk = vec![i; 10];
            sw.append_bytes(&chunk);
            full_data.extend_from_slice(&chunk);
        }
        assert_eq!(sw.current_hash(), ph.hash(&full_data));
    }
}



