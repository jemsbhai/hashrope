//! Sliding window hash rope for bounded-memory streaming hash computation.
//!
//! Provides a generic sliding window over a hash rope without any
//! compression-format dependencies.

use crate::rope::*;

/// Bounded-memory sliding window over a hash rope.
///
/// Maintains a window of at most `W = d_max + m_max` bytes as a rope,
/// plus a prefix hash covering all evicted bytes.
pub struct SlidingWindow {
    arena: Arena,
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
            arena: Arena::with_hash(prime, base),
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

    /// Access the arena.
    pub fn arena(&self) -> &Arena {
        &self.arena
    }

    /// Mutable access to the arena.
    pub fn arena_mut(&mut self) -> &mut Arena {
        &mut self.arena
    }

    /// Current window size in bytes.
    #[inline]
    pub fn window_len(&self) -> u64 {
        self.arena.len(self.r_window)
    }

    /// Total bytes ingested.
    #[inline]
    pub fn pos(&self) -> u64 {
        self.pos
    }

    /// Compute `H(T[0..pos-1])` from sliding state.
    pub fn current_hash(&mut self) -> u64 {
        let wl = self.window_len();
        let wh = self.arena.hash(self.r_window);
        self.arena.hasher_mut().hash_concat(self.h_prefix, wl, wh)
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
        let leaf = self.arena.from_bytes(data);
        self.r_window = self.arena.concat(self.r_window, leaf);
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
            let (_, tmp) = self.arena.split(self.r_window, start);
            let (source, _) = self.arena.split(tmp, length);
            self.r_window = self.arena.concat(self.r_window, source);
        } else {
            // Overlapping: extract pattern of length `offset`, repeat
            let (_, tmp) = self.arena.split(self.r_window, start);
            let (pattern, _) = self.arena.split(tmp, offset);

            let q = length / offset;
            let r = length % offset;

            let mut rep = if q >= 1 {
                self.arena.repeat(pattern, q)
            } else {
                None
            };

            if r > 0 {
                let (partial, _) = self.arena.split(pattern, r);
                rep = self.arena.concat(rep, partial);
            }

            self.r_window = self.arena.concat(self.r_window, rep);
        }

        self.pos += length;
        self.evict();
    }

    fn evict(&mut self) {
        let win_len = self.arena.len(self.r_window);
        if win_len <= self.w {
            return;
        }

        let excess = win_len - self.d_max;
        let (r_old, r_keep) = self.arena.split(self.r_window, excess);

        let old_len = self.arena.len(r_old);
        let old_hash = self.arena.hash(r_old);
        self.h_prefix = self.arena.hasher_mut().hash_concat(self.h_prefix, old_len, old_hash);

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

    #[test]
    fn test_empty() {
        let mut sw = SlidingWindow::default_window();
        assert_eq!(sw.current_hash(), 0);
        assert_eq!(sw.pos(), 0);
    }

    #[test]
    fn test_append_bytes() {
        let mut sw = SlidingWindow::default_window();
        let expected = sw.arena().hash_bytes(b"hello");
        sw.append_bytes(b"hello");
        assert_eq!(sw.current_hash(), expected);
    }

    #[test]
    fn test_append_incremental() {
        let mut sw = SlidingWindow::default_window();
        let expected = sw.arena().hash_bytes(b"hello world");
        sw.append_bytes(b"hello ");
        sw.append_bytes(b"world");
        assert_eq!(sw.current_hash(), expected);
    }

    #[test]
    fn test_non_overlapping_copy() {
        let mut sw = SlidingWindow::default_window();
        let expected = sw.arena().hash_bytes(b"hello worldworld");
        sw.append_bytes(b"hello world");
        sw.append_copy(5, 5);
        assert_eq!(sw.current_hash(), expected);
    }

    #[test]
    fn test_overlapping_copy() {
        let mut sw = SlidingWindow::default_window();
        let expected = sw.arena().hash_bytes(b"abababab");
        sw.append_bytes(b"ab");
        sw.append_copy(2, 6);
        assert_eq!(sw.current_hash(), expected);
    }

    #[test]
    fn test_eviction_preserves_hash() {
        let mut sw = SlidingWindow::new(8, 4, MERSENNE_61, 131);
        let data = b"the quick brown fox jumps over the lazy dog";
        let expected = {
            let mut a = Arena::with_hash(MERSENNE_61, 131);
            a.hash_bytes(data)
        };
        for &byte in data.iter() {
            sw.append_bytes(&[byte]);
        }
        assert_eq!(sw.current_hash(), expected);
    }

    #[test]
    fn test_large_stream() {
        let mut sw = SlidingWindow::new(32, 16, MERSENNE_61, 131);
        let mut full_data = Vec::new();
        for i in 0..50u8 {
            let chunk = vec![i; 10];
            sw.append_bytes(&chunk);
            full_data.extend_from_slice(&chunk);
        }
        let expected = {
            let mut a = Arena::with_hash(MERSENNE_61, 131);
            a.hash_bytes(&full_data)
        };
        assert_eq!(sw.current_hash(), expected);
    }
}
