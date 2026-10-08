//! Sliding window hash rope for streaming hash computation.
//!
//! Provides a generic sliding window over a hash rope without any
//! compression-format dependencies.

use alloc::collections::BTreeMap;

use crate::rope::*;

/// Sliding window over a hash rope.
///
/// Maintains a window of at most `W = d_max + m_max` bytes as a rope,
/// plus a prefix hash covering all evicted bytes.
///
/// [`Self::new`] preserves arena node IDs and retains old arena allocations.
/// Use [`Self::new_bounded`] to reclaim evicted storage while streaming.
pub struct SlidingWindow {
    arena: Arena,
    d_max: u64,
    pub m_max: u64,
    w: u64,
    h_prefix: u64,
    l_prefix: u64,
    r_window: Node,
    pos: u64,
    compact_on_evict: bool,
}

impl SlidingWindow {
    /// Create a sliding window with stable arena node IDs.
    ///
    /// Evicted nodes remain allocated to preserve existing arena handles.
    /// For bounded arena storage, use [`Self::new_bounded`].
    ///
    /// # Panics
    /// Panics if `d_max < 1`, `m_max < 1`, or their sum exceeds `u64::MAX`.
    pub fn new(d_max: u64, m_max: u64, prime: u64, base: u64) -> Self {
        assert!(d_max >= 1, "d_max must be >= 1, got {}", d_max);
        assert!(m_max >= 1, "m_max must be >= 1, got {}", m_max);
        Self {
            arena: Arena::with_hash(prime, base),
            d_max,
            m_max,
            w: d_max
                .checked_add(m_max)
                .expect("window size exceeds u64::MAX"),
            h_prefix: 0,
            l_prefix: 0,
            r_window: None,
            pos: 0,
            compact_on_evict: false,
        }
    }

    /// Create a sliding window that reclaims evicted arena storage.
    ///
    /// Each eviction copies the reachable window graph into a fresh arena,
    /// preserving compressed repeats and shared nodes without expanding them
    /// into bytes. Compaction visits the retained graph and copies leaf data,
    /// releasing old nodes and the old hash power cache.
    /// Arena node IDs observed through [`Self::arena`] may change on append.
    /// Calling [`Self::arena_mut`] permanently disables compaction so that
    /// handles allocated by the caller remain valid.
    ///
    /// # Panics
    /// Panics under the same conditions as [`Self::new`].
    pub fn new_bounded(d_max: u64, m_max: u64, prime: u64, base: u64) -> Self {
        let mut window = Self::new(d_max, m_max, prime, base);
        window.compact_on_evict = true;
        window
    }

    /// Create with default parameters.
    pub fn default_window() -> Self {
        Self::new(32768, 258, crate::polynomial_hash::MERSENNE_61, 131)
    }

    /// Create a reclaiming window with default hash and window parameters.
    pub fn default_bounded_window() -> Self {
        Self::new_bounded(32768, 258, crate::polynomial_hash::MERSENNE_61, 131)
    }

    /// Access the arena.
    pub fn arena(&self) -> &Arena {
        &self.arena
    }

    /// Mutable access to the arena.
    ///
    /// Permanently disables bounded-window compaction to preserve node IDs
    /// allocated by the caller. Subsequent evictions retain arena storage.
    pub fn arena_mut(&mut self) -> &mut Arena {
        self.compact_on_evict = false;
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
    ///
    /// # Panics
    /// Panics if the total stream length exceeds `u64::MAX`.
    pub fn append_bytes(&mut self, data: &[u8]) {
        if data.is_empty() {
            return;
        }
        let new_pos = self
            .pos
            .checked_add(data.len() as u64)
            .expect("stream length exceeds u64::MAX");
        let leaf = self.arena.from_bytes(data);
        self.r_window = self.arena.concat(self.r_window, leaf);
        self.pos = new_pos;
        self.evict();
    }

    /// Append a copy-from-window reference.
    ///
    /// Copies `length` bytes starting `offset` bytes back from the
    /// current end of the window.
    ///
    /// # Panics
    /// Panics if offset exceeds decoded length or window size, offset is zero
    /// for a nonempty copy, or the total stream length exceeds `u64::MAX`.
    pub fn append_copy(&mut self, offset: u64, length: u64) {
        let win_len = self.window_len();

        assert!(
            offset <= self.pos,
            "Invalid copy: offset {} exceeds decoded length {}",
            offset,
            self.pos
        );
        assert!(
            offset <= win_len,
            "Copy offset {} exceeds window size {}. Increase d_max (currently {}) to at least {}.",
            offset,
            win_len,
            self.d_max,
            offset
        );

        // A zero-length copy must not accumulate unreachable split nodes.
        // Keep the offset checks above, including for empty copies.
        if length == 0 {
            return;
        }
        let new_pos = self
            .pos
            .checked_add(length)
            .expect("stream length exceeds u64::MAX");

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

        self.pos = new_pos;
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
        self.h_prefix = self
            .arena
            .hasher_mut()
            .hash_concat(self.h_prefix, old_len, old_hash);

        self.l_prefix = self
            .l_prefix
            .checked_add(excess)
            .expect("evicted length exceeds u64::MAX");
        self.r_window = r_keep;
        if self.compact_on_evict {
            self.compact();
        }
    }

    fn compact(&mut self) {
        fn copy_node(
            old: &Arena,
            new: &mut Arena,
            id: NodeId,
            copied: &mut BTreeMap<NodeId, NodeId>,
        ) -> NodeId {
            if let Some(&new_id) = copied.get(&id) {
                return new_id;
            }
            let node = match old.node(id) {
                NodeInner::Leaf { data, .. } => new.from_bytes(data),
                NodeInner::Internal { left, right, .. } => {
                    let left = copy_node(old, new, *left, copied);
                    let right = copy_node(old, new, *right, copied);
                    new.concat(Some(left), Some(right))
                }
                NodeInner::Repeat { child, reps, .. } => {
                    let child = copy_node(old, new, *child, copied);
                    new.repeat(Some(child), *reps)
                }
            };
            let new_id = node.expect("a nonempty node remains nonempty when copied");
            copied.insert(id, new_id);
            new_id
        }

        let hasher = self.arena.hasher();
        let mut arena = Arena::with_hash(hasher.prime(), hasher.base());
        let root = self
            .r_window
            .map(|id| copy_node(&self.arena, &mut arena, id, &mut BTreeMap::new()));
        self.arena = arena;
        self.r_window = root;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomial_hash::MERSENNE_61;
    use alloc::vec;
    use alloc::vec::Vec;

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
            let a = Arena::with_hash(MERSENNE_61, 131);
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
            let a = Arena::with_hash(MERSENNE_61, 131);
            a.hash_bytes(&full_data)
        };
        assert_eq!(sw.current_hash(), expected);
    }

    #[test]
    fn test_bounded_literal_stream_reclaims_nodes() {
        let mut sw = SlidingWindow::new_bounded(32, 16, MERSENNE_61, 131);
        let mut expected = 0;
        let mut hasher = crate::polynomial_hash::PolynomialHash::default_hash();
        for i in 0..20_000 {
            let byte = (i % 251) as u8;
            sw.append_bytes(&[byte]);
            expected = hasher.hash_concat(expected, 1, hasher.hash(&[byte]));
            assert_eq!(sw.current_hash(), expected);
            assert!(sw.window_len() <= 48);
            assert!(sw.arena().node_count() < 512);
        }
        assert_eq!(sw.pos(), 20_000);
    }

    #[test]
    fn test_bounded_mixed_stream_matches_rebuilt_bytes() {
        let prime = (1 << 31) - 1;
        let mut sw = SlidingWindow::new_bounded(32, 16, prime, 257);
        let mut expected = Vec::new();
        for i in 0..500 {
            let literal = [(i % 251) as u8, 0, 255];
            sw.append_bytes(&literal);
            expected.extend_from_slice(&literal);
            // A non-overlapping copy followed by an overlapping copy.
            for (offset, length) in [(3, 2), (2, 17)] {
                sw.append_copy(offset, length);
                for _ in 0..length {
                    let byte = expected[expected.len() - offset as usize];
                    expected.push(byte);
                }
                assert_eq!(sw.current_hash(), sw.arena().hash_bytes(&expected));
                assert!(sw.arena().node_count() < 256);
            }
        }
        assert_eq!(sw.pos(), expected.len() as u64);
    }

    #[test]
    fn test_bounded_large_literal_releases_evicted_data() {
        let mut sw = SlidingWindow::new_bounded(8, 4, MERSENNE_61, 131);
        let data = vec![123; 256 * 1024];
        sw.append_bytes(&data);
        assert_eq!(sw.current_hash(), sw.arena().hash_bytes(&data));
        assert_eq!(sw.window_len(), 8);
        assert_eq!(sw.arena().node_count(), 1);
        match sw.arena().node(0) {
            NodeInner::Leaf { data, .. } => assert_eq!(data.len(), 8),
            _ => panic!("the retained literal should be a leaf"),
        }
    }

    #[test]
    fn test_bounded_compaction_keeps_large_repeats_compressed() {
        let d_max = 1_u64 << 32;
        let mut sw = SlidingWindow::new_bounded(d_max, 16, MERSENNE_61, 131);
        sw.append_bytes(b"ab");
        sw.append_copy(2, d_max + 30);
        assert_eq!(sw.window_len(), d_max);
        assert!(sw.arena().node_count() < 256);
        let mut hasher = crate::polynomial_hash::PolynomialHash::default_hash();
        let expected = hasher.hash_repeat(hasher.hash(b"ab"), 2, sw.pos() / 2);
        assert_eq!(sw.current_hash(), expected);
    }

    #[test]
    fn test_bounded_compaction_preserves_shared_internal_nodes() {
        let d_max = 1_u64 << 32;
        let mut sw = SlidingWindow::new_bounded(d_max, 1, MERSENNE_61, 131);
        sw.append_bytes(b"a");
        for _ in 0..32 {
            // Non-overlapping doubling builds a DAG whose Internal nodes
            // have identical left/right children, rather than Repeat nodes.
            let length = sw.window_len();
            sw.append_copy(length, length);
        }
        sw.append_copy(1, 2);
        assert_eq!(sw.window_len(), d_max);
        assert!(sw.arena().node_count() < 256);
        let mut hasher = crate::polynomial_hash::PolynomialHash::default_hash();
        let expected = hasher.hash_repeat(hasher.hash(b"a"), 1, sw.pos());
        assert_eq!(sw.current_hash(), expected);
    }

    #[test]
    fn test_mutable_arena_access_preserves_retained_handles() {
        let mut sw = SlidingWindow::new_bounded(8, 4, MERSENNE_61, 131);
        let saved = sw.arena_mut().from_bytes(b"external node");
        let expected = sw.arena().hash(saved);
        for _ in 0..100 {
            sw.append_bytes(b"abc");
            assert_eq!(sw.arena().hash(saved), expected);
            assert_eq!(sw.arena().len(saved), 13);
        }
        assert!(!sw.compact_on_evict);
    }

    #[test]
    fn test_legacy_constructor_preserves_observed_node_ids() {
        let mut sw = SlidingWindow::new(8, 4, MERSENNE_61, 131);
        sw.append_bytes(b"abc");
        let saved = sw.r_window;
        let expected = sw.arena().hash(saved);
        for _ in 0..100 {
            sw.append_bytes(b"def");
        }
        assert_eq!(sw.arena().hash(saved), expected);
        assert!(!sw.compact_on_evict);
        assert!(!SlidingWindow::default_window().compact_on_evict);
        assert!(SlidingWindow::default_bounded_window().compact_on_evict);
    }

    #[test]
    fn test_empty_copies_do_not_accumulate_split_nodes() {
        let mut sw = SlidingWindow::new_bounded(8, 4, MERSENNE_61, 131);
        sw.append_bytes(b"abcdef");
        let allocated = sw.arena().node_count();
        for _ in 0..1_000 {
            sw.append_copy(3, 0);
            sw.append_copy(0, 0);
        }
        assert_eq!(sw.arena().node_count(), allocated);
        assert_eq!(sw.pos(), 6);
        assert_eq!(sw.current_hash(), sw.arena().hash_bytes(b"abcdef"));
    }

    #[test]
    #[should_panic(expected = "exceeds window size")]
    fn test_empty_copy_still_validates_offset() {
        let mut sw = SlidingWindow::new_bounded(8, 4, MERSENNE_61, 131);
        sw.append_bytes(b"abcdefghijklmnop");
        sw.append_copy(10, 0);
    }

    #[test]
    #[should_panic(expected = "window size exceeds u64::MAX")]
    fn test_window_size_overflow() {
        SlidingWindow::new(u64::MAX, 1, MERSENNE_61, 131);
    }

    #[test]
    #[should_panic(expected = "stream length exceeds u64::MAX")]
    fn test_literal_position_overflow() {
        let mut sw = SlidingWindow::default_window();
        sw.pos = u64::MAX;
        sw.append_bytes(b"a");
    }

    #[test]
    #[should_panic(expected = "stream length exceeds u64::MAX")]
    fn test_copy_position_overflow() {
        let mut sw = SlidingWindow::default_window();
        sw.append_bytes(b"a");
        sw.pos = u64::MAX;
        sw.append_copy(1, 1);
    }
}
