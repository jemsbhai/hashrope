//! Exact comparisons preserve arbitrary-precision lengths without trusting hashes.
#![cfg(feature = "bigint")]

use hashrope::bigint::{
    rope_equal, rope_lcp, rope_range_equal, rope_range_lcp, Arena, BigUint, PolynomialHash,
    SlidingWindow, MERSENNE_127, MERSENNE_61,
};
use std::panic::{catch_unwind, AssertUnwindSafe};

fn n(value: usize) -> BigUint {
    value.into()
}

fn byte_lcp(left: &[u8], right: &[u8]) -> usize {
    left.iter().zip(right).take_while(|(a, b)| a == b).count()
}

#[test]
fn default_byte_and_token_collisions_do_not_decide_equality() {
    let pairs = [
        (vec![0, 131], vec![1, 0], n(263)),
        (vec![0, 131, 42], vec![1, 0, 42], n(34_496)),
        (
            33536u32.to_le_bytes().to_vec(),
            1u32.to_le_bytes().to_vec(),
            n(4_513_475),
        ),
    ];
    for (left, right, fingerprint) in pairs {
        let mut arena = Arena::new();
        let a = arena.from_bytes(&left);
        let b = arena.from_bytes_chunked(&right, 1);
        assert_eq!(arena.hash(a), fingerprint);
        assert_eq!(arena.hash(b), fingerprint);
        let count = arena.node_count();
        assert!(!rope_equal(&arena, a, b));
        assert_eq!(rope_lcp(&arena, a, b), n(0));
        assert!(!rope_range_equal(
            &arena,
            a,
            &n(0),
            b,
            &n(0),
            &n(left.len())
        ));
        assert_eq!(
            rope_range_lcp(&arena, a, &n(0), b, &n(0), &n(left.len())),
            n(0)
        );
        assert_eq!(arena.node_count(), count);
        assert_eq!(arena.hash(a), fingerprint);
    }
}

#[test]
fn cross_profile_and_cross_arena_ids_are_not_identity() {
    let mut a = Arena::new();
    let mut b = Arena::with_hash(MERSENNE_127.into(), n(257));
    let x = a.from_bytes(b"abcdefgh");
    let y = b.from_bytes(b"abcdefgX");
    assert_eq!(x, y); // Numeric IDs overlap, but belong to different arenas.
    assert!(!a.equal_across(x, &b, y));
    assert_eq!(a.lcp_across(x, &b, y), n(7));
    let y = b.from_bytes_chunked(b"abcdefgh", 3);
    assert!(a.equal_across(x, &b, y));
    assert_eq!(a.lcp_across(x, &b, y), n(8));
    assert_ne!(a.hash(x), b.hash(y));
    assert!(a.equal_range_across(x, &n(2), &b, y, &n(2), &n(6)));
    assert_eq!(a.lcp_range_across(x, &n(1), &b, y, &n(2), &n(4)), n(0));
}

#[test]
fn different_shapes_and_all_small_ranges_match_plain_bytes() {
    let bytes: Vec<u8> = (0..97).map(|i| (i * 37 + 11) as u8).collect();
    let mut left = Arena::new();
    let mut right = Arena::with_hash(n(127), n(31));
    let a = left.from_bytes_chunked(&bytes, 11);
    let b = right.from_bytes_chunked(&bytes, 7);
    assert!(left.equal_across(a, &right, b));
    let counts = (left.node_count(), right.node_count());
    for start_a in (0..=bytes.len()).step_by(3) {
        for start_b in (0..=bytes.len()).step_by(5) {
            let length = (bytes.len() - start_a).min(bytes.len() - start_b);
            let expected = byte_lcp(
                &bytes[start_a..start_a + length],
                &bytes[start_b..start_b + length],
            );
            assert_eq!(
                left.lcp_range_across(a, &n(start_a), &right, b, &n(start_b), &n(length)),
                n(expected)
            );
            assert_eq!(
                left.equal_range_across(a, &n(start_a), &right, b, &n(start_b), &n(length)),
                expected == length
            );
        }
    }
    assert_eq!((left.node_count(), right.node_count()), counts);
}

#[test]
fn unequal_leaf_word_and_tail_positions_and_offsets_are_exact() {
    // Every position in word-boundary/remainder lengths, with unaligned ranges.
    for length in [
        1, 2, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129, 511, 512,
    ] {
        for prefix in [0, 1, 7] {
            let bytes = vec![42u8; prefix + length + 3];
            let mut arena = Arena::with_hash(n(127), n(31));
            let a = arena.from_bytes_chunked(&bytes, bytes.len());
            for mismatch in 0..length {
                let mut changed = bytes.clone();
                changed[prefix + mismatch] = 43;
                let b = arena.from_bytes_chunked(&changed, changed.len());
                assert_eq!(
                    arena.lcp_range(a, &n(prefix), b, &n(prefix), &n(length)),
                    n(mismatch)
                );
            }
        }
    }
}

#[test]
fn empty_and_identity_shortcuts_still_check_wide_bounds() {
    let mut arena = Arena::new();
    let root = arena.from_bytes(b"abc");
    assert!(arena.equal(None, None));
    assert_eq!(arena.lcp(None, root), n(0));
    assert!(!arena.equal(None, root));
    assert!(arena.equal_range(None, &n(0), root, &n(3), &n(0)));
    let huge = n(1) << 1024usize;
    for (start, length) in [
        (n(4), n(0)),
        (n(3), n(1)),
        (huge.clone(), n(0)),
        (n(0), huge),
    ] {
        assert!(catch_unwind(|| arena.lcp_range(root, &start, root, &start, &length)).is_err());
    }
    assert!(catch_unwind(|| arena.equal_range(None, &n(1), None, &n(1), &n(0))).is_err());
    assert!(catch_unwind(|| arena.equal_range(root, &n(0), None, &n(0), &n(1))).is_err());
}

#[test]
fn thousand_bit_repeats_and_shifted_ranges_remain_symbolic() {
    let copies = n(1) << 1300usize;
    let mut a = Arena::with_hash(n(127), n(31));
    let mut b = Arena::with_hash(MERSENNE_127.into(), n(257));
    let ab = a.from_bytes(b"ab");
    let abab = b.from_bytes_chunked(b"abab", 1);
    let x = a.repeat(ab, &copies);
    let y = b.repeat(abab, &(&copies >> 1usize));
    let expected_len = &copies * 2u8;
    let before = (a.node_count(), b.node_count());
    assert!(a.equal_across(x, &b, y));
    assert_eq!(a.lcp_across(x, &b, y), expected_len);
    let start = (n(1) << 1100usize) + 1u8;
    let length = &copies - 3u8;
    assert!(a.equal_range_across(x, &start, &b, y, &(&start + 2u8), &length));
    assert_eq!(
        a.lcp_range_across(x, &start, &b, y, &(&start + 1u8), &length),
        n(0)
    );
    assert_eq!((a.node_count(), b.node_count()), before);
}

#[test]
fn periodic_witness_checks_its_last_byte_and_independent_phases() {
    let mut arena = Arena::new();
    let a = arena.from_bytes(b"aab");
    let b = arena.from_bytes(b"aaba");
    let a = arena.repeat(a, &(n(1) << 100usize));
    let b = arena.repeat(b, &(n(1) << 100usize));
    assert_eq!(arena.lcp(a, b), n(5)); // p + q - gcd(p,q) == 6.
    let patterns: Vec<Vec<u8>> = (1..=3)
        .flat_map(|len| {
            (0..1usize << len)
                .map(move |bits| (0..len).map(|i| b'a' + ((bits >> i) & 1) as u8).collect())
        })
        .collect();
    for x in &patterns {
        for y in &patterns {
            let raw_a = x.repeat(24);
            let raw_b = y.repeat(24);
            let a = arena.from_bytes_chunked(x, 1);
            let b = arena.from_bytes_chunked(y, 2);
            let a = arena.repeat(a, &n(24));
            let b = arena.repeat(b, &n(24));
            for (offset_a, offset_b) in [(0, 0), (1, 0), (0, 1), (1, 2)] {
                let length = (raw_a.len() - offset_a).min(raw_b.len() - offset_b);
                let expected = byte_lcp(
                    &raw_a[offset_a..offset_a + length],
                    &raw_b[offset_b..offset_b + length],
                );
                assert_eq!(
                    arena.lcp_range(a, &n(offset_a), b, &n(offset_b), &n(length)),
                    n(expected)
                );
            }
        }
    }
}

#[test]
fn repeat_against_internal_and_huge_late_edit_do_not_expand() {
    let copies = n(1) << 512usize;
    let mut arena = Arena::with_hash(n(127), n(31));
    let ab = arena.from_bytes(b"ab");
    let x = arena.repeat(ab, &copies);
    let half = arena.repeat(ab, &(&copies >> 1usize));
    let y = arena.concat(half, half);
    assert!(arena.equal(x, y));
    let length = arena.len(x);
    let prefix = arena.split(x, &(&length - 1u8)).0;
    let final_byte = arena.from_bytes(b"X");
    let changed = arena.concat(prefix, final_byte);
    let before = arena.node_count();
    assert_eq!(arena.lcp(x, changed), &length - 1u8);
    assert!(!arena.equal(x, changed));
    assert_eq!(arena.node_count(), before);
}

#[test]
fn independent_self_concat_dags_and_partial_ranges_keep_wide_lengths() {
    let mut a = Arena::with_hash(n(127), n(31));
    let mut b = Arena::with_hash(n(127), n(3));
    let mut x = a.from_bytes(b"abc");
    let mut y = b.from_bytes_chunked(b"abc", 1);
    for _ in 0..1100 {
        x = a.concat(x, x);
        y = b.concat(y, y);
    }
    let length = a.len(x);
    let start = (n(1) << 1000usize) + 2u8;
    let range_length = &length - &start - 3u8;
    let before = (a.node_count(), b.node_count());
    assert!(a.equal_across(x, &b, y));
    assert_eq!(a.lcp_across(x, &b, y), length);
    assert!(a.equal_range_across(x, &start, &b, y, &(&start + 3u8), &range_length));
    assert_eq!(
        a.lcp_range_across(x, &start, &b, y, &(&start + 1u8), &range_length),
        n(0)
    );
    assert_eq!((a.node_count(), b.node_count()), before);
}

#[test]
fn deep_split_dags_use_heap_frames_and_verified_ranges() {
    fn build() -> (Arena, hashrope::bigint::Node) {
        let mut arena = Arena::with_hash(n(127), n(31));
        let pattern = arena.from_bytes_chunked(b"abcde", 1);
        let count = (n(1) << 1024usize) + 17u8;
        let repeated = arena.repeat(pattern, &count);
        let suffix = arena.split(repeated, &n(1)).1;
        assert!(arena.height(suffix) > 1000);
        (arena, suffix)
    }
    let (a, x) = build();
    let (b, y) = build();
    let before = (a.node_count(), b.node_count());
    assert_eq!(a.lcp_across(x, &b, y), a.len(x));
    assert!(a.equal_across(x, &b, y));
    assert_eq!((a.node_count(), b.node_count()), before);
}

#[test]
fn profile_replacement_fails_before_metadata_or_structure_changes() {
    let mut arena = Arena::new();
    let a = arena.from_bytes(b"immutable bytes");
    let b = arena.from_bytes_chunked(b"immutable bytes", 3);
    let old_hash = arena.hash(a);
    let count = arena.node_count();
    *arena.hasher_mut() = PolynomialHash::new(MERSENNE_127.into(), n(257));
    assert!(arena.equal(a, b));
    assert_eq!(arena.lcp(a, b), n(15));
    assert_eq!(arena.to_bytes(a), b"immutable bytes");
    assert!(catch_unwind(|| arena.hasher()).is_err());
    assert!(catch_unwind(|| arena.hash(a)).is_err());
    assert!(catch_unwind(|| arena.hash(None)).is_err());
    assert!(catch_unwind(|| arena.hash_bytes(b"literal")).is_err());
    assert!(catch_unwind(|| arena.validate(a)).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| arena.from_bytes(b""))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| arena.concat(None, a))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| arena.split(a, &n(0)))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| arena.repeat(a, &n(0)))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| arena.substr_hash(a, &n(0), &n(0)))).is_err());
    assert_eq!(arena.node_count(), count);
    *arena.hasher_mut() = PolynomialHash::default_hash();
    assert_eq!(arena.hash(a), old_hash);
    arena.validate(a);
    let restored = arena.concat(a, b);
    assert_eq!(arena.len(restored), n(30));
}

#[test]
fn sliding_profile_replacement_cannot_corrupt_prefix_or_window() {
    let mut window = SlidingWindow::new_bounded(n(4), n(2), MERSENNE_61.into(), n(131));
    window.append_bytes(b"abcdefghij"); // Includes an evicted prefix.
    let expected_hash = window.current_hash();
    let state = (
        window.pos(),
        window.window_len(),
        window.arena().node_count(),
    );
    *window.arena_mut().hasher_mut() = PolynomialHash::new(MERSENNE_61.into(), n(257));
    assert!(catch_unwind(AssertUnwindSafe(|| window.append_bytes(b""))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| window.append_bytes(b"k"))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| window.append_copy(&n(0), &n(0)))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| window.append_copy(&n(1), &n(4)))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| window.current_hash())).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| window.final_hash())).is_err());
    assert_eq!(
        (
            window.pos(),
            window.window_len(),
            window.arena().node_count()
        ),
        state
    );
    *window.arena_mut().hasher_mut() = PolynomialHash::default_hash();
    assert_eq!(window.current_hash(), expected_hash);
    window.append_bytes(b"k");
    assert_eq!(
        window.current_hash(),
        PolynomialHash::default_hash().hash(b"abcdefghijk")
    );
}
