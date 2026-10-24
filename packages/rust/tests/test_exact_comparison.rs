//! Exact comparisons use bytes and immutable structure, never fingerprints.
use hashrope::{rope_equal, rope_lcp, Arena, Node, PolynomialHash, LAZY_SENTINEL, MERSENNE_61};

fn byte_lcp(a: &[u8], b: &[u8]) -> u64 {
    a.iter().zip(b).take_while(|(x, y)| x == y).count() as u64
}

#[test]
fn exact_ranges_validate_bounds_and_compare_phases_without_nodes() {
    let mut arena = Arena::new_lazy();
    let pattern = arena.from_bytes(b"ab");
    let repeated = arena.repeat(pattern, 1 << 40);
    let count = arena.node_count();
    assert!(arena.equal_range(repeated, 0, repeated, 2, (1 << 40) - 1));
    assert_eq!(arena.lcp_range(repeated, 0, repeated, 1, 100), 0);
    assert!(arena.equal_range(None, 0, repeated, arena.len(repeated), 0));
    assert_eq!(arena.node_count(), count);
    assert_eq!(arena.hash(repeated), LAZY_SENTINEL);
    use std::panic::{catch_unwind, AssertUnwindSafe};
    assert!(catch_unwind(AssertUnwindSafe(|| arena.equal_range(
        repeated,
        u64::MAX,
        repeated,
        0,
        2
    )))
    .is_err());
    assert!(catch_unwind(AssertUnwindSafe(
        || arena.equal_range(repeated, 0, None, 0, 1)
    ))
    .is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| arena.equal_range(None, 1, None, 0, 0))).is_err());
    let mut other = Arena::with_hash(MERSENNE_61, 257);
    let shifted = other.from_bytes(b"zababx");
    assert!(arena.equal_range_across(repeated, 0, &other, shifted, 1, 4));
    assert_eq!(
        arena.lcp_range_across(repeated, 0, &other, shifted, 1, 5),
        4
    );
}

fn chunked(arena: &mut Arena, bytes: &[u8], chunk_size: usize) -> Node {
    let mut root = None;
    for chunk in bytes.chunks(chunk_size) {
        let leaf = arena.from_bytes(chunk);
        root = arena.concat(root, leaf);
    }
    root
}

#[test]
fn collisions_and_encoded_tokens_are_not_equal() {
    let pairs = [
        (vec![0, 131], vec![1, 0], 263),
        (vec![0, 131, 42], vec![1, 0, 42], 34_496),
        (
            33536u32.to_le_bytes().to_vec(),
            1u32.to_le_bytes().to_vec(),
            4_513_475,
        ),
    ];
    for lazy in [false, true] {
        let mut arena = if lazy {
            Arena::new_lazy()
        } else {
            Arena::new()
        };
        for (left_bytes, right_bytes, fingerprint) in &pairs {
            let left = arena.from_bytes(left_bytes);
            let right = arena.from_bytes(right_bytes);
            if !lazy {
                assert_eq!(arena.hash(left), *fingerprint);
                assert_eq!(arena.hash(right), *fingerprint);
            }
            let count = arena.node_count();
            assert!(!rope_equal(&arena, left, right));
            assert_eq!(rope_lcp(&arena, left, right), 0);
            assert_eq!(arena.node_count(), count);
            if lazy {
                assert_eq!(arena.hash(left), LAZY_SENTINEL);
                assert_eq!(arena.hash(right), LAZY_SENTINEL);
            }
        }
    }
}

#[test]
fn independent_shapes_and_prefix_lengths_match_bytes() {
    let mut arena = Arena::new();
    let bytes: Vec<u8> = (0..2049).map(|i| (i * 37) as u8).collect();
    let whole = arena.from_bytes(&bytes);
    let chunks = chunked(&mut arena, &bytes, 17);
    assert!(arena.equal(whole, chunks));
    assert_eq!(arena.lcp(whole, chunks), bytes.len() as u64);
    assert!(arena.equal(None, None));
    assert!(!arena.equal(None, whole));
    assert_eq!(arena.lcp(None, whole), 0);

    for position in [0, 1, 16, 17, 511, 512, 1024, 2048] {
        let mut changed = bytes.clone();
        changed[position] ^= 1;
        let different = chunked(&mut arena, &changed, 31);
        assert!(!arena.equal(whole, different));
        assert_eq!(arena.lcp(whole, different), position as u64);
        let shorter = arena.from_bytes(&bytes[..position]);
        assert_eq!(arena.lcp(whole, shorter), position as u64);
        assert!(!arena.equal(whole, shorter));
    }
}

#[test]
fn shared_edited_prefix_skips_immutable_ranges() {
    let mut arena = Arena::new_lazy();
    let bytes: Vec<u8> = (0..8192).map(|i| (i * 17) as u8).collect();
    let original = arena.from_bytes(&bytes);
    let split_at = 6151;
    let (prefix, rest) = arena.split(original, split_at);
    let (_, tail) = arena.split(rest, 1);
    let replacement = arena.from_bytes(&[bytes[split_at as usize] ^ 1]);
    let edited = arena.concat(prefix, replacement);
    let edited = arena.concat(edited, tail);
    let count = arena.node_count();
    assert_eq!(arena.lcp(original, edited), split_at);
    assert!(!arena.equal(original, edited));
    assert_eq!(arena.node_count(), count);
    assert_eq!(arena.hash(original), LAZY_SENTINEL);
    assert_eq!(arena.hash(edited), LAZY_SENTINEL);
}

#[test]
fn cross_arena_identity_and_profile_do_not_determine_content() {
    let mut left_arena = Arena::new();
    let mut right_arena = Arena::with_hash(MERSENNE_61, 257);
    let left = left_arena.from_bytes(b"abc");
    let wrong = right_arena.from_bytes(b"abd");
    // Numeric IDs are equal, but denote different content in different arenas.
    assert_eq!(left, wrong);
    assert!(!left_arena.equal_across(left, &right_arena, wrong));
    assert_eq!(left_arena.lcp_across(left, &right_arena, wrong), 2);

    let right = chunked(&mut right_arena, b"abc", 1);
    assert_ne!(left_arena.hash(left), right_arena.hash(right));
    assert!(left_arena.equal_across(left, &right_arena, right));
    assert_eq!(left_arena.lcp_across(left, &right_arena, right), 3);

    let mut lazy = Arena::new_lazy();
    let same = chunked(&mut lazy, b"abc", 2);
    assert!(left_arena.equal_across(left, &lazy, same));
    assert_eq!(lazy.hash(same), LAZY_SENTINEL);
}

#[test]
fn huge_repeats_compare_symbolically_across_patterns_shapes_and_profiles() {
    let mut arena = Arena::new_lazy();
    let ab = arena.from_bytes(b"ab");
    let repetitions = 1 << 48;
    let left = arena.repeat(ab, repetitions);
    let independently_allocated_ab = chunked(&mut arena, b"ab", 1);
    let same = arena.repeat(independently_allocated_ab, repetitions);
    let abab = chunked(&mut arena, b"abab", 1);
    let longer_period = arena.repeat(abab, repetitions / 2);
    assert!(arena.equal(left, same));
    assert!(arena.equal(left, longer_period));
    assert_eq!(arena.lcp(left, longer_period), 2 * repetitions);

    let almost = arena.repeat(ab, repetitions - 1);
    let suffix = arena.from_bytes(b"ac");
    let edited = arena.concat(almost, suffix);
    let count = arena.node_count();
    assert_eq!(arena.lcp(left, edited), 2 * repetitions - 1);
    assert!(!arena.equal(left, edited));
    assert_eq!(arena.node_count(), count);

    let mut other = Arena::with_hash(MERSENNE_61, 257);
    let other_pattern = other.from_bytes(b"abab");
    let other_repeat = other.repeat(other_pattern, repetitions / 2);
    assert!(arena.equal_across(left, &other, other_repeat));
    assert_eq!(arena.hash(left), LAZY_SENTINEL);
}

#[test]
fn huge_repeats_split_across_internal_nodes_and_shifted_phases() {
    let mut arena = Arena::new_lazy();
    let repetitions = 1 << 40;
    let ab = arena.from_bytes(b"ab");
    let left = arena.repeat(ab, repetitions + 1);

    let prefix = arena.from_bytes(b"a");
    let ba = arena.from_bytes(b"ba");
    let middle = arena.repeat(ba, repetitions);
    let suffix = arena.from_bytes(b"b");
    let shifted = arena.concat(prefix, middle);
    let shifted = arena.concat(shifted, suffix);
    assert!(arena.equal(left, shifted));
    assert_eq!(arena.lcp(left, shifted), 2 * (repetitions + 1));

    let fresh_ab = arena.from_bytes(b"ab");
    let first_half = arena.repeat(fresh_ab, repetitions / 2);
    let second_half = arena.repeat(fresh_ab, repetitions / 2 + 1);
    let halves = arena.concat(first_half, second_half);
    assert!(arena.equal(left, halves));
}

#[test]
fn independently_allocated_shared_dags_use_verified_range_cache() {
    let mut arena = Arena::new_lazy();
    let mut left = arena.from_bytes(b"ab");
    let mut right = arena.from_bytes(b"ab");
    for _ in 0..40 {
        left = arena.concat(left, left);
        right = arena.concat(right, right);
    }
    let count = arena.node_count();
    assert!(arena.equal(left, right));
    assert_eq!(arena.lcp(left, right), 2 << 40);
    assert_eq!(arena.node_count(), count);
}

#[test]
fn periodic_comparison_checks_the_last_byte_of_its_witness() {
    let mut arena = Arena::new_lazy();
    let p = arena.from_bytes(b"aab");
    let q = arena.from_bytes(b"aaba");
    let left = arena.repeat(p, 1 << 40);
    let right = arena.repeat(q, 1 << 40);
    // Periods 3 and 4 require 3+4-gcd(3,4)=6 matching bytes.
    assert_eq!(arena.lcp(left, right), 5);
    assert!(!arena.equal(left, right));
}

#[test]
fn nested_repeats_and_period_offsets_match_an_independent_byte_model() {
    // Every binary pattern of length 1..=4, many phase offsets, independent
    // factorizations, and lengths around the Fine-Wilf threshold.
    let mut patterns = Vec::new();
    for length in 1..=4 {
        for bits in 0..1 << length {
            patterns.push(
                (0..length)
                    .map(|i| b'a' + ((bits >> i) & 1) as u8)
                    .collect::<Vec<_>>(),
            );
        }
    }
    for (i, p) in patterns.iter().enumerate() {
        let mut arena = Arena::new_lazy();
        let p_leaf = arena.from_bytes(p);
        let p_twice = arena.repeat(p_leaf, 2);
        let p_many = arena.repeat(p_twice, 9);
        let p_bytes = p.repeat(18);
        for q in patterns.iter().skip(i) {
            let q_leaf = chunked(&mut arena, q, 1);
            let q_many = arena.repeat(q_leaf, 17);
            let q_bytes = q.repeat(17);
            for offset in 0..p.len() {
                let (_, suffix) = arena.split(p_many, offset as u64);
                let expected = byte_lcp(&p_bytes[offset..], &q_bytes);
                assert_eq!(
                    arena.lcp(suffix, q_many),
                    expected,
                    "p={p:?} q={q:?} phase={offset}"
                );
                assert_eq!(arena.equal(suffix, q_many), p_bytes[offset..] == q_bytes);
            }
        }
    }
}

#[test]
fn hash_profile_replacement_is_detected_before_composition() {
    let mut arena = Arena::new();
    let old = arena.from_bytes(b"AB");
    let count = arena.node_count();
    *arena.hasher_mut() = PolynomialHash::new(MERSENNE_61, 257);
    use std::panic::{catch_unwind, AssertUnwindSafe};
    assert!(catch_unwind(AssertUnwindSafe(|| arena.from_bytes(b"CD"))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| arena.concat(old, old))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| arena.split(old, 1))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| arena.repeat(old, 2))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| arena.substr_hash(old, 0, 1))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| arena.ensure_hash(old.unwrap()))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| arena.hash(old))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| arena.hash_bytes(b"AB"))).is_err());
    assert!(catch_unwind(AssertUnwindSafe(|| arena.validate(old))).is_err());
    assert_eq!(arena.node_count(), count);
    assert_eq!(arena.to_bytes(old), b"AB");
    assert_eq!(arena.len(old), 2);

    let mut other = Arena::with_hash(MERSENNE_61, 257);
    let matching = other.from_bytes(b"AB");
    assert!(arena.equal_across(old, &other, matching));
    assert_eq!(arena.lcp_across(old, &other, matching), 2);

    // Replacing the cache while preserving parameters is compatible.
    *arena.hasher_mut() = PolynomialHash::default_hash();
    let tail = arena.from_bytes(b"CD");
    let joined = arena.concat(old, tail);
    assert_eq!(arena.hash(joined), arena.hash_bytes(b"ABCD"));
}

#[test]
fn maximum_u64_lengths_remain_exact_without_materialization() {
    let mut a = Arena::new_lazy();
    let mut b = Arena::with_hash(MERSENNE_61, 257);
    let a_unit = a.from_bytes(b"a");
    let b_unit = b.from_bytes(b"a");
    let left = a.repeat(a_unit, u64::MAX);
    let right = b.repeat(b_unit, u64::MAX);
    assert!(a.equal_across(left, &b, right));
    assert_eq!(a.lcp_across(left, &b, right), u64::MAX);
    assert!(a.equal_range_across(left, u64::MAX, &b, right, u64::MAX, 0));
    assert_eq!(a.lcp_range_across(left, u64::MAX - 4, &b, right, 0, 4), 4);

    let prefix = b.repeat(b_unit, u64::MAX - 1);
    let last = b.from_bytes(b"b");
    let changed = b.concat(prefix, last);
    assert!(!a.equal_across(left, &b, changed));
    assert_eq!(a.lcp_across(left, &b, changed), u64::MAX - 1);
    assert_eq!(
        a.lcp_range_across(left, u64::MAX - 4, &b, changed, u64::MAX - 4, 4),
        3
    );
    assert_eq!(a.hash(left), LAZY_SENTINEL);
}

#[test]
fn sliding_hash_profile_replacement_is_rejected_before_state_changes() {
    use hashrope::SlidingWindow;
    use std::panic::{catch_unwind, AssertUnwindSafe};
    for bounded in [false, true] {
        let mut window = if bounded {
            SlidingWindow::new_bounded(4, 3, MERSENNE_61, 131)
        } else {
            SlidingWindow::new(4, 3, MERSENNE_61, 131)
        };
        window.append_bytes(b"abcdefgh");
        let pos = window.pos();
        let length = window.window_len();
        let count = window.arena().node_count();
        *window.arena_mut().hasher_mut() = PolynomialHash::new(MERSENNE_61, 257);
        assert!(catch_unwind(AssertUnwindSafe(|| window.append_bytes(b""))).is_err());
        assert!(catch_unwind(AssertUnwindSafe(|| window.append_bytes(b"x"))).is_err());
        assert!(catch_unwind(AssertUnwindSafe(|| window.append_copy(1, 0))).is_err());
        assert!(catch_unwind(AssertUnwindSafe(|| window.append_copy(1, 2))).is_err());
        assert!(catch_unwind(AssertUnwindSafe(|| window.current_hash())).is_err());
        assert_eq!(window.pos(), pos);
        assert_eq!(window.window_len(), length);
        assert_eq!(window.arena().node_count(), count);
        *window.arena_mut().hasher_mut() = PolynomialHash::default_hash();
        assert_eq!(
            window.current_hash(),
            PolynomialHash::default_hash().hash(b"abcdefgh")
        );
        window.append_copy(1, 2);
        assert_eq!(
            window.current_hash(),
            PolynomialHash::default_hash().hash(b"abcdefghhh")
        );
    }
}
