//! Regression oracles for optimized comparison traversal.
//!
//! The large cases have tiny stored DAGs but enormous expanded content. They
//! must retain structural/repeat shortcuts even when ordinary scans avoid the
//! general comparison cache. Expected results never expand those DAGs.
use hashrope::{Arena, Node, LAZY_SENTINEL, MERSENNE_61};

fn byte_lcp(left: &[u8], right: &[u8]) -> u64 {
    left.iter().zip(right).take_while(|(a, b)| a == b).count() as u64
}

fn chunked(arena: &mut Arena, bytes: &[u8], chunk_size: usize) -> Node {
    let mut root = None;
    for chunk in bytes.chunks(chunk_size) {
        let next = arena.from_bytes(chunk);
        root = arena.concat(root, next);
    }
    root
}

fn next_random(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state
}

fn edited(arena: &mut Arena, root: Node, at: u64, replacement: u8) -> Node {
    let (prefix, rest) = arena.split(root, at);
    let (_, suffix) = arena.split(rest, 1);
    let replacement = arena.from_bytes(&[replacement]);
    let prefix = arena.concat(prefix, replacement);
    arena.concat(prefix, suffix)
}

#[test]
fn seeded_independent_shared_and_range_comparisons_match_bytes() {
    let mut random = 0x98f2_781a_624c_d005;
    for case in 0..32 {
        let size = 513 + (next_random(&mut random) % 3072) as usize;
        let bytes: Vec<u8> = (0..size)
            .map(|_| (next_random(&mut random) >> 32) as u8)
            .collect();
        let mut arena = if case % 2 == 0 {
            Arena::new_lazy()
        } else {
            Arena::new()
        };
        let original = arena.from_bytes(&bytes);
        let independent = chunked(
            &mut arena,
            &bytes,
            1 + (next_random(&mut random) % 97) as usize,
        );
        let mut other = if case % 3 == 0 {
            Arena::new_lazy()
        } else {
            Arena::with_hash(MERSENNE_61, 257)
        };
        let cross_arena = chunked(
            &mut other,
            &bytes,
            1 + (next_random(&mut random) % 113) as usize,
        );
        assert!(arena.equal(original, independent));
        assert_eq!(arena.lcp(original, independent), size as u64);
        assert!(arena.equal_across(original, &other, cross_arena));

        // Exercise first/last bytes, the normal construction leaf boundary,
        // arbitrary middle bytes, and shared edited trees.
        for position in [0, 1, 511, 512, size / 2, size - 1] {
            let mut changed_bytes = bytes.clone();
            changed_bytes[position] ^= 0x80;
            let shared = edited(
                &mut arena,
                original,
                position as u64,
                changed_bytes[position],
            );
            let independent_change = chunked(&mut other, &changed_bytes, 37);
            assert_eq!(arena.lcp(original, shared), position as u64);
            assert!(!arena.equal(original, shared));
            assert_eq!(
                arena.lcp_across(independent, &other, independent_change),
                position as u64
            );
            assert!(!arena.equal_across(independent, &other, independent_change));

            for _ in 0..8 {
                let left_start = (next_random(&mut random) % (size as u64 + 1)) as usize;
                let right_start = (next_random(&mut random) % (size as u64 + 1)) as usize;
                let limit = (size - left_start).min(size - right_start);
                let length = (next_random(&mut random) % (limit as u64 + 1)) as usize;
                let left_bytes = &bytes[left_start..left_start + length];
                let right_bytes = &changed_bytes[right_start..right_start + length];
                let expected = byte_lcp(left_bytes, right_bytes);
                assert_eq!(
                    arena.lcp_range(
                        original,
                        left_start as u64,
                        shared,
                        right_start as u64,
                        length as u64
                    ),
                    expected
                );
                assert_eq!(
                    arena.equal_range(
                        original,
                        left_start as u64,
                        shared,
                        right_start as u64,
                        length as u64
                    ),
                    left_bytes == right_bytes
                );
                assert_eq!(
                    arena.lcp_range_across(
                        independent,
                        left_start as u64,
                        &other,
                        independent_change,
                        right_start as u64,
                        length as u64
                    ),
                    expected
                );
                assert_eq!(
                    arena.equal_range_across(
                        independent,
                        left_start as u64,
                        &other,
                        independent_change,
                        right_start as u64,
                        length as u64
                    ),
                    left_bytes == right_bytes
                );
            }
        }

        let count = arena.node_count();
        let other_count = other.node_count();
        // Equal shifted windows force traversal beyond the first leaf.
        for start in [1, 17, 511, size / 2, size] {
            let length = size - start;
            assert!(arena.equal_range(
                original,
                start as u64,
                independent,
                start as u64,
                length as u64
            ));
            assert_eq!(
                arena.lcp_range_across(
                    original,
                    start as u64,
                    &other,
                    cross_arena,
                    start as u64,
                    length as u64
                ),
                length as u64
            );
        }
        assert_eq!(arena.node_count(), count);
        assert_eq!(other.node_count(), other_count);
        if arena.is_lazy() {
            assert_eq!(arena.hash(original), LAZY_SENTINEL);
            assert_eq!(arena.hash(independent), LAZY_SENTINEL);
        }
        if other.is_lazy() {
            assert_eq!(other.hash(cross_arena), LAZY_SENTINEL);
        }
    }
}

fn ternary_dag(arena: &mut Arena, depth: usize, right_associative: bool) -> Node {
    let mut u = arena.from_bytes(b"ab");
    let mut v = arena.from_bytes(b"ba");
    for _ in 0..depth {
        let next_u;
        let next_v;
        if right_associative {
            let vu = arena.concat(v, u);
            next_u = arena.concat(u, vu);
            let uv = arena.concat(u, v);
            next_v = arena.concat(v, uv);
        } else {
            let uv = arena.concat(u, v);
            next_u = arena.concat(uv, u);
            let vu = arena.concat(v, u);
            next_v = arena.concat(vu, v);
        }
        u = next_u;
        v = next_v;
    }
    u
}

#[test]
fn huge_nonadjacent_shared_dags_and_shifted_ranges_remain_compressed() {
    let mut arena = Arena::new_lazy();
    let left = ternary_dag(&mut arena, 28, false);
    let right = ternary_dag(&mut arena, 28, true);
    let length = 2 * 3u64.pow(28);
    assert_eq!(arena.len(left), length);
    assert_eq!(arena.len(right), length);
    assert!(
        arena.node_count() < 1024,
        "test construction must stay compressed"
    );
    assert!(arena.equal(left, right));
    assert_eq!(arena.lcp(left, right), length);
    for offset in [1, 17, 512, length / 3 - 1] {
        assert!(arena.equal_range(left, offset, right, offset, length - offset - 1));
    }
    let changed = edited(&mut arena, right, length - 1, b'?');
    let count = arena.node_count();
    assert_eq!(arena.lcp(left, changed), length - 1);
    assert!(!arena.equal(left, changed));
    assert_eq!(
        arena.lcp_range(left, 17, changed, 17, length - 17),
        length - 18
    );
    assert_eq!(arena.node_count(), count);
    assert_eq!(arena.hash(left), LAZY_SENTINEL);
    assert_eq!(arena.hash(changed), LAZY_SENTINEL);
}

#[test]
fn huge_repeat_witnesses_can_themselves_be_independent_shared_dags() {
    let mut left_arena = Arena::new_lazy();
    let mut right_arena = Arena::with_hash(MERSENNE_61, 257);
    let mut left_pattern = left_arena.from_bytes(b"abc");
    let mut right_pattern = chunked(&mut right_arena, b"abc", 1);
    for _ in 0..38 {
        left_pattern = left_arena.concat(left_pattern, left_pattern);
        right_pattern = right_arena.concat(right_pattern, right_pattern);
    }
    right_pattern = right_arena.concat(right_pattern, right_pattern);
    let left = left_arena.repeat(left_pattern, 14);
    let right = right_arena.repeat(right_pattern, 7);
    let length = 42u64 << 38;
    assert_eq!(left_arena.len(left), length);
    assert_eq!(right_arena.len(right), length);
    assert!(left_arena.equal_across(left, &right_arena, right));
    assert_eq!(
        left_arena.lcp_range_across(left, 5, &right_arena, right, 8, length - 8),
        length - 8
    );
    assert_eq!(
        left_arena.lcp_range_across(left, 5, &right_arena, right, 7, length - 7),
        0
    );
    assert_eq!(left_arena.hash(left), LAZY_SENTINEL);
}

#[test]
fn same_identity_fast_paths_still_reject_invalid_ranges() {
    use std::panic::{catch_unwind, AssertUnwindSafe};
    let mut arena = Arena::new_lazy();
    let root = arena.from_bytes(b"abc");
    for (start, length) in [(0, 4), (3, 1), (4, 0), (u64::MAX, 0), (1, u64::MAX)] {
        assert!(catch_unwind(AssertUnwindSafe(
            || arena.equal_range(root, start, root, start, length)
        ))
        .is_err());
        assert!(catch_unwind(AssertUnwindSafe(
            || arena.lcp_range(root, start, root, start, length)
        ))
        .is_err());
    }
    assert!(catch_unwind(AssertUnwindSafe(|| arena.equal_range(None, 1, None, 1, 0))).is_err());
    assert!(arena.equal_range(root, 3, root, 3, 0));
    assert!(arena.equal_range(None, 0, None, 0, 0));
}

#[test]
fn symbolic_ranges_near_u64_limit_do_not_wrap_offsets_or_bounds() {
    let mut left_arena = Arena::new_lazy();
    let mut right_arena = Arena::new_lazy();
    let left_pattern = left_arena.from_bytes(b"abc");
    let right_pattern = chunked(&mut right_arena, b"abc", 1);
    // u64::MAX is divisible by three, so this uses the largest valid length.
    let left = left_arena.repeat(left_pattern, u64::MAX / 3);
    let right = right_arena.repeat(right_pattern, u64::MAX / 3);
    assert_eq!(left_arena.len(left), u64::MAX);
    assert!(left_arena.equal_across(left, &right_arena, right));
    assert_eq!(left_arena.lcp_across(left, &right_arena, right), u64::MAX);
    assert!(left_arena.equal_range_across(
        left,
        u64::MAX - 7,
        &right_arena,
        right,
        u64::MAX - 10,
        7
    ));
    assert_eq!(
        left_arena.lcp_range_across(left, u64::MAX - 7, &right_arena, right, u64::MAX - 9, 7),
        0
    );
    assert!(left_arena.equal_range(left, u64::MAX, left, u64::MAX, 0));
    use std::panic::{catch_unwind, AssertUnwindSafe};
    assert!(catch_unwind(AssertUnwindSafe(|| left_arena.equal_range(
        left,
        u64::MAX - 7,
        left,
        0,
        8
    )))
    .is_err());
    assert_eq!(left_arena.hash(left), LAZY_SENTINEL);
    assert_eq!(right_arena.hash(right), LAZY_SENTINEL);
}
