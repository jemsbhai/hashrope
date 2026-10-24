//! Byte-oracle coverage for period-root proofs and aligned branch traversal.
use hashrope::{Arena, Node, LAZY_SENTINEL, MERSENNE_61};

fn chunked(arena: &mut Arena, bytes: &[u8], width: usize) -> Node {
    let mut root = None;
    for part in bytes.chunks(width) {
        let leaf = arena.from_bytes(part);
        root = arena.concat(root, leaf);
    }
    root
}

fn doubled(arena: &mut Arena, mut root: Node, levels: usize) -> Node {
    for _ in 0..levels {
        root = arena.concat(root, root);
    }
    root
}

fn byte_lcp(a: &[u8], b: &[u8]) -> u64 {
    a.iter()
        .zip(b)
        .position(|(x, y)| x != y)
        .unwrap_or(a.len().min(b.len())) as u64
}

fn same_period_oracle(a: &[u8], a_phase: usize, b: &[u8], b_phase: usize, length: u64) -> u64 {
    assert_eq!(a.len(), b.len());
    for i in 0..length.min(a.len() as u64) as usize {
        if a[(a_phase + i) % a.len()] != b[(b_phase + i) % b.len()] {
            return i as u64;
        }
    }
    length
}

#[test]
fn equal_period_proofs_preserve_all_phases_and_short_ranges() {
    for period in [1usize, 2, 3, 7, 16] {
        let pattern: Vec<u8> = (0..period).map(|i| (i * i + 3 * i) as u8).collect();
        for changed in [false, true] {
            let mut other_pattern = pattern.clone();
            if changed {
                other_pattern[period - 1] ^= 0x80;
            }
            let mut a = Arena::new_lazy();
            let mut b = Arena::with_hash(MERSENNE_61, 257);
            let unit_a = chunked(&mut a, &pattern, 3);
            let unit_b = chunked(&mut b, &other_pattern, 5);
            let left = doubled(&mut a, unit_a, 40);
            let right = b.repeat(unit_b, 1 << 40);
            let total = (period as u64) << 40;
            let counts = (a.node_count(), b.node_count());
            for phase_a in 0..period {
                for phase_b in 0..period {
                    let start_a = ((period as u64) << 20) + phase_a as u64;
                    let start_b = ((period as u64) << 21) + phase_b as u64;
                    for span in [
                        0,
                        1,
                        period as u64 - 1,
                        period as u64,
                        period as u64 + 1,
                        2 * period as u64 + 3,
                        total - start_a.max(start_b),
                    ] {
                        let expected =
                            same_period_oracle(&pattern, phase_a, &other_pattern, phase_b, span);
                        assert_eq!(
                            a.lcp_range_across(left, start_a, &b, right, start_b, span),
                            expected
                        );
                        assert_eq!(
                            a.equal_range_across(left, start_a, &b, right, start_b, span),
                            expected == span
                        );
                    }
                }
            }
            assert_eq!((a.node_count(), b.node_count()), counts);
            assert_eq!(a.hash(left), LAZY_SENTINEL);
        }
    }
}

#[test]
fn equal_period_proof_checks_the_last_witness_byte_after_wrapping() {
    let period = 17;
    let pattern = vec![b'a'; period];
    let mut changed = pattern.clone();
    changed[0] = b'b';
    let mut arena = Arena::new_lazy();
    let unit_a = chunked(&mut arena, &pattern, 3);
    let unit_b = chunked(&mut arena, &changed, 4);
    let left = doubled(&mut arena, unit_a, 40);
    let right = doubled(&mut arena, unit_b, 40);
    let span = (period as u64) << 39;
    assert_eq!(arena.lcp_range(left, 1, right, 1, span), 16);
    assert!(!arena.equal_range(left, 1, right, 1, span));
    assert_eq!(arena.lcp_range(left, 1, right, 1, 16), 16);
    assert_eq!(arena.lcp_range(left, 1, right, 1, 17), 16);
}

#[test]
fn three_phase_pieces_keep_large_nonperiodic_units_compressed() {
    let n = 1u64 << 35;
    let mut a = Arena::new_lazy();
    let mut b = Arena::with_hash(MERSENNE_61, 257);
    let aa = a.from_bytes(b"a");
    let ab = a.from_bytes(b"b");
    let ba = b.from_bytes(b"a");
    let bb = b.from_bytes(b"b");
    let a_first = a.repeat(aa, n);
    let a_second = a.repeat(ab, n);
    let b_first = b.repeat(bb, n);
    let b_second = b.repeat(ba, n);
    let unit_a = a.concat(a_first, a_second);
    let unit_b = b.concat(b_first, b_second);
    let left = a.repeat(unit_a, 7);
    let right = b.repeat(unit_b, 7);
    let start_a = n / 2;
    let start_b = n + n / 2;
    let span = 10 * n;
    let counts = (a.node_count(), b.node_count());
    assert_eq!(
        a.lcp_range_across(left, start_a, &b, right, start_b, span),
        span
    );
    assert!(a.equal_range_across(left, start_a, &b, right, start_b, span));
    assert_eq!((a.node_count(), b.node_count()), counts);

    // The only mutation lies immediately before the starting phase, so the
    // last byte of the complete 2N-byte witness must actually be compared.
    let before = b.repeat(ba, n / 2 - 1);
    let changed = b.from_bytes(b"z");
    let after = b.repeat(ba, n / 2);
    let prefix = b.concat(b_first, before);
    let prefix = b.concat(prefix, changed);
    let changed_unit = b.concat(prefix, after);
    let changed_root = b.repeat(changed_unit, 7);
    assert_eq!(b.len(changed_unit), 2 * n);
    assert_eq!(
        a.lcp_range_across(left, start_a, &b, changed_root, start_b, span),
        2 * n - 1
    );
    assert!(!a.equal_range_across(left, start_a, &b, changed_root, start_b, span));
    assert_eq!(a.hash(left), LAZY_SENTINEL);
}

#[test]
fn aligned_and_mismatched_tree_shapes_match_plain_range_oracles() {
    let bytes: Vec<u8> = (0..32_771)
        .map(|i| ((i * 37 + i / 17) % 251) as u8)
        .collect();
    for lazy in [false, true] {
        for (width_a, width_b) in [(31, 31), (256, 256), (511, 512), (2, 3)] {
            for mismatch in [
                None,
                Some(0),
                Some(511),
                Some(16_383),
                Some(bytes.len() - 1),
            ] {
                let mut other = bytes.clone();
                if let Some(i) = mismatch {
                    other[i] ^= 0x80;
                }
                let mut arena = if lazy {
                    Arena::new_lazy()
                } else {
                    Arena::new()
                };
                let left = chunked(&mut arena, &bytes, width_a);
                let right = chunked(&mut arena, &other, width_b);
                let count = arena.node_count();
                assert_eq!(arena.lcp(left, right), byte_lcp(&bytes, &other));
                assert_eq!(arena.equal(left, right), bytes == other);
                for (start_a, start_b, span) in [
                    (0, 0, 513),
                    (1, 1, 8193),
                    (7, 7, bytes.len() - 7),
                    (5, 9, bytes.len() - 9),
                ] {
                    let expected = byte_lcp(
                        &bytes[start_a..start_a + span],
                        &other[start_b..start_b + span],
                    );
                    assert_eq!(
                        arena.lcp_range(left, start_a as u64, right, start_b as u64, span as u64),
                        expected
                    );
                    assert_eq!(
                        arena.equal_range(left, start_a as u64, right, start_b as u64, span as u64),
                        expected == span as u64
                    );
                }
                assert_eq!(arena.node_count(), count);
                if lazy {
                    assert_eq!(arena.hash(left), LAZY_SENTINEL);
                }
            }
        }
    }
}

#[test]
fn aligned_shared_prefix_preserves_repeat_fallback_and_distinct_periods() {
    let mut arena = Arena::new_lazy();
    let prefix = arena.from_bytes(&vec![b'p'; 8192]);
    let unit_a = arena.from_bytes(b"aab");
    let unit_b = arena.from_bytes(b"aaba");
    let repeat_a = arena.repeat(unit_a, 1 << 40);
    let repeat_b = arena.repeat(unit_b, 1 << 40);
    assert_eq!(arena.lcp(repeat_a, repeat_b), 5);
    let left = arena.concat(prefix, repeat_a);
    let right = arena.concat(prefix, repeat_b);
    assert_eq!(arena.lcp(left, right), 8192 + 5);
    assert_eq!(arena.lcp_range(left, 8190, right, 8190, 100), 7);
}
