//! Independent oracles for structural periods and direct child traversal.
use hashrope::{Arena, Node, LAZY_SENTINEL, MERSENNE_61};

fn chunked(arena: &mut Arena, bytes: &[u8], width: usize) -> Node {
    let mut root = None;
    for chunk in bytes.chunks(width) {
        let leaf = arena.from_bytes(chunk);
        root = arena.concat(root, leaf);
    }
    root
}

fn double(arena: &mut Arena, mut node: Node, levels: usize) -> Node {
    for _ in 0..levels {
        node = arena.concat(node, node);
    }
    node
}

fn edit(arena: &mut Arena, node: Node, offset: u64, byte: u8) -> Node {
    let (prefix, rest) = arena.split(node, offset);
    let (_, suffix) = arena.split(rest, 1);
    let replacement = arena.from_bytes(&[byte]);
    let prefix = arena.concat(prefix, replacement);
    arena.concat(prefix, suffix)
}

fn periodic_lcp(a: &[u8], phase_a: usize, b: &[u8], phase_b: usize, limit: u64) -> u64 {
    // This is only used with small authored periods. Inspect their product:
    // every pair of phases recurs by then, so agreement proves the rest.
    for offset in 0..(a.len() * b.len()) {
        if offset as u64 == limit {
            return limit;
        }
        if a[(phase_a + offset) % a.len()] != b[(phase_b + offset) % b.len()] {
            return offset as u64;
        }
    }
    limit
}

#[test]
fn doubled_dags_keep_nonzero_phases_and_cross_profile_ranges_exact() {
    let pattern = b"aabaa";
    let mut a = Arena::new_lazy();
    let mut b = Arena::with_hash(MERSENNE_61, 257);
    let unit_a = a.from_bytes(pattern);
    let unit_b = b.from_bytes(b"aabaaaabaa");
    let left = double(&mut a, unit_a, 40);
    let right = double(&mut b, unit_b, 39);
    let length = 5u64 << 40;
    assert_eq!(a.len(left), length);
    assert_eq!(b.len(right), length);
    let counts = (a.node_count(), b.node_count());
    for left_phase in 0..5 {
        for right_phase in 0..10 {
            // Nonzero large offsets test phase rather than unit-prefix equality.
            let left_start = (5u64 << 18) + left_phase as u64;
            let right_start = (10u64 << 17) + right_phase as u64;
            let span = length - left_start.max(right_start) - 17;
            let expected = periodic_lcp(pattern, left_phase, b"aabaaaabaa", right_phase, span);
            assert_eq!(
                a.lcp_range_across(left, left_start, &b, right, right_start, span),
                expected
            );
            assert_eq!(
                a.equal_range_across(left, left_start, &b, right, right_start, span),
                expected == span
            );
        }
    }
    assert_eq!((a.node_count(), b.node_count()), counts);
    assert_eq!(a.hash(left), LAZY_SENTINEL);
}

#[test]
fn structurally_periodic_witness_checks_its_last_byte() {
    let mut a = Arena::new_lazy();
    let mut b = Arena::new_lazy();
    let unit_a = a.from_bytes(b"aab");
    let unit_b = b.from_bytes(b"aaba");
    let unit_a = a.repeat(unit_a, 7);
    let unit_b = b.repeat(unit_b, 5);
    let left = double(&mut a, unit_a, 35);
    let right = double(&mut b, unit_b, 35);
    // Periods 3 and 4 agree through byte 4, but disagree at byte 5,
    // the last byte of the p+q-gcd(p,q)=6 witness.
    assert_eq!(a.lcp_across(left, &b, right), 5);
    assert!(!a.equal_across(left, &b, right));
    for left_phase in 0..3 {
        for right_phase in 0..4 {
            let expected = periodic_lcp(b"aab", left_phase, b"aaba", right_phase, 10000);
            assert_eq!(
                a.lcp_range_across(
                    left,
                    left_phase as u64,
                    &b,
                    right,
                    right_phase as u64,
                    10000
                ),
                expected
            );
        }
    }
}

#[test]
fn late_edit_breaks_only_the_edited_dag_range() {
    let mut arena = Arena::new_lazy();
    let unit_a = arena.from_bytes(b"abc");
    let unit_b = chunked(&mut arena, b"abc", 1);
    let left = double(&mut arena, unit_a, 40);
    let right = double(&mut arena, unit_b, 40);
    let length = 3u64 << 40;
    let changed = edit(&mut arena, right, length - 1, b'?');
    assert_eq!(arena.lcp(left, changed), length - 1);
    assert!(!arena.equal(left, changed));
    assert!(arena.equal_range(left, 7, changed, 7, length - 8));
    assert_eq!(arena.lcp_range(left, 7, changed, 7, length - 7), length - 8);
    assert_eq!(arena.lcp_range(left, length - 2, changed, length - 2, 2), 1);
    // A previous query must not leak cached equality into another query.
    assert!(arena.equal(left, right));
    assert!(!arena.equal(left, changed));
    assert_eq!(arena.hash(changed), LAZY_SENTINEL);
}

#[test]
fn equal_fingerprints_do_not_define_structural_periods() {
    let mut a = Arena::new();
    let mut b = Arena::with_hash(MERSENNE_61, 131);
    let unit_a = a.from_bytes(&[0, 131]);
    let unit_b = b.from_bytes(&[1, 0]);
    assert_eq!(a.hash(unit_a), b.hash(unit_b));
    let left = double(&mut a, unit_a, 39);
    let right = double(&mut b, unit_b, 39);
    // Numeric node IDs coincide between these independent arenas.
    assert_eq!(left, right);
    assert_eq!(a.lcp_across(left, &b, right), 0);
    assert!(!a.equal_across(left, &b, right));
}

#[test]
fn accepted_custom_profile_hash_inequality_does_not_reject_equal_bytes() {
    // Legacy construction accepts this non-Mersenne modulus even though its
    // reducer is not a general modular arithmetic implementation. Exact APIs
    // promise byte comparisons independently of these cached fingerprints.
    let mut arena = Arena::with_hash(17, 2);
    let flat = arena.from_bytes(&[0, 16]);
    let first = arena.from_bytes(&[0]);
    let second = arena.from_bytes(&[16]);
    let joined = arena.concat(first, second);
    assert_eq!(arena.hash(flat), 0);
    assert_eq!(arena.hash(joined), 2);
    assert_eq!(arena.to_bytes(flat), arena.to_bytes(joined));
    assert!(arena.equal(flat, joined));
    assert_eq!(arena.lcp(flat, joined), 2);
}

#[test]
fn uneven_child_boundaries_and_clipped_ranges_match_byte_oracles() {
    let bytes: Vec<u8> = (0..2053).map(|i| (i * 73 + i / 19) as u8).collect();
    for (left_width, right_width) in [(7, 3), (64, 63), (257, 511), (512, 255)] {
        for position in [0, 1, 63, 64, 255, 256, 511, 512, 513, 1023, 1024, 2052] {
            let mut left_bytes = b"abc".to_vec();
            left_bytes.extend_from_slice(&bytes);
            let mut changed = bytes.clone();
            changed[position] ^= 0x80;
            let mut right_bytes = b"prefix".to_vec();
            right_bytes.extend_from_slice(&changed);
            let mut a = Arena::new_lazy();
            let mut b = Arena::new_lazy();
            let left = chunked(&mut a, &left_bytes, left_width);
            let right = chunked(&mut b, &right_bytes, right_width);
            for start in [0, 1, left_width - 1, left_width, right_width, 1024] {
                for end in [bytes.len(), position.max(start)] {
                    if end < start {
                        continue;
                    }
                    let length = end - start;
                    let expected = if (start..end).contains(&position) {
                        position - start
                    } else {
                        length
                    };
                    assert_eq!(
                        a.lcp_range_across(
                            left,
                            3 + start as u64,
                            &b,
                            right,
                            6 + start as u64,
                            length as u64
                        ),
                        expected as u64
                    );
                }
            }
        }
    }
}

#[test]
fn direct_child_traversal_can_fall_back_after_a_nonzero_prefix() {
    let mut a = Arena::new_lazy();
    let mut b = Arena::new_lazy();
    let prefix = vec![b'x'; 1025];
    let left_prefix = chunked(&mut a, &prefix, 257);
    let right_prefix = chunked(&mut b, &prefix, 63);
    let left_unit = a.from_bytes(b"aab");
    let right_unit = b.from_bytes(b"aaba");
    let left_repeat = a.repeat(left_unit, 1 << 35);
    let right_repeat = b.repeat(right_unit, 1 << 35);
    let left = a.concat(left_prefix, left_repeat);
    let right = b.concat(right_prefix, right_repeat);
    assert_eq!(a.lcp_across(left, &b, right), 1025 + 5);
    assert_eq!(
        a.lcp_range_across(left, 17, &b, right, 17, 10000),
        1025 + 5 - 17
    );
    assert_eq!(a.hash(left), LAZY_SENTINEL);
    assert_eq!(b.hash(right), LAZY_SENTINEL);
}

#[test]
fn every_leaf_mismatch_position_and_bit_has_the_correct_byte_index() {
    for length in 1usize..=512 {
        let bytes: Vec<u8> = (0..length).map(|i| (i * 37 + length) as u8).collect();
        for position in 0..length {
            // Vary the source-slice alignment independently on both sides.
            // The differing bit cycles through a byte, including its high bit.
            let left_offset = (position + length) % 16;
            let right_offset = (position * 3 + length + 1) % 16;
            let mut left_bytes = vec![b'L'; left_offset];
            left_bytes.extend_from_slice(&bytes);
            let mut right_bytes = vec![b'R'; right_offset];
            right_bytes.extend_from_slice(&bytes);
            right_bytes[right_offset + position] ^= 1 << (position % 8);
            let mut arena = Arena::new_lazy();
            let left = arena.from_bytes(&left_bytes);
            let right = arena.from_bytes(&right_bytes);
            assert_eq!(
                arena.lcp_range(
                    left,
                    left_offset as u64,
                    right,
                    right_offset as u64,
                    length as u64
                ),
                position as u64,
                "length={length}, position={position}"
            );
        }
    }
}

#[test]
fn compressed_leaf_mismatch_localization_matches_the_plain_byte_oracle() {
    for length in [1usize, 7, 8, 9, 15, 16, 17, 63, 64, 65, 511, 512] {
        for position in 0..length {
            let bytes: Vec<u8> = (0..length).map(|i| (i * 19 + 3) as u8).collect();
            let mut changed = bytes.clone();
            changed[position] ^= 0x80;
            let mut arena = Arena::new_lazy();
            let left_unit = arena.from_bytes(&bytes);
            let right_unit = arena.from_bytes(&changed);
            // The logical leaf weight forces the compressed comparison engine.
            let left = arena.repeat(left_unit, 1 << 30);
            let right = arena.repeat(right_unit, 1 << 30);
            assert_eq!(arena.lcp(left, right), position as u64);
        }
    }
}

#[test]
fn long_uneven_frontiers_keep_clipped_ranges_and_late_mismatches_exact() {
    let bytes: Vec<u8> = (0..32777).map(|i| (i * 31 + i / 43) as u8).collect();
    for (left_width, right_width) in [(7, 13), (509, 256), (512, 511)] {
        for position in [0, 7, 8, 511, 512, 513, 8191, 8192, 16384, 32776] {
            let mut left_bytes = vec![b'L'; 11];
            left_bytes.extend_from_slice(&bytes);
            let mut right_bytes = vec![b'R'; 19];
            right_bytes.extend_from_slice(&bytes);
            right_bytes[19 + position] ^= 0x40;
            let mut a = Arena::new_lazy();
            let mut b = Arena::new_lazy();
            let left = chunked(&mut a, &left_bytes, left_width);
            let right_same = chunked(&mut a, &right_bytes, right_width);
            let right_across = chunked(&mut b, &right_bytes, right_width);
            for start in [0, 3, 511, 4099] {
                let length = bytes.len() - start;
                let expected = if position >= start {
                    position - start
                } else {
                    length
                };
                assert_eq!(
                    a.lcp_range(
                        left,
                        11 + start as u64,
                        right_same,
                        19 + start as u64,
                        length as u64
                    ),
                    expected as u64
                );
                assert_eq!(
                    a.lcp_range_across(
                        left,
                        11 + start as u64,
                        &b,
                        right_across,
                        19 + start as u64,
                        length as u64
                    ),
                    expected as u64
                );
            }
        }
    }
}
