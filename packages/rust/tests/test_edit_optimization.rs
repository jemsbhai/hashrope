//! Compatibility oracles for edit and metadata optimizations.
//! The scalar hash below freezes published Rust 0.3.2 arithmetic, including
//! accepted non-Mersenne profiles (outside the documented modular contract).
use hashrope::{mersenne_mod, Arena, Node, NodeInner, PolynomialHash, LAZY_SENTINEL, MERSENNE_61};

fn published_reduce(mut value: u128, modulus: u64) -> u64 {
    assert_ne!(modulus, 0);
    let shift = 64 - modulus.leading_zeros();
    while value > modulus as u128 {
        value = (value >> shift) + (value & modulus as u128);
    }
    if value >= modulus as u128 {
        (value - modulus as u128) as u64
    } else {
        value as u64
    }
}

fn published_hash(bytes: &[u8], modulus: u64, base: u64) -> u64 {
    let mut hash = 0;
    for &byte in bytes {
        hash = published_reduce(hash as u128 * base as u128 + byte as u128 + 1, modulus);
    }
    hash
}

fn random_next(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state
}

#[test]
fn reduction_matches_frozen_scalar_for_full_u128_domain_boundaries() {
    // Public reduction accepts u128 inputs, not just the smaller intermediates
    // of hash(). Preserve its exact behavior for non-Mersenne parameters too.
    let mut values = vec![0, 1, 2, u64::MAX as u128, u128::MAX - 1, u128::MAX];
    for bit in 0..128 {
        let power = 1u128 << bit;
        values.extend([power - 1, power, power.saturating_add(1)]);
    }
    let mut random = 0x1f4b_a75c_8273_d610;
    for _ in 0..1024 {
        let high = random_next(&mut random) as u128;
        values.push((high << 64) | random_next(&mut random) as u128);
    }
    for modulus in [
        1,
        2,
        3,
        17,
        31,
        127,
        (1 << 31) - 1,
        MERSENNE_61,
        u64::MAX - 1,
        u64::MAX,
    ] {
        let mut around_modulus = Vec::new();
        for exponent in [0, 1, 2, 31, 61, 63, 64] {
            if let Some(value) = (modulus as u128).checked_mul(1u128 << exponent) {
                around_modulus.extend([value.saturating_sub(1), value, value.saturating_add(1)]);
            }
        }
        let squared = modulus as u128 * modulus as u128;
        around_modulus.extend([
            squared.saturating_sub(1),
            squared,
            squared.saturating_add(1),
        ]);
        for &value in values.iter().chain(around_modulus.iter()) {
            assert_eq!(
                mersenne_mod(value, modulus),
                published_reduce(value, modulus),
                "modulus={modulus}, value={value}"
            );
        }
    }
}

#[test]
fn hash_matches_frozen_scalar_arithmetic_at_all_block_boundaries() {
    let profiles = [
        (MERSENNE_61, 131),
        (MERSENNE_61, 257),
        (MERSENNE_61, MERSENNE_61 - 1),
        (17, 2),
        (3, 2),
        (31, 3),
        (127, 7),
        (1_000_000_007, 131),
    ];
    let mut random = 0x9a34_1ce8_20b7_fdd1;
    let bytes: Vec<u8> = (0..2048)
        .map(|_| (random_next(&mut random) >> 32) as u8)
        .collect();
    for (modulus, base) in profiles {
        let hasher = PolynomialHash::new(modulus, base);
        for length in 0..=1024 {
            for offset in [0, 1, 7, 13] {
                let slice = &bytes[offset..offset + length];
                assert_eq!(
                    hasher.hash(slice),
                    published_hash(slice, modulus, base),
                    "profile=({modulus},{base}), offset={offset}, length={length}"
                );
            }
        }
        // Exercise every byte, every tail residue, and carries at full blocks.
        for byte in 0..=255 {
            let repeated = [byte; 265];
            for length in [
                1, 7, 8, 9, 15, 16, 17, 63, 64, 65, 255, 256, 257, 258, 259, 260, 261, 262, 263,
                264, 265,
            ] {
                assert_eq!(
                    hasher.hash(&repeated[..length]),
                    published_hash(&repeated[..length], modulus, base)
                );
            }
        }
    }
}

#[test]
fn default_hash_long_buffers_and_maximum_block_carries_match_scalar() {
    let hasher = PolynomialHash::default_hash();
    let mut random = 0x7da6_38f2_f910_b3c5;
    let random_bytes: Vec<u8> = (0..1_048_584)
        .map(|_| (random_next(&mut random) >> 32) as u8)
        .collect();
    for pattern in [random_bytes, vec![0; 1_048_584], vec![255; 1_048_584]] {
        for length in [
            4095, 4096, 4097, 65_535, 65_536, 65_537, 1_048_576, 1_048_577, 1_048_583,
        ] {
            assert_eq!(
                hasher.hash(&pattern[1..1 + length]),
                published_hash(&pattern[1..1 + length], MERSENNE_61, 131)
            );
        }
    }
    // Before a block h < p. This exact upper bound covers eight updates with
    // maximal byte contribution 256, without relying on overflow wrapping.
    let mut maximum = (MERSENNE_61 - 1) as u128;
    for _ in 0..8 {
        maximum = maximum.checked_mul(131).unwrap().checked_add(256).unwrap();
    }
    assert!(maximum < (1u128 << 118));
}

fn mix(digest: &mut u64, value: u64) {
    for byte in value.to_le_bytes() {
        *digest = (*digest ^ byte as u64).wrapping_mul(0x100000001b3);
    }
}

fn arena_layout_digest(arena: &Arena) -> u64 {
    let mut digest = 0xcbf29ce484222325;
    for id in 0..arena.node_count() as u32 {
        mix(&mut digest, id as u64);
        match arena.node(id) {
            NodeInner::Leaf {
                data,
                hash_val,
                len,
            } => {
                mix(&mut digest, 0);
                mix(&mut digest, *len);
                mix(&mut digest, *hash_val);
                for &byte in data {
                    mix(&mut digest, byte as u64);
                }
            }
            NodeInner::Internal {
                left,
                right,
                hash_val,
                len,
                weight,
            } => {
                for value in [1, *left as u64, *right as u64, *hash_val, *len, *weight] {
                    mix(&mut digest, value);
                }
            }
            NodeInner::Repeat {
                child,
                reps,
                hash_val,
                len,
                weight,
            } => {
                for value in [2, *child as u64, *reps, *hash_val, *len, *weight] {
                    mix(&mut digest, value);
                }
            }
        }
    }
    digest
}

fn record(
    arena: &mut Arena,
    root: Node,
    expected: Option<&[u8]>,
    digest: &mut u64,
    roots: &mut Vec<u64>,
) {
    arena.validate(root);
    if let Some(expected) = expected {
        assert_eq!(arena.to_bytes(root), expected);
    }
    for value in [
        root.map_or(u64::MAX, |id| id as u64),
        arena.len(root),
        arena.weight(root),
        arena.hash(root),
        arena.node_count() as u64,
        arena_layout_digest(arena),
    ] {
        mix(digest, value);
    }
    let length = arena.len(root);
    for (start, size) in [
        (u64::MAX, 0),
        (length / 3, length / 2),
        (length.saturating_sub(17), length.min(17)),
        (0, length),
    ] {
        let count = arena.node_count();
        mix(digest, arena.substr_hash(root, start, size));
        assert_eq!(arena.node_count(), count);
        mix(digest, arena.hash(root));
    }
    roots.push(arena.hash(root));
    mix(digest, arena_layout_digest(arena));
}

fn edit_trace(modulus: u64, base: u64, lazy: bool) -> (u64, Vec<u64>) {
    let mut arena = if lazy {
        Arena::new_lazy()
    } else {
        Arena::with_hash(modulus, base)
    };
    let mut digest = 0xcbf29ce484222325;
    let mut roots = Vec::new();
    let mut bytes: Vec<u8> = (0..1537).map(|i| (i * 37 + 13) as u8).collect();
    let mut root = arena.from_bytes(&bytes);
    record(&mut arena, root, Some(&bytes), &mut digest, &mut roots);
    let persistent = root;
    let persistent_bytes = bytes.clone();
    for round in 0..12usize {
        let delete_at = (round * 193 + 511) % bytes.len();
        let delete_size = (1 + round % 17).min(bytes.len() - delete_at);
        let (left, rest) = arena.split(root, delete_at as u64);
        let (_, right) = arena.split(rest, delete_size as u64);
        root = arena.concat(left, right);
        bytes.drain(delete_at..delete_at + delete_size);
        let insert_at = (round * 257 + 17) % (bytes.len() + 1);
        let addition: Vec<u8> = (0..1 + round % 19)
            .map(|i| (i * 131 + round * 7) as u8)
            .collect();
        let (left, right) = arena.split(root, insert_at as u64);
        let added = arena.from_bytes(&addition);
        let prefix = arena.concat(left, added);
        root = arena.concat(prefix, right);
        bytes.splice(insert_at..insert_at, addition);
        record(&mut arena, root, Some(&bytes), &mut digest, &mut roots);
        assert_eq!(arena.to_bytes(persistent), persistent_bytes);
    }
    let p1 = arena.from_bytes(&[0, 13]);
    let p2 = arena.from_bytes(&[255, 131, 0]);
    let pattern = arena.concat(p1, p2);
    let repeated = arena.repeat(pattern, 17);
    record(
        &mut arena,
        repeated,
        Some(&[0, 13, 255, 131, 0].repeat(17)),
        &mut digest,
        &mut roots,
    );
    let joined = arena.concat(root, repeated);
    bytes.extend_from_slice(&[0, 13, 255, 131, 0].repeat(17));
    record(&mut arena, joined, Some(&bytes), &mut digest, &mut roots);
    let huge = arena.repeat(pattern, (1 << 32) + 3);
    record(&mut arena, huge, None, &mut digest, &mut roots);
    let (head, tail) = arena.split(huge, 5 * (1 << 31) + 3);
    record(&mut arena, head, None, &mut digest, &mut roots);
    record(&mut arena, tail, None, &mut digest, &mut roots);
    let restored = arena.concat(head, tail);
    record(&mut arena, restored, None, &mut digest, &mut roots);
    assert!(arena.equal(huge, restored));
    (digest, roots)
}

#[test]
fn editing_preserves_baseline_fingerprints_nodes_and_lazy_materialization() {
    // Generated independently from published Rust 0.3.2 source at
    // bd8f9ca5c73e3c80507089cdda204f922869cbc0, without candidate optimizations.
    // Full-node layout digests include old/new nodes, cached hashes, leaf bytes,
    // child IDs, repeat counts and weights. Default/base257/M31 goldens also
    // match the original pre-optimization 0.3.1 baseline.
    let default_roots = [
        1293468461321973878,
        1432555342706740516,
        994328214638166104,
        422262590900319798,
        603509486079471938,
        1730508189894962939,
        1110517824445276321,
        2213302874092148224,
        2297499999215496947,
        1003893823740351525,
        464725185510976067,
        441786105360094352,
        664781496570652600,
        489978388935684972,
        855030794218754678,
        1575601001731378735,
        205919119102156243,
        93404412822878156,
        1575601001731378735,
    ];
    let base257_roots = [
        1775454010760495670,
        1124347367441257714,
        2305728859378639707,
        189402239593377288,
        884149751448301242,
        1037361669114495400,
        320147408247522668,
        1374602024672717724,
        1054500180612572907,
        38228562765633039,
        1138560548323166460,
        1590134785942207830,
        934324226197283608,
        47077818713455230,
        1714905931040353409,
        1299633187386640394,
        2179091956674577550,
        1769979924797985220,
        1299633187386640394,
    ];
    let profile17_roots = [5, 0, 0, 0, 5, 3, 2, 5, 4, 8, 16, 5, 16, 0, 16, 16, 1, 8, 12];
    let profile31_roots = [
        3, 30, 23, 16, 30, 12, 28, 30, 21, 2, 20, 18, 15, 12, 9, 29, 27, 23, 29,
    ];
    let profiles = [
        (MERSENNE_61, 131, false, 2390809490780052657, &default_roots),
        (
            MERSENNE_61,
            257,
            false,
            12360643512613280362,
            &base257_roots,
        ),
        (17, 2, false, 2491766406201968812, &profile17_roots),
        (31, 3, false, 3895836206056079445, &profile31_roots),
        (MERSENNE_61, 131, true, 5507448442885825566, &default_roots),
    ];
    for (modulus, base, lazy, expected_digest, expected_roots) in profiles {
        let (digest, roots) = edit_trace(modulus, base, lazy);
        assert_eq!(
            roots.as_slice(),
            expected_roots,
            "root hashes: profile=({modulus},{base}), lazy={lazy}"
        );
        assert_eq!(
            digest, expected_digest,
            "allocated nodes and materialization: profile=({modulus},{base}), lazy={lazy}"
        );
    }
}

#[test]
fn accepted_custom_profile_retains_structure_sensitive_published_hashes() {
    let mut arena = Arena::with_hash(17, 2);
    let original = arena.from_bytes(&[0, 16]);
    assert_eq!(arena.hash(original), 0);
    let (left, right) = arena.split(original, 1);
    let joined = arena.concat(left, right);
    assert_eq!(arena.hash(joined), 2);
    assert_eq!(arena.substr_hash(original, 0, 2), 0);
    assert_eq!(arena.substr_hash(joined, 0, 2), 2);
    assert!(arena.equal(original, joined));
    assert_eq!(arena.to_bytes(original), arena.to_bytes(joined));
}

#[test]
fn identity_empty_and_clamped_edit_edges_keep_their_behavior() {
    for lazy in [false, true] {
        let mut arena = if lazy {
            Arena::new_lazy()
        } else {
            Arena::new()
        };
        let root = arena.from_bytes(b"abc");
        let count = arena.node_count();
        assert_eq!(arena.from_bytes(b""), None);
        assert_eq!(arena.concat(None, root), root);
        assert_eq!(arena.concat(root, None), root);
        assert_eq!(arena.repeat(root, 0), None);
        assert_eq!(arena.repeat(root, 1), root);
        assert_eq!(arena.repeat(None, u64::MAX), None);
        assert_eq!(arena.split(root, 0), (None, root));
        assert_eq!(arena.split(root, 3), (root, None));
        assert_eq!(arena.split(root, u64::MAX), (root, None));
        assert_eq!(arena.split(None, u64::MAX), (None, None));
        assert_eq!(arena.substr_hash(None, u64::MAX, u64::MAX), 0);
        assert_eq!(arena.substr_hash(root, u64::MAX, 0), 0);
        assert_eq!(arena.node_count(), count);
        if lazy {
            assert_eq!(arena.hash(root), LAZY_SENTINEL);
        }
        let (left, right) = arena.split(root, 1);
        assert_eq!(arena.to_bytes(left), b"a");
        assert_eq!(arena.to_bytes(right), b"bc");
        assert_eq!(arena.to_bytes(root), b"abc");
        if lazy {
            assert_eq!(arena.hash(left), LAZY_SENTINEL);
            assert_eq!(arena.hash(right), LAZY_SENTINEL);
            assert_eq!(arena.hash(root), LAZY_SENTINEL);
            arena.ensure_hash(left.unwrap());
            assert_ne!(arena.hash(left), LAZY_SENTINEL);
            assert_eq!(arena.hash(right), LAZY_SENTINEL);
            assert_eq!(arena.hash(root), LAZY_SENTINEL);
        }
    }
}

#[test]
#[should_panic(expected = "Modulus must be nonzero")]
fn zero_modulus_preserves_published_rejection() {
    mersenne_mod(0, 0);
}
