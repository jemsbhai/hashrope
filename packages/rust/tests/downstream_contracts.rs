//! External-consumer contracts for applications that inspect public rope nodes.
//!
//! These synthetic fixtures use only the published API. The exhaustive match
//! intentionally fails to compile if a variant or field is added, removed, or
//! changed: downstream consumers can depend on that exact public shape.

use hashrope::{
    mersenne_mod, mersenne_mul, phi, rope_concat, rope_from_bytes, rope_hash, rope_height,
    rope_len, rope_repeat, rope_split, rope_substr_hash, rope_to_bytes, validate_rope, Arena, Node,
    NodeId, NodeInner, PolynomialHash, SlidingWindow, LAZY_SENTINEL, MERSENNE_61,
};

#[derive(Debug, PartialEq, Eq)]
enum PublicShape {
    Leaf(Vec<u8>, u64, u64),
    Internal(u32, u32, u64, u64, u64),
    Repeat(u32, u64, u64, u64, u64),
}

fn snapshot(node: &NodeInner) -> PublicShape {
    // No wildcard or `..`: a real external exhaustive match is the contract.
    // The constructor arguments also require the exact public field types.
    match node {
        NodeInner::Leaf {
            data,
            hash_val,
            len,
        } => PublicShape::Leaf(data.clone(), *hash_val, *len),
        NodeInner::Internal {
            left,
            right,
            hash_val,
            len,
            weight,
        } => PublicShape::Internal(*left, *right, *hash_val, *len, *weight),
        NodeInner::Repeat {
            child,
            reps,
            hash_val,
            len,
            weight,
        } => PublicShape::Repeat(*child, *reps, *hash_val, *len, *weight),
    }
}

#[test]
fn public_types_and_free_functions_remain_source_compatible() {
    // Exact aliases and callable signatures, checked from outside the crate.
    let id: NodeId = 17u32;
    let raw_id: u32 = id;
    let node: Node = Some(raw_id);
    let _: Option<u32> = node;
    let _: fn(u128, u64) -> u64 = mersenne_mod;
    let _: fn(u64, u64, u64) -> u64 = mersenne_mul;
    let _: fn(u64, u64, u64) -> u64 = phi;
    let _: fn(&Arena, Node) -> u64 = rope_len;
    let _: fn(&Arena, Node) -> u64 = rope_hash;
    let _: fn(&Arena, Node) -> u64 = rope_height;
    let _: fn(&mut Arena, Node, Node) -> Node = rope_concat;
    let _: fn(&mut Arena, Node, u64) -> (Node, Node) = rope_split;
    let _: fn(&mut Arena, Node, u64) -> Node = rope_repeat;
    let _: fn(&mut Arena, Node, u64, u64) -> u64 = rope_substr_hash;
    let _: fn(&mut Arena, &[u8]) -> Node = rope_from_bytes;
    let _: fn(&Arena, Node) -> Vec<u8> = rope_to_bytes;
    let _: fn(&Arena, Node) = validate_rope;
    let _: fn(u64, u64) -> Arena = Arena::with_hash;

    let mut arena = Arena::new();
    assert_eq!(arena.hasher().prime(), MERSENNE_61);
    assert_eq!(arena.hasher().base(), 131);
    let oracle = PolynomialHash::default_hash();
    let left_bytes = b"\0A\xc3\xa9";
    let right_bytes = b"\xff\r\n";
    let mut pattern = left_bytes.to_vec();
    pattern.extend_from_slice(right_bytes);
    let left = rope_from_bytes(&mut arena, left_bytes);
    let right = rope_from_bytes(&mut arena, right_bytes);
    let joined = rope_concat(&mut arena, left, right);
    let repeated = rope_repeat(&mut arena, joined, 3);
    let expected = pattern.repeat(3);

    assert_eq!(
        snapshot(arena.node(left.unwrap())),
        PublicShape::Leaf(left_bytes.to_vec(), oracle.hash(left_bytes), 4)
    );
    assert_eq!(
        snapshot(arena.node(joined.unwrap())),
        PublicShape::Internal(left.unwrap(), right.unwrap(), oracle.hash(&pattern), 7, 2)
    );
    assert_eq!(
        snapshot(arena.node(repeated.unwrap())),
        PublicShape::Repeat(joined.unwrap(), 3, oracle.hash(&expected), 21, 6)
    );
    let weight: u64 = arena.weight(repeated);
    assert_eq!(weight, 6);
    assert_eq!(rope_len(&arena, repeated), 21);
    assert_eq!(rope_height(&arena, repeated), 2);
    assert_eq!(rope_hash(&arena, repeated), oracle.hash(&expected));
    assert_eq!(rope_to_bytes(&arena, repeated), expected);
    assert_eq!(
        rope_substr_hash(&mut arena, repeated, 2, 11),
        oracle.hash(&expected[2..13])
    );

    let (prefix, suffix) = rope_split(&mut arena, repeated, 5);
    assert_eq!(rope_to_bytes(&arena, prefix), expected[..5]);
    assert_eq!(rope_to_bytes(&arena, suffix), expected[5..]);
    let restored = rope_concat(&mut arena, prefix, suffix);
    validate_rope(&arena, restored);
    assert_eq!(rope_hash(&arena, restored), oracle.hash(&expected));
    assert_eq!(rope_from_bytes(&mut arena, b""), None);
    assert_eq!(rope_hash(&arena, None), 0);
    assert_eq!(rope_len(&arena, None), 0);
    assert_eq!(rope_to_bytes(&arena, None), b"");
}

#[test]
fn lazy_nodes_remain_inspectable_before_hash_materialization() {
    assert_eq!(LAZY_SENTINEL, u64::MAX);
    let mut arena = Arena::new_lazy();
    let left = rope_from_bytes(&mut arena, b"\0\xff");
    let right = rope_from_bytes(&mut arena, b"\xc3\xa9");
    let joined = rope_concat(&mut arena, left, right);
    let repeated = rope_repeat(&mut arena, joined, 3);
    let expected = b"\0\xff\xc3\xa9".repeat(3);
    let ids = [
        left.unwrap(),
        right.unwrap(),
        joined.unwrap(),
        repeated.unwrap(),
    ];

    assert!(arena.is_lazy());
    assert_eq!(
        snapshot(arena.node(left.unwrap())),
        PublicShape::Leaf(b"\0\xff".to_vec(), LAZY_SENTINEL, 2)
    );
    assert_eq!(
        snapshot(arena.node(joined.unwrap())),
        PublicShape::Internal(left.unwrap(), right.unwrap(), LAZY_SENTINEL, 4, 2)
    );
    assert_eq!(
        snapshot(arena.node(repeated.unwrap())),
        PublicShape::Repeat(joined.unwrap(), 3, LAZY_SENTINEL, 12, 6)
    );
    assert_eq!(rope_to_bytes(&arena, repeated), expected);
    validate_rope(&arena, repeated);
    for id in ids {
        assert_eq!(arena.node(id).hash_val(), LAZY_SENTINEL);
        assert_eq!(rope_hash(&arena, Some(id)), LAZY_SENTINEL);
    }

    let oracle = PolynomialHash::default_hash();
    assert_eq!(
        rope_substr_hash(&mut arena, repeated, 0, expected.len() as u64),
        oracle.hash(&expected)
    );
    for id in ids {
        let bytes = rope_to_bytes(&arena, Some(id));
        assert_eq!(arena.node(id).hash_val(), oracle.hash(&bytes));
        assert_eq!(rope_hash(&arena, Some(id)), oracle.hash(&bytes));
    }
}

#[test]
fn legacy_sliding_window_preserves_observed_and_caller_allocated_handles() {
    let _: fn(u64, u64, u64, u64) -> SlidingWindow = SlidingWindow::new;
    let _: fn(&SlidingWindow) -> &Arena = SlidingWindow::arena;
    let _: fn(&mut SlidingWindow) -> &mut Arena = SlidingWindow::arena_mut;
    let _: fn(&mut SlidingWindow) -> u64 = SlidingWindow::current_hash;
    let _: fn(&mut SlidingWindow) -> u64 = SlidingWindow::final_hash;

    let mut window = SlidingWindow::new(4, 3, MERSENNE_61, 131);
    let oracle = PolynomialHash::default_hash();
    let mut decoded = b"seed".to_vec();
    window.append_bytes(&decoded);
    let saved: Vec<_> = (0..window.arena().node_count())
        .map(|index| {
            let id = NodeId::try_from(index).unwrap();
            (id, snapshot(window.arena().node(id)))
        })
        .collect();
    assert!(!saved.is_empty());

    // Observe handles through arena() only. Using arena_mut() here would mask
    // an accidental change to a compacting constructor by disabling compaction.
    for index in 0..12u8 {
        let bytes = [index, 255 - index, b'!'];
        window.append_bytes(&bytes);
        decoded.extend_from_slice(&bytes);
        assert_eq!(window.current_hash(), oracle.hash(&decoded));
        assert_eq!(window.pos(), decoded.len() as u64);
        assert!(window.window_len() <= 7);
        for (id, original) in &saved {
            assert_eq!(&snapshot(window.arena().node(*id)), original);
        }
    }

    window.append_copy(3, 9);
    for _ in 0..9 {
        decoded.push(decoded[decoded.len() - 3]);
    }
    assert_eq!(window.final_hash(), oracle.hash(&decoded));
    for (id, original) in &saved {
        assert_eq!(&snapshot(window.arena().node(*id)), original);
    }

    let retained = window.arena_mut().from_bytes(b"caller-owned bytes");
    let retained_shape = snapshot(window.arena().node(retained.unwrap()));
    window.append_bytes(b"another eviction");
    decoded.extend_from_slice(b"another eviction");
    assert_eq!(window.arena().to_bytes(retained), b"caller-owned bytes");
    assert_eq!(
        snapshot(window.arena().node(retained.unwrap())),
        retained_shape
    );
    assert_eq!(window.final_hash(), oracle.hash(&decoded));
    let maximum_match: u64 = window.m_max;
    assert_eq!(maximum_match, 3);
}
