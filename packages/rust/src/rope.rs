//! Hash rope data structure.
//!
//! Implements Part III of the framework:
//! - Definition 6: Leaf, Internal, RepeatNode
//! - Definition 7: Invariants I1-I9 (BB[2/7] weight balance)
//! - Theorems 6-10: Join, Split, Repeat, SubstrHash, Concat
//!
//! All nodes are immutable. Operations return new trees with structural
//! sharing via `Arc`.

use alloc::sync::Arc;
use alloc::vec;
use alloc::vec::Vec;

use crate::polynomial_hash::{mersenne_mod, mersenne_mul, phi, PolynomialHash};

// ---------------------------------------------------------------------------
// Balance parameter (Definition 7, Invariant I8)
// ---------------------------------------------------------------------------

const ALPHA_NUM: u64 = 2;
const ALPHA_DEN: u64 = 7;

// ---------------------------------------------------------------------------
// Definition 6: Node types
// ---------------------------------------------------------------------------

/// A node in the hash rope. `None` represents the empty rope.
pub type Node = Option<Arc<NodeInner>>;

/// Inner node variants.
#[derive(Debug, Clone)]
pub enum NodeInner {
    Leaf {
        data: Vec<u8>,
        hash_val: u64,
        len: u64,
    },
    Internal {
        left: Arc<NodeInner>,
        right: Arc<NodeInner>,
        hash_val: u64,
        len: u64,
        weight: u64,
    },
    Repeat {
        child: Arc<NodeInner>,
        reps: u64,
        hash_val: u64,
        len: u64,
        weight: u64,
    },
}

impl NodeInner {
    #[inline]
    pub fn hash_val(&self) -> u64 {
        match self {
            NodeInner::Leaf { hash_val, .. } => *hash_val,
            NodeInner::Internal { hash_val, .. } => *hash_val,
            NodeInner::Repeat { hash_val, .. } => *hash_val,
        }
    }

    #[inline]
    pub fn len(&self) -> u64 {
        match self {
            NodeInner::Leaf { len, .. } => *len,
            NodeInner::Internal { len, .. } => *len,
            NodeInner::Repeat { len, .. } => *len,
        }
    }

    #[inline]
    pub fn weight(&self) -> u64 {
        match self {
            NodeInner::Leaf { .. } => 1,
            NodeInner::Internal { weight, .. } => *weight,
            NodeInner::Repeat { weight, .. } => *weight,
        }
    }
}

/// Constructors that enforce invariants.
fn make_leaf(data: Vec<u8>, h: &PolynomialHash) -> Arc<NodeInner> {
    assert!(!data.is_empty(), "Leaf cannot be empty");
    let hash_val = h.hash(&data);
    let len = data.len() as u64;
    Arc::new(NodeInner::Leaf { data, hash_val, len })
}

fn make_internal(left: Arc<NodeInner>, right: Arc<NodeInner>, h: &mut PolynomialHash) -> Arc<NodeInner> {
    let len = left.len() + right.len();
    let weight = left.weight() + right.weight();
    let hash_val = h.hash_concat(left.hash_val(), right.len(), right.hash_val());
    Arc::new(NodeInner::Internal { left, right, hash_val, len, weight })
}

fn make_repeat_node(child: Arc<NodeInner>, reps: u64, h: &mut PolynomialHash) -> Arc<NodeInner> {
    assert!(reps >= 2, "RepeatNode reps must be >= 2, got {}", reps);
    let len = child.len() * reps;
    let weight = child.weight() * reps;
    let x_d = h.power(child.len());
    let phi_val = phi(reps, x_d, h.prime());
    let hash_val = mersenne_mul(child.hash_val(), phi_val, h.prime());
    Arc::new(NodeInner::Repeat { child, reps, hash_val, len, weight })
}

// ---------------------------------------------------------------------------
// Accessors
// ---------------------------------------------------------------------------

/// Total byte length of the string represented by the node.
#[inline]
pub fn rope_len(node: &Node) -> u64 {
    node.as_ref().map_or(0, |n| n.len())
}

/// Hash of the string represented by the node.
#[inline]
pub fn rope_hash(node: &Node) -> u64 {
    node.as_ref().map_or(0, |n| n.hash_val())
}

#[inline]
fn weight(node: &Node) -> u64 {
    node.as_ref().map_or(0, |n| n.weight())
}

// ---------------------------------------------------------------------------
// Balance helpers (Invariant I8)
// ---------------------------------------------------------------------------

#[inline]
fn is_balanced(left: &Arc<NodeInner>, right: &Arc<NodeInner>) -> bool {
    let total = left.weight() + right.weight();
    if total <= 2 {
        return true;
    }
    let wl = left.weight();
    ALPHA_NUM * total <= ALPHA_DEN * wl && ALPHA_DEN * wl <= (ALPHA_DEN - ALPHA_NUM) * total
}

fn rotate_right(left: &Arc<NodeInner>, right: &Arc<NodeInner>, h: &mut PolynomialHash) -> Arc<NodeInner> {
    if let NodeInner::Internal { left: a, right: b, .. } = left.as_ref() {
        let new_right = make_internal(b.clone(), right.clone(), h);
        make_internal(a.clone(), new_right, h)
    } else {
        make_internal(left.clone(), right.clone(), h)
    }
}

fn rotate_left(left: &Arc<NodeInner>, right: &Arc<NodeInner>, h: &mut PolynomialHash) -> Arc<NodeInner> {
    if let NodeInner::Internal { left: b, right: c, .. } = right.as_ref() {
        let new_left = make_internal(left.clone(), b.clone(), h);
        make_internal(new_left, c.clone(), h)
    } else {
        make_internal(left.clone(), right.clone(), h)
    }
}

fn rebalance(left: Arc<NodeInner>, right: Arc<NodeInner>, h: &mut PolynomialHash) -> Arc<NodeInner> {
    let wl = left.weight();
    let wr = right.weight();
    let total = wl + wr;

    if total <= 2 || is_balanced(&left, &right) {
        return make_internal(left, right, h);
    }

    if wl * ALPHA_DEN < ALPHA_NUM * total {
        // Left too light -> rotate left
        if let NodeInner::Internal { left: rl, .. } = right.as_ref() {
            if rl.weight() * ALPHA_DEN > (ALPHA_DEN - ALPHA_NUM) * wr {
                // Double rotation
                let new_right = rotate_right(&left, &right, h);
                if let NodeInner::Internal { left: nr_l, right: nr_r, .. } = new_right.as_ref() {
                    return rotate_left(nr_l, nr_r, h);
                }
            }
            return rotate_left(&left, &right, h);
        }
        return make_internal(left, right, h);
    }

    if wr * ALPHA_DEN < ALPHA_NUM * total {
        // Right too light -> rotate right
        if let NodeInner::Internal { right: lr, .. } = left.as_ref() {
            if lr.weight() * ALPHA_DEN > (ALPHA_DEN - ALPHA_NUM) * wl {
                // Double rotation
                let new_left = rotate_left(&left, &right, h);
                if let NodeInner::Internal { left: nl_l, right: nl_r, .. } = new_left.as_ref() {
                    return rotate_right(nl_l, nl_r, h);
                }
            }
            return rotate_right(&left, &right, h);
        }
        return make_internal(left, right, h);
    }

    make_internal(left, right, h)
}

// ---------------------------------------------------------------------------
// Theorem 6 / Lemma 8: Recursive Join
// ---------------------------------------------------------------------------

fn join(left: Arc<NodeInner>, right: Arc<NodeInner>, h: &mut PolynomialHash) -> Arc<NodeInner> {
    if is_balanced(&left, &right) {
        return make_internal(left, right, h);
    }

    let wl = left.weight();
    let wr = right.weight();

    if wl > wr {
        match left.as_ref() {
            NodeInner::Internal { left: ll, right: lr, .. } => {
                let new_right = join(lr.clone(), right, h);
                rebalance(ll.clone(), new_right, h)
            }
            NodeInner::Repeat { .. } => {
                let mid = left.len() / 2;
                let (ll, lr) = split_inner(&left, mid, h);
                let lr = lr.expect("split of repeat should produce right half");
                let ll = ll.expect("split of repeat should produce left half");
                let new_right = join(lr, right, h);
                join(ll, new_right, h)
            }
            NodeInner::Leaf { .. } => make_internal(left, right, h),
        }
    } else {
        match right.as_ref() {
            NodeInner::Internal { left: rl, right: rr, .. } => {
                let new_left = join(left, rl.clone(), h);
                rebalance(new_left, rr.clone(), h)
            }
            NodeInner::Repeat { .. } => {
                let mid = right.len() / 2;
                let (rl, rr) = split_inner(&right, mid, h);
                let rl = rl.expect("split of repeat should produce left half");
                let rr = rr.expect("split of repeat should produce right half");
                let new_left = join(left, rl, h);
                join(new_left, rr, h)
            }
            NodeInner::Leaf { .. } => make_internal(left, right, h),
        }
    }
}

// ---------------------------------------------------------------------------
// Theorem 10: Concat = Join
// ---------------------------------------------------------------------------

/// Concatenate two ropes. Returns a balanced rope representing `left || right`.
pub fn rope_concat(left: &Node, right: &Node, h: &mut PolynomialHash) -> Node {
    match (left, right) {
        (None, _) => right.clone(),
        (_, None) => left.clone(),
        (Some(l), Some(r)) => Some(join(l.clone(), r.clone(), h)),
    }
}

// ---------------------------------------------------------------------------
// Helper: make_repeat
// ---------------------------------------------------------------------------

fn make_repeat(child: Arc<NodeInner>, reps: u64, h: &mut PolynomialHash) -> Node {
    match reps {
        0 => None,
        1 => Some(child),
        _ => Some(make_repeat_node(child, reps, h)),
    }
}

// ---------------------------------------------------------------------------
// Theorem 7: Split
// ---------------------------------------------------------------------------

/// Split a rope at byte position `pos`.
///
/// Returns `(left, right)` where left has `pos` bytes, right has the rest.
pub fn rope_split(node: &Node, pos: u64, h: &mut PolynomialHash) -> (Node, Node) {
    match node {
        None => (None, None),
        Some(n) => {
            if pos == 0 {
                return (None, Some(n.clone()));
            }
            if pos >= n.len() {
                return (Some(n.clone()), None);
            }
            let (l, r) = split_inner(n, pos, h);
            (l, r)
        }
    }
}

fn split_inner(node: &Arc<NodeInner>, pos: u64, h: &mut PolynomialHash) -> (Node, Node) {
    match node.as_ref() {
        NodeInner::Leaf { data, .. } => {
            let pos = pos as usize;
            let left = make_leaf(data[..pos].to_vec(), h);
            let right = make_leaf(data[pos..].to_vec(), h);
            (Some(left), Some(right))
        }
        NodeInner::Internal { left, right, .. } => {
            let ll = left.len();
            if pos == ll {
                return (Some(left.clone()), Some(right.clone()));
            }
            if pos < ll {
                let (l1, l2) = split_inner(left, pos, h);
                let rejoined = rope_concat(&l2, &Some(right.clone()), h);
                (l1, rejoined)
            } else {
                let (r1, r2) = split_inner(right, pos - ll, h);
                let rejoined = rope_concat(&Some(left.clone()), &r1, h);
                (rejoined, r2)
            }
        }
        NodeInner::Repeat { child, reps, .. } => {
            split_repeat(child, *reps, pos, h)
        }
    }
}

fn split_repeat(child: &Arc<NodeInner>, reps: u64, pos: u64, h: &mut PolynomialHash) -> (Node, Node) {
    let d = child.len();
    let m = pos / d;
    let r = pos % d;

    if r == 0 {
        let left = make_repeat(child.clone(), m, h);
        let right = make_repeat(child.clone(), reps - m, h);
        (left, right)
    } else {
        let (child_left, child_right) = split_inner(child, r, h);

        let left_rep = make_repeat(child.clone(), m, h);
        let left = rope_concat(&left_rep, &child_left, h);

        let right_rep = make_repeat(child.clone(), reps - m - 1, h);
        let right = rope_concat(&child_right, &right_rep, h);

        (left, right)
    }
}

// ---------------------------------------------------------------------------
// Theorem 8: Repeat
// ---------------------------------------------------------------------------

/// Create a RepeatNode representing `node^q`.
pub fn rope_repeat(node: &Node, q: u64, h: &mut PolynomialHash) -> Node {
    match node {
        None => None,
        Some(n) => make_repeat(n.clone(), q, h),
    }
}

// ---------------------------------------------------------------------------
// Theorem 9: SubstrHash (allocation-free)
// ---------------------------------------------------------------------------

/// Compute `H(S[start..start+length-1])` without allocating new nodes.
pub fn rope_substr_hash(node: &Node, start: u64, length: u64, h: &mut PolynomialHash) -> u64 {
    match node {
        None => 0,
        Some(n) if length == 0 => 0,
        Some(n) => hash_range(n, start, length, h),
    }
}

fn hash_range(node: &Arc<NodeInner>, start: u64, length: u64, h: &mut PolynomialHash) -> u64 {
    let p = h.prime();

    match node.as_ref() {
        NodeInner::Leaf { data, .. } => {
            let s = start as usize;
            let e = s + length as usize;
            h.hash(&data[s..e])
        }
        NodeInner::Internal { left, right, .. } => {
            let ll = left.len();
            if start + length <= ll {
                return hash_range(left, start, length, h);
            }
            if start >= ll {
                return hash_range(right, start - ll, length, h);
            }
            // Spanning
            let l_len = ll - start;
            let r_len = length - l_len;
            let h_l = hash_range(left, start, l_len, h);
            let h_r = hash_range(right, 0, r_len, h);
            h.hash_concat(h_l, r_len, h_r)
        }
        NodeInner::Repeat { child, .. } => {
            let d = child.len();
            let first_copy = start / d;
            let start_in_copy = start % d;
            let end = start + length - 1;
            let last_copy = end / d;
            let end_in_copy = end % d;

            if first_copy == last_copy {
                return hash_range(child, start_in_copy, length, h);
            }

            // Tail of first copy
            let tail_len = d - start_in_copy;
            let h_tail = hash_range(child, start_in_copy, tail_len, h);

            // Full copies in the middle
            let full_copies = last_copy - first_copy - 1;
            let (h_full, full_len) = if full_copies > 0 {
                let x_d = h.power(d);
                let phi_val = phi(full_copies, x_d, p);
                (mersenne_mul(child.hash_val(), phi_val, p), full_copies * d)
            } else {
                (0u64, 0u64)
            };

            // Head of last copy
            let head_len = end_in_copy + 1;
            let h_head = hash_range(child, 0, head_len, h);

            // Combine
            let x_full_head = h.power(full_len + head_len);
            let x_head = h.power(head_len);
            mersenne_mod(
                x_full_head as u128 * h_tail as u128
                    + x_head as u128 * h_full as u128
                    + h_head as u128,
                p,
            )
        }
    }
}

// ---------------------------------------------------------------------------
// Reconstruction and construction
// ---------------------------------------------------------------------------

/// Create a rope from a byte slice. Returns `None` for empty.
pub fn rope_from_bytes(data: &[u8], h: &PolynomialHash) -> Node {
    if data.is_empty() {
        None
    } else {
        Some(make_leaf(data.to_vec(), h))
    }
}

/// Reconstruct the byte string from a rope.
pub fn rope_to_bytes(node: &Node) -> Vec<u8> {
    let mut parts = Vec::new();
    if let Some(n) = node {
        collect_bytes(n, &mut parts);
    }
    parts
}

fn collect_bytes(node: &Arc<NodeInner>, parts: &mut Vec<u8>) {
    match node.as_ref() {
        NodeInner::Leaf { data, .. } => parts.extend_from_slice(data),
        NodeInner::Internal { left, right, .. } => {
            collect_bytes(left, parts);
            collect_bytes(right, parts);
        }
        NodeInner::Repeat { child, reps, .. } => {
            for _ in 0..*reps {
                collect_bytes(child, parts);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

/// Validate all rope invariants (Definition 7). Panics on violation.
pub fn validate_rope(node: &Node) {
    if let Some(n) = node {
        validate_inner(n);
    }
}

fn validate_inner(node: &Arc<NodeInner>) {
    match node.as_ref() {
        NodeInner::Leaf { data, len, .. } => {
            assert_eq!(*len, data.len() as u64, "Leaf len mismatch");
        }
        NodeInner::Internal { left, right, len, weight, .. } => {
            assert_eq!(*len, left.len() + right.len(), "Internal len mismatch");
            assert_eq!(*weight, left.weight() + right.weight(), "Internal weight mismatch");
            let total = *weight;
            let wl = left.weight();
            if total > 2 {
                assert!(
                    ALPHA_NUM * total <= ALPHA_DEN * wl,
                    "Left child too light: {}/{}", wl, total
                );
                assert!(
                    ALPHA_DEN * wl <= (ALPHA_DEN - ALPHA_NUM) * total,
                    "Left child too heavy: {}/{}", wl, total
                );
            }
            validate_inner(left);
            validate_inner(right);
        }
        NodeInner::Repeat { child, reps, len, weight, .. } => {
            assert_eq!(*len, child.len() * reps, "RepeatNode len mismatch");
            assert_eq!(*weight, child.weight() * reps, "RepeatNode weight mismatch");
            assert!(*reps >= 2, "RepeatNode reps must be >= 2, got {}", reps);
            validate_inner(child);
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn h() -> PolynomialHash {
        PolynomialHash::default_hash()
    }

    #[test]
    fn test_leaf_basic() {
        let ph = h();
        let leaf = make_leaf(b"hello".to_vec(), &ph);
        assert_eq!(leaf.len(), 5);
        assert_eq!(leaf.weight(), 1);
        assert_eq!(leaf.hash_val(), ph.hash(b"hello"));
    }

    #[test]
    fn test_concat_two_leaves() {
        let mut ph = h();
        let a: Node = Some(make_leaf(b"hello ".to_vec(), &ph));
        let b: Node = Some(make_leaf(b"world".to_vec(), &ph));
        let ab = rope_concat(&a, &b, &mut ph);
        assert_eq!(rope_len(&ab), 11);
        assert_eq!(rope_hash(&ab), ph.hash(b"hello world"));
    }

    #[test]
    fn test_concat_none() {
        let mut ph = h();
        let a: Node = Some(make_leaf(b"test".to_vec(), &ph));
        assert!(rope_concat(&None, &a, &mut ph).is_some());
        assert!(rope_concat(&a, &None, &mut ph).is_some());
        assert!(rope_concat(&None, &None, &mut ph).is_none());
    }

    #[test]
    fn test_split_preserves_hash() {
        let mut ph = h();
        let data = b"the quick brown fox";
        let node = rope_from_bytes(data, &ph);
        for pos in 0..=data.len() as u64 {
            let (left, right) = rope_split(&node, pos, &mut ph);
            let rejoined = rope_concat(&left, &right, &mut ph);
            assert_eq!(rope_hash(&rejoined), ph.hash(data), "Failed at pos={}", pos);
        }
    }

    #[test]
    fn test_repeat_hash() {
        let mut ph = h();
        let node = rope_from_bytes(b"ab", &ph);
        let repeated = rope_repeat(&node, 100, &mut ph);
        let expected: Vec<u8> = b"ab".iter().copied().cycle().take(200).collect();
        assert_eq!(rope_hash(&repeated), ph.hash(&expected));
        assert_eq!(rope_len(&repeated), 200);
    }

    #[test]
    fn test_substr_hash() {
        let mut ph = h();
        let node = rope_from_bytes(b"hello world", &ph);
        assert_eq!(rope_substr_hash(&node, 0, 5, &mut ph), ph.hash(b"hello"));
        assert_eq!(rope_substr_hash(&node, 6, 5, &mut ph), ph.hash(b"world"));
        assert_eq!(rope_substr_hash(&node, 4, 4, &mut ph), ph.hash(b"o wo"));
    }

    #[test]
    fn test_to_bytes_roundtrip() {
        let ph = h();
        let data = b"round trip test";
        let node = rope_from_bytes(data, &ph);
        assert_eq!(rope_to_bytes(&node), data);
    }

    #[test]
    fn test_many_concats_balanced() {
        let mut ph = h();
        let mut node: Node = None;
        for i in 0..100u8 {
            let leaf: Node = Some(make_leaf(vec![i], &ph));
            node = rope_concat(&node, &leaf, &mut ph);
        }
        validate_rope(&node);
        assert_eq!(rope_len(&node), 100);
    }

    #[test]
    fn test_split_repeat_on_boundary() {
        let mut ph = h();
        let leaf: Node = Some(make_leaf(b"abc".to_vec(), &ph));
        let rep = rope_repeat(&leaf, 4, &mut ph); // "abcabcabcabc"
        let (left, right) = rope_split(&rep, 6, &mut ph);
        assert_eq!(rope_to_bytes(&left), b"abcabc");
        assert_eq!(rope_to_bytes(&right), b"abcabc");
    }

    #[test]
    fn test_split_repeat_within() {
        let mut ph = h();
        let leaf: Node = Some(make_leaf(b"abcd".to_vec(), &ph));
        let rep = rope_repeat(&leaf, 3, &mut ph); // "abcdabcdabcd"
        let (left, right) = rope_split(&rep, 5, &mut ph);
        assert_eq!(rope_to_bytes(&left), b"abcda");
        assert_eq!(rope_to_bytes(&right), b"bcdabcd");
    }

    #[test]
    fn test_associativity() {
        let mut ph = h();
        let a: Node = Some(make_leaf(b"aaa".to_vec(), &ph));
        let b: Node = Some(make_leaf(b"bbb".to_vec(), &ph));
        let c: Node = Some(make_leaf(b"ccc".to_vec(), &ph));
        let ab = rope_concat(&a, &b, &mut ph);
        let ab_c = rope_concat(&ab, &c, &mut ph);
        let bc = rope_concat(&b, &c, &mut ph);
        let a_bc = rope_concat(&a, &bc, &mut ph);
        assert_eq!(rope_hash(&ab_c), rope_hash(&a_bc));
        assert_eq!(rope_hash(&ab_c), ph.hash(b"aaabbbccc"));
    }
}
