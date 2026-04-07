//! Hash rope data structure with arena-based allocation.
//!
//! Implements Part III of the framework:
//! - Definition 6: Leaf, Internal, RepeatNode
//! - Definition 7: Invariants I1-I9 (BB[2/7] weight balance)
//! - Theorems 6-10: Join, Split, Repeat, SubstrHash, Concat
//!
//! Nodes are stored in a contiguous arena (`Vec<NodeInner>`) and referenced
//! by `u32` indices. No reference counting, no atomic operations.
//! Structural sharing works via shared indices. The arena is dropped as a unit.

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

/// Index into the arena. Max ~4 billion nodes.
pub type NodeId = u32;

/// A node reference. `None` represents the empty rope.
pub type Node = Option<NodeId>;

/// Inner node variants stored in the arena.
#[derive(Debug, Clone)]
pub enum NodeInner {
    Leaf {
        data: Vec<u8>,
        hash_val: u64,
        len: u64,
    },
    Internal {
        left: NodeId,
        right: NodeId,
        hash_val: u64,
        len: u64,
        weight: u64,
    },
    Repeat {
        child: NodeId,
        reps: u64,
        hash_val: u64,
        len: u64,
        weight: u64,
    },
}

// ---------------------------------------------------------------------------
// Arena: owns all nodes and the polynomial hash state
// ---------------------------------------------------------------------------

/// Arena-based rope storage. All nodes live here.
pub struct Arena {
    nodes: Vec<NodeInner>,
    h: PolynomialHash,
}

impl Arena {
    /// Create a new arena with default polynomial hash parameters.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            h: PolynomialHash::default_hash(),
        }
    }

    /// Create an arena with custom hash parameters.
    pub fn with_hash(prime: u64, base: u64) -> Self {
        Self {
            nodes: Vec::new(),
            h: PolynomialHash::new(prime, base),
        }
    }

    /// Create an arena with a specific PolynomialHash instance.
    pub fn with_hasher(h: PolynomialHash) -> Self {
        Self {
            nodes: Vec::new(),
            h,
        }
    }

    /// Access the polynomial hash state.
    #[inline]
    pub fn hasher(&self) -> &PolynomialHash {
        &self.h
    }

    /// Mutable access to the polynomial hash state.
    #[inline]
    pub fn hasher_mut(&mut self) -> &mut PolynomialHash {
        &mut self.h
    }

    /// Number of nodes allocated in the arena.
    #[inline]
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Access a node by id.
    #[inline]
    pub fn node(&self, id: NodeId) -> &NodeInner {
        &self.nodes[id as usize]
    }

    // -----------------------------------------------------------------------
    // Node accessors
    // -----------------------------------------------------------------------

    #[inline]
    fn hash_val(&self, id: NodeId) -> u64 {
        match &self.nodes[id as usize] {
            NodeInner::Leaf { hash_val, .. } => *hash_val,
            NodeInner::Internal { hash_val, .. } => *hash_val,
            NodeInner::Repeat { hash_val, .. } => *hash_val,
        }
    }

    #[inline]
    fn node_len(&self, id: NodeId) -> u64 {
        match &self.nodes[id as usize] {
            NodeInner::Leaf { len, .. } => *len,
            NodeInner::Internal { len, .. } => *len,
            NodeInner::Repeat { len, .. } => *len,
        }
    }

    #[inline]
    fn node_weight(&self, id: NodeId) -> u64 {
        match &self.nodes[id as usize] {
            NodeInner::Leaf { .. } => 1,
            NodeInner::Internal { weight, .. } => *weight,
            NodeInner::Repeat { weight, .. } => *weight,
        }
    }

    // -----------------------------------------------------------------------
    // Allocation
    // -----------------------------------------------------------------------

    #[inline]
    fn alloc(&mut self, inner: NodeInner) -> NodeId {
        let id = self.nodes.len() as u32;
        self.nodes.push(inner);
        id
    }

    // -----------------------------------------------------------------------
    // Node constructors (enforce invariants)
    // -----------------------------------------------------------------------

    fn make_leaf(&mut self, data: Vec<u8>) -> NodeId {
        assert!(!data.is_empty(), "Leaf cannot be empty");
        let hash_val = self.h.hash(&data);
        let len = data.len() as u64;
        self.alloc(NodeInner::Leaf { data, hash_val, len })
    }

    #[inline]
    fn make_internal(&mut self, left: NodeId, right: NodeId) -> NodeId {
        let len = self.node_len(left) + self.node_len(right);
        let weight = self.node_weight(left) + self.node_weight(right);
        let hash_val = self.h.hash_concat(
            self.hash_val(left),
            self.node_len(right),
            self.hash_val(right),
        );
        self.alloc(NodeInner::Internal { left, right, hash_val, len, weight })
    }

    fn make_repeat_node(&mut self, child: NodeId, reps: u64) -> NodeId {
        assert!(reps >= 2, "RepeatNode reps must be >= 2, got {}", reps);
        let child_len = self.node_len(child);
        let child_weight = self.node_weight(child);
        let child_hash = self.hash_val(child);
        let len = child_len * reps;
        let weight = child_weight * reps;
        let x_d = self.h.power(child_len);
        let phi_val = phi(reps, x_d, self.h.prime());
        let hash_val = mersenne_mul(child_hash, phi_val, self.h.prime());
        self.alloc(NodeInner::Repeat { child, reps, hash_val, len, weight })
    }

    fn make_repeat(&mut self, child: NodeId, reps: u64) -> Node {
        match reps {
            0 => None,
            1 => Some(child),
            _ => Some(self.make_repeat_node(child, reps)),
        }
    }

    // -----------------------------------------------------------------------
    // Public accessors
    // -----------------------------------------------------------------------

    /// Total byte length of the string represented by the node.
    #[inline]
    pub fn len(&self, node: Node) -> u64 {
        node.map_or(0, |id| self.node_len(id))
    }

    /// Hash of the string represented by the node.
    #[inline]
    pub fn hash(&self, node: Node) -> u64 {
        node.map_or(0, |id| self.hash_val(id))
    }

    /// Weight (number of leaves) of the rope.
    #[inline]
    pub fn weight(&self, node: Node) -> u64 {
        node.map_or(0, |id| self.node_weight(id))
    }

    /// Height of the rope tree. Leaves have height 0, `None` has height 0.
    pub fn height(&self, node: Node) -> u64 {
        node.map_or(0, |id| self.node_height(id))
    }

    fn node_height(&self, id: NodeId) -> u64 {
        match &self.nodes[id as usize] {
            NodeInner::Leaf { .. } => 0,
            NodeInner::Internal { left, right, .. } => {
                1 + core::cmp::max(self.node_height(*left), self.node_height(*right))
            }
            NodeInner::Repeat { child, .. } => 1 + self.node_height(*child),
        }
    }

    // -----------------------------------------------------------------------
    // Balance helpers (Invariant I8)
    // -----------------------------------------------------------------------

    #[inline]
    fn is_balanced(&self, left: NodeId, right: NodeId) -> bool {
        let total = self.node_weight(left) + self.node_weight(right);
        if total <= 2 {
            return true;
        }
        let wl = self.node_weight(left);
        ALPHA_NUM * total <= ALPHA_DEN * wl && ALPHA_DEN * wl <= (ALPHA_DEN - ALPHA_NUM) * total
    }

    fn rotate_right(&mut self, left: NodeId, right: NodeId) -> NodeId {
        if let NodeInner::Internal { left: a, right: b, .. } = self.nodes[left as usize].clone() {
            let new_right = self.make_internal(b, right);
            self.make_internal(a, new_right)
        } else {
            self.make_internal(left, right)
        }
    }

    fn rotate_left(&mut self, left: NodeId, right: NodeId) -> NodeId {
        if let NodeInner::Internal { left: b, right: c, .. } = self.nodes[right as usize].clone() {
            let new_left = self.make_internal(left, b);
            self.make_internal(new_left, c)
        } else {
            self.make_internal(left, right)
        }
    }

    fn rebalance(&mut self, left: NodeId, right: NodeId) -> NodeId {
        let wl = self.node_weight(left);
        let wr = self.node_weight(right);
        let total = wl + wr;

        if total <= 2 || self.is_balanced(left, right) {
            return self.make_internal(left, right);
        }

        if wl * ALPHA_DEN < ALPHA_NUM * total {
            // Left too light -> rotate left
            if let NodeInner::Internal { left: rl, .. } = self.nodes[right as usize].clone() {
                if self.node_weight(rl) * ALPHA_DEN > (ALPHA_DEN - ALPHA_NUM) * wr {
                    // Double rotation
                    let new_right = self.rotate_right(left, right);
                    if let NodeInner::Internal { left: nr_l, right: nr_r, .. } = self.nodes[new_right as usize].clone() {
                        return self.rotate_left(nr_l, nr_r);
                    }
                }
                return self.rotate_left(left, right);
            }
            return self.make_internal(left, right);
        }

        if wr * ALPHA_DEN < ALPHA_NUM * total {
            // Right too light -> rotate right
            if let NodeInner::Internal { right: lr, .. } = self.nodes[left as usize].clone() {
                if self.node_weight(lr) * ALPHA_DEN > (ALPHA_DEN - ALPHA_NUM) * wl {
                    // Double rotation
                    let new_left = self.rotate_left(left, right);
                    if let NodeInner::Internal { left: nl_l, right: nl_r, .. } = self.nodes[new_left as usize].clone() {
                        return self.rotate_right(nl_l, nl_r);
                    }
                }
                return self.rotate_right(left, right);
            }
            return self.make_internal(left, right);
        }

        self.make_internal(left, right)
    }

    // -----------------------------------------------------------------------
    // Theorem 6 / Lemma 8: Recursive Join
    // -----------------------------------------------------------------------

    fn join(&mut self, left: NodeId, right: NodeId) -> NodeId {
        if self.is_balanced(left, right) {
            return self.make_internal(left, right);
        }

        let wl = self.node_weight(left);
        let wr = self.node_weight(right);

        if wl > wr {
            match self.nodes[left as usize].clone() {
                NodeInner::Internal { left: ll, right: lr, .. } => {
                    let new_right = self.join(lr, right);
                    self.rebalance(ll, new_right)
                }
                NodeInner::Repeat { .. } => {
                    let mid = self.node_len(left) / 2;
                    let (ll, lr) = self.split_inner(left, mid);
                    let lr = lr.expect("split of repeat should produce right half");
                    let ll = ll.expect("split of repeat should produce left half");
                    let new_right = self.join(lr, right);
                    self.join(ll, new_right)
                }
                NodeInner::Leaf { .. } => self.make_internal(left, right),
            }
        } else {
            match self.nodes[right as usize].clone() {
                NodeInner::Internal { left: rl, right: rr, .. } => {
                    let new_left = self.join(left, rl);
                    self.rebalance(new_left, rr)
                }
                NodeInner::Repeat { .. } => {
                    let mid = self.node_len(right) / 2;
                    let (rl, rr) = self.split_inner(right, mid);
                    let rl = rl.expect("split of repeat should produce left half");
                    let rr = rr.expect("split of repeat should produce right half");
                    let new_left = self.join(left, rl);
                    self.join(new_left, rr)
                }
                NodeInner::Leaf { .. } => self.make_internal(left, right),
            }
        }
    }

    // -----------------------------------------------------------------------
    // Theorem 10: Concat = Join
    // -----------------------------------------------------------------------

    /// Concatenate two ropes. Returns a balanced rope representing `left || right`.
    pub fn concat(&mut self, left: Node, right: Node) -> Node {
        match (left, right) {
            (None, _) => right,
            (_, None) => left,
            (Some(l), Some(r)) => Some(self.join(l, r)),
        }
    }

    // -----------------------------------------------------------------------
    // Theorem 7: Split
    // -----------------------------------------------------------------------

    /// Split a rope at byte position `pos`.
    ///
    /// Returns `(left, right)` where left has `pos` bytes, right has the rest.
    pub fn split(&mut self, node: Node, pos: u64) -> (Node, Node) {
        match node {
            None => (None, None),
            Some(id) => {
                if pos == 0 {
                    return (None, Some(id));
                }
                if pos >= self.node_len(id) {
                    return (Some(id), None);
                }
                self.split_inner(id, pos)
            }
        }
    }

    fn split_inner(&mut self, id: NodeId, pos: u64) -> (Node, Node) {
        match self.nodes[id as usize].clone() {
            NodeInner::Leaf { data, .. } => {
                let pos = pos as usize;
                let left = self.make_leaf(data[..pos].to_vec());
                let right = self.make_leaf(data[pos..].to_vec());
                (Some(left), Some(right))
            }
            NodeInner::Internal { left, right, .. } => {
                let ll = self.node_len(left);
                if pos == ll {
                    return (Some(left), Some(right));
                }
                if pos < ll {
                    let (l1, l2) = self.split_inner(left, pos);
                    let rejoined = self.concat(l2, Some(right));
                    (l1, rejoined)
                } else {
                    let (r1, r2) = self.split_inner(right, pos - ll);
                    let rejoined = self.concat(Some(left), r1);
                    (rejoined, r2)
                }
            }
            NodeInner::Repeat { child, reps, .. } => {
                self.split_repeat(child, reps, pos)
            }
        }
    }

    fn split_repeat(&mut self, child: NodeId, reps: u64, pos: u64) -> (Node, Node) {
        let d = self.node_len(child);
        let m = pos / d;
        let r = pos % d;

        if r == 0 {
            let left = self.make_repeat(child, m);
            let right = self.make_repeat(child, reps - m);
            (left, right)
        } else {
            let (child_left, child_right) = self.split_inner(child, r);

            let left_rep = self.make_repeat(child, m);
            let left = self.concat(left_rep, child_left);

            let right_rep = self.make_repeat(child, reps - m - 1);
            let right = self.concat(child_right, right_rep);

            (left, right)
        }
    }

    // -----------------------------------------------------------------------
    // Theorem 8: Repeat
    // -----------------------------------------------------------------------

    /// Create a RepeatNode representing `node^q`.
    pub fn repeat(&mut self, node: Node, q: u64) -> Node {
        match node {
            None => None,
            Some(id) => self.make_repeat(id, q),
        }
    }

    // -----------------------------------------------------------------------
    // Theorem 9: SubstrHash (allocation-free)
    // -----------------------------------------------------------------------

    /// Compute `H(S[start..start+length-1])` without allocating new nodes.
    pub fn substr_hash(&mut self, node: Node, start: u64, length: u64) -> u64 {
        match node {
            None => 0,
            Some(_) if length == 0 => 0,
            Some(id) => self.hash_range(id, start, length),
        }
    }

    fn hash_range(&mut self, id: NodeId, start: u64, length: u64) -> u64 {
        let p = self.h.prime();

        match self.nodes[id as usize].clone() {
            NodeInner::Leaf { data, .. } => {
                let s = start as usize;
                let e = s + length as usize;
                self.h.hash(&data[s..e])
            }
            NodeInner::Internal { left, right, .. } => {
                let ll = self.node_len(left);
                if start + length <= ll {
                    return self.hash_range(left, start, length);
                }
                if start >= ll {
                    return self.hash_range(right, start - ll, length);
                }
                // Spanning
                let l_len = ll - start;
                let r_len = length - l_len;
                let h_l = self.hash_range(left, start, l_len);
                let h_r = self.hash_range(right, 0, r_len);
                self.h.hash_concat(h_l, r_len, h_r)
            }
            NodeInner::Repeat { child, .. } => {
                let d = self.node_len(child);
                let child_hash = self.hash_val(child);
                let first_copy = start / d;
                let start_in_copy = start % d;
                let end = start + length - 1;
                let last_copy = end / d;
                let end_in_copy = end % d;

                if first_copy == last_copy {
                    return self.hash_range(child, start_in_copy, length);
                }

                // Tail of first copy
                let tail_len = d - start_in_copy;
                let h_tail = self.hash_range(child, start_in_copy, tail_len);

                // Full copies in the middle
                let full_copies = last_copy - first_copy - 1;
                let (h_full, full_len) = if full_copies > 0 {
                    let x_d = self.h.power(d);
                    let phi_val = phi(full_copies, x_d, p);
                    (mersenne_mul(child_hash, phi_val, p), full_copies * d)
                } else {
                    (0u64, 0u64)
                };

                // Head of last copy
                let head_len = end_in_copy + 1;
                let h_head = self.hash_range(child, 0, head_len);

                // Combine
                let x_full_head = self.h.power(full_len + head_len);
                let x_head = self.h.power(head_len);
                mersenne_mod(
                    x_full_head as u128 * h_tail as u128
                        + x_head as u128 * h_full as u128
                        + h_head as u128,
                    p,
                )
            }
        }
    }

    // -----------------------------------------------------------------------
    // Construction and reconstruction
    // -----------------------------------------------------------------------

    /// Create a rope from a byte slice. Returns `None` for empty.
    pub fn from_bytes(&mut self, data: &[u8]) -> Node {
        if data.is_empty() {
            None
        } else {
            Some(self.make_leaf(data.to_vec()))
        }
    }

    /// Reconstruct the byte string from a rope.
    pub fn to_bytes(&self, node: Node) -> Vec<u8> {
        let mut parts = Vec::new();
        if let Some(id) = node {
            self.collect_bytes(id, &mut parts);
        }
        parts
    }

    fn collect_bytes(&self, id: NodeId, parts: &mut Vec<u8>) {
        match &self.nodes[id as usize] {
            NodeInner::Leaf { data, .. } => parts.extend_from_slice(data),
            NodeInner::Internal { left, right, .. } => {
                self.collect_bytes(*left, parts);
                self.collect_bytes(*right, parts);
            }
            NodeInner::Repeat { child, reps, .. } => {
                for _ in 0..*reps {
                    self.collect_bytes(*child, parts);
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Validation
    // -----------------------------------------------------------------------

    /// Validate all rope invariants (Definition 7). Panics on violation.
    pub fn validate(&self, node: Node) {
        if let Some(id) = node {
            self.validate_inner(id);
        }
    }

    fn validate_inner(&self, id: NodeId) {
        match &self.nodes[id as usize] {
            NodeInner::Leaf { data, len, .. } => {
                assert_eq!(*len, data.len() as u64, "Leaf len mismatch");
            }
            NodeInner::Internal { left, right, len, weight, .. } => {
                assert_eq!(*len, self.node_len(*left) + self.node_len(*right), "Internal len mismatch");
                assert_eq!(*weight, self.node_weight(*left) + self.node_weight(*right), "Internal weight mismatch");
                let total = *weight;
                let wl = self.node_weight(*left);
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
                self.validate_inner(*left);
                self.validate_inner(*right);
            }
            NodeInner::Repeat { child, reps, len, weight, .. } => {
                assert_eq!(*len, self.node_len(*child) * reps, "RepeatNode len mismatch");
                assert_eq!(*weight, self.node_weight(*child) * reps, "RepeatNode weight mismatch");
                assert!(*reps >= 2, "RepeatNode reps must be >= 2, got {}", reps);
                self.validate_inner(*child);
            }
        }
    }

    /// Hash a raw byte slice using this arena's polynomial hash.
    pub fn hash_bytes(&self, data: &[u8]) -> u64 {
        self.h.hash(data)
    }
}

// ---------------------------------------------------------------------------
// Free-function wrappers for backward compatibility
// ---------------------------------------------------------------------------

/// Total byte length of the string represented by the node.
#[inline]
pub fn rope_len(arena: &Arena, node: Node) -> u64 {
    arena.len(node)
}

/// Hash of the string represented by the node.
#[inline]
pub fn rope_hash(arena: &Arena, node: Node) -> u64 {
    arena.hash(node)
}

/// Height of the rope tree.
#[inline]
pub fn rope_height(arena: &Arena, node: Node) -> u64 {
    arena.height(node)
}

/// Concatenate two ropes.
#[inline]
pub fn rope_concat(arena: &mut Arena, left: Node, right: Node) -> Node {
    arena.concat(left, right)
}

/// Split a rope at byte position `pos`.
#[inline]
pub fn rope_split(arena: &mut Arena, node: Node, pos: u64) -> (Node, Node) {
    arena.split(node, pos)
}

/// Create a RepeatNode representing `node^q`.
#[inline]
pub fn rope_repeat(arena: &mut Arena, node: Node, q: u64) -> Node {
    arena.repeat(node, q)
}

/// Compute substring hash without allocating new nodes.
#[inline]
pub fn rope_substr_hash(arena: &mut Arena, node: Node, start: u64, length: u64) -> u64 {
    arena.substr_hash(node, start, length)
}

/// Create a rope from a byte slice.
#[inline]
pub fn rope_from_bytes(arena: &mut Arena, data: &[u8]) -> Node {
    arena.from_bytes(data)
}

/// Reconstruct the byte string from a rope.
#[inline]
pub fn rope_to_bytes(arena: &Arena, node: Node) -> Vec<u8> {
    arena.to_bytes(node)
}

/// Validate all rope invariants.
#[inline]
pub fn validate_rope(arena: &Arena, node: Node) {
    arena.validate(node)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn arena() -> Arena {
        Arena::new()
    }

    #[test]
    fn test_leaf_basic() {
        let mut a = arena();
        let leaf = a.from_bytes(b"hello").unwrap();
        assert_eq!(a.node_len(leaf), 5);
        assert_eq!(a.node_weight(leaf), 1);
        assert_eq!(a.hash_val(leaf), a.h.hash(b"hello"));
    }

    #[test]
    fn test_concat_two_leaves() {
        let mut a = arena();
        let l = a.from_bytes(b"hello ");
        let r = a.from_bytes(b"world");
        let lr = a.concat(l, r);
        assert_eq!(a.len(lr), 11);
        assert_eq!(a.hash(lr), a.hash_bytes(b"hello world"));
    }

    #[test]
    fn test_concat_none() {
        let mut a = arena();
        let n = a.from_bytes(b"test");
        assert!(a.concat(None, n).is_some());
        assert!(a.concat(n, None).is_some());
        assert!(a.concat(None, None).is_none());
    }

    #[test]
    fn test_split_preserves_hash() {
        let mut a = arena();
        let data = b"the quick brown fox";
        let expected = a.hash_bytes(data);
        let node = a.from_bytes(data);
        for pos in 0..=data.len() as u64 {
            let (left, right) = a.split(node, pos);
            let rejoined = a.concat(left, right);
            assert_eq!(a.hash(rejoined), expected, "Failed at pos={}", pos);
        }
    }

    #[test]
    fn test_repeat_hash() {
        let mut a = arena();
        let node = a.from_bytes(b"ab");
        let repeated = a.repeat(node, 100);
        let expected: Vec<u8> = b"ab".iter().copied().cycle().take(200).collect();
        assert_eq!(a.hash(repeated), a.hash_bytes(&expected));
        assert_eq!(a.len(repeated), 200);
    }

    #[test]
    fn test_substr_hash() {
        let mut a = arena();
        let node = a.from_bytes(b"hello world");
        assert_eq!(a.substr_hash(node, 0, 5), a.hash_bytes(b"hello"));
        assert_eq!(a.substr_hash(node, 6, 5), a.hash_bytes(b"world"));
        assert_eq!(a.substr_hash(node, 4, 4), a.hash_bytes(b"o wo"));
    }

    #[test]
    fn test_to_bytes_roundtrip() {
        let mut a = arena();
        let data = b"round trip test";
        let node = a.from_bytes(data);
        assert_eq!(a.to_bytes(node), data);
    }

    #[test]
    fn test_many_concats_balanced() {
        let mut a = arena();
        let mut node: Node = None;
        for i in 0..100u8 {
            let leaf = a.from_bytes(&[i]);
            node = a.concat(node, leaf);
        }
        a.validate(node);
        assert_eq!(a.len(node), 100);
    }

    #[test]
    fn test_split_repeat_on_boundary() {
        let mut a = arena();
        let leaf = a.from_bytes(b"abc");
        let rep = a.repeat(leaf, 4); // "abcabcabcabc"
        let (left, right) = a.split(rep, 6);
        assert_eq!(a.to_bytes(left), b"abcabc");
        assert_eq!(a.to_bytes(right), b"abcabc");
    }

    #[test]
    fn test_split_repeat_within() {
        let mut a = arena();
        let leaf = a.from_bytes(b"abcd");
        let rep = a.repeat(leaf, 3); // "abcdabcdabcd"
        let (left, right) = a.split(rep, 5);
        assert_eq!(a.to_bytes(left), b"abcda");
        assert_eq!(a.to_bytes(right), b"bcdabcd");
    }

    #[test]
    fn test_associativity() {
        let mut a = arena();
        let x = a.from_bytes(b"aaa");
        let y = a.from_bytes(b"bbb");
        let z = a.from_bytes(b"ccc");
        let xy = a.concat(x, y);
        let xy_z = a.concat(xy, z);
        let yz = a.concat(y, z);
        let x_yz = a.concat(x, yz);
        assert_eq!(a.hash(xy_z), a.hash(x_yz));
        assert_eq!(a.hash(xy_z), a.hash_bytes(b"aaabbbccc"));
    }

    #[test]
    fn test_cross_validation_hello_world() {
        let mut a = arena();
        let node = a.from_bytes(b"hello world");
        // Known value from Python/Rust cross-validation
        assert_eq!(a.hash(node), 430229793999670395);
    }

    #[test]
    fn test_height_leaf() {
        let mut a = arena();
        assert_eq!(a.height(None), 0);
        let leaf = a.from_bytes(b"x");
        assert_eq!(a.height(leaf), 0);
    }

    #[test]
    fn test_height_grows_log() {
        let mut a = arena();
        let mut node: Node = None;
        for i in 0..1000u16 {
            let leaf = a.from_bytes(&[(i % 256) as u8]);
            node = a.concat(node, leaf);
        }
        let h = a.height(node);
        // BB[2/7] bound: h <= log(1000) / log(7/5) ≈ 20.5
        assert!(h <= 21, "height {} too large for w=1000", h);
    }
}
