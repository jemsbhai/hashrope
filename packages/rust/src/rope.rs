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
        Self::is_balanced_wt(self.node_weight(left), self.node_weight(right))
    }

    /// Weight-only balance check (no node lookup needed).
    #[inline]
    fn is_balanced_wt(wl: u64, wr: u64) -> bool {
        let total = wl + wr;
        if total <= 2 {
            return true;
        }
        ALPHA_NUM * total <= ALPHA_DEN * wl && ALPHA_DEN * wl <= (ALPHA_DEN - ALPHA_NUM) * total
    }

    /// Decompose a node into two children for rotation purposes.
    /// Internal → its children. Repeat → split by reps. Leaf → panic.
    fn decompose(&mut self, id: NodeId) -> (NodeId, NodeId) {
        match self.nodes[id as usize].clone() {
            NodeInner::Internal { left, right, .. } => (left, right),
            NodeInner::Repeat { child, reps, .. } => {
                let half = reps / 2;
                let left = self.make_repeat(child, half).unwrap();
                let right = self.make_repeat(child, reps - half).unwrap();
                (left, right)
            }
            NodeInner::Leaf { .. } => unreachable!("Cannot decompose Leaf in balance"),
        }
    }

    /// Adams-style balance with checked rotations.
    /// Handles RepeatNode via decompose. Single rotation is used only when
    /// BOTH the inner node AND the outer pairing are balanced; otherwise
    /// double rotation decomposes further (bounded by subtree height).
    fn balance(&mut self, left: NodeId, right: NodeId) -> NodeId {
        let wl = self.node_weight(left);
        let wr = self.node_weight(right);
        let total = wl + wr;

        if total <= 2 {
            return self.make_internal(left, right);
        }

        if Self::is_balanced_wt(wl, wr) {
            return self.make_internal(left, right);
        }

        if ALPHA_NUM * total > ALPHA_DEN * wl {
            // Left too light (right too heavy) → rotate left
            let (rl, rr) = self.decompose(right);
            let wrl = self.node_weight(rl);
            let wrr = self.node_weight(rr);
            // Single rotation: inner=(left,rl), outer=(inner, rr)
            // Use single only if BOTH inner and outer would be balanced.
            if Self::is_balanced_wt(wl, wrl)
                && Self::is_balanced_wt(wl + wrl, wrr)
            {
                let new_left = self.make_internal(left, rl);
                self.make_internal(new_left, rr)
            } else {
                // Double rotation: decompose rl, balance all parts
                let (rll, rlr) = self.decompose(rl);
                let new_left = self.balance(left, rll);
                let new_right = self.balance(rlr, rr);
                self.balance(new_left, new_right)
            }
        } else {
            // Left too heavy (right too light) → rotate right
            let (ll, lr) = self.decompose(left);
            let wll = self.node_weight(ll);
            let wlr = self.node_weight(lr);
            // Single rotation: inner=(lr,right), outer=(ll, inner)
            if Self::is_balanced_wt(wlr, wr)
                && Self::is_balanced_wt(wll, wlr + wr)
            {
                let new_right = self.make_internal(lr, right);
                self.make_internal(ll, new_right)
            } else {
                // Double rotation: decompose lr, balance all parts
                let (lrl, lrr) = self.decompose(lr);
                let new_left = self.balance(ll, lrl);
                let new_right = self.balance(lrr, right);
                self.balance(new_left, new_right)
            }
        }
    }

    fn rebalance(&mut self, left: NodeId, right: NodeId) -> NodeId {
        self.balance(left, right)
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
                    self.balance(ll, new_right)
                }
                NodeInner::Repeat { child, reps, .. } => {
                    // Fix 2: split by reps count, not byte position.
                    // Halving reps gives O(log q) recursion depth.
                    let half = reps / 2;
                    let ll = self.make_repeat(child, half).unwrap();
                    let lr = self.make_repeat(child, reps - half).unwrap();
                    let new_right = self.join(lr, right);
                    self.join(ll, new_right)
                }
                NodeInner::Leaf { .. } => self.make_internal(left, right),
            }
        } else {
            match self.nodes[right as usize].clone() {
                NodeInner::Internal { left: rl, right: rr, .. } => {
                    let new_left = self.join(left, rl);
                    self.balance(new_left, rr)
                }
                NodeInner::Repeat { child, reps, .. } => {
                    // Fix 2: split by reps count, not byte position.
                    let half = reps / 2;
                    let rl = self.make_repeat(child, half).unwrap();
                    let rr = self.make_repeat(child, reps - half).unwrap();
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
        // Full-node shortcut: if the entire node is within the requested range,
        // return the precomputed hash in O(1). Without this, the traversal
        // descends to every leaf in the range, degrading to O(range_size).
        if start == 0 && length == self.node_len(id) {
            return self.hash_val(id);
        }

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
    fn test_substr_hash_full_node() {
        // Tests the early-exit path: substr_hash over the entire node
        // must equal the node's stored hash.
        let mut a = arena();
        let node = a.from_bytes(b"hello world");
        assert_eq!(a.substr_hash(node, 0, 11), a.hash(node));
    }

    #[test]
    fn test_substr_hash_multi_leaf_exhaustive() {
        // Build a rope from many small leaves, then verify substr_hash
        // for EVERY possible (start, length) against materialized hash.
        let data = b"abcdefghijklmnop"; // 16 bytes
        let mut a = arena();
        // Build as 4 leaves of 4 bytes each
        let l1 = a.from_bytes(b"abcd");
        let l2 = a.from_bytes(b"efgh");
        let l3 = a.from_bytes(b"ijkl");
        let l4 = a.from_bytes(b"mnop");
        let left = a.concat(l1, l2);
        let right = a.concat(l3, l4);
        let root = a.concat(left, right);

        for start in 0..16u64 {
            for length in 1..=(16 - start) {
                let expected = a.hash_bytes(&data[start as usize..(start + length) as usize]);
                let got = a.substr_hash(root, start, length);
                assert_eq!(got, expected,
                    "substr_hash mismatch at start={}, length={}", start, length);
            }
        }
        // Also test zero-length
        assert_eq!(a.substr_hash(root, 5, 0), 0);
        // Full node
        assert_eq!(a.substr_hash(root, 0, 16), a.hash(root));
    }

    #[test]
    fn test_substr_hash_single_byte_leaves() {
        // Worst-case tree: one byte per leaf. Exercises the early-exit
        // on interior subtrees that are fully within the requested range.
        let data: Vec<u8> = (0..100u8).collect();
        let mut a = arena();
        let mut node: Node = None;
        for &b in &data {
            let leaf = a.from_bytes(&[b]);
            node = a.concat(node, leaf);
        }
        // Full range
        assert_eq!(a.substr_hash(node, 0, 100), a.hash(node));
        // Various sub-ranges
        for &(s, l) in &[(0, 50), (25, 50), (50, 50), (10, 80), (0, 1), (99, 1), (0, 100)] {
            let expected = a.hash_bytes(&data[s..s + l]);
            let got = a.substr_hash(node, s as u64, l as u64);
            assert_eq!(got, expected, "mismatch at start={}, length={}", s, l);
        }
    }

    #[test]
    fn test_substr_hash_repeat_node() {
        // RepeatNode: "abcd" × 100 = 400 bytes
        let mut a = arena();
        let pattern = a.from_bytes(b"abcd");
        let rep = a.repeat(pattern, 100);
        let materialized: Vec<u8> = b"abcd".iter().copied().cycle().take(400).collect();

        // Full range
        assert_eq!(a.substr_hash(rep, 0, 400), a.hash_bytes(&materialized));
        // Within one copy
        assert_eq!(a.substr_hash(rep, 0, 3), a.hash_bytes(b"abc"));
        // Spanning two copies
        assert_eq!(a.substr_hash(rep, 2, 4), a.hash_bytes(b"cdab"));
        // Spanning many copies
        assert_eq!(a.substr_hash(rep, 1, 13), a.hash_bytes(&materialized[1..14]));
        // Last byte
        assert_eq!(a.substr_hash(rep, 399, 1), a.hash_bytes(b"d"));
    }

    #[test]
    fn test_substr_hash_after_split_rejoin() {
        // Verify substr_hash still works correctly on a rope that has been
        // split and rejoined (different tree structure, same logical string).
        let data = b"the quick brown fox jumps over the lazy dog";
        let mut a = arena();
        let original = a.from_bytes(data);
        let (left, right) = a.split(original, 20);
        let rejoined = a.concat(left, right);

        // The rejoined rope has a different tree shape but represents
        // the same string. All substr_hash queries must match.
        for start in (0..data.len()).step_by(5) {
            for length in [1, 3, 7, 10, data.len() - start] {
                let length = length.min(data.len() - start);
                if length == 0 { continue; }
                let expected = a.hash_bytes(&data[start..start + length]);
                let got = a.substr_hash(rejoined, start as u64, length as u64);
                assert_eq!(got, expected,
                    "mismatch on rejoined rope at start={}, length={}", start, length);
            }
        }
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

    // ===================================================================
    // Bug verification tests (from balance bug analysis)
    // ===================================================================

    #[test]
    fn test_bug_a_leaf_concat_repeat() {
        // Bug A: rebalance cannot handle RepeatNode.
        // concat(Leaf_w1, Repeat(child_w1, 9)) → join sees right-heavy,
        // enters RepeatNode branch, after reassembly rebalance is called
        // with a RepeatNode child. The `if let Internal` pattern fails
        // → falls through to make_internal → unbalanced node.
        let mut a = arena();
        let leaf = a.from_bytes(b"x");
        let child = a.from_bytes(b"y");
        let rep = a.repeat(child, 9);
        let result = a.concat(leaf, rep);
        a.validate(result);
    }

    #[test]
    fn test_bug_a_repeat_concat_leaf() {
        // Symmetric direction: Repeat on left, leaf on right.
        let mut a = arena();
        let child = a.from_bytes(b"y");
        let rep = a.repeat(child, 9);
        let leaf = a.from_bytes(b"x");
        let result = a.concat(rep, leaf);
        a.validate(result);
    }

    #[test]
    fn test_bug_b_rotation_unbalanced_children() {
        // Bug B: rotations create unbalanced inner nodes.
        // Build a tree purely from Internal nodes (no RepeatNode)
        // where join triggers a rotation that produces an inner node
        // with child ratio outside [2/7, 5/7].
        // Strategy: build a left-heavy spine of weight ~7, then concat
        // a right subtree of weight ~1. This forces join to walk the
        // spine and rebalance via rotation.
        let mut a = arena();
        // Build left subtree: 7 leaves concatenated sequentially
        // This creates a balanced internal tree of weight 7.
        let mut left: Node = None;
        for i in 0..7u8 {
            let leaf = a.from_bytes(&[b'a' + i]);
            left = a.concat(left, leaf);
        }
        a.validate(left); // Sanity: left is balanced

        // Now concat a single leaf on the right.
        // join(Internal_w7, Leaf_w1) → walks left spine,
        // rebalance may rotate, creating make_internal(w1, w3) = w4
        // with ratio 1/4 < 2/7 → violation.
        let right = a.from_bytes(b"z");
        let result = a.concat(left, right);
        a.validate(result);
    }

    #[test]
    fn test_bug_b_varied_weight_ratios() {
        // Try many weight ratios to find rotation-induced violations.
        // Concat a tree of weight N with a single leaf.
        let mut a = arena();
        for n in 3..=50u8 {
            let mut left: Node = None;
            for i in 0..n {
                let leaf = a.from_bytes(&[i]);
                left = a.concat(left, leaf);
            }
            let right = a.from_bytes(&[255]);
            let result = a.concat(left, right);
            a.validate(result); // Will panic if any ratio violates BB[2/7]
        }
    }

    #[test]
    fn test_stack_overflow_large_repeat_join() {
        // Test the O(reps) vs O(log reps) recursion depth claim.
        // If join's RepeatNode handling is O(reps), this will
        // stack overflow. If O(log reps), it completes fine.
        let mut a = arena();
        let child = a.from_bytes(b"a");
        let rep = a.repeat(child, 100_000);
        let leaf = a.from_bytes(b"z");
        let result = a.concat(leaf, rep);
        // If we get here without stack overflow, the O(reps) claim is wrong.
        assert_eq!(a.len(result), 100_001);
        // Also verify hash correctness
        let mut expected: Vec<u8> = Vec::new();
        expected.push(b'z');
        expected.resize(100_001, b'a');
        assert_eq!(a.hash(result), a.hash_bytes(&expected));
    }

    // ===================================================================
    // Diagnostic helpers
    // ===================================================================

    /// Returns a string describing the node type and weight.
    fn describe(a: &Arena, id: NodeId) -> alloc::string::String {
        use alloc::format;
        let w = a.node_weight(id);
        match &a.nodes[id as usize] {
            NodeInner::Leaf { .. } => format!("Leaf(w={})", w),
            NodeInner::Internal { left, right, .. } => {
                format!("Internal(w={}, L=w{}, R=w{})", w,
                    a.node_weight(*left), a.node_weight(*right))
            }
            NodeInner::Repeat { reps, .. } => format!("Repeat(w={}, reps={})", w, reps),
        }
    }

    /// Check BB[2/7] for a single Internal node. Returns Ok or Err with details.
    fn check_balance(a: &Arena, id: NodeId) -> Result<(), alloc::string::String> {
        use alloc::format;
        if let NodeInner::Internal { left, weight, .. } = &a.nodes[id as usize] {
            let total = *weight;
            let wl = a.node_weight(*left);
            if total > 2 {
                if ALPHA_NUM * total > ALPHA_DEN * wl {
                    return Err(format!("Left too light: {}/{} (min {}/{})",
                        wl, total, ALPHA_NUM, ALPHA_DEN));
                }
                if ALPHA_DEN * wl > (ALPHA_DEN - ALPHA_NUM) * total {
                    return Err(format!("Left too heavy: {}/{} (max {}/{})",
                        wl, total, ALPHA_DEN - ALPHA_NUM, ALPHA_DEN));
                }
            }
        }
        Ok(())
    }

    /// Recursively check all Internal nodes. Returns list of violations.
    fn find_all_violations(a: &Arena, id: NodeId) -> Vec<alloc::string::String> {
        use alloc::format;
        let mut violations = Vec::new();
        match &a.nodes[id as usize].clone() {
            NodeInner::Internal { left, right, .. } => {
                if let Err(msg) = check_balance(a, id) {
                    violations.push(format!("Node {}: {} -> {}", id, describe(a, id), msg));
                }
                violations.extend(find_all_violations(a, *left));
                violations.extend(find_all_violations(a, *right));
            }
            NodeInner::Repeat { child, .. } => {
                violations.extend(find_all_violations(a, *child));
            }
            _ => {}
        }
        violations
    }

    // ===================================================================
    // TERMINATION PROOF: Computational verification
    // ===================================================================

    #[test]
    fn test_termination_max_weight_decreases() {
        // Exhaustively verify that for ALL weight pairs (wl, wr) up to
        // MAX_W, and ALL valid decomposition splits, the ranking function
        // M = max(wl, wr) strictly decreases in every recursive call
        // of the double rotation path.
        //
        // This is the computational core of the termination proof.
        use alloc::format;
        let max_w: u64 = 200;

        fn valid_splits(w: u64) -> Vec<(u64, u64)> {
            let mut s = Vec::new();
            // Internal: children satisfy BB[2/7]
            for a in 1..w {
                let b = w - a;
                let t = a + b;
                if t <= 2 || (ALPHA_NUM * t <= ALPHA_DEN * a
                    && ALPHA_DEN * a <= (ALPHA_DEN - ALPHA_NUM) * t) {
                    s.push((a, b));
                }
            }
            // Repeat: reps-halving
            for c in 1..=w/2 {
                if w % c == 0 {
                    let q = w / c;
                    if q >= 2 {
                        let a = c * (q / 2);
                        let b = c * (q - q / 2);
                        if a >= 1 && b >= 1 { s.push((a, b)); }
                    }
                }
            }
            s.sort();
            s.dedup();
            s
        }

        let mut verified = 0u64;
        let mut violations: Vec<alloc::string::String> = Vec::new();
        let mut leaf_decompose_count = 0u64;
        let mut worst_ratio = 0.0f64;

        for wr in 1..=max_w {
            for wl in 1..wr {
                let total = wl + wr;
                if total <= 2 { continue; }
                // Check: is this pair unbalanced with left too light?
                let bal = ALPHA_NUM * total <= ALPHA_DEN * wl
                    && ALPHA_DEN * wl <= (ALPHA_DEN - ALPHA_NUM) * total;
                if bal { continue; }
                if !(ALPHA_NUM * total > ALPHA_DEN * wl) { continue; }

                for &(wrl, wrr) in &valid_splits(wr) {
                    // Check single rotation
                    let inner_bal = {
                        let t = wl + wrl;
                        t <= 2 || (ALPHA_NUM * t <= ALPHA_DEN * wl
                            && ALPHA_DEN * wl <= (ALPHA_DEN - ALPHA_NUM) * t)
                    };
                    let outer_bal = {
                        let t = wl + wrl + wrr;
                        let wl2 = wl + wrl;
                        t <= 2 || (ALPHA_NUM * t <= ALPHA_DEN * wl2
                            && ALPHA_DEN * wl2 <= (ALPHA_DEN - ALPHA_NUM) * t)
                    };
                    if inner_bal && outer_bal { continue; }

                    if wrl < 2 { leaf_decompose_count += 1; continue; }

                    for &(wrll, wrlr) in &valid_splits(wrl) {
                        let max1 = core::cmp::max(wl, wrll);
                        let max2 = core::cmp::max(wrlr, wrr);
                        let nl = wl + wrll;
                        let nr = wrlr + wrr;
                        let max3 = core::cmp::max(nl, nr);

                        if max1 >= wr {
                            violations.push(format!("CALL1: bal({},{}) max {} >= {}", wl, wr, max1, wr));
                        }
                        if max2 >= wr {
                            violations.push(format!("CALL2: bal({},{}) max {} >= {}", wl, wr, max2, wr));
                        }
                        if max3 >= wr {
                            violations.push(format!("CALL3: bal({},{})→bal({},{}) max {} >= {}", wl, wr, nl, nr, max3, wr));
                        }

                        let ratio = max3 as f64 / wr as f64;
                        if ratio > worst_ratio { worst_ratio = ratio; }
                        verified += 1;
                    }
                }
            }
        }

        assert!(violations.is_empty(),
            "Termination violated in {} of {} cases. First: {}",
            violations.len(), verified,
            violations.first().map(|s| s.as_str()).unwrap_or("none"));
        assert_eq!(leaf_decompose_count, 0,
            "decompose would be called on {} Leaf nodes", leaf_decompose_count);
        assert!(worst_ratio < 1.0,
            "Contraction ratio {} >= 1.0", worst_ratio);
        // Verify the analytical bound: worst case is 223/245 ≈ 0.9184
        assert!(worst_ratio < 0.92,
            "Contraction ratio {} exceeds 0.92", worst_ratio);
    }

    // ===================================================================
    // FAULT ISOLATION: Direct rebalance probes by node-type combination
    // ===================================================================

    #[test]
    fn test_diag_rebalance_internal_internal() {
        // Rebalance with both sides Internal. Tests if rotations work
        // correctly when the code path fully matches.
        use alloc::format;
        let mut failures: Vec<alloc::string::String> = Vec::new();

        for wl in 1u64..=15 {
            for wr in 1u64..=15 {
                let mut a = arena();
                let l = build_tree(&mut a, wl);
                let r = build_tree(&mut a, wr);
                let result = a.rebalance(l, r);
                let violations = find_all_violations(&a, result);
                if !violations.is_empty() {
                    failures.push(format!(
                        "rebalance(Internal_w{}, Internal_w{}): {:?}",
                        wl, wr, violations
                    ));
                }
            }
        }

        if !failures.is_empty() {
            panic!(
                "rebalance(Internal, Internal) produced {} violations:\n{}",
                failures.len(),
                failures.join("\n")
            );
        }
    }

    #[test]
    fn test_diag_rebalance_repeat_leaf() {
        // Rebalance with Repeat on left, Leaf on right.
        // This is Bug A's core path: the `if let Internal` will fail.
        use alloc::format;
        let mut failures: Vec<alloc::string::String> = Vec::new();

        for q in 2u64..=30 {
            let mut a = arena();
            let child = a.make_leaf(Vec::from(b"c" as &[u8]));
            let rep = a.make_repeat_node(child, q);
            let leaf = a.make_leaf(Vec::from(b"z" as &[u8]));

            let result = a.rebalance(rep, leaf);
            let violations = find_all_violations(&a, result);
            if !violations.is_empty() {
                failures.push(format!(
                    "rebalance(Repeat(w={}), Leaf(w=1)): {:?}",
                    q, violations
                ));
            }

            // Also test reverse direction
            let mut a2 = arena();
            let child2 = a2.make_leaf(Vec::from(b"c" as &[u8]));
            let rep2 = a2.make_repeat_node(child2, q);
            let leaf2 = a2.make_leaf(Vec::from(b"z" as &[u8]));

            let result2 = a2.rebalance(leaf2, rep2);
            let violations2 = find_all_violations(&a2, result2);
            if !violations2.is_empty() {
                failures.push(format!(
                    "rebalance(Leaf(w=1), Repeat(w={})): {:?}",
                    q, violations2
                ));
            }
        }

        if !failures.is_empty() {
            panic!(
                "rebalance(Repeat, Leaf) produced {} violations:\n{}",
                failures.len(),
                failures.join("\n")
            );
        }
    }

    #[test]
    fn test_diag_rebalance_repeat_internal() {
        // Rebalance with Repeat on one side, Internal on the other.
        use alloc::format;
        let mut failures: Vec<alloc::string::String> = Vec::new();

        for q in 2u64..=20 {
            for wi in 1u64..=20 {
                // Repeat left, Internal right
                let mut a = arena();
                let child = a.make_leaf(Vec::from(b"c" as &[u8]));
                let rep = a.make_repeat_node(child, q);
                let internal = build_tree(&mut a, wi);

                let result = a.rebalance(rep, internal);
                let violations = find_all_violations(&a, result);
                if !violations.is_empty() {
                    failures.push(format!(
                        "rebalance(Repeat(w={}), Internal(w={})): {:?}",
                        q, wi, violations
                    ));
                }

                // Internal left, Repeat right
                let mut a2 = arena();
                let internal2 = build_tree(&mut a2, wi);
                let child2 = a2.make_leaf(Vec::from(b"c" as &[u8]));
                let rep2 = a2.make_repeat_node(child2, q);

                let result2 = a2.rebalance(internal2, rep2);
                let violations2 = find_all_violations(&a2, result2);
                if !violations2.is_empty() {
                    failures.push(format!(
                        "rebalance(Internal(w={}), Repeat(w={})): {:?}",
                        wi, q, violations2
                    ));
                }
            }
        }

        if !failures.is_empty() {
            panic!(
                "rebalance(Repeat, Internal) produced {} violations:\n{}",
                failures.len(),
                failures.join("\n")
            );
        }
    }

    #[test]
    fn test_diag_rebalance_repeat_repeat() {
        // Rebalance with RepeatNode on both sides.
        use alloc::format;
        let mut failures: Vec<alloc::string::String> = Vec::new();

        for ql in 2u64..=20 {
            for qr in 2u64..=20 {
                let mut a = arena();
                let cl = a.make_leaf(Vec::from(b"L" as &[u8]));
                let repl = a.make_repeat_node(cl, ql);
                let cr = a.make_leaf(Vec::from(b"R" as &[u8]));
                let repr = a.make_repeat_node(cr, qr);

                let result = a.rebalance(repl, repr);
                let violations = find_all_violations(&a, result);
                if !violations.is_empty() {
                    failures.push(format!(
                        "rebalance(Repeat(w={}), Repeat(w={})): {:?}",
                        ql, qr, violations
                    ));
                }
            }
        }

        if !failures.is_empty() {
            panic!(
                "rebalance(Repeat, Repeat) produced {} violations:\n{}",
                failures.len(),
                failures.join("\n")
            );
        }
    }

    // ===================================================================
    // CAUSAL TRACE: Step-by-step join for minimal failure
    // ===================================================================

    #[test]
    fn test_diag_trace_join_repeat9_leaf1() {
        // Manually trace join(Repeat(child_w1, 9), Leaf_w1) step by step.
        // At each step, record the node type and weights of arguments
        // and the result. Identify the exact step that creates a violation.
        use alloc::format;
        let mut a = arena();

        let child = a.make_leaf(Vec::from(b"y" as &[u8]));  // id=0, w=1
        let rep9 = a.make_repeat_node(child, 9);             // id=1, w=9
        let leaf = a.make_leaf(Vec::from(b"x" as &[u8]));   // id=2, w=1

        // join(rep9=w9, leaf=w1): wl=9 > wr=1, left is Repeat
        // Step 1: split rep9 at byte midpoint
        let mid = a.node_len(rep9) / 2; // = 4
        let (ll_opt, lr_opt) = a.split_inner(rep9, mid);
        let ll = ll_opt.unwrap();
        let lr = lr_opt.unwrap();

        let trace1 = format!(
            "Step 1: split Repeat(w=9) at byte {} -> ll={}, lr={}",
            mid, describe(&a, ll), describe(&a, lr)
        );

        // Step 2: new_right = join(lr, leaf)
        let new_right = a.join(lr, leaf);
        let trace2 = format!(
            "Step 2: join({}, {}) -> {}",
            describe(&a, lr), describe(&a, leaf), describe(&a, new_right)
        );
        let violations_nr = find_all_violations(&a, new_right);
        let trace2v = format!("  violations in new_right: {:?}", violations_nr);

        // Step 3: join(ll, new_right)
        let final_result = a.join(ll, new_right);
        let trace3 = format!(
            "Step 3: join({}, {}) -> {}",
            describe(&a, ll), describe(&a, new_right), describe(&a, final_result)
        );
        let violations_final = find_all_violations(&a, final_result);
        let trace3v = format!("  violations in final: {:?}", violations_final);

        let full_trace = format!(
            "\n=== TRACE: join(Repeat(w=9), Leaf(w=1)) ===\n{}\n{}\n{}\n{}\n{}\n",
            trace1, trace2, trace2v, trace3, trace3v
        );

        // Print the trace, panic only if violations found
        if !violations_final.is_empty() {
            panic!("{}", full_trace);
        }
    }

    // ===================================================================
    // FAULT BOUNDARY: Bug B — find minimal split+concat that triggers
    // rotation violation without RepeatNode
    // ===================================================================

    #[test]
    fn test_diag_bug_b_minimal_reproduction() {
        // Build Internal-only trees via sequential concat (known balanced),
        // then split at every position and concat halves in reverse.
        // Find the smallest total weight where a violation occurs.
        use alloc::format;
        let mut found: Vec<alloc::string::String> = Vec::new();

        for total in 3u64..=30 {
            let mut a = arena();
            let mut rope: Node = None;
            for i in 0..total as u8 {
                let leaf = a.from_bytes(&[i]);
                rope = a.concat(rope, leaf);
            }
            // Verify the sequential build is clean
            a.validate(rope);

            for split_pos in 1..total {
                let mut a2 = arena();
                let mut rope2: Node = None;
                for i in 0..total as u8 {
                    let leaf = a2.from_bytes(&[i]);
                    rope2 = a2.concat(rope2, leaf);
                }
                let (left, right) = a2.split(rope2, split_pos);
                // Concat in reverse order
                let reversed = a2.concat(right, left);
                if let Some(rid) = reversed {
                    let violations = find_all_violations(&a2, rid);
                    if !violations.is_empty() {
                        found.push(format!(
                            "total={}, split_at={}: {:?}",
                            total, split_pos, violations
                        ));
                    }
                }
            }

            if !found.is_empty() {
                panic!(
                    "Bug B minimal reproduction found at total={}.\nFirst violations:\n{}",
                    total,
                    found.iter().take(10).cloned().collect::<Vec<_>>().join("\n")
                );
            }
        }
    }

    /// Build a balanced tree of exactly the given weight using only leaves.
    /// For w=1, returns a Leaf. For w>=2, recursively builds balanced halves.
    fn build_tree(a: &mut Arena, w: u64) -> NodeId {
        if w == 1 {
            a.make_leaf(Vec::from(&[b'a'][..]))
        } else {
            let lw = w / 2;
            let rw = w - lw;
            let l = build_tree(a, lw);
            let r = build_tree(a, rw);
            a.make_internal(l, r)
        }
    }

    // ===================================================================
    // Aggressive bug-hunting tests
    // ===================================================================

    #[test]
    fn test_bug_a_sweep_repeat_weights() {
        // Sweep: concat(Repeat(child_w1, q), Leaf_w1) for q in 2..=200
        // and the reverse direction. Validate after every concat.
        for q in 2u64..=200 {
            let mut a = arena();
            let child = a.from_bytes(b"c");
            let rep = a.repeat(child, q);
            let leaf = a.from_bytes(b"z");

            // Repeat on left
            let r1 = a.concat(rep, leaf);
            a.validate(r1);

            // Repeat on right
            let leaf2 = a.from_bytes(b"z");
            let rep2 = a.repeat(child, q);
            let r2 = a.concat(leaf2, rep2);
            a.validate(r2);
        }
    }

    #[test]
    fn test_bug_a_repeat_concat_repeat() {
        // Two RepeatNodes concatenated. Both sides trigger RepeatNode
        // branches in join.
        for ql in [2, 3, 5, 10, 50, 100].iter() {
            for qr in [2, 3, 5, 10, 50, 100].iter() {
                let mut a = arena();
                let cl = a.from_bytes(b"L");
                let cr = a.from_bytes(b"R");
                let repl = a.repeat(cl, *ql);
                let repr = a.repeat(cr, *qr);
                let result = a.concat(repl, repr);
                a.validate(result);

                // Verify content correctness
                let bytes = a.to_bytes(result);
                assert_eq!(bytes.len(), (*ql + *qr) as usize);
            }
        }
    }

    #[test]
    fn test_bug_a_multi_leaf_concat_repeat() {
        // Build an Internal tree of weight W, then concat with Repeat of weight Q.
        // This tests rebalance receiving RepeatNode children deeper in the tree.
        for w in [2, 3, 5, 7, 10, 20].iter() {
            for q in [2, 3, 5, 9, 20, 100].iter() {
                let mut a = arena();
                let mut left: Node = None;
                for i in 0..*w as u8 {
                    let leaf = a.from_bytes(&[b'a' + (i % 26)]);
                    left = a.concat(left, leaf);
                }
                let child = a.from_bytes(b"R");
                let rep = a.repeat(child, *q);
                let result = a.concat(left, rep);
                a.validate(result);

                // Also reverse direction
                let mut a2 = arena();
                let child2 = a2.from_bytes(b"L");
                let rep2 = a2.repeat(child2, *q);
                let mut right: Node = None;
                for i in 0..*w as u8 {
                    let leaf = a2.from_bytes(&[b'a' + (i % 26)]);
                    right = a2.concat(right, leaf);
                }
                let result2 = a2.concat(rep2, right);
                a2.validate(result2);
            }
        }
    }

    #[test]
    fn test_cdh_sort_lz77_pattern() {
        // Simulate cdh-sort's rope_builder LZ77 back-reference pattern:
        // 1. Build initial rope from literal bytes
        // 2. Split at some offset, take a prefix
        // 3. Repeat that prefix (the back-reference)
        // 4. Concat it back
        // This is the exact pattern that triggers the bug in production.
        let mut a = arena();

        // Step 1: initial literal "abcdefgh"
        let mut rope = a.from_bytes(b"abcdefgh");

        // Simulate several LZ77 back-references
        for &(offset, length, distance) in &[
            (0u64, 4u64, 1u64),   // copy 4 bytes starting at offset 0, distance 1
            (2, 6, 2),
            (0, 3, 3),
            (1, 10, 1),           // distance=1 is the run-length trigger
        ] {
            // Extract the source pattern
            let (_left, _) = a.split(rope, offset + distance);
            let (_, pattern) = a.split(rope, offset);

            // Take just `distance` bytes from pattern as the repeat unit
            let (unit, _) = a.split(pattern, distance);

            if let Some(unit_id) = unit {
                // Repeat to cover `length` bytes
                let reps_needed = (length + distance - 1) / distance;
                let repeated = a.repeat(Some(unit_id), reps_needed);

                // Trim to exact length
                let (exact, _) = a.split(repeated, length);

                // Concat onto the rope
                rope = a.concat(rope, exact);
                a.validate(rope);
            }
        }

        // Final validation
        a.validate(rope);
        assert!(a.len(rope) > 8, "Rope should have grown from back-references");
    }

    #[test]
    fn test_cdh_sort_distance_1_intensive() {
        // distance=1 (run-length) is the exact cdh-sort trigger mentioned
        // in the analysis. Build increasingly long run-length chains.
        let mut a = arena();
        let mut rope = a.from_bytes(b"a");

        for len in 1..=50u64 {
            let child = a.from_bytes(b"a");
            let rep = a.repeat(child, len);
            rope = a.concat(rope, rep);
            a.validate(rope);

            // Verify hash
            let expected_len = a.len(rope);
            let bytes = a.to_bytes(rope);
            assert_eq!(bytes.len() as u64, expected_len);
            assert_eq!(a.hash(rope), a.hash_bytes(&bytes));
        }
    }

    #[test]
    fn test_split_repeat_rejoin_sweep() {
        // Split a RepeatNode at every possible byte position, rejoin,
        // validate balance and hash.
        let mut a = arena();
        let child = a.from_bytes(b"abc");
        let rep = a.repeat(child, 10); // "abc" x 10 = 30 bytes
        let full_hash = a.hash(rep);

        for pos in 1..30u64 {
            let (left, right) = a.split(rep, pos);
            let rejoined = a.concat(left, right);
            a.validate(rejoined);
            assert_eq!(a.hash(rejoined), full_hash,
                "Hash mismatch after split at pos={}", pos);
        }
    }

    #[test]
    fn test_many_small_repeats_accumulated() {
        // Accumulate many small repeat nodes via concat.
        // This stresses rebalance with RepeatNode children throughout.
        let mut a = arena();
        let mut rope: Node = None;
        for i in 0..100u8 {
            let child = a.from_bytes(&[b'a' + (i % 26)]);
            let rep = a.repeat(child, ((i % 5) as u64) + 2);
            rope = a.concat(rope, rep);
            a.validate(rope);
        }
        // Verify content
        let bytes = a.to_bytes(rope);
        assert_eq!(bytes.len() as u64, a.len(rope));
        assert_eq!(a.hash(rope), a.hash_bytes(&bytes));
    }

    #[test]
    fn test_nested_repeat_of_concat() {
        // repeat(concat(A, B), q) — the child of the repeat is an
        // Internal node, not a Leaf. Tests deeper interaction.
        let mut a = arena();
        let la = a.from_bytes(b"ab");
        let lb = a.from_bytes(b"cd");
        let inner = a.concat(la, lb); // Internal("ab", "cd"), weight=2
        let rep = a.repeat(inner, 20); // weight=40
        a.validate(rep);

        // Now concat with various things
        let leaf = a.from_bytes(b"z");
        let r1 = a.concat(rep, leaf);
        a.validate(r1);

        let r2 = a.concat(leaf, rep);
        a.validate(r2);

        // Repeat of repeat-of-concat
        let child2 = a.from_bytes(b"xy");
        let rep2 = a.repeat(child2, 15);
        let r3 = a.concat(rep, rep2);
        a.validate(r3);

        // Verify hash correctness
        let bytes = a.to_bytes(r3);
        assert_eq!(a.hash(r3), a.hash_bytes(&bytes));
    }

    #[test]
    fn test_bug_b_split_then_concat_internals() {
        // Try to trigger Bug B by creating non-standard Internal tree
        // shapes via split, then concatenating them.
        // Split can produce trees with unusual weight distributions
        // that sequential concat cannot.
        let mut a = arena();
        for total in [10, 20, 50, 100].iter() {
            let mut rope: Node = None;
            for i in 0..*total as u8 {
                let leaf = a.from_bytes(&[i]);
                rope = a.concat(rope, leaf);
            }
            // Split at various points, then cross-concat the halves
            for split_at in [1, 2, 3, *total / 4, *total / 3, *total / 2, *total - 1].iter() {
                let (left, right) = a.split(rope, *split_at as u64);
                // Concat in reverse order
                let reversed = a.concat(right, left);
                a.validate(reversed);
                // Concat each half with a fresh leaf
                let extra = a.from_bytes(b"!");
                let r1 = a.concat(left, extra);
                a.validate(r1);
                let r2 = a.concat(extra, right);
                a.validate(r2);
            }
        }
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
