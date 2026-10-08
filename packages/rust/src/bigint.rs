//! Arbitrary-precision hashes and compressed ropes, enabled by `bigint`.
//!
//! This is an additive API: the crate's original `u64` API and public node
//! layout are unchanged. Hashes, lengths, repetition counts and positions use
//! [`BigUint`]. Defaults remain the 61-bit modulus and base 131; use
//! `BigUint::from(MERSENNE_127)` for the Python 127-bit profile. Bytes are hashed
//! as supplied: pass `text.as_bytes()` for UTF-8, with no normalization.
//!
//! Concatenation and splitting retain compressed repeats. Substring hashing
//! allocates no rope nodes, although arbitrary-precision arithmetic and the
//! power cache allocate memory. Like the original arena, nodes are immutable
//! and handles belong to the arena that allocated them. This module is eager;
//! the original API retains its separate lazy-hash functionality.
//!
//! ```
//! use hashrope::bigint::{Arena, BigUint, MERSENNE_127};
//!
//! let mut arena = Arena::with_hash(MERSENNE_127.into(), 131u8.into());
//! let pattern = arena.from_bytes(b"ab");
//! let copies = BigUint::from(1u8) << 100usize;
//! let rope = arena.repeat(pattern, &copies);
//! assert_eq!(arena.len(rope), &copies * 2u8);
//! // Query a byte range without materializing the enormous repeated string.
//! let start = &copies + 1u8;
//! assert_eq!(arena.substr_hash(rope, &start, &3u8.into()), arena.hash_bytes(b"bab"));
//! ```

use alloc::collections::{BTreeMap, BTreeSet};
use alloc::vec::Vec;
pub use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};

pub use crate::{MERSENNE_127, MERSENNE_61};

/// Reduce a nonnegative integer modulo a positive modulus.
///
/// In particular this supports Mersenne moduli of arbitrary width. Panics for
/// a zero modulus. It does not certify primality.
pub fn mersenne_mod(a: &BigUint, p: &BigUint) -> BigUint {
    assert!(!p.is_zero(), "Modulus must be positive");
    a % p
}

/// Compute `(a * b) mod p` without a fixed-width intermediate.
pub fn mersenne_mul(a: &BigUint, b: &BigUint, p: &BigUint) -> BigUint {
    mersenne_mod(&(a * b), p)
}

/// Compute `sum(alpha^i, i=0..q)` using inverse-free repeated doubling.
///
/// Uses O(log q) modular multiplications, including when `alpha == 1`.
pub fn phi(q: &BigUint, alpha: &BigUint, p: &BigUint) -> BigUint {
    assert!(!p.is_zero(), "Modulus must be positive");
    if q.is_zero() {
        return BigUint::zero();
    }
    let mut sum = BigUint::one() % p;
    let mut power = alpha % p;
    for bit in (0..q.bits() - 1).rev() {
        sum = (&sum * (&power + 1u8)) % p;
        power = (&power * &power) % p;
        if q.bit(bit) {
            sum = (sum * alpha + 1u8) % p;
            power = (power * alpha) % p;
        }
    }
    sum
}

/// Python-compatible polynomial hash with arbitrary-precision parameters.
///
/// `H(s) = sum((s[i] + 1) * base^(len(s)-1-i)) mod prime`.
/// The power cache has one entry per distinct requested exponent.
#[derive(Clone, Debug)]
pub struct PolynomialHash {
    p: BigUint,
    x: BigUint,
    powers: BTreeMap<BigUint, BigUint>,
}

impl PolynomialHash {
    /// Create a hasher for a Mersenne modulus and a base in `[2, prime)`.
    ///
    /// Panics if the modulus is not positive and of the form `2^k - 1`, or
    /// the base is outside its range. Primality is not checked; composite
    /// Mersenne moduli retain the same arithmetic behavior as Python.
    pub fn new(prime: BigUint, base: BigUint) -> Self {
        assert!(
            !prime.is_zero() && (&prime & (&prime + 1u8)).is_zero(),
            "Modulus must be a positive Mersenne number"
        );
        assert!(
            base >= BigUint::from(2u8) && base < prime,
            "Base must be in [2, p-1]"
        );
        Self {
            p: prime,
            x: base,
            powers: BTreeMap::new(),
        }
    }

    /// The unchanged default profile: `2^61 - 1`, base 131.
    pub fn default_hash() -> Self {
        Self::new(MERSENNE_61.into(), 131u8.into())
    }

    /// Modulus (not a primality certificate).
    pub fn prime(&self) -> &BigUint {
        &self.p
    }

    /// Polynomial base.
    pub fn base(&self) -> &BigUint {
        &self.x
    }

    /// Compute `base^n mod prime`, caching by the full exponent.
    pub fn power(&mut self, n: &BigUint) -> BigUint {
        if let Some(value) = self.powers.get(n) {
            return value.clone();
        }
        let value = self.x.modpow(n, &self.p);
        self.powers.insert(n.clone(), value.clone());
        value
    }

    /// Compute a byte hash; empty input hashes to zero.
    pub fn hash(&self, data: &[u8]) -> BigUint {
        let mut hash = BigUint::zero();
        for &byte in data {
            hash = (hash * &self.x + (u16::from(byte) + 1)) % &self.p;
        }
        hash
    }

    /// Compute `H(A || B)` from `H(A)`, the byte length of B, and `H(B)`.
    pub fn hash_concat(&mut self, h_a: &BigUint, len_b: &BigUint, h_b: &BigUint) -> BigUint {
        (h_a * self.power(len_b) + h_b) % &self.p
    }

    /// Compute `H(S^q)` without expanding the repeated bytes.
    pub fn hash_repeat(&mut self, h_s: &BigUint, d: &BigUint, q: &BigUint) -> BigUint {
        let x_d = self.power(d);
        (h_s * phi(q, &x_d, &self.p)) % &self.p
    }

    /// Compute `H(P^(length/d) || P[..length%d])` for a nonempty pattern.
    ///
    /// `h_prefix` must be the hash of the final partial pattern, or zero if
    /// there is no remainder. Panics when `d` is zero.
    pub fn hash_overlap(
        &mut self,
        h_p: &BigUint,
        d: &BigUint,
        length: &BigUint,
        h_prefix: &BigUint,
    ) -> BigUint {
        assert!(!d.is_zero(), "Pattern length must be positive");
        let q = length / d;
        let r = length % d;
        let h_repeated = self.hash_repeat(h_p, d, &q);
        self.hash_concat(&h_repeated, &r, h_prefix)
    }
}

impl Default for PolynomialHash {
    fn default() -> Self {
        Self::default_hash()
    }
}

/// Arena-local node index. The number of allocated nodes is at most `2^32`.
pub type NodeId = u32;
/// An arena-local rope handle; `None` is the empty rope.
pub type Node = Option<NodeId>;

/// Immutable arbitrary-precision node metadata.
#[derive(Clone, Debug)]
pub enum NodeInner {
    Leaf {
        data: Vec<u8>,
        hash_val: BigUint,
        len: BigUint,
    },
    Internal {
        left: NodeId,
        right: NodeId,
        hash_val: BigUint,
        len: BigUint,
        weight: BigUint,
    },
    Repeat {
        child: NodeId,
        reps: BigUint,
        hash_val: BigUint,
        len: BigUint,
        weight: BigUint,
    },
}

impl NodeInner {
    /// Stored eager hash.
    pub fn hash_val(&self) -> &BigUint {
        match self {
            Self::Leaf { hash_val, .. }
            | Self::Internal { hash_val, .. }
            | Self::Repeat { hash_val, .. } => hash_val,
        }
    }
}

/// BB[2/7] balanced, compressed rope arena with arbitrary-size metadata.
///
/// All handles are stable for the arena's lifetime. Nodes and cached powers
/// remain allocated until the arena is dropped; streaming callers that need
/// reclamation can use [`SlidingWindow::new_bounded`]. Materializing bytes is
/// naturally limited by addressable memory; structural operations are not.
pub struct Arena {
    nodes: Vec<NodeInner>,
    h: PolynomialHash,
}

impl Arena {
    /// Construct an arena using the unchanged default hash profile.
    pub fn new() -> Self {
        Self::with_hasher(PolynomialHash::default_hash())
    }

    /// Construct an arena with explicit arbitrary-precision hash parameters.
    pub fn with_hash(prime: BigUint, base: BigUint) -> Self {
        Self::with_hasher(PolynomialHash::new(prime, base))
    }

    /// Construct an arena owning the supplied hasher.
    pub fn with_hasher(h: PolynomialHash) -> Self {
        Self {
            nodes: Vec::new(),
            h,
        }
    }

    /// Hash parameters and cache.
    pub fn hasher(&self) -> &PolynomialHash {
        &self.h
    }

    /// Mutable hasher access, e.g. for compositional hash calculations.
    /// Replacing the hasher with different parameters invalidates existing
    /// node hashes; use a fresh arena to change the hash profile.
    pub fn hasher_mut(&mut self) -> &mut PolynomialHash {
        &mut self.h
    }

    /// Number of allocated nodes, including old persistent versions.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Inspect a node belonging to this arena.
    pub fn node(&self, id: NodeId) -> &NodeInner {
        &self.nodes[id as usize]
    }

    fn node_len(&self, id: NodeId) -> &BigUint {
        match self.node(id) {
            NodeInner::Leaf { len, .. }
            | NodeInner::Internal { len, .. }
            | NodeInner::Repeat { len, .. } => len,
        }
    }

    fn node_weight(&self, id: NodeId) -> BigUint {
        match self.node(id) {
            NodeInner::Leaf { .. } => BigUint::one(),
            NodeInner::Internal { weight, .. } | NodeInner::Repeat { weight, .. } => weight.clone(),
        }
    }

    fn alloc(&mut self, node: NodeInner) -> NodeId {
        let id =
            NodeId::try_from(self.nodes.len()).expect("Arena node count exceeds NodeId capacity");
        self.nodes.push(node);
        id
    }

    fn make_leaf(&mut self, data: Vec<u8>) -> NodeId {
        assert!(!data.is_empty(), "Leaf cannot be empty");
        let hash_val = self.h.hash(&data);
        let len = data.len().into();
        self.alloc(NodeInner::Leaf {
            data,
            hash_val,
            len,
        })
    }

    fn make_internal(&mut self, left: NodeId, right: NodeId) -> NodeId {
        let len = self.node_len(left) + self.node_len(right);
        let weight = self.node_weight(left) + self.node_weight(right);
        let h_left = self.node(left).hash_val().clone();
        let h_right = self.node(right).hash_val().clone();
        let len_right = self.node_len(right).clone();
        let hash_val = self.h.hash_concat(&h_left, &len_right, &h_right);
        self.alloc(NodeInner::Internal {
            left,
            right,
            hash_val,
            len,
            weight,
        })
    }

    fn make_repeat(&mut self, child: NodeId, reps: &BigUint) -> Node {
        if reps.is_zero() {
            return None;
        }
        if reps.is_one() {
            return Some(child);
        }
        let child_len = self.node_len(child).clone();
        let len = &child_len * reps;
        let weight = self.node_weight(child) * reps;
        let child_hash = self.node(child).hash_val().clone();
        let hash_val = self.h.hash_repeat(&child_hash, &child_len, reps);
        Some(self.alloc(NodeInner::Repeat {
            child,
            reps: reps.clone(),
            hash_val,
            len,
            weight,
        }))
    }

    /// Total logical byte length, including compressed repeats.
    pub fn len(&self, node: Node) -> BigUint {
        node.map_or_else(BigUint::zero, |id| self.node_len(id).clone())
    }

    /// Hash of the logical bytes; empty ropes hash to zero.
    pub fn hash(&self, node: Node) -> BigUint {
        node.map_or_else(BigUint::zero, |id| self.node(id).hash_val().clone())
    }

    /// Number of logical leaves (repeats multiply the child weight).
    pub fn weight(&self, node: Node) -> BigUint {
        node.map_or_else(BigUint::zero, |id| self.node_weight(id))
    }

    /// Structural height; leaves and the empty rope have height zero.
    /// Shared subtrees are visited only once.
    pub fn height(&self, node: Node) -> u64 {
        let Some(root) = node else {
            return 0;
        };
        let mut heights = BTreeMap::<NodeId, u64>::new();
        let mut pending = alloc::vec![(root, false)];
        while let Some((id, expanded)) = pending.pop() {
            if heights.contains_key(&id) {
                continue;
            }
            if expanded {
                let height = match self.node(id) {
                    NodeInner::Leaf { .. } => 0,
                    NodeInner::Internal { left, right, .. } => {
                        1 + core::cmp::max(heights[left], heights[right])
                    }
                    NodeInner::Repeat { child, .. } => 1 + heights[child],
                };
                heights.insert(id, height);
            } else {
                pending.push((id, true));
                match self.node(id) {
                    NodeInner::Leaf { .. } => {}
                    NodeInner::Internal { left, right, .. } => {
                        pending.push((*right, false));
                        pending.push((*left, false));
                    }
                    NodeInner::Repeat { child, .. } => pending.push((*child, false)),
                }
            }
        }
        heights[&root]
    }

    fn balanced(wl: &BigUint, wr: &BigUint) -> bool {
        let total = wl + wr;
        total <= BigUint::from(2u8) || (&total * 2u8 <= wl * 7u8 && wl * 7u8 <= &total * 5u8)
    }

    fn decompose(&mut self, id: NodeId) -> (NodeId, NodeId) {
        match self.node(id).clone() {
            NodeInner::Internal { left, right, .. } => (left, right),
            NodeInner::Repeat { child, reps, .. } => {
                let half = &reps >> 1usize;
                let left = self.make_repeat(child, &half).unwrap();
                let right = self.make_repeat(child, &(reps - half)).unwrap();
                (left, right)
            }
            NodeInner::Leaf { .. } => unreachable!("Cannot decompose a leaf"),
        }
    }

    fn join(&mut self, left: NodeId, right: NodeId) -> NodeId {
        // Join and rotation continuations live on the heap: a compressed
        // repetition count can have thousands of bits, exceeding a native
        // thread's call-stack capacity even though the graph fits in memory.
        enum Task {
            Join(NodeId, NodeId),
            Balance(NodeId, NodeId),
            WithLeft(NodeId, bool),
            WithRight(NodeId, bool),
            BalanceResults,
        }
        let mut pending = alloc::vec![Task::Join(left, right)];
        let mut results = Vec::new();
        while let Some(task) = pending.pop() {
            match task {
                Task::WithLeft(left, balance) => {
                    let right = results.pop().unwrap();
                    pending.push(if balance {
                        Task::Balance(left, right)
                    } else {
                        Task::Join(left, right)
                    });
                }
                Task::WithRight(right, balance) => {
                    let left = results.pop().unwrap();
                    pending.push(if balance {
                        Task::Balance(left, right)
                    } else {
                        Task::Join(left, right)
                    });
                }
                Task::BalanceResults => {
                    let right = results.pop().unwrap();
                    let left = results.pop().unwrap();
                    pending.push(Task::Balance(left, right));
                }
                Task::Join(left, right) => {
                    let wl = self.node_weight(left);
                    let wr = self.node_weight(right);
                    if Self::balanced(&wl, &wr) {
                        results.push(self.make_internal(left, right));
                    } else if wl > wr {
                        let balance = matches!(self.node(left), NodeInner::Internal { .. });
                        let (ll, lr) = self.decompose(left);
                        pending.push(Task::WithLeft(ll, balance));
                        pending.push(Task::Join(lr, right));
                    } else {
                        let balance = matches!(self.node(right), NodeInner::Internal { .. });
                        let (rl, rr) = self.decompose(right);
                        pending.push(Task::WithRight(rr, balance));
                        pending.push(Task::Join(left, rl));
                    }
                }
                Task::Balance(left, right) => {
                    let wl = self.node_weight(left);
                    let wr = self.node_weight(right);
                    if Self::balanced(&wl, &wr) {
                        results.push(self.make_internal(left, right));
                    } else if (&wl + &wr) * 2u8 > &wl * 7u8 {
                        let (rl, rr) = self.decompose(right);
                        let wrl = self.node_weight(rl);
                        let wrr = self.node_weight(rr);
                        if Self::balanced(&wl, &wrl) && Self::balanced(&(&wl + &wrl), &wrr) {
                            let new_left = self.make_internal(left, rl);
                            results.push(self.make_internal(new_left, rr));
                        } else {
                            let (rll, rlr) = self.decompose(rl);
                            pending.push(Task::BalanceResults);
                            pending.push(Task::Balance(rlr, rr));
                            pending.push(Task::Balance(left, rll));
                        }
                    } else {
                        let (ll, lr) = self.decompose(left);
                        let wll = self.node_weight(ll);
                        let wlr = self.node_weight(lr);
                        if Self::balanced(&wlr, &wr) && Self::balanced(&wll, &(&wlr + &wr)) {
                            let new_right = self.make_internal(lr, right);
                            results.push(self.make_internal(ll, new_right));
                        } else {
                            let (lrl, lrr) = self.decompose(lr);
                            pending.push(Task::BalanceResults);
                            pending.push(Task::Balance(lrr, right));
                            pending.push(Task::Balance(ll, lrl));
                        }
                    }
                }
            }
        }
        results.pop().unwrap()
    }

    /// Concatenate without expanding compressed repeats.
    pub fn concat(&mut self, left: Node, right: Node) -> Node {
        match (left, right) {
            (None, _) => right,
            (_, None) => left,
            (Some(left), Some(right)) => Some(self.join(left, right)),
        }
    }

    /// Split at a byte position in `[0, len]`, preserving compressed repeats.
    /// Panics for an out-of-range position.
    pub fn split(&mut self, node: Node, pos: &BigUint) -> (Node, Node) {
        assert!(*pos <= self.len(node), "Split position exceeds rope length");
        if pos.is_zero() {
            return (None, node);
        }
        if *pos == self.len(node) {
            return (node, None);
        }
        self.split_inner(
            node.expect("Nonzero position requires a nonempty rope"),
            pos,
        )
    }

    fn split_inner(&mut self, id: NodeId, pos: &BigUint) -> (Node, Node) {
        let mut id = id;
        let mut pos = pos.clone();
        let mut ancestors = Vec::new();
        let (mut left_part, mut right_part) = loop {
            if pos.is_zero() {
                break (None, Some(id));
            }
            if pos == *self.node_len(id) {
                break (Some(id), None);
            }
            match self.node(id).clone() {
                NodeInner::Leaf { data, .. } => {
                    let pos = pos
                        .to_usize()
                        .expect("Leaf position exceeds addressable memory");
                    break (
                        Some(self.make_leaf(data[..pos].to_vec())),
                        Some(self.make_leaf(data[pos..].to_vec())),
                    );
                }
                NodeInner::Internal { left, right, .. } => {
                    let left_len = self.node_len(left).clone();
                    if pos == left_len {
                        break (Some(left), Some(right));
                    }
                    if pos < left_len {
                        ancestors.push((None, Some(right)));
                        id = left;
                    } else {
                        ancestors.push((Some(left), None));
                        id = right;
                        pos -= left_len;
                    }
                }
                NodeInner::Repeat { child, reps, .. } => {
                    let d = self.node_len(child);
                    let q = &pos / d;
                    let r = &pos % d;
                    if r.is_zero() {
                        break (
                            self.make_repeat(child, &q),
                            self.make_repeat(child, &(reps - q)),
                        );
                    }
                    let before = self.make_repeat(child, &q);
                    let after = self.make_repeat(child, &(reps - q - 1u8));
                    ancestors.push((before, after));
                    id = child;
                    pos = r;
                }
            }
        };
        while let Some((before, after)) = ancestors.pop() {
            left_part = self.concat(before, left_part);
            right_part = self.concat(right_part, after);
        }
        (left_part, right_part)
    }

    /// Represent `node` repeated `q` times with arbitrary-size metadata.
    pub fn repeat(&mut self, node: Node, q: &BigUint) -> Node {
        node.and_then(|id| self.make_repeat(id, q))
    }

    /// Hash a valid byte range without allocating rope nodes or expanding
    /// repetitions. Panics unless `start + length <= len`.
    pub fn substr_hash(&mut self, node: Node, start: &BigUint, length: &BigUint) -> BigUint {
        assert!(
            start + length <= self.len(node),
            "Substring range exceeds rope length"
        );
        if length.is_zero() {
            return BigUint::zero();
        }
        self.hash_range(
            node.expect("Nonempty substring requires a rope"),
            start,
            length,
        )
    }

    fn hash_range(&mut self, id: NodeId, start: &BigUint, length: &BigUint) -> BigUint {
        if length.is_zero() {
            return BigUint::zero();
        }
        if start.is_zero() && length == self.node_len(id) {
            return self.node(id).hash_val().clone();
        }
        enum Part {
            Range(NodeId, BigUint, BigUint),
            Hashed(BigUint, BigUint),
        }
        let mut pending = alloc::vec![Part::Range(id, start.clone(), length.clone())];
        let mut result = BigUint::zero();
        // Visit the selected segments left-to-right. Whole subtrees and full
        // repetition runs contribute one hash regardless of logical size.
        while let Some(part) = pending.pop() {
            match part {
                Part::Hashed(hash, length) => {
                    result = self.h.hash_concat(&result, &length, &hash);
                }
                Part::Range(id, start, length) => {
                    if length.is_zero() {
                        continue;
                    }
                    if start.is_zero() && length == *self.node_len(id) {
                        let hash = self.node(id).hash_val().clone();
                        result = self.h.hash_concat(&result, &length, &hash);
                        continue;
                    }
                    match self.node(id).clone() {
                        NodeInner::Leaf { data, .. } => {
                            let start = start
                                .to_usize()
                                .expect("Leaf position exceeds addressable memory");
                            let len = length
                                .to_usize()
                                .expect("Leaf length exceeds addressable memory");
                            let hash = self.h.hash(&data[start..start + len]);
                            result = self.h.hash_concat(&result, &length, &hash);
                        }
                        NodeInner::Internal { left, right, .. } => {
                            let left_len = self.node_len(left).clone();
                            if &start + &length <= left_len {
                                pending.push(Part::Range(left, start, length));
                            } else if start >= left_len {
                                pending.push(Part::Range(right, start - left_len, length));
                            } else {
                                let left_part = left_len - &start;
                                let right_part = length - &left_part;
                                pending.push(Part::Range(right, BigUint::zero(), right_part));
                                pending.push(Part::Range(left, start, left_part));
                            }
                        }
                        NodeInner::Repeat { child, .. } => {
                            let d = self.node_len(child).clone();
                            let offset = start % &d;
                            let tail_len = core::cmp::min(&d - &offset, length.clone());
                            let remaining = length - &tail_len;
                            let copies = &remaining / &d;
                            let head_len = &remaining % &d;
                            if !head_len.is_zero() {
                                pending.push(Part::Range(child, BigUint::zero(), head_len));
                            }
                            if !copies.is_zero() {
                                let child_hash = self.node(child).hash_val().clone();
                                let full_hash = self.h.hash_repeat(&child_hash, &d, &copies);
                                pending.push(Part::Hashed(full_hash, copies * &d));
                            }
                            pending.push(Part::Range(child, offset, tail_len));
                        }
                    }
                }
            }
        }
        result
    }

    /// Build a balanced rope from bytes, using leaves of at most 512 bytes.
    pub fn from_bytes(&mut self, data: &[u8]) -> Node {
        self.from_bytes_chunked(data, 512)
    }

    /// Build from bytes with a caller-selected positive leaf size.
    pub fn from_bytes_chunked(&mut self, data: &[u8], chunk_size: usize) -> Node {
        assert!(chunk_size > 0, "Chunk size must be positive");
        let mut root = None;
        for chunk in data.chunks(chunk_size) {
            let leaf = Some(self.make_leaf(chunk.to_vec()));
            root = self.concat(root, leaf);
        }
        root
    }

    /// Materialize logical bytes. Panics before traversal when the logical
    /// length exceeds the platform's `Vec` address limit. Allocation can still
    /// fail for an otherwise addressable length; use substring hashing for
    /// enormous compressed ropes.
    pub fn to_bytes(&self, node: Node) -> Vec<u8> {
        let len = self
            .len(node)
            .to_usize()
            .expect("Rope byte length exceeds usize::MAX");
        assert!(
            len <= isize::MAX as usize,
            "Rope byte length exceeds isize::MAX"
        );
        let mut bytes = Vec::with_capacity(len);
        if let Some(id) = node {
            self.collect_bytes(id, &mut bytes);
        }
        bytes
    }

    fn collect_bytes(&self, id: NodeId, bytes: &mut Vec<u8>) {
        let mut pending = alloc::vec![(id, 1usize)];
        while let Some((id, copies)) = pending.pop() {
            if copies > 1 {
                pending.push((id, copies - 1));
            }
            match self.node(id) {
                NodeInner::Leaf { data, .. } => bytes.extend_from_slice(data),
                NodeInner::Internal { left, right, .. } => {
                    pending.push((*right, 1));
                    pending.push((*left, 1));
                }
                NodeInner::Repeat { child, reps, .. } => {
                    // to_bytes checks the full logical size before traversal.
                    pending.push((
                        *child,
                        reps.to_usize()
                            .expect("Repeat count exceeds addressable memory"),
                    ));
                }
            }
        }
    }

    /// Validate lengths, weights, hashes, and BB[2/7] balance. Panics on any
    /// invariant violation. Shared subtrees are visited once.
    pub fn validate(&self, node: Node) {
        let mut seen = BTreeSet::new();
        let mut pending: Vec<NodeId> = node.into_iter().collect();
        while let Some(id) = pending.pop() {
            if !seen.insert(id) {
                continue;
            }
            match self.node(id) {
                NodeInner::Leaf {
                    data,
                    len,
                    hash_val,
                } => {
                    assert!(!data.is_empty(), "Empty leaf");
                    assert_eq!(*len, BigUint::from(data.len()), "Leaf length mismatch");
                    assert_eq!(*hash_val, self.h.hash(data), "Leaf hash mismatch");
                }
                NodeInner::Internal {
                    left,
                    right,
                    hash_val,
                    len,
                    weight,
                } => {
                    assert_eq!(
                        *len,
                        self.node_len(*left) + self.node_len(*right),
                        "Internal length mismatch"
                    );
                    let wl = self.node_weight(*left);
                    let wr = self.node_weight(*right);
                    assert_eq!(*weight, &wl + &wr, "Internal weight mismatch");
                    assert!(Self::balanced(&wl, &wr), "BB[2/7] balance violation");
                    let power = self.h.x.modpow(self.node_len(*right), &self.h.p);
                    assert_eq!(
                        *hash_val,
                        (self.node(*left).hash_val() * power + self.node(*right).hash_val())
                            % &self.h.p,
                        "Internal hash mismatch"
                    );
                    pending.push(*right);
                    pending.push(*left);
                }
                NodeInner::Repeat {
                    child,
                    reps,
                    hash_val,
                    len,
                    weight,
                } => {
                    assert!(
                        *reps >= BigUint::from(2u8),
                        "Repeat count must be at least 2"
                    );
                    assert_eq!(*len, self.node_len(*child) * reps, "Repeat length mismatch");
                    assert_eq!(
                        *weight,
                        self.node_weight(*child) * reps,
                        "Repeat weight mismatch"
                    );
                    let power = self.h.x.modpow(self.node_len(*child), &self.h.p);
                    assert_eq!(
                        *hash_val,
                        (self.node(*child).hash_val() * phi(reps, &power, &self.h.p)) % &self.h.p,
                        "Repeat hash mismatch"
                    );
                    pending.push(*child);
                }
            }
        }
    }

    /// Hash raw bytes using this arena's parameters.
    pub fn hash_bytes(&self, data: &[u8]) -> BigUint {
        self.h.hash(data)
    }
}

impl Default for Arena {
    fn default() -> Self {
        Self::new()
    }
}

/// Total logical byte length.
pub fn rope_len(arena: &Arena, node: Node) -> BigUint {
    arena.len(node)
}
/// Hash of logical bytes.
pub fn rope_hash(arena: &Arena, node: Node) -> BigUint {
    arena.hash(node)
}
/// Structural height.
pub fn rope_height(arena: &Arena, node: Node) -> u64 {
    arena.height(node)
}
/// Concatenate ropes.
pub fn rope_concat(arena: &mut Arena, left: Node, right: Node) -> Node {
    arena.concat(left, right)
}
/// Split at a valid byte position.
pub fn rope_split(arena: &mut Arena, node: Node, pos: &BigUint) -> (Node, Node) {
    arena.split(node, pos)
}
/// Compressed repetition.
pub fn rope_repeat(arena: &mut Arena, node: Node, q: &BigUint) -> Node {
    arena.repeat(node, q)
}
/// Hash a valid range without allocating rope nodes.
pub fn rope_substr_hash(
    arena: &mut Arena,
    node: Node,
    start: &BigUint,
    length: &BigUint,
) -> BigUint {
    arena.substr_hash(node, start, length)
}
/// Construct from bytes.
pub fn rope_from_bytes(arena: &mut Arena, data: &[u8]) -> Node {
    arena.from_bytes(data)
}
/// Materialize logical bytes when their length is addressable.
pub fn rope_to_bytes(arena: &Arena, node: Node) -> Vec<u8> {
    arena.to_bytes(node)
}
/// Validate all structural and hash invariants.
pub fn validate_rope(arena: &Arena, node: Node) {
    arena.validate(node)
}

/// Sliding byte window with arbitrary-size hash and stream position.
///
/// [`Self::new`] retains stable arena handles. [`Self::new_bounded`] reclaims
/// unreachable nodes after eviction, preserving the compressed retained graph.
pub struct SlidingWindow {
    arena: Arena,
    d_max: BigUint,
    /// Maximum window slack used to determine the eviction threshold.
    pub m_max: BigUint,
    w: BigUint,
    h_prefix: BigUint,
    l_prefix: BigUint,
    r_window: Node,
    pos: BigUint,
    compact_on_evict: bool,
}

impl SlidingWindow {
    /// Construct a stable-handle window; `d_max` and `m_max` must be positive.
    pub fn new(d_max: BigUint, m_max: BigUint, prime: BigUint, base: BigUint) -> Self {
        assert!(
            !d_max.is_zero() && !m_max.is_zero(),
            "Window sizes must be positive"
        );
        let w = &d_max + &m_max;
        Self {
            arena: Arena::with_hash(prime, base),
            d_max,
            m_max,
            w,
            h_prefix: BigUint::zero(),
            l_prefix: BigUint::zero(),
            r_window: None,
            pos: BigUint::zero(),
            compact_on_evict: false,
        }
    }

    /// Reclaim unreachable nodes and cached powers after each eviction.
    ///
    /// Handles observed through `arena()` may change on append. Calling
    /// `arena_mut()` permanently disables compaction to preserve caller-owned
    /// handles. Compaction costs O(retained graph size plus leaf data) and
    /// never expands repeats.
    pub fn new_bounded(d_max: BigUint, m_max: BigUint, prime: BigUint, base: BigUint) -> Self {
        let mut window = Self::new(d_max, m_max, prime, base);
        window.compact_on_evict = true;
        window
    }

    /// Defaults: distance 32768, slack 258, modulus `2^61-1`, base 131.
    pub fn default_window() -> Self {
        Self::new(
            32768u32.into(),
            258u16.into(),
            MERSENNE_61.into(),
            131u8.into(),
        )
    }

    /// Default parameters with eviction-time reclamation.
    pub fn default_bounded_window() -> Self {
        Self::new_bounded(
            32768u32.into(),
            258u16.into(),
            MERSENNE_61.into(),
            131u8.into(),
        )
    }

    /// Inspect the backing arena; bounded-mode handles may change on append.
    pub fn arena(&self) -> &Arena {
        &self.arena
    }

    /// Access the backing arena, permanently disabling compaction.
    pub fn arena_mut(&mut self) -> &mut Arena {
        self.compact_on_evict = false;
        &mut self.arena
    }

    /// Current retained logical byte length.
    pub fn window_len(&self) -> BigUint {
        self.arena.len(self.r_window)
    }

    /// Total bytes ingested, including bytes evicted from the window.
    pub fn pos(&self) -> BigUint {
        self.pos.clone()
    }

    /// Hash of the entire stream so far.
    pub fn current_hash(&mut self) -> BigUint {
        let len = self.window_len();
        let hash = self.arena.hash(self.r_window);
        self.arena.h.hash_concat(&self.h_prefix, &len, &hash)
    }

    /// Alias for [`Self::current_hash`].
    pub fn final_hash(&mut self) -> BigUint {
        self.current_hash()
    }

    /// Append literal bytes and evict excess window contents.
    pub fn append_bytes(&mut self, data: &[u8]) {
        if data.is_empty() {
            return;
        }
        let leaf = self.arena.from_bytes(data);
        self.r_window = self.arena.concat(self.r_window, leaf);
        self.pos += data.len();
        self.evict();
    }

    /// Copy `length` bytes starting `offset` bytes back, allowing overlap.
    ///
    /// Panics if offset exceeds retained or decoded length, or if a nonempty
    /// copy has offset zero. A zero-length copy does not allocate rope nodes.
    pub fn append_copy(&mut self, offset: &BigUint, length: &BigUint) {
        let win_len = self.window_len();
        assert!(
            *offset <= self.pos && *offset <= win_len,
            "Copy offset exceeds available window"
        );
        if length.is_zero() {
            return;
        }
        assert!(
            !offset.is_zero(),
            "A nonempty copy requires a positive offset"
        );
        let (_, source) = self.arena.split(self.r_window, &(&win_len - offset));
        let copy = if offset >= length {
            self.arena.split(source, length).0
        } else {
            let q = length / offset;
            let r = length % offset;
            let repeat = self.arena.repeat(source, &q);
            let partial = self.arena.split(source, &r).0;
            self.arena.concat(repeat, partial)
        };
        self.r_window = self.arena.concat(self.r_window, copy);
        self.pos += length;
        self.evict();
    }

    fn evict(&mut self) {
        let len = self.window_len();
        if len <= self.w {
            return;
        }
        let excess = len - &self.d_max;
        let (old, keep) = self.arena.split(self.r_window, &excess);
        let hash = self.arena.hash(old);
        self.h_prefix = self.arena.h.hash_concat(&self.h_prefix, &excess, &hash);
        self.l_prefix += excess;
        self.r_window = keep;
        if self.compact_on_evict {
            self.compact();
        }
    }

    fn compact(&mut self) {
        let mut arena = Arena::with_hash(self.arena.h.p.clone(), self.arena.h.x.clone());
        let mut copied = BTreeMap::<NodeId, NodeId>::new();
        let mut pending: Vec<(NodeId, bool)> =
            self.r_window.into_iter().map(|id| (id, false)).collect();
        while let Some((id, expanded)) = pending.pop() {
            if copied.contains_key(&id) {
                continue;
            }
            if expanded {
                let mut node = self.arena.node(id).clone();
                match &mut node {
                    NodeInner::Leaf { .. } => {}
                    NodeInner::Internal { left, right, .. } => {
                        *left = copied[left];
                        *right = copied[right];
                    }
                    NodeInner::Repeat { child, .. } => {
                        *child = copied[child];
                    }
                }
                copied.insert(id, arena.alloc(node));
            } else {
                pending.push((id, true));
                match self.arena.node(id) {
                    NodeInner::Leaf { .. } => {}
                    NodeInner::Internal { left, right, .. } => {
                        pending.push((*right, false));
                        pending.push((*left, false));
                    }
                    NodeInner::Repeat { child, .. } => pending.push((*child, false)),
                }
            }
        }
        let root = self.r_window.map(|id| copied[&id]);
        self.arena = arena;
        self.r_window = root;
    }
}

impl Default for SlidingWindow {
    fn default() -> Self {
        Self::default_window()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn n(value: u64) -> BigUint {
        value.into()
    }
    fn decimal(value: &str) -> BigUint {
        BigUint::parse_bytes(value.as_bytes(), 10).unwrap()
    }
    fn wide_arena() -> Arena {
        Arena::with_hash(MERSENNE_127.into(), n(131))
    }

    #[test]
    fn python_127_and_521_bit_vectors() {
        // Values generated by Python hashrope 0.2.2's PolynomialHash.
        let bytes: Vec<u8> = (0..=255).collect();
        let h = PolynomialHash::new(MERSENNE_127.into(), n(131));
        assert_eq!(h.hash(b"hello world"), decimal("157448002584923450032577"));
        assert_eq!(
            h.hash(&bytes),
            decimal("7653535697974194259553846731380633687")
        );
        let p = (BigUint::one() << 521usize) - 1u8;
        let mut h = PolynomialHash::new(p, n(131));
        assert_eq!(h.hash(&bytes), decimal("316186461957924461420614969997580723271260356953460080109619744107785345125878648903495405730214033287464229196557173316066057309411852489971865261098995235"));
        let exponent = (BigUint::one() << 80usize) + 3u8;
        assert_eq!(h.power(&exponent), decimal("2927773399229949557588963459749824909320332130862810652948295505970238977053625730812354219135310649755846964280756659433838798605367230816207040720198003025"));
    }

    #[test]
    fn default_profile_preserves_legacy_bytes() {
        let data: Vec<u8> = (0..=255).cycle().take(4097).collect();
        let old = crate::PolynomialHash::default_hash();
        let wide = PolynomialHash::default_hash();
        assert_eq!(wide.hash(&data), BigUint::from(old.hash(&data)));
        assert_eq!(wide.hash(b"hello world"), n(430229793999670395));
        assert_eq!(wide.hash(&[]), BigUint::zero());
        assert_eq!(wide.hash(&[0]), BigUint::one());
    }

    #[test]
    fn modular_arithmetic_and_inverse_free_phi() {
        let huge = (BigUint::one() << 1024usize) + 256u16;
        assert_eq!(mersenne_mod(&huge, &n(3)), n(2));
        assert_eq!(mersenne_mod(&n(256), &n(3)), n(1));
        assert_eq!(mersenne_mul(&huge, &huge, &n(3)), n(1));
        for modulus in [3, 7, 31, 127] {
            let p = n(modulus);
            for alpha in [0, 1, 2, 131] {
                let alpha = n(alpha);
                let mut sum = BigUint::zero();
                let mut power = BigUint::one();
                for q in 0..50 {
                    assert_eq!(phi(&n(q), &alpha, &p), sum);
                    sum = (sum + &power) % &p;
                    power = (power * &alpha) % &p;
                }
            }
        }
        let q = (BigUint::one() << 130usize) + 7u8;
        assert_eq!(phi(&q, &BigUint::one(), &n(127)), q % 127u8);
    }

    #[test]
    fn unicode_remains_utf8_bytes_and_byte_offsets() {
        let data = "café 🦀 e\u{301}".as_bytes();
        let mut a = wide_arena();
        let node = a.from_bytes_chunked(data, 2);
        assert_eq!(a.len(node), BigUint::from(data.len()));
        assert_eq!(a.to_bytes(node), data);
        assert_eq!(a.hash(node), a.hash_bytes(data));
        // Every valid byte position, including the middle of UTF-8 sequences.
        for start in 0..=data.len() {
            for length in 0..=data.len() - start {
                assert_eq!(
                    a.substr_hash(node, &start.into(), &length.into()),
                    a.hash_bytes(&data[start..start + length])
                );
            }
        }
        a.validate(node);
    }

    #[test]
    fn exhaustive_repeat_split_substrings_and_rejoin() {
        let mut a = wide_arena();
        let seed = a.from_bytes_chunked(b"abcde", 2);
        let repeated = a.repeat(seed, &n(9));
        let prefix = a.from_bytes(b"prefix:");
        let suffix = a.from_bytes(b":suffix");
        let node = a.concat(prefix, repeated);
        let node = a.concat(node, suffix);
        let bytes = a.to_bytes(node);
        let node_count = a.node_count();
        for start in 0..=bytes.len() {
            for length in 0..=bytes.len() - start {
                assert_eq!(
                    a.substr_hash(node, &start.into(), &length.into()),
                    a.hash_bytes(&bytes[start..start + length])
                );
            }
        }
        assert_eq!(
            a.node_count(),
            node_count,
            "Substring hashing allocated nodes"
        );
        for pos in 0..=bytes.len() {
            let (left, right) = a.split(node, &pos.into());
            assert_eq!(a.to_bytes(left), bytes[..pos]);
            assert_eq!(a.to_bytes(right), bytes[pos..]);
            a.validate(left);
            a.validate(right);
            let joined = a.concat(left, right);
            assert_eq!(a.hash(joined), a.hash(node));
            a.validate(joined);
        }
    }

    #[test]
    fn repetitions_positions_and_lengths_above_u64_stay_compressed() {
        let mut a = wide_arena();
        let seed = a.from_bytes(b"abcd");
        let q = (BigUint::one() << 130usize) + 7u8;
        let repeated = a.repeat(seed, &q);
        assert_eq!(a.len(repeated), &q * 4u8);
        assert_eq!(a.weight(repeated), q);
        assert_eq!(a.node_count(), 2);
        assert_eq!(
            a.hash(repeated),
            decimal("94497087302789760135495345347933475721")
        );
        let start = (BigUint::one() << 130usize) + 1u8;
        let allocated = a.node_count();
        assert_eq!(
            a.substr_hash(repeated, &start, &n(11)),
            a.hash_bytes(b"bcdabcdabcd")
        );
        assert_eq!(a.node_count(), allocated);
        let (left, right) = a.split(repeated, &start);
        assert_eq!(a.len(left), start);
        assert_eq!(a.len(right), a.len(repeated) - a.len(left));
        a.validate(left);
        a.validate(right);
        let joined = a.concat(left, right);
        assert_eq!(a.hash(joined), a.hash(repeated));
        a.validate(joined);
        assert!(
            a.node_count() < 10000,
            "A compressed operation expanded repeats"
        );
    }

    #[test]
    fn thousand_bit_repeat_traversals_use_heap_frames() {
        // Previously validate() and partial substr_hash() exhausted the
        // default Windows stack after this split. The resulting graph is
        // small in memory but over 1000 levels deep.
        let mut a = Arena::with_hash(n(127), n(31));
        let pattern = a.from_bytes_chunked(b"abcde", 1);
        let count = (BigUint::one() << 1024usize) + 17u8;
        let repeated = a.repeat(pattern, &count);
        let (left, right) = a.split(repeated, &n(1));
        assert_eq!(a.to_bytes(left), b"a");
        assert!(a.height(right) > 1000);
        a.validate(right);
        let allocated = a.node_count();
        assert_eq!(a.substr_hash(right, &n(0), &n(4)), a.hash_bytes(b"bcde"));
        assert_eq!(a.node_count(), allocated);
        let (first, remainder) = a.split(right, &n(4));
        assert_eq!(a.to_bytes(first), b"bcde");
        let rejoined = a.concat(first, remainder);
        let rejoined = a.concat(left, rejoined);
        assert_eq!(a.hash(rejoined), a.hash(repeated));
        a.validate(rejoined);

        // Copy exactly the reachable deep graph without traversing the call
        // stack. Graph copying must preserve sharing and node metadata.
        let mut window = SlidingWindow::new_bounded(n(1), n(1), n(127), n(31));
        window.arena = a;
        window.r_window = right;
        let expected = window.arena.hash(right);
        window.compact();
        assert_eq!(window.arena.hash(window.r_window), expected);
        assert!(window.arena.height(window.r_window) > 1000);
        window.arena.validate(window.r_window);
    }

    #[test]
    fn extreme_weight_ratios_and_nested_repeats_remain_balanced() {
        let mut a = wide_arena();
        let left = a.from_bytes(b"L");
        let pattern = a.from_bytes_chunked(b"abc", 1);
        let q = (BigUint::one() << 80usize) + 1u8;
        let repeat = a.repeat(pattern, &q);
        let right = a.from_bytes(b"R");
        let prefixed = a.concat(left, repeat);
        let both = a.concat(prefixed, right);
        a.validate(both);
        let nested = a.repeat(both, &q);
        let edge = a.len(both) - 1u8;
        assert_eq!(a.substr_hash(nested, &edge, &n(5)), a.hash_bytes(b"RLabc"));
        a.validate(nested);
        assert!(a.height(both) < 200);
        assert!(a.node_count() < 10000);
    }

    #[test]
    fn repeated_concats_and_persistent_versions() {
        let mut a = wide_arena();
        let mut root = None;
        let mut bytes = Vec::new();
        for i in 0..300u16 {
            let byte = (i % 256) as u8;
            let leaf = a.from_bytes(&[byte]);
            let old = root;
            root = a.concat(root, leaf);
            assert_eq!(a.hash(old), a.hash_bytes(&bytes));
            bytes.push(byte);
            a.validate(root);
            assert_eq!(a.hash(root), a.hash_bytes(&bytes));
        }
        assert!(a.height(root) < 20);
    }

    #[test]
    fn overlap_matches_literal_bytes() {
        let mut h = PolynomialHash::new(MERSENNE_127.into(), n(131));
        let pattern = b"abc";
        let pattern_hash = h.hash(pattern);
        for length in 0..100usize {
            let bytes: Vec<u8> = pattern.iter().copied().cycle().take(length).collect();
            let prefix_hash = h.hash(&pattern[..length % pattern.len()]);
            assert_eq!(
                h.hash_overlap(&pattern_hash, &n(3), &length.into(), &prefix_hash),
                h.hash(&bytes)
            );
        }
    }

    #[test]
    fn sliding_matches_rebuilt_stream_across_evictions() {
        for bounded in [false, true] {
            let mut window = if bounded {
                SlidingWindow::new_bounded(n(8), n(3), MERSENNE_127.into(), n(131))
            } else {
                SlidingWindow::new(n(8), n(3), MERSENNE_127.into(), n(131))
            };
            let mut bytes = b"abcde".to_vec();
            window.append_bytes(&bytes);
            for round in 0..100usize {
                let offset = 1 + round % 5;
                let length = round % 13;
                let start = bytes.len() - offset;
                for index in 0..length {
                    bytes.push(bytes[start + index]);
                }
                window.append_copy(&offset.into(), &length.into());
                assert_eq!(window.pos(), BigUint::from(bytes.len()));
                assert_eq!(window.current_hash(), window.arena.hash_bytes(&bytes));
                assert!(window.window_len() <= n(11));
                window.arena.validate(window.r_window);
                let literal = [round as u8];
                window.append_bytes(&literal);
                bytes.extend_from_slice(&literal);
                assert_eq!(window.final_hash(), window.arena.hash_bytes(&bytes));
                if bounded {
                    assert!(window.arena.node_count() < 100);
                }
            }
        }
    }

    #[test]
    fn sliding_huge_copy_is_compressed_and_reclaimed() {
        let mut window = SlidingWindow::new_bounded(n(8), n(3), MERSENNE_127.into(), n(131));
        window.append_bytes(b"ab");
        let q = (BigUint::one() << 80usize) + 3u8;
        let length = &q * 2u8;
        window.append_copy(&n(2), &length);
        let mut h = PolynomialHash::new(MERSENNE_127.into(), n(131));
        let expected = h.hash_repeat(&h.hash(b"ab"), &n(2), &(&q + 1u8));
        assert_eq!(window.current_hash(), expected);
        assert_eq!(window.pos(), length + 2u8);
        assert_eq!(window.window_len(), n(8));
        assert_eq!(window.arena.to_bytes(window.r_window), b"abababab");
        assert!(window.arena.node_count() < 20);
        window.append_bytes(b"!");
        assert_eq!(
            window.current_hash(),
            h.hash_concat(&expected, &n(1), &h.hash(b"!"))
        );
    }

    #[test]
    fn mutable_arena_disables_reclamation_to_preserve_handles() {
        let mut window = SlidingWindow::new_bounded(n(3), n(2), MERSENNE_127.into(), n(131));
        let saved = window.arena_mut().from_bytes(b"saved");
        for _ in 0..20 {
            window.append_bytes(b"stream");
        }
        assert_eq!(window.arena().to_bytes(saved), b"saved");
        assert!(!window.compact_on_evict);
    }

    #[test]
    #[should_panic(expected = "Substring range exceeds rope length")]
    fn invalid_substring_does_not_synthesize_repeat_bytes() {
        let mut a = wide_arena();
        let seed = a.from_bytes(b"ab");
        let node = a.repeat(seed, &n(3));
        a.substr_hash(node, &n(5), &n(2));
    }

    #[test]
    #[should_panic(expected = "Split position exceeds rope length")]
    fn invalid_split_is_rejected() {
        let mut a = wide_arena();
        let node = a.from_bytes(b"ab");
        a.split(node, &n(3));
    }

    #[test]
    #[should_panic(expected = "Rope byte length exceeds usize::MAX")]
    fn impossible_materialization_fails_before_expansion() {
        let mut a = wide_arena();
        let seed = a.from_bytes(b"a");
        let node = a.repeat(seed, &(BigUint::one() << 100usize));
        a.to_bytes(node);
    }

    #[test]
    #[should_panic(expected = "positive Mersenne")]
    fn non_mersenne_profile_is_rejected_in_new_api() {
        PolynomialHash::new(n(257), n(131));
    }
}
