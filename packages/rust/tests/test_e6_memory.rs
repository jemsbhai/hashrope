//! E6: Memory profiling — rope overhead vs flat buffer
//!
//! Claim C6: Memory overhead of the rope tree is bounded and predictable.
//!
//! Method: Walk the rope tree, count nodes and sum memory usage analytically.
//! Each node type has a known heap size (struct + Arc + Vec allocation).
//! Compare total rope memory to n (raw payload bytes).

use hashrope::{
    rope_concat, rope_from_bytes, rope_repeat, rope_len, validate_rope,
    Node, PolynomialHash,
};
use std::sync::Arc;

/// Sizes derived from Rust's memory layout.
/// Arc<NodeInner> = pointer (8 bytes) + ArcInner overhead (strong + weak counts = 16 bytes)
/// NodeInner::Leaf = enum tag(8) + Vec<u8>(24) + hash_val(8) + len(8) = 48 + Vec heap alloc
/// NodeInner::Internal = enum tag(8) + 2×Arc(16) + hash_val(8) + len(8) + weight(8) = 48
/// NodeInner::Repeat = enum tag(8) + Arc(8) + reps(8) + hash_val(8) + len(8) + weight(8) = 48
const ARC_INNER_OVERHEAD: usize = 16; // strong_count + weak_count
const NODE_INNER_SIZE: usize = 48;    // approximate, all variants similar
const VEC_HEADER_SIZE: usize = 24;    // ptr + len + capacity on stack (inside NodeInner)

struct MemStats {
    leaf_count: usize,
    internal_count: usize,
    repeat_count: usize,
    payload_bytes: usize, // actual data bytes in leaves
}

impl MemStats {
    fn new() -> Self {
        Self { leaf_count: 0, internal_count: 0, repeat_count: 0, payload_bytes: 0 }
    }

    fn walk(&mut self, node: &Node) {
        if let Some(n) = node {
            self.walk_inner(n);
        }
    }

    fn walk_inner(&mut self, node: &Arc<hashrope::rope::NodeInner>) {
        match node.as_ref() {
            hashrope::rope::NodeInner::Leaf { data, .. } => {
                self.leaf_count += 1;
                self.payload_bytes += data.len();
            }
            hashrope::rope::NodeInner::Internal { left, right, .. } => {
                self.internal_count += 1;
                self.walk_inner(left);
                self.walk_inner(right);
            }
            hashrope::rope::NodeInner::Repeat { child, .. } => {
                self.repeat_count += 1;
                self.walk_inner(child);
            }
        }
    }

    fn total_nodes(&self) -> usize {
        self.leaf_count + self.internal_count + self.repeat_count
    }

    /// Estimated heap bytes for the rope structure.
    /// Each node: ArcInner(16) + NodeInner(48) = 64 bytes
    /// Each leaf also has a Vec heap allocation of data.len() bytes
    fn estimated_heap_bytes(&self) -> usize {
        let per_node = ARC_INNER_OVERHEAD + NODE_INNER_SIZE;
        let node_bytes = self.total_nodes() * per_node;
        let data_bytes = self.payload_bytes; // Vec heap allocations
        node_bytes + data_bytes
    }
}

fn build_sequential_rope(n: usize, ph: &mut PolynomialHash) -> Node {
    let mut rope: Node = None;
    for i in 0..n {
        let leaf = rope_from_bytes(&[(i % 256) as u8], ph);
        rope = rope_concat(&rope, &leaf, ph);
    }
    rope
}

#[test]
fn e6_memory_sequential_rope() {
    let sizes: Vec<usize> = vec![100, 1_000, 10_000, 100_000, 1_000_000];

    println!("\n=== Sequential rope (one byte per leaf) ===");
    println!("{:<10} {:>8} {:>8} {:>8} {:>14} {:>10}",
        "n", "leaves", "intern", "total", "est_bytes", "overhead");
    println!("{}", "-".repeat(66));

    for &n in &sizes {
        let mut ph = PolynomialHash::default_hash();
        let rope = build_sequential_rope(n, &mut ph);
        validate_rope(&rope);

        let mut stats = MemStats::new();
        stats.walk(&rope);

        let est = stats.estimated_heap_bytes();
        let overhead = est as f64 / n as f64;

        println!("{:<10} {:>8} {:>8} {:>8} {:>14} {:>10.1}x",
            n, stats.leaf_count, stats.internal_count,
            stats.total_nodes(), est, overhead);
    }
}

#[test]
fn e6_memory_chunked_rope() {
    let sizes: Vec<usize> = vec![1_000, 10_000, 100_000, 1_000_000];
    let chunk_sizes = vec![64, 256, 1024, 4096];

    println!("\n=== Chunked rope (multi-byte leaves) ===");
    println!("{:<10} {:<8} {:>8} {:>8} {:>14} {:>10}",
        "n", "chunk", "leaves", "total", "est_bytes", "overhead");
    println!("{}", "-".repeat(66));

    for &n in &sizes {
        for &chunk in &chunk_sizes {
            if chunk > n { continue; }

            let mut ph = PolynomialHash::default_hash();
            let mut rope: Node = None;
            let mut offset = 0usize;
            while offset < n {
                let sz = chunk.min(n - offset);
                let data: Vec<u8> = (offset..offset + sz).map(|i| (i % 256) as u8).collect();
                let leaf = rope_from_bytes(&data, &ph);
                rope = rope_concat(&rope, &leaf, &mut ph);
                offset += sz;
            }

            let mut stats = MemStats::new();
            stats.walk(&rope);

            let est = stats.estimated_heap_bytes();
            let overhead = est as f64 / n as f64;

            println!("{:<10} {:<8} {:>8} {:>8} {:>14} {:>10.2}x",
                n, chunk, stats.leaf_count, stats.total_nodes(), est, overhead);
        }
    }
}

#[test]
fn e6_memory_repeat_node() {
    let q_values: Vec<u64> = vec![10, 100, 1_000, 10_000, 100_000, 1_000_000, 100_000_000];
    let base_size = 100usize;

    println!("\n=== RepeatNode memory vs materialized ===");
    println!("{:<12} {:>8} {:>14} {:>14} {:>10}",
        "q", "nodes", "rope_bytes", "material_bytes", "savings");
    println!("{}", "-".repeat(64));

    for &q in &q_values {
        let mut ph = PolynomialHash::default_hash();
        let base: Vec<u8> = (0..base_size).map(|i| (i % 256) as u8).collect();
        let base_rope = rope_from_bytes(&base, &ph);
        let repeated = rope_repeat(&base_rope, q, &mut ph);

        assert_eq!(rope_len(&repeated), q * base_size as u64);

        let mut stats = MemStats::new();
        stats.walk(&repeated);

        let est = stats.estimated_heap_bytes();
        let materialized = q as usize * base_size;
        let savings = materialized as f64 / est as f64;

        println!("{:<12} {:>8} {:>14} {:>14} {:>10.0}x",
            q, stats.total_nodes(), est, materialized, savings);
    }
}
