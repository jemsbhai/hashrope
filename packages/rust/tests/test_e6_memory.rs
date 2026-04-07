//! E6: Memory profiling — rope overhead vs flat buffer
//!
//! With arena allocation, memory model changes:
//! - No Arc overhead (no strong/weak counts)
//! - Nodes stored contiguously in Vec<NodeInner>
//! - Per node: ~48 bytes (enum tag + fields)
//! - Leaf nodes additionally store payload in Vec<u8>

use hashrope::{Arena, Node};
use hashrope::rope::NodeInner;

const NODE_INNER_SIZE: usize = 48; // approximate per-node size in the arena Vec

struct MemStats {
    leaf_count: usize,
    internal_count: usize,
    repeat_count: usize,
    payload_bytes: usize,
}

impl MemStats {
    fn new() -> Self {
        Self { leaf_count: 0, internal_count: 0, repeat_count: 0, payload_bytes: 0 }
    }

    fn walk(&mut self, arena: &Arena, node: Node) {
        if let Some(id) = node {
            self.walk_inner(arena, id);
        }
    }

    fn walk_inner(&mut self, arena: &Arena, id: hashrope::NodeId) {
        match arena.node(id) {
            NodeInner::Leaf { data, .. } => {
                self.leaf_count += 1;
                self.payload_bytes += data.len();
            }
            NodeInner::Internal { left, right, .. } => {
                self.internal_count += 1;
                self.walk_inner(arena, *left);
                self.walk_inner(arena, *right);
            }
            NodeInner::Repeat { child, .. } => {
                self.repeat_count += 1;
                self.walk_inner(arena, *child);
            }
        }
    }

    fn total_nodes(&self) -> usize {
        self.leaf_count + self.internal_count + self.repeat_count
    }

    /// Estimated heap bytes: arena Vec storage + leaf data Vecs.
    /// No Arc overhead anymore!
    fn estimated_heap_bytes(&self) -> usize {
        let node_bytes = self.total_nodes() * NODE_INNER_SIZE;
        node_bytes + self.payload_bytes
    }
}

#[test]
fn e6_memory_sequential_rope() {
    let sizes: Vec<usize> = vec![100, 1_000, 10_000, 100_000, 1_000_000];

    println!("\n=== Sequential rope (one byte per leaf, arena) ===");
    println!("{:<10} {:>8} {:>8} {:>8} {:>14} {:>10}",
        "n", "leaves", "intern", "total", "est_bytes", "overhead");
    println!("{}", "-".repeat(66));

    for &n in &sizes {
        let mut a = Arena::new();
        let mut rope: Node = None;
        for i in 0..n {
            let leaf = a.from_bytes(&[(i % 256) as u8]);
            rope = a.concat(rope, leaf);
        }
        a.validate(rope);

        let mut stats = MemStats::new();
        stats.walk(&a, rope);

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

    println!("\n=== Chunked rope (multi-byte leaves, arena) ===");
    println!("{:<10} {:<8} {:>8} {:>8} {:>14} {:>10}",
        "n", "chunk", "leaves", "total", "est_bytes", "overhead");
    println!("{}", "-".repeat(66));

    for &n in &sizes {
        for &chunk in &chunk_sizes {
            if chunk > n { continue; }

            let mut a = Arena::new();
            let mut rope: Node = None;
            let mut offset = 0usize;
            while offset < n {
                let sz = chunk.min(n - offset);
                let data: Vec<u8> = (offset..offset + sz).map(|i| (i % 256) as u8).collect();
                let leaf = a.from_bytes(&data);
                rope = a.concat(rope, leaf);
                offset += sz;
            }

            let mut stats = MemStats::new();
            stats.walk(&a, rope);

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

    println!("\n=== RepeatNode memory vs materialized (arena) ===");
    println!("{:<12} {:>8} {:>14} {:>14} {:>10}",
        "q", "nodes", "rope_bytes", "material_bytes", "savings");
    println!("{}", "-".repeat(64));

    for &q in &q_values {
        let mut a = Arena::new();
        let base: Vec<u8> = (0..base_size).map(|i| (i % 256) as u8).collect();
        let base_rope = a.from_bytes(&base);
        let repeated = a.repeat(base_rope, q);

        assert_eq!(a.len(repeated), q * base_size as u64);

        let mut stats = MemStats::new();
        stats.walk(&a, repeated);

        let est = stats.estimated_heap_bytes();
        let materialized = q as usize * base_size;
        let savings = materialized as f64 / est as f64;

        println!("{:<12} {:>8} {:>14} {:>14} {:>10.0}x",
            q, stats.total_nodes(), est, materialized, savings);
    }
}
