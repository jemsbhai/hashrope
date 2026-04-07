//! E2: Rope concat/split vs naive rehash
//!
//! Claim C2: Rope concat/split is asymptotically faster than rebuilding
//! the hash from scratch.
//!
//! Method:
//!   - Build a rope of n bytes via sequential single-byte concat
//!   - Concat benchmark: split into halves, time the rejoin. Compare to
//!     hashing all n bytes from scratch.
//!   - Split benchmark: time split_at(n/2). Compare to hashing each half
//!     from scratch.
//!
//! Expected: O(log n) for rope ops vs O(n) for naive rehash.

use hashrope::{
    rope_concat, rope_from_bytes, rope_hash, rope_split, validate_rope,
    Node, PolynomialHash,
};
use std::time::Instant;

/// Build a rope by sequential single-byte concats (creates a deep tree).
fn build_sequential_rope(n: usize, ph: &mut PolynomialHash) -> Node {
    let mut rope: Node = None;
    for i in 0..n {
        let leaf = rope_from_bytes(&[(i % 256) as u8], ph);
        rope = rope_concat(&rope, &leaf, ph);
    }
    rope
}

/// Build the flat byte buffer matching the sequential rope.
fn build_flat_buffer(n: usize) -> Vec<u8> {
    (0..n).map(|i| (i % 256) as u8).collect()
}

#[test]
fn e2_concat_vs_rehash() {
    let sizes: Vec<usize> = vec![100, 1_000, 10_000, 100_000];

    println!("\n=== CONCAT: rejoin two halves vs hash-from-scratch ===");
    println!("{:<10} {:>14} {:>14} {:>10}",
        "n", "join_ns", "rehash_ns", "speedup");
    println!("{}", "-".repeat(54));

    for &n in &sizes {
        let mut ph = PolynomialHash::default_hash();
        let rope = build_sequential_rope(n, &mut ph);
        validate_rope(&rope);
        let flat = build_flat_buffer(n);

        // Split into halves
        let (left, right) = rope_split(&rope, n as u64 / 2, &mut ph);

        // Scale iterations: join is fast, rehash is slow for large n
        let iters = (1_000_000u64 / n as u64).max(10);

        // -- Rope join --
        for _ in 0..5 {
            let _ = std::hint::black_box(rope_concat(&left, &right, &mut ph));
        }
        let start = Instant::now();
        for _ in 0..iters {
            let _ = std::hint::black_box(rope_concat(&left, &right, &mut ph));
        }
        let join_ns = start.elapsed().as_nanos() as f64 / iters as f64;

        // -- Naive rehash --
        for _ in 0..5 {
            let _ = std::hint::black_box(ph.hash(&flat));
        }
        let start = Instant::now();
        for _ in 0..iters {
            let _ = std::hint::black_box(ph.hash(&flat));
        }
        let rehash_ns = start.elapsed().as_nanos() as f64 / iters as f64;

        // Correctness
        let rejoined = rope_concat(&left, &right, &mut ph);
        assert_eq!(rope_hash(&rejoined), ph.hash(&flat), "Hash mismatch at n={}", n);

        println!("{:<10} {:>14.1} {:>14.1} {:>10.0}x",
            n, join_ns, rehash_ns, rehash_ns / join_ns);
    }

    println!("\n=== SPLIT: split_at(n/2) vs rehash each half ===");
    println!("{:<10} {:>14} {:>14} {:>10}",
        "n", "split_ns", "rehash_ns", "speedup");
    println!("{}", "-".repeat(54));

    for &n in &sizes {
        let mut ph = PolynomialHash::default_hash();
        let rope = build_sequential_rope(n, &mut ph);
        let flat = build_flat_buffer(n);
        let mid = n / 2;

        let iters = (1_000_000u64 / n as u64).max(10);

        // -- Rope split --
        for _ in 0..5 {
            let _ = std::hint::black_box(rope_split(&rope, mid as u64, &mut ph));
        }
        let start = Instant::now();
        for _ in 0..iters {
            let _ = std::hint::black_box(rope_split(&rope, mid as u64, &mut ph));
        }
        let split_ns = start.elapsed().as_nanos() as f64 / iters as f64;

        // -- Naive: hash each half from scratch --
        let left_flat = &flat[..mid];
        let right_flat = &flat[mid..];
        for _ in 0..5 {
            let _ = std::hint::black_box(ph.hash(left_flat));
            let _ = std::hint::black_box(ph.hash(right_flat));
        }
        let start = Instant::now();
        for _ in 0..iters {
            let _ = std::hint::black_box(ph.hash(left_flat));
            let _ = std::hint::black_box(ph.hash(right_flat));
        }
        let rehash_ns = start.elapsed().as_nanos() as f64 / iters as f64;

        // Correctness
        let (l, r) = rope_split(&rope, mid as u64, &mut ph);
        assert_eq!(rope_hash(&l), ph.hash(left_flat), "Left hash mismatch at n={}", n);
        assert_eq!(rope_hash(&r), ph.hash(right_flat), "Right hash mismatch at n={}", n);

        println!("{:<10} {:>14.1} {:>14.1} {:>10.0}x",
            n, split_ns, rehash_ns, rehash_ns / split_ns);
    }
}
