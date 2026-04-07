//! E2: Rope concat/split vs naive rehash

use hashrope::{Arena, PolynomialHash};
use std::time::Instant;

fn build_flat_buffer(n: usize) -> Vec<u8> {
    (0..n).map(|i| (i % 256) as u8).collect()
}

#[test]
fn e2_concat_vs_rehash() {
    let sizes: Vec<usize> = vec![100, 1_000, 10_000, 100_000];
    let ph = PolynomialHash::default_hash();

    println!("\n=== CONCAT: rejoin two halves vs hash-from-scratch ===");
    println!("{:<10} {:>14} {:>14} {:>10}", "n", "join_ns", "rehash_ns", "speedup");
    println!("{}", "-".repeat(54));

    for &n in &sizes {
        let mut a = Arena::new();
        let mut rope = None;
        for i in 0..n {
            let leaf = a.from_bytes(&[(i % 256) as u8]);
            rope = a.concat(rope, leaf);
        }
        a.validate(rope);
        let flat = build_flat_buffer(n);

        let (left, right) = a.split(rope, n as u64 / 2);
        let iters = (1_000_000u64 / n as u64).max(10);

        // -- Rope join --
        for _ in 0..5 {
            let _ = std::hint::black_box(a.concat(left, right));
        }
        let start = Instant::now();
        for _ in 0..iters {
            let _ = std::hint::black_box(a.concat(left, right));
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
        let rejoined = a.concat(left, right);
        assert_eq!(a.hash(rejoined), ph.hash(&flat), "Hash mismatch at n={}", n);

        println!("{:<10} {:>14.1} {:>14.1} {:>10.0}x", n, join_ns, rehash_ns, rehash_ns / join_ns);
    }

    println!("\n=== SPLIT: split_at(n/2) vs rehash each half ===");
    println!("{:<10} {:>14} {:>14} {:>10}", "n", "split_ns", "rehash_ns", "speedup");
    println!("{}", "-".repeat(54));

    for &n in &sizes {
        let mut a = Arena::new();
        let mut rope = None;
        for i in 0..n {
            let leaf = a.from_bytes(&[(i % 256) as u8]);
            rope = a.concat(rope, leaf);
        }
        let flat = build_flat_buffer(n);
        let mid = n / 2;
        let iters = (1_000_000u64 / n as u64).max(10);

        for _ in 0..5 {
            let _ = std::hint::black_box(a.split(rope, mid as u64));
        }
        let start = Instant::now();
        for _ in 0..iters {
            let _ = std::hint::black_box(a.split(rope, mid as u64));
        }
        let split_ns = start.elapsed().as_nanos() as f64 / iters as f64;

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

        let (l, r) = a.split(rope, mid as u64);
        assert_eq!(a.hash(l), ph.hash(left_flat), "Left hash mismatch at n={}", n);
        assert_eq!(a.hash(r), ph.hash(right_flat), "Right hash mismatch at n={}", n);

        println!("{:<10} {:>14.1} {:>14.1} {:>10.0}x", n, split_ns, rehash_ns, rehash_ns / split_ns);
    }
}
