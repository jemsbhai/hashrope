//! Diagnostic: substr_hash scaling test
//!
//! Measures substr_hash time across different tree sizes to verify
//! O(log w) complexity after the full-node shortcut fix.
//! Uses std::time::Instant for direct measurement — no Criterion.

use hashrope::Arena;
use std::time::Instant;

/// Build a rope of `n` single-byte leaves (worst case for tree depth).
fn build_rope(n: usize) -> (Arena, Option<u32>) {
    let mut arena = Arena::new();
    let mut node = None;
    for i in 0..n {
        let leaf = arena.from_bytes(&[(i % 256) as u8]);
        node = arena.concat(node, leaf);
    }
    (arena, node)
}

#[test]
fn diag_substr_hash_scaling() {
    let sizes = [100, 1_000, 10_000, 100_000, 1_000_000];
    let iterations = 10_000;

    println!("\n=== substr_hash scaling (middle-half query, {} iters) ===", iterations);
    println!("{:<12} {:>8} {:>12} {:>12} {:>10}",
        "n", "height", "total_ms", "per_call_ns", "log2(n)");
    println!("{}", "-".repeat(60));

    for &n in &sizes {
        let (mut arena, node) = build_rope(n);

        // Query the middle half: start = n/4, length = n/2
        let start = n as u64 / 4;
        let length = n as u64 / 2;

        // Verify correctness first
        let expected = {
            let bytes = arena.to_bytes(node);
            arena.hash_bytes(&bytes[start as usize..(start + length) as usize])
        };
        let got = arena.substr_hash(node, start, length);
        assert_eq!(got, expected, "Correctness check failed for n={}", n);

        // Warm up the power cache
        let _ = arena.substr_hash(node, start, length);

        // Time it
        let t = Instant::now();
        for _ in 0..iterations {
            let _ = arena.substr_hash(node, start, length);
        }
        let elapsed = t.elapsed();

        let h = arena.height(node);
        let total_ms = elapsed.as_secs_f64() * 1000.0;
        let per_call_ns = elapsed.as_nanos() as f64 / iterations as f64;
        let log2_n = (n as f64).log2();

        println!("{:<12} {:>8} {:>12.3} {:>12.1} {:>10.1}",
            n, h, total_ms, per_call_ns, log2_n);
    }

    println!();
    println!("=== substr_hash full-node query (early exit path) ===");
    println!("{:<12} {:>12}", "n", "per_call_ns");
    println!("{}", "-".repeat(30));

    for &n in &sizes {
        let (mut arena, node) = build_rope(n);

        // Full-node query — should hit early exit immediately
        let _ = arena.substr_hash(node, 0, n as u64);

        let t = Instant::now();
        for _ in 0..iterations {
            let _ = arena.substr_hash(node, 0, n as u64);
        }
        let elapsed = t.elapsed();
        let per_call_ns = elapsed.as_nanos() as f64 / iterations as f64;

        println!("{:<12} {:>12.1}", n, per_call_ns);
    }

    println!();
    println!("=== split_rejoin for comparison ===");
    println!("{:<12} {:>12}", "n", "per_call_ns");
    println!("{}", "-".repeat(30));

    // Only test smaller sizes since split allocates new nodes
    for &n in &[100, 1_000, 10_000] {
        let iters = if n <= 1_000 { 10_000 } else { 1_000 };
        let (mut arena, node) = build_rope(n);

        let t = Instant::now();
        for _ in 0..iters {
            let (left, right) = arena.split(node, n as u64 / 2);
            let rejoined = arena.concat(left, right);
            std::hint::black_box(arena.hash(rejoined));
        }
        let elapsed = t.elapsed();
        let per_call_ns = elapsed.as_nanos() as f64 / iters as f64;

        println!("{:<12} {:>12.1}", n, per_call_ns);
    }
}
