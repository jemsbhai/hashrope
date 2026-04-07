//! E3: RepeatNode scaling — O(log q) vs O(q) materialized
//!
//! Claim C3: RepeatNode computes the hash of q repetitions in O(log q) time.
//!
//! Method: Create a base rope of 100 bytes, then:
//!   a) rope_repeat(B, q) — should be O(log q)
//!   b) Naive: build q*100 bytes in a Vec, hash from scratch — O(q·|B|)

use hashrope::{rope_from_bytes, rope_hash, rope_repeat, PolynomialHash};
use std::time::Instant;

/// Build a deterministic 100-byte base payload.
fn base_bytes() -> Vec<u8> {
    (0..100u8).collect()
}

#[test]
fn e3_repeat_node_scaling() {
    let mut ph = PolynomialHash::default_hash();
    let base = base_bytes();
    let base_rope = rope_from_bytes(&base, &ph);
    let base_len = base.len();

    let q_values: Vec<u64> = vec![
        10, 100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000, 100_000_000,
    ];

    let iters_rope = 10_000u64;

    println!("\n{:<12} {:>14} {:>14} {:>10}",
        "q", "repeat_ns", "naive_ns", "speedup");
    println!("{}", "-".repeat(56));

    for &q in &q_values {
        // -- Rope repeat --
        // Warm up
        for _ in 0..100 {
            let _ = std::hint::black_box(rope_repeat(&base_rope, q, &mut ph));
        }
        let start = Instant::now();
        for _ in 0..iters_rope {
            let _ = std::hint::black_box(rope_repeat(&base_rope, q, &mut ph));
        }
        let rope_ns = start.elapsed().as_nanos() as f64 / iters_rope as f64;

        // -- Naive baseline --
        let total_bytes = q * base_len as u64;
        // Build materialized buffer and hash it
        let materialized: Vec<u8> = base.iter().copied().cycle().take(total_bytes as usize).collect();
        // Scale iterations down for large buffers
        let iters_naive = if total_bytes > 1_000_000_000 {
            1
        } else if total_bytes > 100_000_000 {
            2
        } else if total_bytes > 1_000_000 {
            5
        } else if total_bytes > 10_000 {
            100
        } else {
            1_000
        };
        // Warm up (skip for very large buffers)
        if total_bytes <= 1_000_000_000 {
            for _ in 0..2 {
                let _ = std::hint::black_box(ph.hash(&materialized));
            }
        }
        let start = Instant::now();
        for _ in 0..iters_naive {
            let _ = std::hint::black_box(ph.hash(&materialized));
        }
        let naive_ns = start.elapsed().as_nanos() as f64 / iters_naive as f64;

        let speedup = format!("{:.0}x", naive_ns / rope_ns);

        println!("{:<12} {:>14.1} {:>14.1} {:>10}",
            q, rope_ns, naive_ns, speedup);

        // Correctness: verify rope hash matches naive hash
        let repeated = rope_repeat(&base_rope, q, &mut ph);
        assert_eq!(
            rope_hash(&repeated),
            ph.hash(&materialized),
            "Hash mismatch at q={}",
            q
        );
    }
}
