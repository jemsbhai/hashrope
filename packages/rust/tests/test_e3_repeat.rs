//! E3: RepeatNode scaling — O(log q) vs O(q) materialized

use hashrope::{Arena, PolynomialHash};
use std::time::Instant;

fn base_bytes() -> Vec<u8> {
    (0..100u8).collect()
}

#[test]
fn e3_repeat_node_scaling() {
    let base = base_bytes();
    let base_len = base.len();
    let ph = PolynomialHash::default_hash();

    let q_values: Vec<u64> = vec![
        10, 100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000, 100_000_000,
    ];

    let iters_rope = 10_000u64;

    println!("\n{:<12} {:>14} {:>14} {:>10}",
        "q", "repeat_ns", "naive_ns", "speedup");
    println!("{}", "-".repeat(56));

    for &q in &q_values {
        // -- Rope repeat --
        let mut a = Arena::new();
        let base_rope = a.from_bytes(&base);

        for _ in 0..100 {
            let _ = std::hint::black_box(a.repeat(base_rope, q));
        }
        let start = Instant::now();
        for _ in 0..iters_rope {
            let _ = std::hint::black_box(a.repeat(base_rope, q));
        }
        let rope_ns = start.elapsed().as_nanos() as f64 / iters_rope as f64;

        // -- Naive baseline --
        let total_bytes = q * base_len as u64;
        let materialized: Vec<u8> = base.iter().copied().cycle().take(total_bytes as usize).collect();
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
        println!("{:<12} {:>14.1} {:>14.1} {:>10}", q, rope_ns, naive_ns, speedup);

        // Correctness
        let repeated = a.repeat(base_rope, q);
        assert_eq!(a.hash(repeated), ph.hash(&materialized), "Hash mismatch at q={}", q);
    }
}
