//! E4: Phi accumulator — verify O(log q) time scaling
//!
//! Claim C4: The phi accumulator (geometric series via repeated doubling)
//! runs in O(log q) time.
//!
//! Method: Time phi(q, alpha, p) for q from 10 to 10^7.
//! On a log-log plot, O(log q) growth should show a very gentle slope
//! (sub-linear in log q — essentially near-constant for practical q).

use hashrope::{phi, PolynomialHash, MERSENNE_61};
use std::time::Instant;

/// Number of iterations per measurement to get stable timings.
const ITERS: u64 = 100_000;

#[test]
fn e4_phi_scaling() {
    let ph = PolynomialHash::default_hash();
    let alpha = 257u64; // typical base
    let p = MERSENNE_61;

    let q_values: Vec<u64> = vec![
        10, 100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000,
    ];

    println!("\n{:<12} {:>12} {:>12} {:>8}", "q", "time_ns", "log2(q)", "bits(q)");
    println!("{}", "-".repeat(50));

    for &q in &q_values {
        // Warm up
        for _ in 0..1000 {
            let _ = std::hint::black_box(phi(q, alpha, p));
        }

        let start = Instant::now();
        for _ in 0..ITERS {
            let _ = std::hint::black_box(phi(q, alpha, p));
        }
        let elapsed = start.elapsed();
        let ns_per_call = elapsed.as_nanos() as f64 / ITERS as f64;
        let log2_q = (q as f64).log2();
        let bits = 64 - q.leading_zeros();

        println!("{:<12} {:>12.1} {:>12.2} {:>8}",
            q, ns_per_call, log2_q, bits);
    }

    // Sanity: phi should produce correct values
    assert_eq!(phi(1, alpha, p), 1);
    assert_eq!(phi(2, alpha, p), hashrope::mersenne_mod((1 + alpha) as u128, p));
}
