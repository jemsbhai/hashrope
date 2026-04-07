//! E1: Sequential hashing throughput (MB/s)
//!
//! Claim C1: Polynomial hash throughput is competitive with general-purpose hashes.
//!
//! Method: Feed identical byte slices to each hasher, measure wall-clock time,
//! compute throughput in MB/s.
//!
//! Expected: hashrope will be slower than BLAKE3/xxhash3 (SIMD-optimized).
//! The purpose is to quantify the gap honestly.
//!
//! Fairness: Each hashrope iteration uses a fresh PolynomialHash (no pre-warmed
//! power cache). All hashers process the same byte slices.

use hashrope::PolynomialHash;
use sha2::{Sha256, Digest};
use std::time::Instant;

/// Generate deterministic pseudo-random bytes (simple LCG, not for crypto).
fn gen_random_bytes(n: usize, seed: u64) -> Vec<u8> {
    let mut buf = vec![0u8; n];
    let mut state = seed;
    for byte in buf.iter_mut() {
        // LCG parameters from Numerical Recipes
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *byte = (state >> 56) as u8;
    }
    buf
}

fn bench_hashrope(data: &[u8], iters: u64) -> f64 {
    // Warm up
    for _ in 0..3 {
        let ph = PolynomialHash::default_hash();
        let _ = std::hint::black_box(ph.hash(data));
    }
    let start = Instant::now();
    for _ in 0..iters {
        let ph = PolynomialHash::default_hash();
        let _ = std::hint::black_box(ph.hash(data));
    }
    start.elapsed().as_secs_f64() / iters as f64
}

fn bench_blake3(data: &[u8], iters: u64) -> f64 {
    // Warm up
    for _ in 0..3 {
        let _ = std::hint::black_box(blake3::hash(data));
    }
    let start = Instant::now();
    for _ in 0..iters {
        let _ = std::hint::black_box(blake3::hash(data));
    }
    start.elapsed().as_secs_f64() / iters as f64
}

fn bench_sha256(data: &[u8], iters: u64) -> f64 {
    // Warm up
    for _ in 0..3 {
        let mut hasher = Sha256::new();
        hasher.update(data);
        let _ = std::hint::black_box(hasher.finalize());
    }
    let start = Instant::now();
    for _ in 0..iters {
        let mut hasher = Sha256::new();
        hasher.update(data);
        let _ = std::hint::black_box(hasher.finalize());
    }
    start.elapsed().as_secs_f64() / iters as f64
}

fn bench_xxh3(data: &[u8], iters: u64) -> f64 {
    // Warm up
    for _ in 0..3 {
        let _ = std::hint::black_box(xxhash_rust::xxh3::xxh3_64(data));
    }
    let start = Instant::now();
    for _ in 0..iters {
        let _ = std::hint::black_box(xxhash_rust::xxh3::xxh3_64(data));
    }
    start.elapsed().as_secs_f64() / iters as f64
}

fn throughput_mbps(bytes: usize, secs: f64) -> f64 {
    (bytes as f64 / (1024.0 * 1024.0)) / secs
}

#[test]
fn e1_sequential_throughput() {
    let sizes: Vec<usize> = vec![
        1_024,          // 1 KB
        10_240,         // 10 KB
        102_400,        // 100 KB
        1_048_576,      // 1 MB
        10_485_760,     // 10 MB
        104_857_600,    // 100 MB
    ];

    println!("\n{:<12} {:>12} {:>12} {:>12} {:>12}",
        "size", "hashrope", "blake3", "sha256", "xxh3");
    println!("{:<12} {:>12} {:>12} {:>12} {:>12}",
        "", "(MB/s)", "(MB/s)", "(MB/s)", "(MB/s)");
    println!("{}", "-".repeat(64));

    for &size in &sizes {
        let data = gen_random_bytes(size, 42);

        // Scale iterations inversely with size
        let iters = (100_000_000u64 / size as u64).max(3);

        let t_hr = bench_hashrope(&data, iters);
        let t_b3 = bench_blake3(&data, iters);
        let t_sha = bench_sha256(&data, iters);
        let t_xx = bench_xxh3(&data, iters);

        let label = if size >= 1_048_576 {
            format!("{} MB", size / 1_048_576)
        } else {
            format!("{} KB", size / 1_024)
        };

        println!("{:<12} {:>12.1} {:>12.1} {:>12.1} {:>12.1}",
            label,
            throughput_mbps(size, t_hr),
            throughput_mbps(size, t_b3),
            throughput_mbps(size, t_sha),
            throughput_mbps(size, t_xx));
    }
}
