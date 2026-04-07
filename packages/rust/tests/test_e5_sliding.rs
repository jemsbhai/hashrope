//! E5: Sliding window throughput — hashrope vs Rabin-Karp
//!
//! Claim C5: Sliding window throughput is competitive with Rabin-Karp.
//!
//! Method: Feed data one byte at a time through each sliding window,
//! measure total throughput (MB/s).
//!
//! Rabin-Karp baseline uses the same Mersenne-prime modular arithmetic
//! and base as hashrope to isolate structural overhead (tree maintenance
//! vs simple ring buffer).
//!
//! Note: Rabin-Karp rolling hash is O(1) per byte. hashrope's sliding
//! window is O(log W) per byte due to rope operations. The question is
//! how large the constant-factor gap is in practice.

use hashrope::{mersenne_mod, mersenne_mul, PolynomialHash, MERSENNE_61};
use hashrope::SlidingWindow;
use std::time::Instant;

/// Simple Rabin-Karp rolling hash using same arithmetic as hashrope.
struct RabinKarp {
    hash: u64,
    base: u64,
    prime: u64,
    window: Vec<u8>,
    w: usize,
    pos: usize,
    /// base^w mod p — precomputed for removing the oldest byte
    base_pow_w: u64,
}

impl RabinKarp {
    fn new(window_size: usize, base: u64, prime: u64) -> Self {
        // Precompute base^w mod p
        let mut base_pow_w = 1u64;
        for _ in 0..window_size {
            base_pow_w = mersenne_mul(base_pow_w, base, prime);
        }
        Self {
            hash: 0,
            base,
            prime,
            window: vec![0u8; window_size],
            w: window_size,
            pos: 0,
            base_pow_w,
        }
    }

    #[inline]
    fn slide(&mut self, byte: u8) {
        let p = self.prime;
        let idx = self.pos % self.w;
        let old_byte = self.window[idx];
        self.window[idx] = byte;

        if self.pos < self.w {
            // Still filling the window
            self.hash = mersenne_mod(
                self.hash as u128 * self.base as u128 + byte as u128,
                p,
            );
        } else {
            // Remove old_byte * base^w, shift, add new byte
            let remove = mersenne_mul(old_byte as u64, self.base_pow_w, p);
            let shifted = if self.hash >= remove {
                self.hash - remove
            } else {
                self.hash + p - remove
            };
            self.hash = mersenne_mod(
                shifted as u128 * self.base as u128 + byte as u128,
                p,
            );
        }
        self.pos += 1;
    }
}

/// Generate deterministic pseudo-random bytes.
fn gen_random_bytes(n: usize, seed: u64) -> Vec<u8> {
    let mut buf = vec![0u8; n];
    let mut state = seed;
    for byte in buf.iter_mut() {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *byte = (state >> 56) as u8;
    }
    buf
}

#[test]
fn e5_sliding_window_throughput() {
    let window_sizes: Vec<usize> = vec![32, 64, 128, 256, 1024];
    let data_sizes: Vec<usize> = vec![1_024, 10_240, 102_400, 1_048_576];

    println!("\n{:<10} {:<8} {:>14} {:>14} {:>10}",
        "data", "window", "hashrope", "rabin_karp", "ratio");
    println!("{:<10} {:<8} {:>14} {:>14} {:>10}",
        "", "", "(MB/s)", "(MB/s)", "(rk/hr)");
    println!("{}", "-".repeat(62));

    for &data_size in &data_sizes {
        let data = gen_random_bytes(data_size, 42);

        for &w in &window_sizes {
            let iters = (10_000_000u64 / data_size as u64).max(1);

            // -- hashrope SlidingWindow --
            // d_max = w, m_max = 1 (keep exactly w bytes in window)
            let start = Instant::now();
            for _ in 0..iters {
                let mut sw = SlidingWindow::new(w as u64, 1, MERSENNE_61, 131);
                for &byte in &data {
                    sw.append_bytes(&[byte]);
                }
                let _ = std::hint::black_box(sw.final_hash());
            }
            let hr_secs = start.elapsed().as_secs_f64() / iters as f64;
            let hr_mbps = (data_size as f64 / (1024.0 * 1024.0)) / hr_secs;

            // -- Rabin-Karp --
            let start = Instant::now();
            for _ in 0..iters {
                let mut rk = RabinKarp::new(w, 131, MERSENNE_61);
                for &byte in &data {
                    rk.slide(byte);
                }
                let _ = std::hint::black_box(rk.hash);
            }
            let rk_secs = start.elapsed().as_secs_f64() / iters as f64;
            let rk_mbps = (data_size as f64 / (1024.0 * 1024.0)) / rk_secs;

            let label = if data_size >= 1_048_576 {
                format!("{} MB", data_size / 1_048_576)
            } else {
                format!("{} KB", data_size / 1_024)
            };

            println!("{:<10} {:<8} {:>14.1} {:>14.1} {:>10.1}x",
                label, w, hr_mbps, rk_mbps, rk_mbps / hr_mbps);
        }
    }
}
