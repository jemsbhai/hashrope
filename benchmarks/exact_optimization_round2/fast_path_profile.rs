//! Targeted, untimed-fixture diagnostics for exact ordinary-tree traversal.
use hashrope::{Arena, Node};
use std::{hint::black_box, time::Instant};

fn chunked(a: &mut Arena, bytes: &[u8], chunk: usize) -> Node {
    let mut root = None;
    for bytes in bytes.chunks(chunk) {
        let leaf = a.from_bytes(bytes);
        root = a.concat(root, leaf);
    }
    root
}

fn measure(name: &str, n: usize, iterations: usize, mut f: impl FnMut() -> u64, expected: u64) {
    assert_eq!(f(), expected);
    let mut samples = Vec::new();
    for _ in 0..7 {
        let start = Instant::now();
        for _ in 0..iterations { assert_eq!(black_box(f()), expected); }
        samples.push(start.elapsed().as_secs_f64() * 1e6 / iterations as f64);
    }
    samples.sort_by(f64::total_cmp);
    println!("{{\"case\":\"{name}\",\"bytes\":{n},\"median_us\":{},\"samples_us\":{:?},\"iterations\":{iterations}}}", samples[3], samples);
}

fn main() {
    for n in [10_000usize, 100_000, 1_000_000] {
        let bytes: Vec<u8> = (0..n).map(|i| (i * 37) as u8).collect();
        let mut left_arena = Arena::new_lazy();
        let left = left_arena.from_bytes(&bytes);
        let mut other = Arena::new_lazy();
        let independent = other.from_bytes(&bytes);
        let rechunked = chunked(&mut other, &bytes, 257);
        let mut changed = bytes.clone();
        changed[n - 1] ^= 1;
        let late = other.from_bytes(&changed);
        changed[0] ^= 1;
        let early = other.from_bytes(&changed);
        let rounds = if n < 1_000_000 { 1000 } else { 100 };
        measure("independent_equal", n, rounds, || left_arena.lcp_across(left, &other, independent), n as u64);
        measure("rechunked_equal", n, rounds, || left_arena.lcp_across(left, &other, rechunked), n as u64);
        measure("late_mismatch", n, rounds, || left_arena.lcp_across(left, &other, late), n as u64 - 1);
        measure("early_mismatch", n, rounds, || left_arena.lcp_across(left, &other, early), 0);
        measure("shifted_equal_range", n, rounds, || left_arena.lcp_range_across(left, 13, &other, rechunked, 13, n as u64 - 27), n as u64 - 27);
        let (prefix, _) = left_arena.split(left, n as u64 - 1);
        let replacement = left_arena.from_bytes(&[bytes[n - 1] ^ 1]);
        let shared = left_arena.concat(prefix, replacement);
        measure("shared_late_mismatch", n, rounds, || left_arena.lcp(left, shared), n as u64 - 1);
    }
    for position in [0usize, 1, 7, 8, 63, 127, 255, 510, 511] {
        let bytes = [b'x'; 512];
        let mut changed = bytes;
        changed[position] = b'y';
        let mut arena = Arena::new_lazy();
        let a = arena.from_bytes(&bytes);
        let b = arena.from_bytes(&changed);
        measure(&format!("leaf_mismatch_{position}"), 512, 100_000, || arena.lcp(a, b), position as u64);
    }
}
