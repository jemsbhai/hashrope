//! E7: BB[2/7] tree height verification
//!
//! Claim C7: BB[2/7] balancing maintains O(log w) tree height in practice.
//!
//! Method: Build ropes by sequential single-byte insertions (worst case),
//! measure actual height, compare to theoretical upper bound.
//!
//! Theoretical bound for BB[α] with α = 2/7:
//!   h ≤ log(w) / log(1 / (1 − α))
//!     = log(w) / log(7/5)
//!     ≈ 2.06 · log₂(w)    (using ln(7/5) ≈ 0.3365)
//!
//! We also test random-order insertions via join of random-sized sub-ropes.

use hashrope::{rope_concat, rope_from_bytes, rope_height, rope_len, validate_rope, Node, PolynomialHash};

/// Theoretical height bound for BB[2/7].
/// h_max = ceil(log(w) / log(7/5))
fn bb27_height_bound(w: u64) -> u64 {
    if w <= 1 {
        return 0;
    }
    // log(w) / log(7/5) = ln(w) / ln(1.4)
    let bound = (w as f64).ln() / (7.0_f64 / 5.0).ln();
    bound.ceil() as u64
}

/// Build a rope by appending single bytes sequentially (worst-case for balance).
fn build_sequential(n: usize) -> (Node, PolynomialHash) {
    let mut ph = PolynomialHash::default_hash();
    let mut rope: Node = None;
    for i in 0..n {
        let byte = (i % 256) as u8;
        let leaf = rope_from_bytes(&[byte], &ph);
        rope = rope_concat(&rope, &leaf, &mut ph);
    }
    (rope, ph)
}

/// Build a rope by joining random-sized chunks (more realistic usage).
fn build_chunked(n: usize, chunk_sizes: &[usize]) -> (Node, PolynomialHash) {
    let mut ph = PolynomialHash::default_hash();
    let mut rope: Node = None;
    let mut offset = 0usize;
    let mut chunk_idx = 0;
    while offset < n {
        let sz = chunk_sizes[chunk_idx % chunk_sizes.len()].min(n - offset);
        let data: Vec<u8> = (offset..offset + sz).map(|i| (i % 256) as u8).collect();
        let chunk = rope_from_bytes(&data, &ph);
        rope = rope_concat(&rope, &chunk, &mut ph);
        offset += sz;
        chunk_idx += 1;
    }
    (rope, ph)
}

#[test]
fn e7_sequential_insertion() {
    let sizes: Vec<usize> = vec![100, 1_000, 10_000, 100_000];

    println!("\n{:<10} {:>8} {:>8} {:>10} {:>6}", "n", "weight", "height", "bound", "pass");
    println!("{}", "-".repeat(50));

    for &n in &sizes {
        let (rope, _ph) = build_sequential(n);
        validate_rope(&rope);

        let h = rope_height(&rope);
        let w = n as u64; // one leaf per byte, so weight = n
        let bound = bb27_height_bound(w);

        println!("{:<10} {:>8} {:>8} {:>10} {:>6}",
            n, w, h, bound, if h <= bound { "OK" } else { "FAIL" });

        assert!(
            h <= bound,
            "Height {} exceeds BB[2/7] bound {} for n={} (weight={})",
            h, bound, n, w
        );
    }
}

#[test]
fn e7_chunked_insertion() {
    // Varied chunk sizes to simulate realistic usage
    let chunk_sizes = vec![1, 3, 7, 15, 31, 64, 128, 5, 11, 42];
    let sizes: Vec<usize> = vec![100, 1_000, 10_000, 100_000];

    println!("\n{:<10} {:>8} {:>8} {:>10} {:>6}", "n", "weight", "height", "bound", "pass");
    println!("{}", "-".repeat(50));

    for &n in &sizes {
        let (rope, _ph) = build_chunked(n, &chunk_sizes);
        validate_rope(&rope);

        let h = rope_height(&rope);
        let len = rope_len(&rope);
        assert_eq!(len, n as u64);

        // Weight counts leaves, not bytes — with multi-byte leaves, weight < n.
        // We need weight from the rope itself.
        let w = rope.as_ref().map_or(0, |node| node.weight());
        let bound = bb27_height_bound(w);

        println!("{:<10} {:>8} {:>8} {:>10} {:>6}",
            n, w, h, bound, if h <= bound { "OK" } else { "FAIL" });

        assert!(
            h <= bound,
            "Height {} exceeds BB[2/7] bound {} for n={} (weight={})",
            h, bound, n, w
        );
    }
}
