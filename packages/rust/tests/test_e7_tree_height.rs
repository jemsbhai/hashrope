//! E7: BB[2/7] tree height verification

use hashrope::Arena;

fn bb27_height_bound(w: u64) -> u64 {
    if w <= 1 { return 0; }
    let bound = (w as f64).ln() / (7.0_f64 / 5.0).ln();
    bound.ceil() as u64
}

fn build_sequential(n: usize) -> Arena {
    let mut a = Arena::new();
    let mut rope = None;
    for i in 0..n {
        let leaf = a.from_bytes(&[(i % 256) as u8]);
        rope = a.concat(rope, leaf);
    }
    // Store the final rope id somewhere — we'll return the arena
    // and reconstruct. Actually let's just return (Arena, Node).
    a
}

#[test]
fn e7_sequential_insertion() {
    let sizes: Vec<usize> = vec![100, 1_000, 10_000, 100_000];

    println!("\n{:<10} {:>8} {:>8} {:>10} {:>6}", "n", "weight", "height", "bound", "pass");
    println!("{}", "-".repeat(50));

    for &n in &sizes {
        let mut a = Arena::new();
        let mut rope = None;
        for i in 0..n {
            let leaf = a.from_bytes(&[(i % 256) as u8]);
            rope = a.concat(rope, leaf);
        }
        a.validate(rope);

        let h = a.height(rope);
        let w = n as u64;
        let bound = bb27_height_bound(w);

        println!("{:<10} {:>8} {:>8} {:>10} {:>6}",
            n, w, h, bound, if h <= bound { "OK" } else { "FAIL" });

        assert!(h <= bound, "Height {} exceeds bound {} for n={}", h, bound, n);
    }
}

#[test]
fn e7_chunked_insertion() {
    let chunk_sizes = vec![1, 3, 7, 15, 31, 64, 128, 5, 11, 42];
    let sizes: Vec<usize> = vec![100, 1_000, 10_000, 100_000];

    println!("\n{:<10} {:>8} {:>8} {:>10} {:>6}", "n", "weight", "height", "bound", "pass");
    println!("{}", "-".repeat(50));

    for &n in &sizes {
        let mut a = Arena::new();
        let mut rope = None;
        let mut offset = 0usize;
        let mut chunk_idx = 0;
        while offset < n {
            let sz = chunk_sizes[chunk_idx % chunk_sizes.len()].min(n - offset);
            let data: Vec<u8> = (offset..offset + sz).map(|i| (i % 256) as u8).collect();
            let chunk = a.from_bytes(&data);
            rope = a.concat(rope, chunk);
            offset += sz;
            chunk_idx += 1;
        }
        a.validate(rope);

        let h = a.height(rope);
        assert_eq!(a.len(rope), n as u64);
        let w = a.weight(rope);
        let bound = bb27_height_bound(w);

        println!("{:<10} {:>8} {:>8} {:>10} {:>6}",
            n, w, h, bound, if h <= bound { "OK" } else { "FAIL" });

        assert!(h <= bound, "Height {} exceeds bound {} for n={} (w={})", h, bound, n, w);
    }
}
