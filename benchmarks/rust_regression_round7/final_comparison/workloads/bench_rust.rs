// Standalone release-mode harness, compiled by bench_rust.py against either tree.
use hashrope::{Arena, NodeInner};
#[cfg(feature = "exact")]
use hashrope::Node;
use std::{hint::black_box, mem::size_of, time::Instant};

fn measure<T>(name: &str, iterations: usize, mut f: impl FnMut() -> T) {
    for _ in 0..3 { black_box(f()); }
    for sample in 0..7 {
        let start = Instant::now();
        for _ in 0..iterations { black_box(f()); }
        let ns = start.elapsed().as_nanos() as f64 / iterations as f64;
        println!("time\t{name}\t{sample}\t{iterations}\t{ns}");
    }
}

#[cfg(feature = "exact")]
fn hash_lcp(arena: &mut Arena, a: Node, b: Node) -> u64 {
    let (mut lo, mut hi) = (0, arena.len(a).min(arena.len(b)));
    while lo < hi {
        let mid = lo + (hi - lo + 1) / 2;
        if arena.substr_hash(a, 0, mid) == arena.substr_hash(b, 0, mid) {
            lo = mid;
        } else { hi = mid - 1; }
    }
    lo
}

fn main() {
    let data: Vec<u8> = (0..65536).map(|i| (i % 256) as u8).collect();
    measure("core/from_bytes_64KiB", 20, || {
        let mut a = Arena::new();
        black_box(a.from_bytes(black_box(&data)));
    });
    let mut a = Arena::new();
    let left = a.from_bytes(&data);
    let right = a.from_bytes(&data);
    measure("core/concat_balanced_64KiB", 2000, || a.concat(black_box(left), black_box(right)));
    measure("core/split_chunked_64KiB", 100, || a.split(black_box(left), 32769));
    measure("core/substr_warm_64KiB", 1000, || a.substr_hash(black_box(left), 1, 65534));
    measure("core/root_hash", 1000000, || a.hash(black_box(left)));
    let pattern = a.from_bytes(b"ab");
    let repeated = a.repeat(pattern, 1u64 << 32);
    let suffix = a.from_bytes(b"!");
    measure("core/repeat_join_2pow32", 100, || a.concat(black_box(repeated), suffix));
    let mut memory_arena = Arena::new();
    memory_arena.from_bytes(&data);
    println!("memory\tnode_inner_bytes\t{}", size_of::<NodeInner>());
    println!("memory\tarena_inline_bytes\t{}", size_of::<Arena>());
    println!("memory\tnodes_for_64KiB\t{}", memory_arena.node_count());

    #[cfg(feature = "exact")]
    {
        for size in [65536usize, 1048576] {
            let raw: Vec<u8> = (0..size).map(|i| (i % 256) as u8).collect();
            let mut arena = Arena::new();
            let first = arena.from_bytes(&raw);
            let independent = arena.from_bytes(&raw);
            let mut changed = raw.clone();
            changed[0] ^= 1;
            let early = arena.from_bytes(&changed);
            changed[0] ^= 1;
            changed[size - 1] ^= 1;
            let late = arena.from_bytes(&changed);
            let (prefix, _) = arena.split(first, size as u64 - 1);
            let ending = arena.from_bytes(&changed[size - 1..]);
            let shared_late = arena.concat(prefix, ending);
            for (label, b, expected) in [("same_root", first, size as u64),
                ("independent_equal", independent, size as u64),
                ("early_mismatch", early, 0), ("late_mismatch", late, size as u64 - 1),
                ("shared_late_mismatch", shared_late, size as u64 - 1)] {
                assert_eq!(arena.lcp(first, b), expected);
                assert_eq!(arena.equal(first, b), expected == size as u64);
                assert_eq!(hash_lcp(&mut arena, first, b), expected);
                let name = format!("compare/{size}/{label}");
                let iterations = if label == "same_root" || label == "early_mismatch" { 2000 } else { 100 };
                measure(&format!("{name}/exact_lcp"), iterations, || arena.lcp(black_box(first), black_box(b)));
                measure(&format!("{name}/exact_equal"), iterations, || arena.equal(black_box(first), black_box(b)));
                measure(&format!("{name}/hash_only_lcp"), 100, || hash_lcp(&mut arena, black_box(first), black_box(b)));
            }
        }
        let mut arena = Arena::new();
        let p = arena.from_bytes(b"ab");
        let a = arena.repeat(p, 1u64 << 40);
        let q = arena.from_bytes(b"abab");
        let b = arena.repeat(q, 1u64 << 39);
        assert_eq!(arena.lcp(a, b), 1u64 << 41);
        assert!(arena.equal(a, b));
        measure("compare/repeat_2TiB/exact_lcp", 2000, || arena.lcp(black_box(a), black_box(b)));
        measure("compare/repeat_2TiB/exact_equal", 2000, || arena.equal(black_box(a), black_box(b)));
        measure("compare/repeat_2TiB/hash_only_lcp", 100, || hash_lcp(&mut arena, black_box(a), black_box(b)));
        let mut a = arena.from_bytes(b"ab");
        let mut b = arena.from_bytes(b"ab");
        for _ in 0..40 {
            a = arena.concat(a, a);
            b = arena.concat(b, b);
        }
        assert_eq!(arena.lcp(a, b), 1u64 << 41);
        assert!(arena.equal(a, b));
        measure("compare/self_concat_DAG_2TiB/exact_lcp", 100, || arena.lcp(black_box(a), black_box(b)));
        measure("compare/self_concat_DAG_2TiB/exact_equal", 100, || arena.equal(black_box(a), black_box(b)));
        measure("compare/self_concat_DAG_2TiB/hash_only_lcp", 100, || hash_lcp(&mut arena, black_box(a), black_box(b)));
    }
}
