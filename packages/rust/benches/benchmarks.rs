use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use hashrope::{
    rope_concat, rope_from_bytes, rope_hash, rope_repeat, rope_split, rope_substr_hash,
    Node, PolynomialHash, SlidingWindow,
};
use hashrope::polynomial_hash::{phi, mersenne_mul, MERSENNE_61};

fn bench_phi(c: &mut Criterion) {
    let alpha = 131u64;
    let x_d = mersenne_mul(alpha, alpha, MERSENNE_61); // some alpha^d

    let mut group = c.benchmark_group("phi");
    for q in [10, 100, 1_000, 10_000, 100_000, 1_000_000] {
        group.bench_with_input(BenchmarkId::from_parameter(q), &q, |b, &q| {
            b.iter(|| phi(black_box(q), black_box(x_d), MERSENNE_61));
        });
    }
    group.finish();
}

fn bench_rope_concat_sequential(c: &mut Criterion) {
    let mut group = c.benchmark_group("rope_concat_sequential");
    for n in [100, 1_000, 10_000, 100_000] {
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &n| {
            b.iter(|| {
                let ph = PolynomialHash::default_hash();
                let mut h = PolynomialHash::default_hash();
                let mut node: Node = None;
                for i in 0..n {
                    let leaf = rope_from_bytes(&[i as u8], &ph);
                    node = rope_concat(&node, &leaf, &mut h);
                }
                black_box(rope_hash(&node))
            });
        });
    }
    group.finish();
}

fn bench_repeat_vs_materialized(c: &mut Criterion) {
    let mut group = c.benchmark_group("repeat_vs_materialized");
    let ph = PolynomialHash::default_hash();
    let pattern = rope_from_bytes(b"abcdefgh", &ph);

    for q in [10, 100, 1_000, 10_000, 100_000, 1_000_000] {
        group.bench_with_input(BenchmarkId::new("repeat_node", q), &q, |b, &q| {
            let mut h = PolynomialHash::default_hash();
            b.iter(|| {
                let rep = rope_repeat(&pattern, black_box(q), &mut h);
                black_box(rope_hash(&rep))
            });
        });

        if q <= 10_000 {
            group.bench_with_input(BenchmarkId::new("materialized", q), &q, |b, &q| {
                b.iter(|| {
                    let data: Vec<u8> = b"abcdefgh".iter().copied().cycle().take(8 * q as usize).collect();
                    let h = PolynomialHash::default_hash();
                    black_box(h.hash(&data))
                });
            });
        }
    }
    group.finish();
}

fn bench_split_rejoin(c: &mut Criterion) {
    let ph = PolynomialHash::default_hash();
    let mut h = PolynomialHash::default_hash();

    // Build a rope of 10K leaves
    let mut node: Node = None;
    for i in 0..10_000u16 {
        let leaf = rope_from_bytes(&[i as u8], &ph);
        node = rope_concat(&node, &leaf, &mut h);
    }

    c.bench_function("split_rejoin_10k", |b| {
        let mut h = PolynomialHash::default_hash();
        b.iter(|| {
            let (left, right) = rope_split(&node, black_box(5000), &mut h);
            let rejoined = rope_concat(&left, &right, &mut h);
            black_box(rope_hash(&rejoined))
        });
    });
}

fn bench_substr_hash(c: &mut Criterion) {
    let ph = PolynomialHash::default_hash();
    let mut h = PolynomialHash::default_hash();

    let mut node: Node = None;
    for i in 0..10_000u16 {
        let leaf = rope_from_bytes(&[i as u8], &ph);
        node = rope_concat(&node, &leaf, &mut h);
    }

    c.bench_function("substr_hash_10k", |b| {
        let mut h = PolynomialHash::default_hash();
        b.iter(|| {
            black_box(rope_substr_hash(&node, black_box(2500), black_box(5000), &mut h))
        });
    });
}

fn bench_sliding_window(c: &mut Criterion) {
    let mut group = c.benchmark_group("sliding_window");
    for size in [1_000, 10_000, 100_000] {
        let data: Vec<u8> = (0..size).map(|i| (i % 256) as u8).collect();
        group.bench_with_input(BenchmarkId::from_parameter(size), &data, |b, data| {
            b.iter(|| {
                let mut sw = SlidingWindow::default_window();
                // Feed in 64-byte chunks
                for chunk in data.chunks(64) {
                    sw.append_bytes(chunk);
                }
                black_box(sw.final_hash())
            });
        });
    }
    group.finish();
}

fn bench_mersenne_mul(c: &mut Criterion) {
    c.bench_function("mersenne_mul", |b| {
        b.iter(|| {
            mersenne_mul(black_box(123456789), black_box(987654321), MERSENNE_61)
        });
    });
}

criterion_group!(
    benches,
    bench_mersenne_mul,
    bench_phi,
    bench_rope_concat_sequential,
    bench_repeat_vs_materialized,
    bench_split_rejoin,
    bench_substr_hash,
    bench_sliding_window,
);
criterion_main!(benches);
