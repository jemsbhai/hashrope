use criterion::{black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use hashrope::{
    Arena, Node, PolynomialHash, SlidingWindow,
};
use hashrope::polynomial_hash::{phi, mersenne_mul, MERSENNE_61};

fn bench_phi(c: &mut Criterion) {
    let alpha = 131u64;
    let x_d = mersenne_mul(alpha, alpha, MERSENNE_61);

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
                let mut arena = Arena::new();
                let mut node: Node = None;
                for i in 0..n {
                    let leaf = arena.from_bytes(&[i as u8]);
                    node = arena.concat(node, leaf);
                }
                black_box(arena.hash(node))
            });
        });
    }
    group.finish();
}

fn bench_repeat_vs_materialized(c: &mut Criterion) {
    let mut group = c.benchmark_group("repeat_vs_materialized");

    for q in [10, 100, 1_000, 10_000, 100_000, 1_000_000] {
        group.bench_with_input(BenchmarkId::new("repeat_node", q), &q, |b, &q| {
            b.iter_batched(
                || {
                    let mut arena = Arena::new();
                    let pattern = arena.from_bytes(b"abcdefgh");
                    (arena, pattern)
                },
                |(mut arena, pattern)| {
                    let rep = arena.repeat(pattern, black_box(q));
                    black_box(arena.hash(rep))
                },
                BatchSize::SmallInput,
            );
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

/// Helper: build a 10K-leaf rope for benchmarking split/substr operations.
fn build_10k_rope() -> (Arena, Node) {
    let mut arena = Arena::new();
    let mut node: Node = None;
    for i in 0..10_000u16 {
        let leaf = arena.from_bytes(&[i as u8]);
        node = arena.concat(node, leaf);
    }
    (arena, node)
}

fn bench_split_rejoin(c: &mut Criterion) {
    c.bench_function("split_rejoin_10k", |b| {
        b.iter_batched(
            build_10k_rope,
            |(mut arena, node)| {
                let (left, right) = arena.split(node, black_box(5000));
                let rejoined = arena.concat(left, right);
                black_box(arena.hash(rejoined))
            },
            BatchSize::SmallInput,
        );
    });
}

fn bench_substr_hash(c: &mut Criterion) {
    c.bench_function("substr_hash_10k", |b| {
        b.iter_batched(
            build_10k_rope,
            |(mut arena, node)| {
                black_box(arena.substr_hash(node, black_box(2500), black_box(5000)))
            },
            BatchSize::SmallInput,
        );
    });
}

fn bench_sliding_window(c: &mut Criterion) {
    let mut group = c.benchmark_group("sliding_window");
    for size in [1_000, 10_000, 100_000] {
        let data: Vec<u8> = (0..size).map(|i| (i % 256) as u8).collect();
        group.bench_with_input(BenchmarkId::from_parameter(size), &data, |b, data| {
            b.iter(|| {
                let mut sw = SlidingWindow::default_window();
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
