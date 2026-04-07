# hashrope

A BB[2/7] weight-balanced binary tree augmented with polynomial hash metadata at every node.

**Zero dependencies. Rust 2021 edition.**

## What it does

hashrope is a rope data structure where every node carries a polynomial hash. This enables:

- **O(k log w) concat and split** with strict BB[2/7] rebalancing via recursive join
- **O(k log q) repetition encoding** via RepeatNode
- **O(k log w) substring hashing** without allocating new nodes
- **O(1) whole-sequence hash** at any point during streaming ingestion
- **Bounded-memory sliding window** with eviction

All arithmetic uses Mersenne prime modular reduction (bit-shift, no division). Nodes are immutable with structural sharing via `Arc`.

## Install

```toml
[dependencies]
hashrope = "0.1"
```

## Quick start

```rust
use hashrope::*;

let ph = PolynomialHash::default_hash();
let mut h = PolynomialHash::default_hash();

let a = rope_from_bytes(b"hello ", &ph);
let b = rope_from_bytes(b"world", &ph);
let ab = rope_concat(&a, &b, &mut h);

assert_eq!(rope_hash(&ab), ph.hash(b"hello world"));

// Split at any byte position
let (left, right) = rope_split(&ab, 6, &mut h);

// Substring hash without allocation
assert_eq!(rope_substr_hash(&ab, 0, 5, &mut h), ph.hash(b"hello"));
```

## Benchmarks

```sh
cargo bench
```

Uses [Criterion.rs](https://github.com/bheisler/criterion.rs) with HTML reports.

## Provenance

Translated from the [Python implementation](../python/) which was extracted from [UHC](https://github.com/jemsbhai/uhc).

## License

MIT
