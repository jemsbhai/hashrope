# hashrope

A BB[2/7] weight-balanced binary tree augmented with polynomial hash metadata at every node.

**Zero dependencies. `no_std + alloc` compatible. Rust 2021 edition.**

## What it does

hashrope is a rope data structure where every node carries a polynomial hash over a Mersenne prime field. This enables:

- **O(log w) concat and split** with strict BB[2/7] weight-balance via recursive join
- **O(log q) repetition encoding** via `RepeatNode` — represent `s^q` without materializing
- **O(log w) substring hashing** without allocating new nodes
- **O(1) whole-sequence hash** at any node
- **Bounded-memory sliding window** with eviction

All arithmetic uses Mersenne prime (2⁶¹ − 1) modular reduction via bit-shift — no division. Nodes are stored in a contiguous arena (`Vec<NodeInner>`) for cache-friendly access. No reference counting, no atomic operations.

## Install

```toml
[dependencies]
hashrope = "0.2"
```

## Quick start

```rust
use hashrope::Arena;

let mut a = Arena::new();

// Build ropes from byte slices
let hello = a.from_bytes(b"hello ");
let world = a.from_bytes(b"world");

// Concat: O(log w)
let hw = a.concat(hello, world);
assert_eq!(a.hash(hw), a.hash_bytes(b"hello world"));
assert_eq!(a.len(hw), 11);

// Split at any byte position: O(log w)
let (left, right) = a.split(hw, 6);
assert_eq!(a.to_bytes(left), b"hello ");
assert_eq!(a.to_bytes(right), b"world");

// Rejoin preserves hash
let rejoined = a.concat(left, right);
assert_eq!(a.hash(rejoined), a.hash(hw));

// Repeat: O(log q) — no materialization
let pattern = a.from_bytes(b"ab");
let repeated = a.repeat(pattern, 1_000_000); // "ab" × 10⁶
assert_eq!(a.len(repeated), 2_000_000);

// Substring hash without allocation: O(log w)
let h = a.substr_hash(hw, 0, 5); // hash of "hello"
assert_eq!(h, a.hash_bytes(b"hello"));

// Validate BB[2/7] invariants (for debugging)
a.validate(hw);
```

## Node types

| Type | Description | Weight |
|------|-------------|--------|
| `Leaf(data)` | Raw byte slice | 1 |
| `Internal(left, right)` | Binary node with BB[2/7] balance | w(left) + w(right) |
| `Repeat(child, q)` | Virtual repetition, q ≥ 2 | w(child) × q |

## Arena allocation

All nodes live in a single `Arena`. Node references are `u32` indices. The arena is dropped as a unit — no per-node deallocation. This gives cache-friendly traversal and zero-cost structural sharing.

```rust
let mut a = Arena::new();
let node = a.from_bytes(b"data");
println!("Nodes allocated: {}", a.node_count());
```

## Sliding window

For streaming applications, `SlidingWindow` maintains a bounded-size rope with O(log w) append and eviction:

```rust
use hashrope::{SlidingWindow, MERSENNE_61};

let mut sw = SlidingWindow::new(512, 512, MERSENNE_61, 257);
sw.append_bytes(b"incoming data...");
// Oldest bytes are evicted when the window exceeds d_max + m_max
```

## Polynomial hash properties

The hash function is a polynomial rolling hash over GF(2⁶¹ − 1):

- **Homomorphic under concatenation**: `H(A‖B) = H(A)·x^|B| + H(B)`
- **Homomorphic under repetition**: `H(s^q) = H(s)·Φ(q, x^|s|)` where Φ is computed in O(log q)
- **Collision probability**: ≤ n/p per query for strings of length n, where p = 2⁶¹ − 1

## Benchmarks

```sh
cargo bench
```

Uses [Criterion.rs](https://github.com/bheisler/criterion.rs) with HTML reports.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).

## License

MIT
