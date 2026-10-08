# hashrope

A BB[2/7] weight-balanced binary tree augmented with polynomial hash metadata at every node.

**Zero dependencies by default. `no_std + alloc` compatible. Rust 2021 edition.**

The optional `bigint` feature adds Python-compatible arbitrary-precision hashes
and compressed rope lengths. The additions described below are unreleased;
use this source checkout until a new crate release is published. See
[PARITY.md](PARITY.md) for the Python API matrix, exact fixture provenance,
downstream compatibility checks, and remaining adoption limitations.

## What it does

hashrope is a rope data structure where every node carries a polynomial hash over a Mersenne prime field. This enables:

- **O(log w) concat and split** with strict BB[2/7] weight-balance via recursive join
- **O(log q) repetition encoding** via `RepeatNode` — represent `s^q` without materializing
- **O(log w) substring hashing** without allocating new nodes
- **O(1) whole-sequence hash** at any node
- **Sliding windows** with eviction and opt-in arena reclamation

The default u64 API uses Mersenne modular reduction by bit folding. The optional bigint API uses BigUint modular arithmetic. Nodes live in contiguous arenas without reference counting or atomic operations.

## Install

```toml
[dependencies]
hashrope = "0.3"
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

For streaming applications, `SlidingWindow` maintains a bounded-size logical
window. The existing constructors preserve arena handles and retain evicted
allocations. The new `new_bounded` constructor also reclaims arena storage:

```rust
use hashrope::{SlidingWindow, MERSENNE_61};

let mut sw = SlidingWindow::new_bounded(512, 512, MERSENNE_61, 257);
sw.append_bytes(b"incoming data...");
// Oldest bytes are evicted when the window exceeds d_max + m_max
```

Reclamation copies the retained graph after eviction, preserving shared nodes
and compressed repeats. Its cost depends on the retained graph and leaf data;
it does not expand repeated bytes. Node IDs observed through `arena()` may
change on append in this opt-in mode. Calling `arena_mut()` permanently disables
reclamation so handles allocated by the caller remain valid. Keep the existing
constructors when stable arena handles are required.

## Arbitrary precision

Enable `bigint` on the source dependency:

```toml
[dependencies]
hashrope = { path = "path/to/hashrope/packages/rust", features = ["bigint"] }
```

```rust
# #[cfg(feature = "bigint")] {
use hashrope::bigint::{Arena, BigUint, MERSENNE_127};

let mut arena = Arena::with_hash(MERSENNE_127.into(), 131u8.into());
let pattern = arena.from_bytes("é🙂".as_bytes());
let count = (BigUint::from(1u8) << 80usize) + BigUint::from(17u8);
let repeated = arena.repeat(pattern, &count);
assert_eq!(arena.len(repeated), &count * BigUint::from(6u8));
let first_byte = arena.substr_hash(repeated, &0u8.into(), &1u8.into());
assert_eq!(first_byte, arena.hash_bytes(&[0xc3]));
# }
```

`bigint::PolynomialHash`, `bigint::Arena`, its rope wrappers, and
`bigint::SlidingWindow` use `BigUint` for hashes, moduli, lengths, positions,
and repetition counts. They support the 127-bit Python profile and wider
Mersenne moduli. This feature uses `num-bigint` and `num-traits` without their
`std` features. The original `u64` types, enum layouts, and default profile
remain unchanged. `BigUint` values can be displayed as decimal strings without
precision loss; applications own any external serialization format.

## Byte and compatibility contracts

All inputs and offsets are bytes. Text callers explicitly choose an encoding
(for example UTF-8 with `str::as_bytes()`); Hashrope does not normalize Unicode.
Splitting inside a UTF-8 code point preserves exact bytes, although either half
may not decode as text. Substring queries take **start and length**, not end.
Supply a valid in-bounds range. The new bigint API rejects malformed ranges;
the existing API retains its established range behavior.

Both APIs default to `p = 2^61 - 1`, base `131`, with the `byte + 1` polynomial
convention. The `u64` API checks representational limits before storing node
metadata and panics on overflow. Materializing a compressed rope requires
addressable memory. Hashes are probabilistic fingerprints, not cryptographic
digests or collision-free equality proofs.

Rust lazy arenas remain available. `Arena::hash` and `rope_hash` return
`LAZY_SENTINEL` until materialization; call `ensure_hash(root_id)` or query the
full range with `substr_hash` first. The new bigint API computes hashes eagerly.

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
