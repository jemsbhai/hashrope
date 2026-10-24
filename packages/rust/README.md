# hashrope

A BB[2/7] weight-balanced binary tree augmented with polynomial hash metadata at every node.

**Zero dependencies by default. `no_std + alloc` compatible. Rust 2021 edition.**

The optional `bigint` feature adds Python-compatible arbitrary-precision hashes
and compressed rope lengths in Rust 0.3.2. The compatible 0.3 release line
retains existing consumers' dependency requirements. See
[PARITY.md](PARITY.md) for the Python API matrix, exact fixture provenance,
downstream compatibility checks, and remaining adoption limitations.

## What it does

hashrope is a rope data structure where every node carries a polynomial hash over a Mersenne prime field. This enables:

- **O(log w) concat and split** with strict BB[2/7] weight-balance via recursive join
- **O(log q) repetition encoding** via `RepeatNode` — represent `s^q` without materializing
- **O(log w) substring hashing** without allocating new nodes
- **O(1) whole-sequence hash** at any node
- **Exact equality and longest common prefix** with shared-subtree and periodicity shortcuts
- **Sliding windows** with eviction and opt-in arena reclamation

The logarithmic tree traversal bounds exclude uncached exponentiation and repeat arithmetic. The u64 `from_bytes` constructor bounds leaves to 512 bytes. Substring hashing allocates no rope nodes but may grow the exponent cache. Lazy arenas require hash materialization.

The default u64 API uses Mersenne modular reduction by bit folding. The optional bigint API uses BigUint modular arithmetic. Nodes live in contiguous arenas without reference counting or atomic operations.

## Install

Version 0.3.3 adds exact whole-content and range comparison and hash-profile
checks while preserving the APIs and default fingerprints from 0.3.2.
The optional `bigint` backend provides the same exact-comparison contract
with arbitrary-precision lengths. See the release validation for measured
performance and remaining limitations.

```toml
[dependencies]
hashrope = "0.3.3"
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

All nodes live in a single `Arena`. Node references are `u32` indices. They identify nodes only within their owning arena: passing another arena's numeric ID is invalid even if that ID happens to exist locally. The arena is dropped as a unit — no per-node deallocation. This gives cache-friendly traversal and zero-cost structural sharing.

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

Enable `bigint` on the dependency:

```toml
[dependencies]
hashrope = { version = "0.3.2", features = ["bigint"] }
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
- **Noncryptographic fingerprints**: useful for candidate matching, never proof of content identity

The default base remains **131** for compatibility. `[0,131]` and `[1,0]` both hash to **263**, before modular reduction. A larger modulus or an additional modulus with the same base does not separate this pair. Existing fingerprint values and default parameters are unchanged.

For the default profile, leaf hashing evaluates eight bytes per modular
reduction using `H(prefix) * 131^8 + H(block)`. This is the same polynomial:
the block fits in `u64`, and the combined accumulator is below `2^118`.
The tail uses the original byte recurrence. Other accepted profiles retain
the per-byte arithmetic from published 0.3.2, including its small-modulus
reduction fix. Splitting copies only the two output slices; lazy materialization
and partial substring hashing borrow the original leaf bytes.

A polynomial collision bound requires a probability model. For two distinct equal-length byte strings of length `n`, fixed before a base is chosen uniformly from `{2,...,p-1}`, the usual root bound is at most `(n-1)/(p-2)` for `p > 256`. The fixed public default base has no such random-base guarantee, particularly for inputs chosen after the parameters are known.

## Exact comparison

```rust
use hashrope::Arena;

let mut arena = Arena::new();
let a = arena.from_bytes(&[0, 131, 42]);
let b = arena.from_bytes(&[1, 0, 42]);
assert_eq!(arena.hash(a), arena.hash(b)); // A candidate match only.
assert!(!arena.equal(a, b));
assert_eq!(arena.lcp(a, b), 0);
```

`Arena::equal(a,b)` compares represented bytes regardless of tree shape. `Arena::lcp(a,b)` returns the exact longest common prefix length **in bytes**. Free functions `rope_equal(&arena,a,b)` and `rope_lcp(&arena,a,b)` provide the same operations. Across arenas, use `left_arena.equal_across(a,&right_arena,b)` or `lcp_across`; each handle belongs to the corresponding arena. Different prime/base profiles and eager/lazy modes are supported.

For windows, `equal_range(a,a_start,b,b_start,length)` and `lcp_range` compare equal-length ranges without splitting or hashing. Their `_across` variants accept the other arena before `b`. Free wrappers are `rope_range_equal` and `rope_range_lcp`. Bounds are checked without overflowing: invalid ranges panic, and empty ranges at a rope's end are valid.

The `bigint::Arena` provides the same method names, with `&BigUint` range
arguments and a `BigUint` LCP result. It uses an explicit task stack and a cache
of verified ranges, so deep compressed structures do not require recursive
comparison calls. Its offsets and periodic witnesses remain arbitrary precision;
only proven leaf ranges convert to machine-sized indices. Big-integer and cache
operations allocate. Neither backend flattens content or uses hashes as an exact
comparison decision.

The following traversal details apply to the default fixed-width `Arena`.
These operations read bytes and immutable structure, without materializing lazy hashes, splitting ropes, or allocating flattened content. Identical shared ranges are skipped. Ordinary trees use borrowed leaf slices without heap allocations. Long cross-arena ranges (at least 8 KiB) use two iterative traversal stacks to avoid revisiting branches when the tree boundaries differ. Their 32-frame capacity uses 1,600 bytes of frontier storage on 64-bit targets, excluding other local variables. Smaller and same-arena ranges retain a recursive traversal that enters known children directly. An unequal leaf is localized using safe 8-byte comparisons.

Normalization and frontier node examinations consume a visit budget of eight times the combined number of stored arena nodes. Exhausting that budget, encountering a repeat, or filling a frontier stack transfers the unverified suffix to the cached engine; the verified prefix is not rescanned. Clearly compressed roots go directly to that engine. This bounds uncached traversal before fallback without adding metadata to nodes. Node layout and the 512-byte leaf construction policy are unchanged.

In the cached engine, repeated ranges with periods `p` and `q` are proved equal by checking at most `p+q-gcd(p,q)` bytes before skipping further copies (the Fine-Wilf periodicity bound); that checking can itself skip shared or repeated structure. At query entry, chains of `Repeat` nodes and concatenations of the identical child establish smaller periods structurally. The witness checks the original ranges, including their phase offsets. Reducing these chains once per query avoids repeatedly walking them on late mismatches. A per-call cache records ranges already proved equal, avoiding repeated traversal of other shared DAGs.

Exact comparison is **not a logarithmic worst-case operation**. Independently stored equal content can require reading the entire common prefix plus tree traversal. If fallback is needed, cache storage grows with the number of ranges proved equal, and its ordered lookups add logarithmic overhead. Recursive stack depth follows the compared structure. Sharing and repetition often make comparisons much cheaper, but do not provide a universal bound in compressed size.

## Hash profile consistency

Both arena backends fix the prime/base profile at construction. `hasher_mut()` remains available for operations such as warming the exponent cache. Replacing the hasher with the same parameters is supported. Replacing it with different parameters causes subsequent arena hashing, construction, edits and validation to panic before using inconsistent metadata; create a new arena to change profiles. Byte access, length access and exact comparisons remain available because they do not depend on fingerprint parameters.

`validate` checks tree structure and unchanged arena parameters; it does not recompute every cached fingerprint. Lazy `hash` still returns `LAZY_SENTINEL` before materialization. Exact comparison handles lazy ropes directly.

## Benchmarks

### Release 0.3.3 integration

The [release integration check](../../releases/0.3.3/performance/rust_integration/README.md)
compares the newly integrated default Rust backend with the completed round-four
optimized snapshot, which already had exact APIs. Three alternating process
pairs retain the unchanged 42-cell microbenchmark and 24-cell edit/comparison
workload. This is a targeted integration diagnostic, not the full historical
matrix: outer Python interleaving is omitted and the binaries are rebuilt.

Results are mixed. Independent-equal 1 MiB LCP changes from 57.420 to 62.475 µs
(median of process medians; paired geometric-mean ratio 1.105, slower in 3/3
pairs). The 2 TiB self-concatenated DAG LCP changes from 0.252 to 0.340 µs
(paired ratio 1.308, slower in 3/3 pairs), while DAG equality is approximately
unchanged. Independent 1 MiB equality and late-mismatch LCP have paired ratios
of 0.983 and 0.983. Exact edit/comparison paired changes range from +1.0% to
+8.0%; 100K independent targets are slower in all three pairs (+4.5%), while
the other exact edit cells have mixed directions. Layout remains 48 bytes per
`NodeInner`, 2,136 inline bytes per `Arena`, and 255 nodes for the 64 KiB fixture.

All samples and cells are retained. The three-pair check does not isolate a
cause or justify dismissing recurring slowdowns as noise. It establishes
neither universal efficiency nor a no-regression guarantee. Bigint exact APIs
have dedicated compressed-graph correctness checks but are not timed by this
default-feature diagnostic.

### Historical optimization studies

The [round-four regression investigation](../../benchmarks/exact_optimization_round4/README.md)
reexamines the third round's 12–17% slower comparison-only cells. No further
production change was selected; the equality penalty recurs in the original
short-sample protocol, while the separate diagnostics do not show a consistent
penalty of that size. Comparator source and the normalized emitted
instructions of seven inspected routines are identical across the third-round
optimization; that alone does not establish identical linked layout or runtime.

The separate diagnostic uses 36 fresh processes with longer paired warmups,
balanced execution order, swapped implementation slots, identical-code controls,
and separate pinned/unpinned arms. Neither primary penalty reproduces across
these conditions. Twelve additional paired runs use the exact retained original
executables and the entire unchanged Rust benchmark sequence. Their geometric
mean third-round/pre-third-round ratios are 0.9959 for independent 1 MiB equality and
0.9914 for late-mismatch LCP; exploratory per-cell 95% log-t intervals are
[0.9458, 1.0486] and [0.9326, 1.0540], respectively. These intervals are not
adjusted for the many cells and do not bound every deployment. Individual
process pairs vary substantially. The final full-protocol rerun below again
shows slower equality, so these diagnostic estimates neither establish universal
equivalence nor dismiss that signal as noise. No causal mechanism has been
isolated to justify a production change.

```sh
cargo bench
```

Uses [Criterion.rs](https://github.com/bheisler/criterion.rs) with HTML reports.

The [full-matrix replication report](../../benchmarks/exact_optimization_round4/README.md)
reruns the same workload against Ropey 1.6.1, comparing the pre-third-round core
with the unchanged third-round core. The earlier optimization grouped default-profile
hash evaluation and removed redundant leaf copies. Each batch performs
50 edits followed by comparison with prebuilt targets, half equal and half
different at the final byte. Values are **µs per cycle**, mean ± sample SD
across nine fresh-process medians. Setup is excluded; every boolean is checked.

| Bytes / target history | Pre-third-round core | Unchanged third-round core | Ropey, current cohort |
|---|---:|---:|---:|
| 10,000 / independent | 7.43 ± 0.64 | 4.47 ± 0.24 | 1.17 ± 0.19 |
| 10,000 / shared | 6.64 ± 0.53 | 3.70 ± 0.15 | 0.91 ± 0.09 |
| 100,000 / independent | 23.63 ± 1.08 | 19.18 ± 1.90 | 11.18 ± 1.60 |
| 100,000 / shared | 12.05 ± 0.73 | 7.39 ± 0.58 | 5.36 ± 0.60 |
| 1,000,000 / independent | 153.36 ± 14.66 | 144.17 ± 7.33 | 106.70 ± 7.96 |
| 1,000,000 / shared | 14.68 ± 0.85 | 9.17 ± 0.52 | 34.17 ± 2.69 |

Both 10 KB cases improve in 9/9 paired runs, by 1.66× and 1.79× respectively.
Ropey still wins both regimes at 10 KB and 100 KB, and every independent-target
size. The 1 MB independent gain is 6.0% (7/9 paired improvements),
with a remaining 1.35× Hashrope/Ropey time ratio. Shared history gives Hashrope
a 3.73× advantage at 1 MB in this cohort. These scoped ratios do not establish
universal superiority or statistical significance.

The retained third-round [decomposition diagnostic](../../benchmarks/exact_optimization_round3/edit_profile.md)
measures 10 KB independent edit-only time at 5.065 ± 0.236 → 1.983 ± 0.068 µs,
versus Ropey's 0.437 ± 0.009 µs (three process medians per core). Each edit
still hashes about 709 leaf bytes and creates 40.41 arena nodes on average;
fresh allocation requests fall from 9.88 to 6.92. The original warmup leaves
zero power-cache misses. Persistent path rebuilding and cached metadata remain
costs, while long independent comparisons still scan content. Separate phase
times are not additive; lazy-mode differences do not isolate a precise hashing
percentage. An extra 1 KB diagnostic still favors Ropey.

In that third-round diagnostic, the 512-byte hash primitive improves 12.84×,
but the one-byte probe is about
0.423 ns slower. Both fingerprint-candidate arms use Hashrope's hash kernel,
including the Ropey rehash control; their gains are distinct from exact equality.

The fourth-round structural microbenchmarks retain three alternating rounds of
seven samples, reported as median [minimum, maximum] of round medians. Building
64 KiB changes from 288.730 [288.275, 309.450] to 47.020 [42.640, 51.905] µs;
splitting from 2.553 [2.386, 2.576] to 0.460 [0.422, 0.507] µs; and warm
substring hashing from 4.491 [4.362, 4.590] to 0.492 [0.480, 0.515] µs.
The 2 TiB self-concatenation DAG LCP is 0.280 [0.251, 0.280] →
0.280 [0.256, 0.283] µs. These engines are unchanged in this round.

Independent 1 MiB equality increases from 54.926 [54.543, 58.861] to
65.895 [55.910, 65.990] µs (+20.0%), and late-mismatch LCP from
57.184 [54.632, 57.449] to 60.059 [57.426, 65.809] µs (+5.0%).
Equality is slower in all three pairs, by approximately 20.8%, 12.1% and 1.8%.
Its direction recurs from the original third-round result, while the LCP
penalty is smaller. Secondary 64 KiB late-mismatch cells also worsen: LCP
3.082→3.695 µs (+19.9%) and equality 3.083→3.979 µs (+29.1%), with
nonoverlapping three-round ranges; these ranges are not confidence intervals.
The longer controlled diagnostic and original-executable
replication above do not show the same consistent effect. The full rerun
interleaves Python work and uses rebuilt executables with different hashes;
these differences were not isolated, and no causal explanation is established.
This leaves a context-sensitive equality signal unresolved. The unchanged source
does not establish unchanged runtime, and the short protocol does not establish
a general comparator penalty. All cells and raw samples are retained in the
report. These independent microbench roots share one
arena and use the recursive comparator; the independent-target churn arm
uses a separate target arena and long comparisons use the frontier path.
The data and sizes also differ, so these microbench timings cannot explain the churn timings by
subtraction. The construction probe includes arena setup and destruction.
Rust layout measurements remain 48 bytes per
`NodeInner`, 2,136 inline bytes per `Arena`, and 255 nodes for the 64 KiB
fixture; these are layout measurements, not total memory usage.

The [original third-round study](../../benchmarks/exact_optimization_round3/README.md),
[previous traversal study](../../benchmarks/exact_optimization_round2/README.md)
and [first optimization study](../../benchmarks/exact_optimization/README.md)
remain preserved separately. Fingerprints, default parameters, node layout,
the leaf-size policy and exact-comparison code are unchanged in this round.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).

## License

MIT
