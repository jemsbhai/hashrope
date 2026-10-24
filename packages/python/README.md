# hashrope

A BB[2/7] weight-balanced binary tree augmented with polynomial hash metadata at every node.

[![PyPI version](https://img.shields.io/pypi/v/hashrope)](https://pypi.org/project/hashrope/)
[![Python 3.10+](https://img.shields.io/pypi/pyversions/hashrope)](https://pypi.org/project/hashrope/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

**Zero dependencies. Pure Python. Python 3.10+.**

## What it does

hashrope is a rope data structure where every node carries a polynomial hash. This enables:

- **O(k log w) concat and split** with strict BB[2/7] rebalancing via Adams-style balanced join
- **O(k log q) repetition encoding** via RepeatNode -- represents "repeat this subtree q times" without materializing the data
- **O(k log w) substring hashing** without allocating new nodes
- **O(1) whole-sequence hash** at any point during streaming ingestion
- **Bounded-memory sliding window** that evicts old data while maintaining a running hash of everything seen

All arithmetic uses Mersenne prime modular reduction (bit-shift, no division).

Fingerprints are noncryptographic candidate-match keys. Matching fingerprints
do not prove equal content. Use `rope_equal` or `rope_lcp` for exact decisions;
these compare bytes while taking advantage of shared structure and repetition.

## Install

Version 0.3.3 adds exact comparison and hash-profile checks while preserving
existing fingerprints and defaults. Version numbers are synchronized across
the core packages; language-specific APIs and performance are documented separately.

```
pip install "hashrope>=0.3.3"
```

## Quick start

```python
from hashrope import PolynomialHash, Leaf, rope_concat, rope_split, rope_substr_hash, rope_hash

h = PolynomialHash()

# Build a rope
a = Leaf(b"hello ", h)
b = Leaf(b"world", h)
ab = rope_concat(a, b, h)

# Whole-string hash matches direct computation
assert rope_hash(ab) == h.hash(b"hello world")

# Split at any byte position
left, right = rope_split(ab, 6, h)

# Substring hash without allocating rope nodes
sub_h = rope_substr_hash(ab, 0, 5, h)
assert sub_h == h.hash(b"hello")
```

### Exact comparison and hash profiles

```python
from hashrope import Leaf, PolynomialHash, rope_equal, rope_hash, rope_lcp

h = PolynomialHash()  # Existing defaults remain p=2**61-1, base=131.
a = Leaf(bytes([0, 131, 42]), h)
b = Leaf(bytes([1, 0, 42]), h)
assert rope_hash(a) == rope_hash(b)  # A deliberate legacy collision.
assert not rope_equal(a, b)
assert rope_lcp(a, b) == 0

# Exact comparisons also work across profiles and tree shapes.
c = Leaf(a.data, PolynomialHash(base=257))
assert rope_equal(a, c)
assert rope_lcp(a, c) == len(a.data)
```

`rope_lcp` returns a byte count; for fixed-width four-byte token encodings,
`rope_lcp(a, b) // 4` gives the exact token LCP. Neither exact API depends on
fingerprints. Empty ropes (`None`) are supported. Node `==` retains its existing
structural meaning; use `rope_equal` for content equality across tree shapes.

For slices, `rope_range_equal(left, left_start, right, right_start, length)` and
`rope_range_lcp(left, left_start, right, right_start, length)` compare byte ranges
directly, without splitting or recomputing hashes. The latter returns a count
relative to the two starts. Both reject negative or out-of-bounds ranges with
`ValueError`; an empty range may start at the rope's end. They also work across
hash profiles.

The comparison walks virtual ranges without splitting or flattening ropes.
Identical immutable subtrees are skipped, and a per-call cache reuses previously
verified non-leaf ranges in shared DAGs. Two periodic spans can be certified by
comparing at most `p + q - gcd(p, q)` bytes, where `p` and `q` are their periods
([Fine–Wilf periodicity theorem](https://doi.org/10.1090/S0002-9939-1965-0174934-9)).
This also handles different period lengths and shifted phases. Unshared inputs
can still require reading the entire common prefix, plus traversal, cache and
integer-arithmetic costs: exact LCP has no universal logarithmic bound. Temporary storage includes the traversal
stack, verified-range cache, and bounded leaf-comparison slices (at most
192 KiB of temporary byte payload). Unequal chunks use native exact byte
comparisons to locate late differences; a short scalar prefix check keeps
early differences cheap. The sum of byte work over those shrinking probes is
linear in the chunk size, and the overall exact-comparison bound is unchanged.

Each hasher exposes an immutable `HashProfile(prime, base)` through `.profile`;
nodes retain that profile. Separate hashers with equal parameters are compatible.
`Internal`, `RepeatNode`, `rope_concat`, `rope_split`, `rope_repeat`, and
`rope_substr_hash` raise `ValueError` on incompatible profiles, including calls
that otherwise return an unchanged node or cached hash. Replacing a sliding
window's hasher with incompatible parameters is rejected before hashing or
mutating the window. Create a new hasher and rebuild content to change profiles.
Hasher subclasses must preserve the fingerprint algorithm when overriding
methods for instrumentation.

All existing fingerprint values, defaults, and valid same-profile operation
results are unchanged. No stable node-pickle format is promised: rebuild nodes
serialized by older versions, which lack profile provenance. `validate_rope`
checks structure and profile consistency; it does not recompute fingerprints.
Initial safety measurements are in
[`benchmarks/comparison_safety`](../../benchmarks/comparison_safety/README.md).
The [exact-path optimization rerun](../../benchmarks/exact_optimization/README.md)
profiles leaf mismatch localization and compares archived cores on the same
consumer workloads. The optimization improves late mismatches; it does not
speed up every offset or remove the traversal cost of independent equal trees.
The report includes the narrow regression at the 64-byte probe boundary.

### RepeatNode

Encode repetitions in O(log q) time and O(1) space:

```python
from hashrope import PolynomialHash, Leaf, rope_repeat, rope_hash

h = PolynomialHash()
pattern = Leaf(b"ab", h)
repeated = rope_repeat(pattern, 1000, h)  # "ab" x 1000

# Hash computed in O(log 1000), not O(2000)
assert rope_hash(repeated) == h.hash(b"ab" * 1000)
```

### SlidingWindow

Bounded-memory streaming hash with copy-from-window support:

```python
from hashrope import PolynomialHash, SlidingWindow

h = PolynomialHash()
sw = SlidingWindow(d_max=1024, m_max=256)

# Stream data in
sw.append_bytes(b"hello ")
sw.append_bytes(b"world")

# Hash of everything seen so far
assert sw.current_hash() == h.hash(b"hello world")

# Copy-from-window (like LZ77 back-references)
sw.append_copy(offset=5, length=5)  # copies "world" again
assert sw.final_hash() == h.hash(b"hello worldworld")
```

## API reference

### Byte-index contract

Every index and length in this API (`rope_split` `pos`, `rope_substr_hash`
`start`/`length`, `SlidingWindow` offsets) is a byte offset into the stored
byte string. The library stores, splits, and hashes bytes; it has no notion
of characters or encodings. When the bytes are UTF-8-encoded text, an
arbitrary byte offset may fall inside a multibyte character's encoding. That
split is still exact at the byte level: concatenating the two halves (or
`rope_to_bytes` of each) reproduces the original bytes, and every hash is
computed over bytes, so nothing is lost or corrupted. Only the two halves,
viewed on their own, may not each be valid UTF-8. Callers who need splits on
character boundaries must supply byte offsets that land on encoded character
boundaries (for example `len(text[:k].encode("utf-8"))` for a split after
the first `k` characters).

### Polynomial hash

| Function / Class | Description |
|---|---|
| `PolynomialHash(prime, base)` | Hash function over Z/pZ. Default: p = 2^61 - 1, base = 131 |
| `HashProfile(prime, base)` | Immutable fingerprint parameters; available as `hasher.profile` and `node.profile` |
| `PolynomialHash.hash(data)` | H(data) with +1 byte offset to prevent zero collisions |
| `PolynomialHash.hash_concat(h_a, len_b, h_b)` | H(A \|\| B) from H(A), \|B\|, H(B) in O(1) |
| `PolynomialHash.hash_repeat(h_s, d, q)` | H(S^q) from H(S) and \|S\| in O(log q) |
| `PolynomialHash.hash_overlap(h_p, d, l, h_prefix)` | H(P^q \|\| P[0..r-1]) for overlapping references |
| `phi(q, alpha, p)` | Geometric accumulator: sum of alpha^i for i in [0, q) in O(log q) |
| `mersenne_mod(a, p)` | Bit-shift reduction mod Mersenne prime |
| `mersenne_mul(a, b, p)` | Modular multiply via mersenne_mod |
| `MERSENNE_61` | 2^61 - 1 (default prime) |
| `MERSENNE_127` | 2^127 - 1 (wider noncryptographic fingerprint) |

### Rope nodes

| Class | Description |
|---|---|
| `Leaf(data, h)` | Immutable leaf storing a byte sequence. Weight = 1 |
| `Internal(left, right, h)` | Internal node. Hash = H(left) \* x^len(right) + H(right) |
| `RepeatNode(child, reps, h)` | Represents child^reps. Hash via geometric accumulator |

### Rope operations

| Function | Time | Description |
|---|---|---|
| `rope_concat(left, right, h)` | O(k \|h_L - h_R\|) | Concatenate with Adams-style balanced join and BB[2/7] rebalancing |
| `rope_split(node, pos, h)` | O(k log w) | Split at byte position |
| `rope_repeat(node, q, h)` | O(k log q) | Create RepeatNode |
| `rope_substr_hash(node, start, length, h)` | O(k log w) | Substring hash without allocating rope nodes |
| `rope_len(node)` | O(1) | Total byte length |
| `rope_hash(node)` | O(1) | Hash of represented string |
| `rope_equal(left, right)` | Content-dependent; may scan unshared content | Exact byte equality, including across profiles |
| `rope_lcp(left, right)` | Content-dependent; may scan unshared content | Exact common-prefix length in bytes |
| `rope_range_equal(left, start_a, right, start_b, length)` | Content-dependent | Checked exact range equality without splitting |
| `rope_range_lcp(left, start_a, right, start_b, length)` | Content-dependent | Checked exact range LCP, relative to the starts |
| `rope_from_bytes(data, h)` | O(n) | Create Leaf from bytes |
| `rope_to_bytes(node)` | O(n) | Reconstruct bytes (for testing) |
| `validate_rope(node, h=None)` | O(w) | Check structure and profile consistency, without rehashing |

### SlidingWindow

| Method | Description |
|---|---|
| `SlidingWindow(d_max, m_max, prime, base)` | Create window. W = d_max + m_max |
| `.append_bytes(data)` | Append literal bytes |
| `.append_copy(offset, length)` | Copy from window (handles overlapping) |
| `.current_hash()` | H(everything seen) in O(1) |
| `.final_hash()` | Same as current_hash (alias) |
| `.window_len` | Current window size in bytes |
| `.pos` | Total bytes ingested |

## Mathematical foundation

The data structure and its correctness proofs come from the compressed-domain hashing framework (24 theorems, 14 lemmas, 5 corollaries). Key results used by hashrope:

| Result | Statement |
|---|---|
| Theorem 1 (Concatenation) | H(A \|\| B) = H(A) * x^{\|B\|} + H(B) |
| Theorem 2 (Repetition) | H(S^q) = H(S) * Phi(q, x^d) |
| Theorem 3 (Phi computation) | Phi(q, alpha) computable in O(log q) without modular inverse |
| Theorem 6 (Join) | Adams-style balanced join with BB[2/7] rebalancing |
| Theorem 7 (Split) | Split in O(k log w) with RepeatNode splittability |
| Theorem 8 (Repeat) | RepeatNode in O(k log q), O(1) space |
| Theorem 9 (SubstrHash) | Allocation-free substring hashing in O(k log w) |
| Theorem 11 (Sliding window) | Bounded-memory streaming with eviction |

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history and details on the v0.2.1 balance bugfix.

## Provenance

Extracted from [UHC](https://github.com/jemsbhai/uhc) (Unified Hash-Compression Engine). The rope and polynomial hash modules are verbatim copies with import paths updated. The sliding window is adapted from UHC's SlidingRopeState with the LZ77 dependency removed.

## License

MIT
