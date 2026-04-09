# hashrope

A BB[2/7] weight-balanced binary tree augmented with polynomial hash metadata at every node.

[![PyPI version](https://img.shields.io/pypi/v/hashrope)](https://pypi.org/project/hashrope/)
[![Python 3.10+](https://img.shields.io/pypi/pyversions/hashrope)](https://pypi.org/project/hashrope/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://img.shields.io/badge/tests-141%20passing-brightgreen)](https://github.com/jemsbhai/hashrope)

**Zero dependencies. Pure Python. Python 3.10+.**

## What it does

hashrope is a rope data structure where every node carries a polynomial hash. This enables:

- **O(k log w) concat and split** with strict BB[2/7] rebalancing via Adams-style balanced join
- **O(k log q) repetition encoding** via RepeatNode -- represents "repeat this subtree q times" without materializing the data
- **O(k log w) substring hashing** without allocating new nodes
- **O(1) whole-sequence hash** at any point during streaming ingestion
- **Bounded-memory sliding window** that evicts old data while maintaining a running hash of everything seen

All arithmetic uses Mersenne prime modular reduction (bit-shift, no division).

## Install

```
pip install hashrope
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

# Substring hash without allocation
sub_h = rope_substr_hash(ab, 0, 5, h)
assert sub_h == h.hash(b"hello")
```

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

### Polynomial hash

| Function / Class | Description |
|---|---|
| `PolynomialHash(prime, base)` | Hash function over Z/pZ. Default: p = 2^61 - 1, base = 131 |
| `PolynomialHash.hash(data)` | H(data) with +1 byte offset to prevent zero collisions |
| `PolynomialHash.hash_concat(h_a, len_b, h_b)` | H(A \|\| B) from H(A), \|B\|, H(B) in O(1) |
| `PolynomialHash.hash_repeat(h_s, d, q)` | H(S^q) from H(S) and \|S\| in O(log q) |
| `PolynomialHash.hash_overlap(h_p, d, l, h_prefix)` | H(P^q \|\| P[0..r-1]) for overlapping references |
| `phi(q, alpha, p)` | Geometric accumulator: sum of alpha^i for i in [0, q) in O(log q) |
| `mersenne_mod(a, p)` | Bit-shift reduction mod Mersenne prime |
| `mersenne_mul(a, b, p)` | Modular multiply via mersenne_mod |
| `MERSENNE_61` | 2^61 - 1 (default prime) |
| `MERSENNE_127` | 2^127 - 1 (high-security prime) |

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
| `rope_substr_hash(node, start, length, h)` | O(k log w) | Substring hash without allocation |
| `rope_len(node)` | O(1) | Total byte length |
| `rope_hash(node)` | O(1) | Hash of represented string |
| `rope_from_bytes(data, h)` | O(n) | Create Leaf from bytes |
| `rope_to_bytes(node)` | O(n) | Reconstruct bytes (for testing) |
| `validate_rope(node)` | O(w) | Assert all 9 invariants hold |

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
