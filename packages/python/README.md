# hashrope

A BB[2/7] weight-balanced binary tree augmented with polynomial hash metadata at every node.

**Zero dependencies. Pure Python. Python 3.10+.**

## What it does

hashrope is a rope data structure where every node carries a cryptographically useful polynomial hash. This enables:

- **O(k * log w) concat and split** with automatic rebalancing
- **O(k * log q) repetition encoding** via RepeatNode (represents "repeat this subtree q times" without materializing the data)
- **O(k * log w) substring hashing** without allocating new nodes
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

### RepeatNode (O(log q) repetition)

```python
from hashrope import Leaf, rope_repeat, rope_hash

h = PolynomialHash()
pattern = Leaf(b"ab", h)
repeated = rope_repeat(pattern, 1000, h)  # "ab" x 1000

# Hash computed in O(log 1000), not O(2000)
assert rope_hash(repeated) == h.hash(b"ab" * 1000)
```

### SlidingWindow (bounded memory streaming)

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

## Mathematical foundation

The data structure and its correctness proofs come from the compressed-domain hashing framework (24 theorems, 14 lemmas, 5 corollaries). Key results used by hashrope:

| Result | Statement |
|--------|-----------|
| Theorem 1 | H(A \|\| B) = H(A) * x^{\|B\|} + H(B) |
| Theorem 2 | H(S^q) = H(S) * Phi(q, x^d) |
| Theorem 3 | Phi(q, alpha) computable in O(log q) without modular inverse |
| Theorem 7 | Split in O(k * log w) |
| Theorem 9 | Allocation-free SubstrHash in O(k * log w) |

## License

MIT
