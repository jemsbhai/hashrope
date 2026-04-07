# hashrope — Experiment Findings

**Machine**: Windows, 64 GB RAM, NVIDIA RTX 4090 (CPU-only benchmarks)
**Rust toolchain**: rustc 1.94.0 (4a4ef493e 2026-03-02)
**Compiler flags**: `--release` (opt-level 3)
**Date started**: 2026-04-07

---

## E7: BB[2/7] Tree Height Verification

**Status**: PASS — Claim C7 confirmed.

Theoretical bound: h ≤ log(w) / log(7/5) ≈ 2.06 · log₂(w)

### Sequential insertion (one byte per leaf, weight = n)

| n | weight | height | bound | pass |
|---|--------|--------|-------|------|
| 100 | 100 | 8 | 14 | OK |
| 1,000 | 1,000 | 11 | 21 | OK |
| 10,000 | 10,000 | 15 | 28 | OK |
| 100,000 | 100,000 | 20 | 35 | OK |

Actual heights ≈ 1.1–1.2 · log₂(w), well under the 2.06 · log₂(w) ceiling.

### Chunked insertion (varied chunk sizes, multi-byte leaves)

| n (bytes) | weight | height | bound | pass |
|-----------|--------|--------|-------|------|
| 100 | 6 | 3 | 6 | OK |
| 1,000 | 36 | 6 | 11 | OK |
| 10,000 | 327 | 11 | 18 | OK |
| 100,000 | 3,257 | 15 | 25 | OK |

Multi-byte leaves reduce weight dramatically, giving even lower heights.

---

## E4: Phi Accumulator Scaling

**Status**: PASS — Claim C4 confirmed (O(log q) scaling).

100,000 iterations per measurement, `--release` mode.

| q | time (ns) | log₂(q) | bits(q) |
|---|-----------|---------|---------|
| 10 | 22.1 | 3.32 | 4 |
| 100 | 43.9 | 6.64 | 7 |
| 1,000 | 77.5 | 9.97 | 10 |
| 10,000 | 91.4 | 13.29 | 14 |
| 100,000 | 116.6 | 16.61 | 17 |
| 1,000,000 | 137.3 | 19.93 | 20 |
| 10,000,000 | 154.7 | 23.25 | 24 |

Time scales linearly with bit-length of q (i.e., log₂(q)), confirming O(log q).
From q=10 to q=10⁷ (6 orders of magnitude), time increases only ~7×,
consistent with the ~6× increase in bit-length (4 → 24 bits).

---

## E3: RepeatNode Scaling

**Status**: PASS — Claim C3 confirmed (O(log q) vs O(q)).

Base rope: 100 bytes. Rope repeat timed over 10,000 iterations.
Naive baseline: materialize q×100 bytes into Vec, hash from scratch.
Hash correctness verified at every q.

| q | repeat (ns) | naive (ns) | speedup |
|---|-------------|------------|---------|
| 10 | 60.9 | 4,584 | 75× |
| 100 | 80.8 | 46,453 | 575× |
| 1,000 | 112.4 | 431,307 | 3,838× |
| 10,000 | 121.0 | 4,293,374 | 35,468× |
| 100,000 | 148.8 | 43,895,740 | 294,919× |
| 1,000,000 | 169.8 | 442,455,740 | 2,606,053× |
| 10,000,000 | 196.2 | 4,450,499,350 | 22,685,795× |
| 100,000,000 | 249.4 | 44,415,002,200 | 178,101,701× |

RepeatNode: 61 ns → 249 ns across 7 orders of magnitude (O(log q)).
Naive: 4.6 µs → 44.4 s, scaling linearly (O(q·|B|)).
At q=10⁸, RepeatNode is **178 million× faster** while producing the identical hash.

---
