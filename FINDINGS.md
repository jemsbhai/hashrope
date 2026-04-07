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

## E1: Sequential Hashing Throughput

**Status**: COMPLETE — Claim C1 context established.

Synthetic random data (LCG seeded), `--release` mode.
Iterations scaled inversely with size for stable timings.
Each hashrope iteration uses a fresh `PolynomialHash` (no pre-warmed cache).

| Size | hashrope (MB/s) | BLAKE3 (MB/s) | SHA-256 (MB/s) | xxh3 (MB/s) |
|------|-----------------|---------------|----------------|-------------|
| 1 KB | 193.2 | 1,433.7 | 2,312.8 | 32,643.2 |
| 10 KB | 215.0 | 4,503.5 | 2,476.3 | 37,204.0 |
| 100 KB | 216.8 | 5,686.8 | 2,516.8 | 36,315.1 |
| 1 MB | 217.4 | 5,749.5 | 2,491.4 | 36,147.8 |
| 10 MB | 214.7 | 5,523.5 | 2,463.3 | 38,356.6 |
| 100 MB | 213.5 | 4,653.1 | 2,402.6 | 11,437.7 |

**Analysis**:
- hashrope steady-state: ~215 MB/s (byte-at-a-time polynomial evaluation)
- BLAKE3: ~5,500 MB/s (~26× faster, SIMD-optimized)
- SHA-256: ~2,450 MB/s (~11× faster)
- xxh3: ~36,000 MB/s (~170× faster, drops to ~11,400 at 100 MB — likely L3 cache pressure)

hashrope's raw sequential throughput is not competitive with SIMD-optimized hashes.
This is expected and by design: the polynomial hash's value is structural composability
(O(log n) concat/split, O(log q) repeat) — capabilities E2, E3, and E5 demonstrate.
The ~215 MB/s throughput is the cost of byte-at-a-time Mersenne-prime arithmetic,
a known optimization target (SIMD unrolling planned for Phase 2).

---

## E2: Rope Concat/Split vs Naive Rehash

**Status**: PASS — Claim C2 confirmed.

Ropes built by sequential single-byte concat (worst case). `--release` mode.
Hash correctness verified at every size.

### Concat (rejoin two halves)

| n | join (ns) | rehash (ns) | speedup |
|---|-----------|-------------|---------|
| 100 | 55.0 | 482.0 | 9× |
| 1,000 | 52.2 | 4,599.3 | 88× |
| 10,000 | 62.0 | 47,625.0 | 768× |
| 100,000 | 80.0 | 432,410.0 | 5,405× |

Join time is nearly constant (~55–80 ns), confirming O(log n).
Rehash scales linearly with n. At n=100K, join is **5,405× faster**.

### Split at midpoint

| n | split (ns) | rehash (ns) | speedup |
|---|------------|-------------|---------|
| 100 | 495.8 | 424.1 | 0.9× |
| 1,000 | 976.0 | 4,566.0 | 5× |
| 10,000 | 1,010.0 | 43,783.0 | 43× |
| 100,000 | 1,260.0 | 452,560.0 | 359× |

Split has higher constant overhead than join due to node allocation along the
split path. At n=100, split is actually *slower* than rehashing — the tree
traversal + node creation cost dominates for tiny inputs. The crossover occurs
around n ≈ 200–500. By n=100K, split is **359× faster** and the asymptotic
advantage (O(log n) vs O(n)) is clear.

---

## E5: Sliding Window Throughput

**Status**: COMPLETE — Claim C5 result: not competitive on raw throughput.

Synthetic random data (LCG seeded), `--release` mode.
hashrope: `SlidingWindow` with `d_max=W, m_max=1`, one byte at a time.
Rabin-Karp: same Mersenne-prime base (131) and modulus (2⁶¹−1), ring buffer.

| Data | Window | hashrope (MB/s) | Rabin-Karp (MB/s) | Ratio |
|------|--------|-----------------|-------------------|-------|
| 1 KB | 32 | 1.9 | 282.6 | 150× |
| 1 KB | 64 | 1.6 | 290.5 | 181× |
| 1 KB | 128 | 1.3 | 278.6 | 219× |
| 1 KB | 256 | 1.1 | 254.2 | 239× |
| 1 KB | 1024 | 1.3 | 206.0 | 159× |
| 10 KB | 32 | 1.7 | 297.9 | 174× |
| 10 KB | 64 | 1.5 | 306.6 | 198× |
| 10 KB | 128 | 1.1 | 289.0 | 255× |
| 10 KB | 256 | 0.8 | 268.1 | 317× |
| 10 KB | 1024 | 0.9 | 275.0 | 317× |
| 100 KB | 32 | 1.9 | 298.5 | 160× |
| 100 KB | 64 | 1.5 | 309.6 | 200× |
| 100 KB | 128 | 1.1 | 300.7 | 267× |
| 100 KB | 256 | 0.9 | 297.6 | 348× |
| 100 KB | 1024 | 0.9 | 308.6 | 342× |
| 1 MB | 32 | 1.9 | 315.3 | 169× |
| 1 MB | 64 | 1.6 | 311.5 | 200× |
| 1 MB | 128 | 1.1 | 303.9 | 268× |
| 1 MB | 256 | 0.9 | 300.6 | 350× |
| 1 MB | 1024 | 0.9 | 311.3 | 349× |

**Analysis**:
- hashrope sliding window: ~0.9–1.9 MB/s (O(log W) per byte — rope split/concat/evict at every step)
- Rabin-Karp: ~275–315 MB/s (O(1) per byte — ring buffer with single multiply)
- Gap: 150–350×, increasing with window size (more tree work per eviction)

The hashrope sliding window is **not competitive with Rabin-Karp on raw throughput**.
This is the expected cost of maintaining a full structural rope vs a scalar rolling hash.
The sliding window's value is that it provides the complete rope data structure at every
position — enabling O(log n) substr_hash, split, and repeat within the window — capabilities
Rabin-Karp cannot offer. In the UHC context, the sliding window processes LZ77 copy
references structurally, which is impossible with a scalar rolling hash.

The ~1–2 MB/s throughput is dominated by per-byte allocation overhead (Arc, Vec for each
leaf). This is a primary optimization target for Phase 2 (arena allocator, batched ingestion).

---

## E6: Memory Profiling

**Status**: PASS — Claim C6 confirmed (bounded, predictable overhead).

Memory estimated analytically: walk tree, count nodes.
Per node: ~64 bytes (16 ArcInner overhead + 48 NodeInner struct).
Leaf nodes additionally store payload bytes in a heap-allocated Vec.

### Sequential rope (one byte per leaf — worst case)

| n | leaves | internal | total nodes | est. bytes | overhead |
|---|--------|----------|-------------|------------|----------|
| 100 | 100 | 99 | 199 | 12,836 | 128.4× |
| 1,000 | 1,000 | 999 | 1,999 | 128,936 | 128.9× |
| 10,000 | 10,000 | 9,999 | 19,999 | 1,289,936 | 129.0× |
| 100,000 | 100,000 | 99,999 | 199,999 | 12,899,936 | 129.0× |
| 1,000,000 | 1,000,000 | 999,999 | 1,999,999 | 128,999,936 | 129.0× |

Overhead is a **constant 129×** regardless of n — each 1-byte leaf costs ~64 bytes
of node metadata. This is the degenerate case (no user would store one byte per leaf
in practice).

### Chunked rope (multi-byte leaves — realistic)

| n | chunk size | leaves | total nodes | est. bytes | overhead |
|---|------------|--------|-------------|------------|----------|
| 1,000 | 64 | 16 | 31 | 2,984 | 2.98× |
| 1,000 | 256 | 4 | 7 | 1,448 | 1.45× |
| 10,000 | 64 | 157 | 313 | 30,032 | 3.00× |
| 10,000 | 1024 | 10 | 19 | 11,216 | 1.12× |
| 10,000 | 4096 | 3 | 5 | 10,320 | 1.03× |
| 100,000 | 256 | 391 | 781 | 149,984 | 1.50× |
| 100,000 | 1024 | 98 | 195 | 112,480 | 1.12× |
| 100,000 | 4096 | 25 | 49 | 103,136 | 1.03× |
| 1,000,000 | 1024 | 977 | 1,953 | 1,124,992 | 1.12× |
| 1,000,000 | 4096 | 245 | 489 | 1,031,296 | 1.03× |

With realistic chunk sizes, overhead drops to **1.03–3×** — converging toward
the payload size as chunks grow. At 4 KB chunks, the tree metadata is only 3%
of total memory.

### RepeatNode memory vs materialized

| q | nodes | rope bytes | materialized bytes | savings |
|---|-------|------------|-------------------|---------|
| 10 | 2 | 228 | 1,000 | 4× |
| 100 | 2 | 228 | 10,000 | 44× |
| 1,000 | 2 | 228 | 100,000 | 439× |
| 10,000 | 2 | 228 | 1,000,000 | 4,386× |
| 100,000 | 2 | 228 | 10,000,000 | 43,860× |
| 1,000,000 | 2 | 228 | 100,000,000 | 438,596× |
| 100,000,000 | 2 | 228 | 10,000,000,000 | 43,859,649× |

RepeatNode uses exactly 2 nodes (228 bytes) regardless of q — the child leaf
plus the RepeatNode itself. At q=10⁸, this is **43.8 million× less memory**
than materialization.

---
