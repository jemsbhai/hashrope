# hashrope — Experiment Protocol

**Version**: 0.1 (draft)
**Date**: 2026-04-07
**Authors**: [your name]
**Status**: Design phase — no results collected yet

---

## 1. Claims

The hashrope paper makes the following empirical claims.
Each claim maps to one or more experiments below.

| ID | Claim | Experiments |
|----|-------|-------------|
| C1 | Polynomial hash throughput is competitive with general-purpose hashes on raw byte streams | E1 |
| C2 | Rope concat/split is asymptotically faster than naive rehash-from-scratch | E2 |
| C3 | RepeatNode gives O(log q) hash computation vs O(q) materialized | E3 |
| C4 | Phi accumulator runs in O(log q) time | E4 |
| C5 | Sliding window throughput is competitive with Rabin-Karp rolling hash | E5 |
| C6 | Memory overhead of the tree structure is bounded and predictable | E6 |
| C7 | BB[2/7] balancing maintains O(log w) tree height in practice | E7 |

---

## 2. Baselines

| Baseline | Crate | Version | Why |
|----------|-------|---------|-----|
| BLAKE3 | `blake3` | 1.x | Fast cryptographic hash, SIMD-optimized; strong throughput baseline |
| SHA-256 | `sha2` | 0.10.x | Standard cryptographic hash; widely understood reference point |
| xxhash3 | `xxhash-rust` (feature `xxh3`) | 0.8.x | Fast non-cryptographic hash; closest competitor class |
| Rabin-Karp | hand-rolled | — | Classical rolling hash; direct competitor to sliding window |
| ropey | `ropey` | 1.x | Production rope crate (UTF-8 text); structural operations baseline |
| Vec<u8> + rehash | hand-rolled | — | Naive baseline: append bytes, rehash entire buffer |

**Note**: BLAKE3, SHA-256, and xxhash3 are *not* structurally composable (no split/concat in O(log n)).
The throughput comparison (E1) is apples-to-apples on raw sequential hashing only.
Experiments E2–E5 demonstrate capabilities these baselines cannot provide at all, or only at O(n) cost.

---

## 3. Datasets

### 3.1 Canonical corpora

| Corpus | Source | Purpose |
|--------|--------|---------|
| Canterbury corpus | https://corpus.canterbury.ac.nz/ | Standard compression benchmark; mixed file types |
| Silesia corpus | https://sun.aei.polsl.pl/~sdeor/index.php?page=silesia | Larger (211 MB); covers text, binaries, images, databases |

### 3.2 Synthetic data

| Dataset | Description | Sizes |
|---------|-------------|-------|
| `random_bytes` | Uniform random via `rand` crate (seeded) | 1 KB, 10 KB, 100 KB, 1 MB, 10 MB, 100 MB |
| `english_text` | Repeated English paragraph (controlled) | Same sizes |
| `all_zeros` | Degenerate case: all 0x00 bytes | Same sizes |
| `short_strings` | Random 8–64 byte strings | 10K strings |

### 3.3 Seeds and reproducibility

- All random data generated with `rand::SeedableRng` / `StdRng::seed_from_u64(42)`
- Corpus files identified by SHA-256 checksum logged in results
- Exact corpus file versions recorded

---

## 4. Experiments

### E1: Sequential hashing throughput (MB/s)

**Claim**: C1
**Question**: How does hashrope's polynomial hash compare to BLAKE3, SHA-256, and xxhash3 on raw sequential byte streams?

**Method**:
1. For each dataset d ∈ {Canterbury files, Silesia files, synthetic 1KB–100MB}:
2. For each hasher h ∈ {hashrope polynomial, BLAKE3, SHA-256, xxhash3}:
3. Feed d to h as a single contiguous byte slice
4. Measure wall-clock time via Criterion (minimum 10 iterations, 5s warm-up)
5. Compute throughput = |d| / time (MB/s)

**Output**: Table and line plot of throughput (MB/s) vs input size, one line per hasher.

**Expected outcome**: hashrope's polynomial hash will be slower than BLAKE3/xxhash3 (which use SIMD).
The point is to quantify the gap, not to claim superiority on raw throughput.
The value proposition of hashrope is structural composability (E2–E5), not raw speed.

**Fairness**: All hashers process the same byte slices. No pre-warming of hashrope's power cache between trials (each iteration starts with a fresh `PolynomialHash` instance).

---

### E2: Rope concat/split vs naive rehash

**Claim**: C2
**Question**: How does hashrope's O(log n) concat/split compare to rebuilding the hash from scratch?

**Method**:
1. Build a rope R of n random bytes (one byte per leaf), n ∈ {100, 1K, 10K, 100K, 1M}
2. **Concat benchmark**: Split R into two halves, time the `join` operation. Compare to: hash the full n bytes from scratch.
3. **Split benchmark**: Time `split_at(n/2)`. Compare to: re-hash each half from scratch.
4. Criterion, ≥10 iterations per size.

**Output**: Log-log plot of time vs n for each operation. Expect O(log n) for rope ops, O(n) for naive.

---

### E3: RepeatNode scaling

**Claim**: C3
**Question**: Does RepeatNode compute the hash of q repetitions in O(log q) time?

**Method**:
1. Create a base rope B of 100 random bytes.
2. For q ∈ {10, 100, 1K, 10K, 100K, 1M, 10M, 100M}:
   a. `repeat_node(B, q)` — time the construction (should be O(log q))
   b. Naive baseline: materialize q copies of B into a `Vec<u8>`, hash from scratch — time is O(q·|B|)
3. Criterion, ≥10 iterations per q.

**Output**: Log-log plot of time vs q. Fit slopes; expect slope ≈ 0 for repeat_node (O(log q) on log scale ≈ flat for large q), slope ≈ 1 for naive.

**Memory**: For large q (10M, 100M), the naive baseline materializes 1–10 GB. Skip the naive baseline for q > 1M if memory is insufficient, and note this in results. The RepeatNode itself should use O(log q) memory — verify via allocation counting.

---

### E4: Phi accumulator

**Claim**: C4
**Question**: Does the phi accumulator (geometric series computation) run in O(log q)?

**Method**:
1. For q ∈ {10, 100, 1K, 10K, 100K, 1M, 10M}:
   a. Time `phi(base_hash, base_pow, q)` via Criterion
2. Plot time vs q on log-log axes. Fit slope; expect ≈ log(q) growth.

**Output**: Log-log plot with fitted line. Report slope and R².

---

### E5: Sliding window throughput

**Claim**: C5
**Question**: How does hashrope's sliding window compare to a classical Rabin-Karp rolling hash?

**Method**:
1. For each dataset d ∈ {synthetic random 1KB–100MB, Canterbury, Silesia}:
2. Window size w ∈ {32, 64, 128, 256, 1024} bytes
3. Slide over d one byte at a time, computing the window hash at each position
4. Compare: hashrope `SlidingRopeState` vs hand-rolled Rabin-Karp (same polynomial base, same modulus for fairness)
5. Criterion, ≥10 iterations.

**Output**: Throughput (MB/s) vs input size, one line per method per window size.

**Rabin-Karp baseline implementation**: Use the same Mersenne-prime modular arithmetic and base as hashrope to isolate the structural overhead (tree maintenance vs simple ring buffer).

---

### E6: Memory profiling

**Claim**: C6
**Question**: What is the memory overhead of the rope tree structure vs a flat byte buffer?

**Method**:
1. For n ∈ {1K, 10K, 100K, 1M} bytes:
   a. Build a rope from n bytes (one byte per leaf)
   b. Measure peak RSS (via `jemalloc` stats or OS-level measurement)
   c. Compare to n bytes in a `Vec<u8>`
   d. Compute overhead ratio = (rope RSS) / n
2. For RepeatNode: measure RSS for `repeat_node(B, q)` with q up to 10^8. Compare to q·|B|.

**Output**: Table of overhead ratios. Expected: constant factor overhead for rope (each internal node stores ~64 bytes of metadata).

**Tool**: `dhat` crate (heap profiler for Rust) or `jemalloc-ctl` for allocation stats.

---

### E7: BB[2/7] tree height verification

**Claim**: C7
**Question**: Does the BB[2/7] balancing scheme maintain O(log w) tree height in practice?

**Method**:
1. Build ropes by sequential single-byte insertions (worst-case for balance):
   n ∈ {100, 1K, 10K, 100K, 1M}
2. After construction, measure actual tree height h(n)
3. Compare to theoretical bound: h ≤ c · log₂(n) where c = 1/log₂(1/(1−2/7)) ≈ 2.58 for BB[2/7]
4. Also build ropes by random-order insertions and repeated join of random-sized sub-ropes.

**Output**: Table of (n, h_actual, h_theoretical_bound). Plot h vs log₂(n).

---

## 5. Environment

All benchmarks run on the same machine:
- CPU: [to be filled — run `wmic cpu get name`]
- RAM: 64 GB
- GPU: NVIDIA RTX 4090 (not used for benchmarks — CPU only)
- OS: Windows [version to be filled]
- Rust toolchain: [to be filled — run `rustc --version`]
- Compiler flags: `--release` (default `opt-level = 3`)
- Criterion config: default (5s warm-up, 5s measurement, auto-tuned iterations)

Benchmarks must NOT run with other CPU-intensive processes active.
CPU frequency scaling should be noted (Windows typically uses "Balanced" power plan — switch to "High Performance" for benchmarks).

---

## 6. Reporting standards

- All times reported as median ± MAD (median absolute deviation) from Criterion
- Throughput in MB/s = input_bytes / median_time
- Log-log plots for scaling experiments with least-squares fitted slopes
- R² values for all fitted lines
- Raw Criterion JSON output archived in `benches/results/` for reproducibility
- No cherry-picking: all runs reported, including outliers (Criterion handles this)

---

## 7. Known limitations (to be stated in paper)

- Polynomial hash is not cryptographic — no collision resistance claims
- Mersenne-prime arithmetic is u64; no u128 or bignum benchmarks
- Current implementation is single-threaded; no parallel hashing benchmarks
- `Arc`-based node sharing adds atomic reference counting overhead (optimization planned)
- Power cache has O(1) lookup for exponents < 256, O(log n) for larger — benchmarks should test both regimes

---

## 8. Experiment execution order

1. E7 (tree height) — cheapest, validates data structure correctness
2. E4 (phi accumulator) — micro-benchmark, fast
3. E3 (RepeatNode) — builds on E4
4. E1 (throughput) — requires baseline crate dependencies
5. E2 (concat/split) — requires ropey dependency
6. E5 (sliding window) — requires Rabin-Karp baseline
7. E6 (memory) — requires profiling tooling, run last
