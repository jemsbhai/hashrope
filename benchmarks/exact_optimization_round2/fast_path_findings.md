# Ordinary-tree comparison diagnostics

These experiments isolate the Rust fast path from the compressed-DAG changes.
They use the immutable round-two baseline, copied into separate release-build
projects. Fingerprints, arena layout, node construction and default parameters
are identical across the diagnostic variants. The final unchanged benchmark
matrix is the source for overall performance claims; these smaller runs explain
which changes were selected.

The pinned Ropey 1.6.1 implementation was also inspected locally:
`src/rope.rs` delegates equality to `RopeSlice`, and `src/slice.rs` advances
two chunk iterators while comparing their overlapping byte slices. This gives
a concrete comparison point for Hashrope's repeated recursive alignment of
different tree boundaries. Hashrope's construction still uses its existing
512-byte maximum leaves and binary structure; neither layout nor editing
policy is changed to improve the benchmark.

## Selected changes

- A recursive comparison enters the relevant Internal children directly and
  advances directly to the right child when its boundary is consumed. This
  avoids inspecting the same parent again during normalization.
- Once native slice comparison proves inequality, safe eight-byte loads and a
  trailing-zero count locate the first differing byte. The conversion is
  explicitly little endian, so the byte position is portable. No unsafe loads
  or fingerprint decisions are used.
- Cross-arena ranges of at least 8 KiB use two explicit traversal frontiers.
  Each frontier has 32 pending-range slots: together they occupy 1,600 bytes
  of frontier storage on 64-bit targets, excluding other locals, with no heap
  allocation. Repeats, a full
  frontier stack, or an exhausted node-visit budget return the already verified
  prefix to the general exact engine. Same-arena comparisons retain recursion
  because it performed better on shared subtrees.

The public APIs remain independent of fingerprint metadata. In particular, an
equality-only fingerprint rejection shortcut was not added to favor the mixed
equal/unequal benchmark. The accepted custom-profile regression demonstrates
that identical bytes can have shape-dependent stored fingerprints under a
non-Mersenne profile, making that separation consequential.

## Controlled fast-path ablations

`fast_path_diagnostics.json` records three alternating rounds, each with seven
timed samples per cell. Fixtures are built outside timing, and every returned
LCP is asserted against the known byte result. Values below are medians of the
three per-round medians, in microseconds.

| Case | Round-two baseline | Direct children | Direct children + word locator |
|---|---:|---:|---:|
| 1 MB independently built equal inputs | 66.151 | 61.885 | 61.696 |
| 1 MB equal inputs, different chunk boundaries | 216.296 | 206.293 | 208.945 |
| 1 MB shared input, final-byte mismatch | 0.152 | 0.097 | 0.106 |
| 512-byte leaf, mismatch at byte 511 | 0.12488 | 0.16806 | 0.02917 |

The first-byte mismatch in a 512-byte leaf changed from 8.82 ns to 9.52 ns
with the combined variant, an approximately 0.7 ns increase in this diagnostic.
The word locator alone measured 8.94 ns. The gain is therefore strongest for
later mismatches; it should not be described as a universal improvement in
early-exit latency. The cross-arena frontier's first-byte result below showed
no corresponding increase for the long-input case.

The separate frontier comparison in `fast_path_frontier_diagnostics.json` used
the same three-round/seven-sample protocol. Comparing within this run avoids
attributing between-run timing variation to a code change.

| Case | Direct children + word locator | Cross-arena frontier |
|---|---:|---:|
| 1 MB independently built equal inputs | 57.132 | 57.867 |
| 1 MB equal inputs, different chunk boundaries | 198.401 | 150.966 |
| 1 MB equal ranges with nonzero offsets | 199.413 | 148.913 |
| 1 MB input, first-byte mismatch | 0.082 | 0.079 |
| 1 MB shared input, final-byte mismatch | 0.100 | 0.090 |

Enabling the frontier for same-arena inputs increased the shared late-mismatch
cell from 0.100 to 0.132 microseconds, so that variant was rejected. The
production dispatcher chooses between recursive and frontier traversal once per
fast comparison; recursive prefixes do not repeat the dispatch check. The
diagnostic prototype placed that check inside the recursive entry point. Final
matrix measurements cover the production version.

## Unchanged edit-plus-comparison workload

`fast_path_churn.json` runs the existing `edit_compare.rs` harness unchanged:
50 seeded edits followed by exact comparison, half equal targets and half
final-byte mismatches, with independent and shared histories separated. These
are preliminary single-process results for seed 1, with three samples per cell.
They are not the final nine-process matrix.

| Independent-history edit + exact comparison | Baseline Hashrope | Frontier Hashrope | Frontier-run Ropey |
|---|---:|---:|---:|
| 10 KB | 8.980 us | 7.898 us | 1.240 us |
| 100 KB | 32.580 us | 28.462 us | 12.082 us |
| 1 MB | 180.508 us | 169.396 us | 113.198 us |

The diagnostic controls vary appreciably between processes, especially shared
Ropey comparisons, so the final matrix must determine the size of the gain. The
10 KB Hashrope edit-plus-fingerprint-candidate arm already takes about 4.6
microseconds, versus about 1.2 microseconds for Ropey's edit-plus-exact arm.
This indicates that the small-document gap includes edit and metadata costs;
it is not all exact-comparison traversal.

## Reproduction and verification

`fast_path_profile.py --build-only` prepares the first four variants.
`--run --output NEW_FILE` executes their timing runs. Its `--variants` option
also accepts `children_words_frontier` and `children_words_frontier_cross`.
`fast_path_churn.py` builds/runs the unchanged edit harness against those
diagnostic libraries. Both runners refuse to overwrite result JSON and record
source, harness and compiled-binary hashes. Cached build directories are
disposable; the baseline archive and recorded result files are preserved.

Private regression tests force budget exhaustion at both node inspections,
frontier depth exhaustion before and after a verified prefix, Repeat fallback
after an 8 KiB prefix, and memoized continuation from each fallback position.
The independent review suite additionally checks uneven tree boundaries,
nonzero range offsets, all first-mismatch positions through 512-byte leaves,
lazy and cross-profile comparisons, and enormous compressed DAGs.
