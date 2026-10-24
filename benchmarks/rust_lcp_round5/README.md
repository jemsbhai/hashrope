# Rust LCP round five

The candidate reduces the intended LCP work and improves the primary LCP
fixtures against a freshly measured published-0.3.3 baseline. Across all nine
pairs per cohort, independent-equal LCP improves by 17.0% at 64KiB and 13.0% at
1MiB, and the 2TiB self-concatenation DAG improves by 58.3%, using paired
geometric-mean latency ratios. This is **not** a blanket no-regression result:
independent 1M exact churn retains a 5.1% slower pooled estimate with mixed
directions, and the concat control is 8.4% slower.

The [original three-pair report](REPORT.md) and its data remain intact. A
[separately predeclared confirmation](confirmation/PLAN.md) added six pairs
using the same four archived executables, with no rebuilding, source edits,
dependency changes or excluded samples. Both phases are reported separately
below, alongside all nine pooled pairs.

## Exact comparison results

Times are medians of nine process medians. Ratios are geometric means of the
candidate/baseline **paired process ratios** in each phase, so they need not
equal the ratio of the displayed times. Smaller is faster. Ranges are observed
paired-ratio ranges, **not confidence intervals**. The last column counts slower
candidate pairs among all nine.

| Exact operation | Baseline | Candidate | Original 3 ratio | Supplemental 6 ratio | Pooled 9 ratio (range) | Slower |
|---|---:|---:|---:|---:|---:|---:|
| 64KiB independent-equal LCP | 3.377 µs | 2.926 µs | 0.862 | 0.815 | 0.830 (0.689–1.062) | 1/9 |
| 1MiB independent-equal LCP | 58.282 µs | 51.985 µs | 0.864 | 0.873 | 0.870 (0.743–0.962) | 0/9 |
| 2TiB self-concat DAG LCP | 310 ns | 139 ns | 0.423 | 0.415 | 0.417 (0.328–0.507) | 0/9 |
| 64KiB last-byte-mismatch LCP | 3.458 µs | 2.947 µs | 0.831 | 0.840 | 0.837 (0.683–1.017) | 1/9 |
| 1MiB last-byte-mismatch LCP | 57.255 µs | 52.667 µs | 0.879 | 0.929 | 0.912 (0.855–0.974) | 0/9 |
| 1MiB independent equality | 59.670 µs | 51.588 µs | 0.865 | 0.869 | 0.868 (0.767–0.937) | 0/9 |
| 1MiB last-byte-mismatch equality | 57.584 µs | 52.236 µs | 1.038 | 0.886 | 0.934 (0.848–1.134) | 2/9 |
| 2TiB self-concat DAG equality | 281 ns | 158 ns | 0.561 | 0.537 | 0.545 (0.472–0.627) | 0/9 |
| 2TiB Repeat LCP | 80.45 ns | 80.10 ns | 1.035 | 1.003 | 1.013 (0.900–1.100) | 5/9 |
| 2TiB Repeat equality | 82.80 ns | 82.25 ns | 1.033 | 0.991 | 1.005 (0.878–1.108) | 5/9 |

The primary long independent LCP, DAG LCP and 1MiB late-LCP improvements repeat
in all nine pairs. The initially mixed 64KiB LCP improves in all six supplemental
pairs; its original slower pair remains visible. Early-mismatch and shared-late
exact paths improve in every pair at both sizes. Same-root calls remain around
2 ns for LCP and 4 ns for equality, including loop/dispatch overhead.

Repeat LCP is 1.35% slower pooled, with five slower pairs. Core concat is
16.35→17.50 ns, with original/supplemental/pooled ratios 1.079/1.087/1.084
(range 0.902–1.306; six slower pairs). Split has a mixed 2.5% slower pooled
estimate. These controls remain visible despite their algorithms being
unchanged. All 42 micro cells, including unsafe historical hash-only controls,
are retained in the full data.

## Exact edit-and-compare results and Ropey controls

Each timing is per edit followed by equality comparison, with 50 edits and half
equal/half final-byte-different targets. Sizes are decimal bytes. Both histories
and all four original arms are retained; this table shows the two exact arms.
Baseline/candidate refer to the Hashrope binary build; Ropey's source and
version are unchanged between them.

| Input/history/arm | Baseline | Candidate | Original 3 ratio | Supplemental 6 ratio | Pooled 9 ratio (range) | Slower |
|---|---:|---:|---:|---:|---:|---:|
| 10K independent, Hashrope | 4.836 µs | 4.598 µs | 1.026 | 0.957 | 0.979 (0.853–1.042) | 5/9 |
| 10K independent, Ropey | 1.114 µs | 1.082 µs | 1.097 | 1.052 | 1.067 (0.818–1.567) | 5/9 |
| 10K shared, Hashrope | 3.854 µs | 3.716 µs | 0.973 | 0.923 | 0.939 (0.842–1.008) | 1/9 |
| 10K shared, Ropey | 0.880 µs | 0.882 µs | 1.106 | 0.977 | 1.018 (0.883–1.162) | 5/9 |
| 100K independent, Hashrope | 19.160 µs | 18.650 µs | 0.955 | 0.975 | 0.968 (0.870–1.071) | 3/9 |
| 100K independent, Ropey | 10.752 µs | 11.342 µs | 1.047 | 0.973 | 0.997 (0.719–1.204) | 4/9 |
| 100K shared, Hashrope | 7.906 µs | 7.426 µs | 0.923 | 0.970 | 0.954 (0.781–1.326) | 2/9 |
| 100K shared, Ropey | 5.216 µs | 5.058 µs | 0.946 | 1.044 | 1.010 (0.783–1.211) | 4/9 |
| 1M independent, Hashrope | 144.530 µs | 149.156 µs | 1.115 | 1.020 | 1.051 (0.933–1.307) | 5/9 |
| 1M independent, Ropey | 107.020 µs | 106.148 µs | 1.299 | 1.011 | 1.099 (0.798–2.038) | 5/9 |
| 1M shared, Hashrope | 9.932 µs | 8.420 µs | 0.888 | 0.874 | 0.879 (0.780–1.004) | 1/9 |
| 1M shared, Ropey | 34.466 µs | 35.582 µs | 1.211 | 0.982 | 1.053 (0.752–1.288) | 8/9 |

Independent 1M exact Hashrope churn remains an unresolved tradeoff: the initial
11.5% slower estimate becomes 2.0% slower in the supplement and 5.1% slower
pooled, with five of nine directions slower. The initial 10K independent
slowdown reverses in the supplemental phase and has a slightly favorable pooled
ratio, but five of nine pairs are still slower. Shared 1M has a favorable ratio
in both phases and improves in eight of nine pairs.

Independent-target churn uses the unchanged long cross-arena frontier, whereas
the primary same-arena micro roots use the newly optimized aligned recursion.
The concat algorithm is also unchanged. That distinction does not rule out
effects from generated-code layout or runtime context, establish noise, or
isolate a cause. The unchanged Ropey controls vary substantially between phases;
their variation cannot cancel the observed Hashrope tradeoffs. The small-input
Ropey edit gap remains: these LCP changes do not make Hashrope faster in every
exact-edit workload.

## Mechanism and correctness

The [mechanism audit](mechanism_audit.md) reviews two guarded changes: aligned
Internal children avoid repeated normalization/boundary discovery, and equal
proven periods compare bounded underlying-unit ranges directly. Neither uses
fingerprint equality as proof or changes fingerprint arithmetic.

Separate instrumented copies show 1MiB aligned independent normalizer visits
falling from 8,190 to 4,096, with the same 2,048 leaf comparisons over 1,048,576
input bytes. For the DAG, visits fall from 84 to 2 and a 368-byte requested heap
allocation is removed; 82 period-chain inspections remain. These are source
events and requested allocations, not CPU instructions, actual memory-load
counts, peak memory or timings. All sixteen probe oracles passed. No timed
ablation assigns a separate speedup to each change.

Layout rows are unchanged in every process: 48-byte `NodeInner`, 2,136-byte
inline `Arena`, and 255 stored nodes for the 64KiB fixture. They do not measure
peak or retained memory. Only `rope.rs` changes production content; Cargo.toml's
LF/CRLF-only difference is explicitly recorded, with normalized content required
to match. Version, default hash profile, public API and fingerprint arithmetic
remain unchanged.

The [validation report](VALIDATION.md) and [machine-readable evidence](validation/results.json)
record 219 all-feature release tests including experiments and doctest, default
doctests, format, strict all-feature/all-target Clippy, wasm32 checks in both
feature modes, and warning-free documentation. Collision/range/phase and
published fingerprint compatibility checks are part of these gates. The narrow
default backend is timed; bigint correctness is tested separately.

## Data and reproducibility

Both phases retain **36 fresh processes and 6,588 timing samples**: 5,292 micro
and 1,296 churn samples. Every process passed its existing outside-timing
oracles. The original 12 processes/2,196 samples and supplemental 24
processes/4,392 samples are stored separately. No raw data, process or sample was
replaced. The supplement continues the alternating sequence; churn seeds 1/2/3
are repeated with invocation 0. All four binaries were built before the first
exclusive timing slot and reused byte for byte in the second slot.

- [Original plan](PLAN.md), [manifest](manifest.json), [build record](build.json),
  [raw data](raw/), [results](results.json), [summary](summary.json),
  [all cells](all_cells.csv), and [original report](REPORT.md).
- [Confirmation plan](confirmation/PLAN.md), [manifest](confirmation/manifest.json),
  [raw data](confirmation/raw/), [results](confirmation/results.json),
  [separate six/pooled-nine summaries](confirmation/summary.json), and
  [all supplemental/pooled cells](confirmation/all_cells.csv).
- `baseline.zip` preserves committed published source at `ffba5ee`;
  `candidate.zip` preserves the reviewed implementation. `drivers.zip` and
  `binaries.zip` retain the drivers and original executables. Source, lock,
  binary and raw hashes are checked before and after timing. The original 66
  summaries are reconstructed from raw data again during confirmation analysis.

Workload sources are byte-identical to the release-integration study, using its
exact locks (Ropey 1.6.1, smallvec 1.16.4 and str_indices 0.4.4).

```text
# Original phases, in a fresh directory with frozen baseline/workloads:
python benchmarks/rust_lcp_round5/run.py freeze
python benchmarks/rust_lcp_round5/run.py build
python benchmarks/rust_lcp_round5/run.py run
python benchmarks/rust_lcp_round5/summarize.py

# Supplemental protocol; reuses original archived executables:
python benchmarks/rust_lcp_round5/confirmation/confirm.py freeze
python benchmarks/rust_lcp_round5/confirmation/confirm.py run
python benchmarks/rust_lcp_round5/confirmation/summarize_confirmation.py
```

The scripts refuse to overwrite completed evidence. Timing requires a separately
coordinated exclusive CPU slot; neither run phase builds anything.

## Limits and historical context

This is one unpinned Windows host across two declared phases. CPU frequency,
thermal state and OS background work are uncontrolled. Nine pairs with repeated
seeds are not nine independent workload draws. Ranges are descriptive, and
favorable pooled results do not bound every smaller regression or all future
workloads. The retained mixed results warrant that qualification even though
the principal LCP gains reproduce.

The [release-integration study](../../releases/0.3.3/performance/rust_integration/README.md)
and [round-four evidence](../exact_optimization_round4/controlled_findings.md)
are historical context. No new time is divided by an old time to claim a
controlled gain. Their exact engine matches after uniform formatting, which
does not establish identical machine code or explain the earlier regressions.
The full round-four outer Python interleaving is absent here, and its prior
protocol-sensitive equality observations remain unresolved.
