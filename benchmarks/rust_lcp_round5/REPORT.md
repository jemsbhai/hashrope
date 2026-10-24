# Rust LCP round-five results

The new candidate reduced the intended LCP work and improved the principal long
LCP cells against a freshly measured published-0.3.3 baseline. The fixed study
does not establish that every workload improved: independent 1M exact churn has
an unresolved 11.5% slower paired estimate, and some smaller paths vary between
pairs. No raw sample was excluded or replaced.

## Paired comparison results

The times below are medians of three process medians. Ratios are the geometric
means of the **three process-paired candidate/baseline ratios**; ranges show those
three ratios, not confidence intervals. A ratio below one is faster. Differences
between the ratio of displayed times and the paired estimator are intentional.

| Exact operation | Baseline | Candidate | Paired ratio (range) | Candidate slower pairs |
|---|---:|---:|---:|---:|
| 64KiB independent-equal LCP | 3.377 µs | 2.926 µs | 0.862 (0.695–1.062) | 1/3 |
| 1MiB independent-equal LCP | 59.866 µs | 51.164 µs | 0.864 (0.838–0.902) | 0/3 |
| 2TiB self-concat DAG LCP | 339 ns | 142 ns | 0.423 (0.366–0.498) | 0/3 |
| 64KiB last-byte-mismatch LCP | 3.458 µs | 3.012 µs | 0.831 (0.683–0.958) | 0/3 |
| 1MiB last-byte-mismatch LCP | 58.186 µs | 52.667 µs | 0.879 (0.855–0.917) | 0/3 |
| 1MiB independent equality | 60.670 µs | 51.181 µs | 0.865 (0.837–0.917) | 0/3 |
| 1MiB last-byte-mismatch equality | 57.329 µs | 61.532 µs | 1.038 (0.911–1.134) | 2/3 |
| 2TiB self-concat DAG equality | 281 ns | 158 ns | 0.561 (0.510–0.615) | 0/3 |
| 2TiB Repeat LCP | 76.55 ns | 80.10 ns | 1.035 (0.982–1.079) | 2/3 |
| 2TiB Repeat equality | 77.80 ns | 82.25 ns | 1.033 (0.973–1.070) | 2/3 |

The primary 1MiB independent-equal LCP and DAG LCP estimates improve by 13.6% and
57.7%, respectively, with all three paired directions agreeing. The 64KiB
primary estimate improves by 13.8%, but one pair is 6.2% slower. The long late-LCP
estimate improves by 12.1% in all pairs; the subsequent late-equality cell remains
mixed and cannot be described as an equivalent timing win.

Early-mismatch LCP and equality, and both shared-late-mismatch sizes, improved in
all three pairs. Same-root calls remain approximately 2 ns for LCP and 4 ns for
equality; these tiny cells include loop and dispatch overhead. No microbenchmark
cell exceeded the predeclared threshold of a 10% slower geometric mean with all
three paired directions slower. This flag does not rule out smaller changes:
the Repeat cells above are around 3% slower, and concat/split estimates are about
8% slower with mixed directions. All 42 micro cells, including unsafe historical
hash-only controls, remain in [all_cells.csv](all_cells.csv) and
[summary.json](summary.json).

## Edit-and-compare outcomes

Each timing is per edit followed by comparison; each process repeats the same
50-edit, 50%-equal fixture. Input sizes here are decimal bytes, unlike the KiB/MiB
comparison fixtures above. The target histories and all four arms are unchanged.

| Exact Hashrope churn | Baseline | Candidate | Paired ratio (range) | Candidate slower pairs |
|---|---:|---:|---:|---:|
| 10K independent | 4.836 µs | 4.960 µs | 1.026 (1.015–1.036) | 3/3 |
| 10K shared | 3.824 µs | 3.754 µs | 0.973 (0.938–1.000) | 0/3 |
| 100K independent | 20.538 µs | 19.092 µs | 0.955 (0.930–0.976) | 0/3 |
| 100K shared | 7.906 µs | 7.338 µs | 0.923 (0.867–0.982) | 0/3 |
| 1M independent | 150.730 µs | 171.468 µs | 1.115 (0.933–1.307) | 2/3 |
| 1M shared | 10.474 µs | 9.576 µs | 0.888 (0.780–0.969) | 0/3 |

The independent 1M result is a material unresolved slowdown estimate. Its paired
ratios are 0.933, 1.138 and 1.307; it falls below the all-three-direction flag only
because the first pair is faster. The 10K independent estimate is a smaller 2.6%
slowdown with all three directions agreeing. Neither is hidden by the favorable
shared-history results.

The unchanged Ropey exact controls also vary: independent 1M has a 1.299 paired
ratio with a 0.906–2.038 range. Shared 10K and 1M Ropey exact cells trigger the
predeclared flags (1.106 and 1.211, respectively, all three pairs slower). This
shows variability in the process cohort; it does **not** prove that the Hashrope
changes are noise, identify a cause, or cancel the observed Hashrope slowdowns.
Every one of the 24 churn cells remains in the complete data, including the
hash-only candidate and Ropey-rehash arms.

## Mechanism and correctness evidence

Only `rope.rs` has a production content change; the manifest's raw CRLF/LF
difference is recorded separately and its normalized content is identical.
Fingerprint arithmetic, package versions, default profile and public API are
unchanged. The [mechanism audit](mechanism_audit.md) reviews the two guarded paths:
aligned child boundaries avoid rediscovering the same Internal-node boundaries;
equal proven periods compare bounded underlying-unit ranges directly.

Separate instrumented copies reduce normalizer node visits from 8,190 to 4,096
for the 1MiB aligned independent fixture, while retaining 2,048 leaf comparisons
over 1,048,576 input bytes. For the 2TiB DAG, normalizer visits fall from 84 to 2
and one 368-byte requested heap allocation is removed. Period-chain inspections
remain 82. These are actual source events and allocation requests, not CPU
instructions, memory-load counts, peak memory or timings. The instrumented probe
is separate from the unmodified benchmark binaries, and all sixteen probe
oracles passed. These work reductions support a mechanism for improvement; no
timed ablation isolates each change's individual contribution.

The candidate passed [full validation](validation/results.json): 219 all-feature
release tests including the experiments and doctest, default doctests, format,
strict all-feature/all-target Clippy, wasm32 checks with and without bigint, and
warning-free all-feature documentation. Published-profile/fingerprint fixtures,
collision/range tests and new adversarial phase/boundary cases are retained.
All twelve benchmark processes also passed their untimed answer oracles.

Layout rows are unchanged in every process: `NodeInner` 48 bytes, inline `Arena`
2,136 bytes, and 255 stored nodes for the 64KiB construction. These rows do not
measure peak or retained process memory.

## Provenance and limits

[PLAN.md](PLAN.md) was frozen before timing. The run used three alternating fresh
process pairs for each unchanged workload: baseline/candidate,
candidate/baseline, baseline/candidate; churn seeds were 1/2/3 with invocation 0.
All four binaries were built before the exclusive timing window. The compiler,
byte-identical workload and lock hashes, six-file core snapshots and executable
archive are recorded in [manifest.json](manifest.json) and [build.json](build.json).
The process cohort retains 1,764 micro samples and 432 churn samples; raw stdout,
stderr, arguments and oracle success are preserved under [raw/](raw/).

The immutable baseline comes from committed published-0.3.3 source at `ffba5ee`.
It is independently measured here. Earlier
[release-integration](../../releases/0.3.3/performance/rust_integration/README.md)
and [round-four](../exact_optimization_round4/controlled_findings.md) findings
remain historical context; no new time is divided by an old time to infer a
controlled improvement. The historical exact engine is identical after uniform
formatting, which does not imply identical generated code or explain its earlier
regressions. The round-four protocol-sensitive equality observations remain
unresolved.

This is the specified three-pair diagnostic on one unpinned Windows host.
Frequency, thermal state and OS background activity were not controlled. The
full historical outer Python interleaving was omitted, consistently with the
release-integration protocol. The data support improvements in the intended LCP
fixtures, with mixed outcomes elsewhere; they do not establish universal
improvement or absence of smaller regressions. Any follow-up must predeclare its
fixed cohort separately and retain these results intact.
