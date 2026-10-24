# hashrope

A BB[2/7] weight-balanced binary tree augmented with polynomial hash metadata at every node.

Zero dependencies. Multiple language implementations.

## Packages

| Language | Path | Status |
|----------|------|--------|
| Python | [`packages/python`](packages/python) | Implemented |
| Rust | [`packages/rust`](packages/rust) | Implemented |
| JavaScript / TypeScript | [`packages/npm`](packages/npm) | Node, ESM, CommonJS and browser builds |
| Java | [hashrope-java](https://github.com/jemsbhai/hashrope-java) | Maven Central package |
| C | `packages/c` | Planned |

Version 0.3.4 synchronizes the four core packages and both Python/Rust
`hashrope-bio` packages. This patch includes the accepted Rust LCP and concat
optimizations; the other runtimes retain their 0.3.3 behavior. Exact equality
and longest-common-prefix APIs verify
bytes, including compressed ranges; profile checks prevent combining inconsistent
hash metadata. Existing fingerprints remain candidate keys and can still collide.
Applications using hash equality alone must adopt the exact APIs.

The [release record](releases/0.3.4/README.md) distinguishes current artifact
validation from the historical benchmark studies below. Matching version numbers
do not imply identical language APIs, timings, or a logarithmic exact-comparison
bound.

## What it does

hashrope is a rope data structure where every node carries a polynomial hash, enabling O(log n) concat, split, substring hashing, and O(log q) repetition encoding -- all with strict BB[2/7] balance guarantees and Mersenne prime modular arithmetic.

Fingerprints can collide. Python's `rope_equal` / `rope_lcp` and Rust's
`Arena::equal` / `Arena::lcp` provide exact byte comparisons, with shortcuts for
shared structure and repeated content. Exact comparisons may scan unshared
content; they do not have a general logarithmic bound. Hash profiles are checked
when building ropes and making hash queries; existing fingerprint defaults are unchanged.
See the [comparison benchmarks](benchmarks/comparison_safety/README.md) for the
measured time and memory tradeoffs.

The latest [Rust regression investigation](benchmarks/rust_regression_round7/README.md)
retains the round-six implementation. Twelve fresh pairs using the original
executables estimate +2.0% for independent 100 KB edit/comparison and −0.7%
for independent 1 MiB LCP, with mixed directions and intervals including one.
The earlier 4.4% and 6.8% slowdowns do not repeat at stable magnitudes. However,
controlled pinned LCP retains an approximately 2% penalty in both slot orders,
and a secondary comparison-only estimate is slower. These costs remain
unresolved, not proven noise.

A trial that removes redundant join weight lookups passed 221 tests but did
not establish a useful gain for either target in another twelve-pair full
comparison. Some early-mismatch cases slowed, so the trial was rejected and
the previous core restored exactly. All results and the remaining Ropey gaps
are documented. Correctness, fingerprints, defaults and versions are unchanged;
all 221 tests pass again on the restored core. The accepted Rust optimizations
are included in 0.3.4; the rejected trial is excluded.

The [round-six study](benchmarks/rust_regression_round6/README.md) records the
retained concat optimization's 14.5% and 12.2% improvements in separate cohorts.
Its gain repeats at 15.2% in the latest original-executable replication and
13.6–15.5% in the controlled diagnostic. Those results apply to the measured
concat fixture and do not establish uniform comparison performance.

The earlier [Rust LCP study](benchmarks/rust_lcp_round5/README.md) compares
the round-five optimization with published 0.3.3 using the unchanged workloads.
Across nine fresh process pairs, independent-equal LCP takes 17.0% less time
at 64 KiB and 13.0% less at 1 MiB; self-concatenated 2 TiB DAG LCP takes 58.3%
less. These are paired geometric-mean changes, with improvements in 8/9, 9/9
and 9/9 pairs respectively. Aligned branches avoid redundant inspections, and
equal-period DAGs verify bytes directly at their period roots. Fingerprints,
defaults and node layout were unchanged; its 219 Rust release tests passed.
That cohort's broader results were mixed: independent 1 MB edit/comparison was 5.1% slower
in the pooled cohort (5/9 pairs slower), and the unchanged concat control was
8.4% slower (6/9). The report preserves the original three pairs and six-pair
confirmation separately, including variation and remaining Ropey gaps. These
changes are included in 0.3.4; the measured snapshots retain version 0.3.3.

The [edit and metadata optimization study](benchmarks/exact_optimization_round3/README.md)
groups the existing default hash into eight-byte evaluations and removes
redundant leaf copies. The [fourth-round rerun](benchmarks/exact_optimization_round4/README.md)
repeats that comparison with no further code change. At 10 KB, edit-and-exact-compare time falls from
7.43 to 4.47 µs for independent targets (1.66× faster), and from 6.64 to
3.70 µs for shared-history targets (1.79×). Both improve in all nine paired
runs, although Ropey remains faster at 1.17 and 0.91 µs respectively.

At 1 MB, the independent case changes from 0.153 to 0.144 ms,
still 1.35× Ropey's time. Shared-history time falls from 0.0147 to
0.00917 ms, 3.73× faster than Ropey in this fixture. These are means over
nine fresh processes with prebuilt targets; the report includes dispersion,
remaining smaller-input gaps to Ropey, and unchanged Python controls.
The original third-round comparison-only microbenchmarks include slower equality
and late-mismatch cells, retained in that report. That study kept fingerprints,
defaults, node layout and exact-comparison code unchanged.

The [follow-up regression investigation](benchmarks/exact_optimization_round4/README.md)
finds an unresolved, context-sensitive comparison cost. The original short-sample
protocol again reports slower 1 MiB equality: 54.926→65.895 µs (+20.0%);
late-mismatch LCP changes from 57.184→60.059 µs (+5.0%). These are medians
of three round medians, with substantial variation between rounds.
In contrast, 36 controlled processes and 12 paired runs of the exact original
executables do not show a consistent penalty of the original size. The latter's
paired geometric time ratios are 0.996 for equality and 0.991 for LCP.
The full-protocol rerun uses rebuilt executables and interleaves Python runs;
these context differences were not isolated. The evidence neither establishes
a general comparator slowdown nor dismisses the repeated equality result as
noise. That investigation selected no further production change and retained
the third-round implementation. The later Rust LCP work is recorded separately.

The [earlier traversal study](benchmarks/exact_optimization_round2/README.md)
and [first optimization study](benchmarks/exact_optimization/README.md)
remain preserved. The traversal study's structural-period proof reversed
the compressed-DAG regression without expanding the symbolic 2 TiB inputs.

See [`packages/python/README.md`](packages/python/README.md) for full API documentation.

## License

MIT
