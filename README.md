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

Version 0.3.3 synchronizes the four core packages and both Python/Rust
`hashrope-bio` packages. Exact equality and longest-common-prefix APIs verify
bytes, including compressed ranges; profile checks prevent combining inconsistent
hash metadata. Existing fingerprints remain candidate keys and can still collide.
Applications using hash equality alone must adopt the exact APIs.

The [release record](releases/0.3.3/README.md) distinguishes current artifact
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
and late-mismatch cells, retained in that report. Fingerprints, defaults,
node layout and exact-comparison code are unchanged.

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
noise. No further production change was justified, so correctness and
fingerprints are preserved by keeping the implementation byte-identical to
the third round.

The [earlier traversal study](benchmarks/exact_optimization_round2/README.md)
and [first optimization study](benchmarks/exact_optimization/README.md)
remain preserved. The traversal study's structural-period proof reversed
the compressed-DAG regression without expanding the symbolic 2 TiB inputs.

See [`packages/python/README.md`](packages/python/README.md) for full API documentation.

## License

MIT
