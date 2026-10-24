# Rust release integration measurements

The 0.3.3 integration passes the unchanged benchmark oracles and retains compressed
comparison, but these measurements do **not** support a blanket no-regression
claim. Independent-equal LCP and self-concatenated DAG LCP are consistently slower
in this three-pair diagnostic. Other cells improve or vary. No cause is isolated
and no additional production optimization was selected from this run.

The reference is the **completed round-four optimized source**, whose package
manifest still says 0.3.1 and which already includes the new exact APIs. It is
not the published 0.3.2 package. The candidate adds the published 0.3.2 safety and
parity changes to that earlier work, plus the new bigint exact implementation;
these measurements build the default narrow Rust backend. Both sides use Rust
1.94.0 on Windows x86-64. The original workload sources are unchanged.

## Protocol and evidence

[PLAN.md](PLAN.md) was recorded before timing. Sources and workloads were frozen
in [manifest.json](manifest.json), with candidate/reference/driver archives.
[build.json](build.json) records all four binary hashes, compiler details, exact
Cargo manifests and lock contents; [binaries.zip](binaries.zip) preserves the
executables. Ropey 1.6.1, smallvec 1.16.4 and str_indices 0.4.4 are unchanged.
The path dependency version is the only difference in the two churn lockfiles.

There are three alternating fresh process pairs for each cohort, ordered
reference/candidate, candidate/reference, reference/candidate. The microbenchmark
retains all 42 cells and seven samples per process. Churn retains all 24 cells and
three samples per process for seeds 1, 2 and 3, invocation zero. Its 50 edits,
50% equal/50% final-byte-different targets, independent/shared histories and
10K/100K/1M sizes are unchanged. All 12 processes completed successfully with
their existing exact answer checks. All 2,196 timing samples are retained in
[raw/](raw/) and [results.json](results.json). [summary.json](summary.json) and
[all_cells.csv](all_cells.csv) retain every cell, all paired process medians and
ratios. The summarizer independently verifies raw hashes, sample counts, output
seeds, process order, all medians and the archived binary hashes. The measured
candidate snapshot still matches the working Rust source at analysis time.

The initial driver parser expected a textual churn row label; source inspection
corrected it to the original numeric seed before building or timing. The
preflight manifest and driver archive are retained separately. No workload or
timed result was replaced.

## Exact comparison results

Times are the median of the three per-process medians, in microseconds. The ratio
is the geometric mean of the three **paired** candidate/reference process-median
ratios; it need not equal the ratio of the displayed times. Ratios above one are
slower. Ranges are observed paired ratios, not confidence intervals.

| Cell | Reference µs | Candidate µs | Paired ratio | Paired ratio range | Slower pairs |
|---|---:|---:|---:|---:|---:|
| 64KiB independent equal: LCP | 3.075 | 3.370 | 1.088 | 1.064–1.113 | 3/3 |
| 64KiB independent equal: equality | 3.074 | 3.364 | 1.057 | 0.993–1.104 | 2/3 |
| 1MiB independent equal: LCP | 57.420 | 62.475 | **1.105** | 1.016–1.187 | **3/3** |
| 1MiB independent equal: equality | 58.703 | 56.688 | 0.983 | 0.944–1.009 | 1/3 |
| 1MiB late mismatch: LCP | 57.812 | 57.114 | 0.983 | 0.970–0.993 | 0/3 |
| 1MiB late mismatch: equality | 57.779 | 58.672 | 1.048 | 1.015–1.075 | 3/3 |
| 1MiB shared late mismatch: LCP | 0.125 | 0.115 | 0.954 | 0.920–1.000 | 0/3 |
| 2TiB Repeat: LCP | 0.0816 | 0.0803 | 0.979 | 0.910–1.033 | 1/3 |
| 2TiB Repeat: equality | 0.0845 | 0.0764 | 0.944 | 0.890–0.986 | 0/3 |
| 2TiB self-concat DAG: LCP | 0.252 | 0.340 | **1.308** | 1.219–1.361 | **3/3** |
| 2TiB self-concat DAG: equality | 0.279 | 0.282 | 0.994 | 0.903–1.075 | 2/3 |

The DAG LCP slowdown is about 88 ns in the displayed medians, while the represented
2TiB sequence remains compressed. Its repeatable direction merits disclosure even
though its absolute cost is small. Equality and LCP can share core traversal yet
show different changes in this fixed-sequence microbenchmark; that observation
does not identify an explanation or invalidate either measured cell.

## Edit plus exact comparison

Same aggregation as above; values are microseconds per edit-and-comparison cycle.
These are three paired processes, not the prior nine-pair benchmark cohort.

| Bytes | Target history | Reference µs | Candidate µs | Paired ratio | Paired ratio range | Slower pairs |
|---|---|---:|---:|---:|---:|---:|
| 10K | Independent | 4.648 | 4.734 | 1.059 | 0.974–1.165 | 2/3 |
| 10K | Shared | 3.752 | 3.890 | 1.057 | 0.911–1.198 | 2/3 |
| 100K | Independent | 19.450 | 20.246 | 1.045 | 1.030–1.064 | 3/3 |
| 100K | Shared | 7.632 | 8.000 | 1.010 | 0.892–1.088 | 2/3 |
| 1M | Independent | 150.186 | 157.636 | 1.018 | 0.917–1.095 | 2/3 |
| 1M | Shared | 9.068 | 9.816 | 1.080 | 0.963–1.169 | 2/3 |

All four original churn arms remain in the CSV. The fingerprint-only candidate
arms are unsafe as general equality decisions and are not substitutes for exact
comparison. Their 10K and 1M shared-history cells also slow in all three pairs
(paired ratios 1.142 and 1.142). The unchanged Ropey exact 100K shared-history
control slows in all three pairs as well (1.146); this demonstrates that build/run
context can affect the controls, but does not prove that any Hashrope slowdown is
noise. In the candidate run, Hashrope remains slower than Ropey for these small
independent edits, while shared history gives Hashrope an advantage at 1M bytes.

## Primitive costs and memory

| Cell | Reference | Candidate | Unit | Paired ratio |
|---|---:|---:|---|---:|
| Build 64KiB | 38.730 | 30.780 | µs | 0.906 |
| Balanced concat | 15.050 | 15.900 | ns | 1.086 |
| Split 64KiB | 433.000 | 418.000 | ns | 0.977 |
| Warm substring hash | 487.300 | 454.500 | ns | 0.951 |
| Root hash read | 0.5275 | 0.5288 | ns | 1.019 |
| Join repeated span | 6.755 | 5.244 | µs | 0.784 |

`NodeInner` remains 48 bytes, inline `Arena` remains 2,136 bytes, and constructing
64KiB allocates 255 nodes in both arms and all runs. These are layout/node-count
measurements, not peak or retained process memory. Sub-nanosecond/nanosecond cells
include loop/compiler effects and do not establish a transferable latency floor.

## Limits and reproduction

Only three process pairs were collected, without CPU pinning or measured
frequency/thermal control. The fixed original Rust sequence is preserved, but
the full historical matrix's **outer Python interleaving is omitted**. Binaries
were rebuilt, and both the source integration and compilation context differ.
This check does not resolve the previously reported full-matrix protocol
sensitivity or establish universal efficiency across packages, profiles, shapes
or machines. Bigint exact behavior is validated by dedicated correctness and
compressed-graph tests, not timed here. Independent byte comparisons can still
require a full scan; periodic and shared-structure shortcuts are not a universal
logarithmic bound.

The script refuses to overwrite snapshots or results. For an independent rerun,
use another checkout with the captured candidate sources, put `run.py` and
`PLAN.md` under the same release directory hierarchy, and capture into its fresh
study directory before building and running under an exclusive timing slot.
`python run.py freeze --historical <original-workspace>` is the original capture
command; `python run.py build`, `python run.py run`, and `python summarize.py`
record the subsequent phases. The original snapshots and raw measurements should
remain immutable.
