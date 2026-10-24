# Independent review of the round-seven candidate decision

**Rejecting the direction-flag candidate and restoring the frozen round-six
runtime is a defensible conservative decision.** The fixed full-workload cohort
does not establish improvement in either predeclared target, and it adds measured
secondary costs. This is not proof that the candidate is slower everywhere, nor
does restoration fix or explain the previously observed 100K or 1 MiB costs.

[final_review.py](final_review.py) and [final_review.json](final_review.json)
independently reparse and verify all **48 processes, 8,784 raw samples, 66 paired
statistic vectors and 66 CSV rows**. They check all four archived executable
hashes, source/driver/lock archives, copied build inputs, process status,
arguments, order, seeds, sidecars, raw hashes, medians, paired ratios, ranges,
direction counts, order/seed strata and exploratory log-t intervals. All checks
passed. This audit ran no builds or timings.

The analysis preserves the distinction between source at measurement and source
after selection. The producer summary records the measured working candidate;
the independent audit observed restoration already completed and verified that
**all six production files byte-match the frozen round-six baseline**. Both
candidate and baseline archives, executables and all trial observations remain.

## What the candidate comparison establishes

Ratios below are candidate/baseline geometric means of 12 paired process
medians, with exploratory 95% log-t intervals. They are unadjusted for multiple
cells and are not equivalence bounds.

| Case | Paired ratio [interval] | Candidate slower |
|---|---:|---:|
| Primary: independent 100K exact edit/comparison | 0.9898 [0.9479, 1.0336] | 5/12 |
| Primary: independent-equal 1 MiB exact LCP | 0.9964 [0.9602, 1.0340] | 6/12 |
| Early-mismatch 64 KiB exact LCP | 1.0491 [1.0121, 1.0874] | 10/12 |
| Early-mismatch 1 MiB exact equality | 1.0394 [1.0054, 1.0745] | 9/12 |
| Shared-history 100K exact edit/comparison | 1.0419 [0.9873, 1.0995] | 7/12 |
| Late-mismatch 1 MiB exact LCP | 0.9707 [0.9446, 0.9976] | 1/12 |

The two primary point estimates are only about 1.0% and 0.4% lower, and their
execution-order strata have mixed directions. Their intervals permit both
improvement and increased time. That supports saying **benefit was not
established**, rather than claiming the targets are equivalent or unchanged.

The early-mismatch medians are 38.025 to 40.750 ns for 64 KiB LCP and 58.225
to 59.550 ns for 1 MiB equality. The absolute changes are small, but their
direction counts and intervals remain relevant secondary costs. The shared
100K result is less conclusive and is also retained. Other cells include gains:
late-mismatch 1 MiB exact LCP is about 2.9% lower, and the 1M independent
fingerprint-candidate edit control is about 7.2% lower. Several hash-only controls
are slower. None of these secondary results is removed or substituted for the
predeclared targets; all 66 cells remain in the archived comparison.

## Mechanism and disposition

The source change correctly reused the known `wl > wr` direction and removed
two redundant helper weight lookups. Source review found no profile, overflow,
hash, byte, node-sequence or exactness defect. The inspected Windows call path
retained a tail jump, but the wrapper grew from 82 to 91 static instructions and
added register-save work. Fewer source-level lookups do not guarantee faster
machine execution. The measurements do not isolate that code-generation change
as the cause of any particular cell's increase.

Given no established target benefit, secondary negatives and extra factoring
complexity, restoring the already validated baseline avoids retaining an
unjustified optimization. That selection uses the fixed cohort as declared;
there is no selective rerun, sample exclusion or favorable pooled replacement.
The completed restored-source validation was independently checked below.

The retained baseline still has unresolved context-sensitive performance costs:
the separate controlled study shows a roughly 2% pinned 1 MiB LCP increase and
a positive secondary 100K comparison-only accumulated-time estimate. A/A
variation, identical comparator instruction bodies or a failed optimization
cannot dismiss those observations as noise. This round yields a stronger
diagnosis and a documented rejected candidate, not a new production speedup,
universal no-regression claim, version bump or publication.

## Final documentation and restored-source checks

[documentation_review.py](documentation_review.py) and its
[saved result](documentation_review.json) verify all 66 rendered final-table
cells, all 12 headline trial rows and all six retained-baseline Ropey rows.
The latter ratios use paired baseline-arm process medians, not rejected-candidate
measurements or the ratio of the displayed aggregate times. Numeric formatting,
intervals, slower counts, ranges and execution-order estimates match the saved
summary. The fingerprint-candidate percentage is correctly rounded to -7.2%.

Both validation records independently pass seven checks and contain 221 passing
tests, zero failed tests and zero ignored tests. Every command log, source
archive and archived validation-driver hash matches its record. All 25 restored
source, fixture and lock files match their recorded hashes; all six production
files byte-match the frozen round-six baseline. The candidate and restored
validation archives remain distinct and intact.

The reviewed main report, library documentation, parity notes and changelog
consistently distinguish the frozen round-five/round-six diagnostics, rejected
direction-argument trial, and retained round-six production implementation.
The unresolved pinned LCP and comparison-only costs are preserved rather than
called noise or fixed. No material numerical, provenance or interpretive issue
remains. This final audit ran only lightweight reads and arithmetic; it did not
rerun builds, tests or timings.
