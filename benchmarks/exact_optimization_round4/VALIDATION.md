# Comparison regression investigation validation

No round-four production change was selected. The completed-round-three working
state is preserved in `baseline.zip`; the measured historical reference is
separately staged in `reference.zip`. The frozen measurement plan identifies
both roles explicitly.

## Carried-forward correctness validation

All 20 current production and Rust integration-test files are byte-identical to
the completed-round-three tested sources, verified from both archives and the
working tree in [carried_validation.json](carried_validation.json). The following **658 selected tests
are carried forward from that study**, rather than presented as freshly run
round-four suites:

| Suite | Previously passed |
|---|---:|
| Python core | 203 |
| Rust core and integrations, release mode | 144 |
| Paper helpers and consumers | 255 |
| Biology Python | 41 |
| Biology Rust, including its compile-only doctest | 15 |
| **Total distinct selected tests** | **658** |

The six edit/fingerprint compatibility tests also passed in debug mode and
against the pre-round-three source; those repeated executions are not counted
twice. The release `no_std + alloc` check passed, and the core doctest command
succeeded with zero tests discovered. The 10 GB materialization experiment E3
and sliding-throughput experiment E5 remained explicitly filtered. See the
[retained validation and logs](../exact_optimization_round3/VALIDATION.md).

The fresh source audit in [unchanged_code.json](unchanged_code.json) compares
the selected source with the preserved completed-round-three working state:
98 existing Rust production functions, zero changes or additions, and no
changes to production text outside function bodies (including inline/codegen
attributes). Six complete modules and 34 Python definitions also match.
The measured pre-round-three reference intentionally differs in the previously
reviewed hash kernel and three leaf-copy paths. The final source verification
checks the working production and integration-test files against the measured
current archive.

## Newly executed diagnostic checks

The controlled comparison study validated ordinary byte oracles, compressed
construction oracles and a default-fingerprint collision check. Its 36-process
dataset retains all 3,456 timed blocks. The separate original-binary replication
verified both retained executable hashes against their original manifests,
validated the unchanged benchmark sequence and retained all 7,056 samples from
24 processes. No samples were dropped. The independent audit checked raw file
hashes, answers and paired estimates; see [REVIEW.md](REVIEW.md) and
[controlled_findings.md](controlled_findings.md).

These runtime checks do not replace the broader carried-forward correctness
suites. Separately emitted normalized instruction text for seven comparator
routines matched, with the dispatcher section-ordinal difference retained in
the audit. That check is not a proof of final linked-binary identity or of equal
runtime in every environment.

## Final matrix

The full matrix completed using frozen historical-reference and unchanged-current
sources, the same workloads and process/sample counts as round three, and new
output paths: 18 Rust, 18 paper CPU and 18 real-pair/K workers, two full-corpus
biology workers, plus three alternating Python/Rust micro rounds with seven
samples each. All byte/token/per-edit oracles and all five context audits passed.
`summary.json` records `complete_matrix: true`. The full short-protocol equality
median is higher again; correctness completion does not imply performance
equivalence or dismiss that result.

All 20 measured core/Rust-test files match current working files byte for byte;
376 earlier benchmark source/artifact files remain unchanged. Alternate-output
aggregation reproduced the same summary, and verified biology resume neither
reran timing nor rewrote results. All new Python scripts parsed and root, paper
and biology `git diff --check` passed with line-ending warnings only. The
machine-readable checks are in `final_audit.json`.

The independent final review found no material numeric or interpretive issue
in the full report, paper, claim ledger or findings. The desktop LaTeX compiler
reported success for the current revised `paper/main.tex`. The paper LOGBOOK
completion entry preserves every preceding byte; `logbook_completion.json`
records that check. No measured production source changed after freeze.
