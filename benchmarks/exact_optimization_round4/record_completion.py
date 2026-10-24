"""Append the completed, reviewed study record without altering earlier log bytes."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1].parent / "hashropeictaipaper"
path = PAPER / "LOGBOOK.md"
before = path.read_bytes()
marker = "### Round-four completion: protocol-sensitive comparison cost, no new code"
assert marker.encode() not in before
entry = """

### Round-four completion: protocol-sensitive comparison cost, no new code

**Completed and independently reviewed, 2026-10-24.** The frozen final plan
remained pre-round-three reference versus unchanged completed-third-round
source. No production code, fingerprint default, profile behavior or layout
changed. The separate preservation baseline remains intact.

All planned cohorts completed sequentially under the team CPU timing pause:
18 Rust edit-and-compare workers, 18 paper CPU workers, 18 real-pair/K workers,
two full-corpus HIV workers, and the original three alternating Python/Rust
micro rounds with seven samples each. All byte/token/per-edit oracles passed,
including 7,200 real-pair and 90 K-cell checks across cores. All five context
audits passed. The original nine fresh-process design remains distinct from
biology's nine within-process pairs and micro's three round medians.

The final evidence is protocol-sensitive. The 36-process controlled diagnostic
(3,456 timed blocks) and fixed 12-pair retained-original-executable replication
(7,056 raw samples) did not establish the original stable primary penalty.
Retained-executable paired geometric current/reference ratios were equality
.9959 [.9458–1.0486] and late LCP .9914 [.9326–1.0540], with exploratory per-cell
95% intervals. These do not replace the historical median-of-three statistic.

The full short-protocol rerun again records a slower equality median:
54.926 [54.543–58.861]→65.895 [55.910–65.990] µs (+20.0%), with all three paired
round medians higher. Late LCP is 57.184 [54.632–57.449]→60.059 [57.426–65.809]
µs (+5.0%). Secondary 64 KiB late LCP/equality also rise 19.9%/29.1%, with
nonoverlapping observed round ranges. The full matrix rebuilds the executables
and interleaves the original outer Python/Rust rounds; both linked build context
and preceding execution differ from the retained-binary diagnostic. Neither is
isolated as the cause. No equivalence, universal no-regression, noise dismissal
or validated new production remedy is claimed. Earlier measurements and the
pre-final decision record remain unchanged.

The existing edit/hash improvement repeats: 10K independent/shared churn is
7.428±.639→4.468±.237 µs and 6.643±.525→3.703±.147 µs (1.66×/1.79×, both 9/9),
but still 3.82×/4.06× current-cohort Ropey. At 1M independent targets the remaining
Hashrope/Ropey ratio is 1.351; shared-history Hashrope takes 26.8% of Ropey's
current-cohort time. Unchanged Ropey shared control also varies between cohorts.
The full report retains all cells and uncertainty, rather than treating these
fresh repetitions as a fourth optimization. DAG/Repeat remain symbolic and
inline storage layouts remain unchanged, without a total-memory claim.

Unchanged Python controls retain the earlier limits: EXP005 slopes .963/.998
both fail the historical ≤.3 gate. Exact beats radix in all nine 2M-token runs
but remains slower than NumPy; the dispersion-band crossover is the smallest
grid point, with a vacuous lower condition. Real pairs favor exact over radix,
while 64K/K100 favors indexed radix. R2's 1M exact query excludes 999±35 ms of
query construction. HIV median batches vary with no Rust attribution; adaptive
inner counts match 4/2/1. No GPU serving, TTFT, energy or chromosome-scale work
was rerun.

The 658 selected tests are explicitly carried forward from round three, not
reported as freshly executed. Both immutable archives and the current working
tree are byte-identical for all 20 production/Rust integration-test files.
The fresh source audit checks 98 existing Rust production functions, zero
changes/additions and unchanged outside-function attributes/text. All 376
previous benchmark source/artifact files remain byte-identical. Final source,
context/harness/data, locked dependency and binary identities pass verification.
Alternate-output aggregation reproduces the same summary; verified biology
resume neither reruns timing nor rewrites results. New Python scripts parse;
root, paper and biology diff checks pass with line-ending warnings only.

The independent review checked all 490 aggregate statistic vectors, 510 raw
micro medians, 95 rendered result rows, five context hashes, 57 worker bindings,
31 source/helper hashes and 20 tested source files, then reviewed paper/ledger/
findings and both library READMEs. It found no material numeric or interpretive
issue. The desktop LaTeX compiler reported success for the current revised
`paper/main.tex`. No production or measured-source edits followed the freeze.

Artifacts: `../hashrope/benchmarks/exact_optimization_round4/README.md`,
`summary.json`, `measurement_plan.json`, `controlled_findings.md`, `REVIEW.md`,
`VALIDATION.md`, `final_audit.json`, `carried_validation.json`, source archives,
and all raw results. The generated summary records `complete_matrix: true`.
"""
path.write_bytes(before + entry.encode("utf-8"))
assert path.read_bytes().startswith(before)
record = {"logbook": str(path), "prior_bytes": len(before),
          "prior_sha256": hashlib.sha256(before).hexdigest(),
          "final_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
          "prior_bytes_preserved": True}
(HERE / "logbook_completion.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
audit_path = HERE / "final_audit.json"
audit = json.loads(audit_path.read_text())
audit["logbook_prior_bytes_preserved"] = True
audit["final_paper_compilation"] = "Desktop compile_latex_document succeeded for current revised paper/main.tex; reported by root and recorded in independent review"
audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
print("Completion recorded; all preceding LOGBOOK bytes preserved")
