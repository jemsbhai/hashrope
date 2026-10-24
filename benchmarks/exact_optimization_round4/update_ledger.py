"""Apply the completed round-four result to the paper claim ledger once."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
path = ROOT.parent / "hashropeictaipaper/CLAIMS.md"
text = path.read_text(encoding="utf-8")
old = """unchanged-code controls. T3/S3/B1 retain explicitly identified first-study
results; B2 records the latest third-round cohort. All 658 selected tests and
benchmark oracles pass; final core/test sources match the measured archive and
243 earlier benchmark source/artifact files remain byte-identical."""
new = """unchanged-code controls. The [fourth investigation](../hashrope/benchmarks/exact_optimization_round4/README.md)
selects no production change and repeats pre-third-round versus unchanged current
source. Its full short-protocol equality penalty recurs, despite the different
controlled/retained-binary estimates; the mechanism remains unresolved.
T3/S3/B1 retain explicitly identified earlier results; B2 records the latest
historical-reference/current retest. All fresh benchmark oracles pass. The 658
selected tests are carried forward by identity of all 20 core/Rust-test files;
they are not claimed as newly executed suites. All 376 earlier benchmark
source/artifact files remain byte-identical."""
assert old in text
text = text.replace(old, new, 1)
lines = text.splitlines()
for index, line in enumerate(lines):
    if line.startswith("| T3 |"):
        extra = (" Fourth-round controlled (36 processes) and retained-original-binary (12 pairs) diagnostics do not establish a stable primary penalty: "
                 "retained equality current/reference geometric ratio .9959 [.9458–1.0486], late LCP .9914 [.9326–1.0540] (exploratory per-cell 95% intervals). "
                 "But the complete short-protocol rerun again records higher equality, 54.926 [54.543–58.861]→65.895 [55.910–65.990] µs (+20.0%, all three paired rounds higher), "
                 "and late LCP 57.184 [54.632–57.449]→60.059 [57.426–65.809] µs (+5.0%), using median-of-three-round-medians statistics. "
                 "Rebuilt executable identities and outer Python interleaving differ from the retained-binary diagnostic, without an isolated cause. "
                 "No further production change, equivalence claim or dismissal as noise is justified. Unchanged Python slopes .963/.998 still fail ≤.3. "
                 "[Fourth report](../hashrope/benchmarks/exact_optimization_round4/README.md).")
        lines[index] = line[:-2] + extra + " |"
    if line.startswith("| B2 |"):
        lines[index] = (
            "| B2 | Fingerprint maintenance and exact edit-and-compare are different workloads | REFRAMED; measured shared-ancestry win | "
            "Historical EXP-020 369.7× at 1M and 2901.8× at 10M measured fingerprint retrieval/recomputation. "
            "The fourth-round full-matrix retest compares pre-third-round source against unchanged current source; no new production change is selected. "
            "It retains 50 ASCII edits, half equal targets and nine fresh processes/core. At 10K, independent exact edit-and-compare is "
            "7.428±.639→4.468±.237 µs (1.66×), shared 6.643±.525→3.703±.147 µs (1.79×), both 9/9; "
            "Ropey still wins at 1.169±.186/.912±.087 µs (remaining ratios 3.82×/4.06×). At 1M independent targets, "
            ".1534±.0147→.1442±.0073 ms (7/9), still 1.351× Ropey .1067±.0080. Shared .0147±.0008→.00917±.00052 ms (9/9), "
            "versus current-cohort Ropey .03417±.00269 (3.73× faster). Ropey's unchanged shared control itself moves .04405±.01071→.03417±.00269 ms. "
            "Both 100K cells improve in 9/9 and still lose to Ropey. The retained third-round diagnostics support faster leaf hashing and eliminated source-copy allocations; "
            "phases are not additive and both fingerprint-candidate arms use Hashrope from their measured core version. "
            "The separate same-arena short-protocol equality penalty recurs in the fourth round, with unresolved protocol/build sensitivity; see T3. "
            "Setup excluded, every boolean checked, no universal or statistical-significance claim. "
            "[Fourth-round retest](../hashrope/benchmarks/exact_optimization_round4/README.md); "
            "[third-round optimization and diagnostics](../hashrope/benchmarks/exact_optimization_round3/README.md); "
            "[earlier rounds](../hashrope/benchmarks/exact_optimization_round2/README.md). |")
path.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(path)
