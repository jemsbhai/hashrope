# Round-seven runtime and generated-code investigation

**Disposition:** the direction-argument trial below was measured and rejected.
Production retains the frozen round-six implementation. See the
[full report](../README.md) and [decision](../decision.json); fewer metadata
lookups did not establish a useful gain in the two primary workloads.

No timing was performed by these scripts. `inspect_final_binaries.py` reads the
four retained round-six final-comparison executables and their actual rlibs,
checks executable hashes against the original build manifest, and records the
tool/input hashes. `retained_round6/manifest.json` contains the results.

The final micro and churn builds both retain identical normalized instruction
bodies for `exact_lcp` (507 lines), the recursive comparator (546), frontier
expansion (69), and `split_inner`. Symbol hashes and anonymous data identifiers
are normalized; instruction offsets, branch targets and registers remain.
This is not proof of identical linked data targets or execution costs.

Exact 64-byte prologue signatures were located uniquely in each PE `.text`
section and cross-checked against x64 runtime-function entries. The linked
locations are:

| Cohort / function | Round-five baseline RVA | Round-six candidate RVA | Change |
|---|---:|---:|---:|
| Micro dispatcher `exact_lcp` | `0xA9B0` | `0xAB00` | +`0x150` |
| Micro recursive comparator | `0x6840` | `0x6840` | none |
| Micro frontier expansion | `0x6710` | `0x6710` | none |
| Churn dispatcher `exact_lcp` | `0x228B0` | `0x22A00` | +`0x150` |
| Churn recursive comparator | `0x1F0A0` | `0x1F0A0` | none |
| Churn frontier expansion | `0x1EF70` | `0x1EF70` | none |

Dispatcher alignment changes from 48 modulo 64 to zero; the other two function
locations do not change within each cohort. The 1 MiB same-arena micro traverses
the recursive comparator, while long independent cross-arena churn uses the
frontier in the dispatcher. These placement observations do not establish that
either timing difference was caused by placement. Runtime data placement,
process state, relocation targets and the full workload are not controlled by
this inspection. An initial extractor counted adjacent unwind funclets as part
of some bodies; that output is preserved under `retained_round6_pre_parser/`.
The corrected extractor stops at every unindented label.

## Narrow candidate: pass the known heavier side

Round-six `join` computes both weights to test balance, then its unbalanced
helper looks both weights up again solely to decide which side is heavier.
The original pre-round-six implementation also reloaded the two weight values,
but could retain node pointers and tags in its single function; the factored
helper has to inspect both nodes again. Thus repeated metadata work is real,
although this does not identify the cause of the earlier timing estimates.

Both callers already know the heavier side. The candidate passes one boolean
as the fourth Windows x64 register argument. The helper branches on that value
before inspecting only the selected node. It no longer loads two weights to
rediscover the direction. The full-width balance formula, public profile guard,
length-then-weight overflow preflight, construction/hash order and node layout
remain unchanged. Neither caller mutates the arena between computing weights
and passing their comparison to the helper.

`build_direction_probe.py` builds frozen-current and direction-only copies with
the same unchanged churn workload and locked dependencies; it never executes
them. `direction_build_manifest.json` records hashes. Normalized assembly and
diffs show that the four-argument helper still receives a tail jump from `join`.
The primary helper body drops from 247 to 241 instruction lines; the wrapper
grows from 82 to 91 and saves/restores two additional registers to keep the
comparison operands available across the full-width balance calculation.
Public concat grows from 226 to 235 lines. These are explicit tradeoffs, not a
claim that fewer metadata lookups guarantee lower latency. Comparator bodies
remain identical in this scratch build.

The exact inspected source was drafted as the production candidate, with SHA-256
`5931fcb08942a8ea887bb38ec25f57e651665583a39a42142af090f687f231d6`.
`direction.patch` isolates it from the frozen round-six source; the complete
source snapshot is `direction.rope.rs.txt`. It passed 27 selected release tests
(two private constructor/overflow tests, seven fingerprint/edit-layout goldens,
13 exact-comparison tests and five period/LCP regressions), plus the two private
tests in debug mode. Their complete logs are retained here. The parent task
owns full validation, measurement and the final acceptance decision. No exact
comparison, arithmetic, fingerprint default, public API or version change is
part of this candidate.
