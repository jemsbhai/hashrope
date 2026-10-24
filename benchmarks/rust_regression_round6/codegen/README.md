# Retained round-five code generation

`manifest.json` records the executable, rlib, tool and script hashes. All four
executable hashes were checked against the retained round-five build manifest.
The dumps come directly from those rlibs, without rebuilding or executing a
benchmark. The accompanying script refuses to overwrite its output directory.

For both the micro and churn builds, normalized instructions are identical for
`concat`, `join`, `make_internal`, `split_inner` and `ExactFrontier::expand`.
Normalization removes only Rust/LLVM symbol hashes and anonymous data symbol
identifiers. Relative offsets, registers and branch targets remain. This is not
a comparison of final linked placement or proof that anonymous data targets
have identical contents. The differing anonymous references in these unchanged
methods occur on cold panic paths.

Balanced concat still executes redundant work in both versions. Its public
preflight checks length and weight, calls `join`, which reloads weights and
checks balance, then tail-calls `make_internal`, which reloads and checks both
sums again. These are separate emitted functions, with separate prologues.
Removing that work can improve the absolute path; its presence in both versions
does not explain the historical difference by itself.

The long cross-arena frontier is inlined into `exact_lcp`. That emitted function
changes from 490 to 507 instruction lines, despite unchanged frontier source.
The stack reservation remains `0x6F8`. Register allocation and spills differ:
the candidate saves and reloads an additional live phase value around the leaf
`memcmp` call. The standalone frontier expansion helper is unchanged. These
facts justify a focused code-generation experiment, but do not establish the
cause or cost of the observed timing difference. The recursive exact comparator
also changes, as expected from the round-five aligned-child optimization.

To repeat with the retained local rlibs, run `inspect_retained.py --root <repo>
--output <new-directory> --dumpbin <MSVC-dumpbin.exe>`. Raw dumps, normalized
functions and diffs are preserved here. An initial scratch parser also matched
exception destructor symbols and stopped; the retained successful script
explicitly selects the primary Hashrope function symbols.

## Separate scratch frontier outline experiment

`build_outline_probe.py` copies the frozen round-five core twice, changes only
the `frontier_lcp` inline attribute to `#[inline(never)]` in one copy, and builds
the unchanged churn workload with locked dependencies. It never runs the
executables. `outline_build_manifest.json` binds those source, executable,
library and tool hashes. Both complete dumps and extracted functions remain.

The unchanged copy reproduces the 507-instruction dispatcher. Outlining reduces
that dispatcher to 135 instructions with `0x98` stack reservation and creates a
480-instruction frontier with `0x698` reservation. The total is 615 rather than
507 static instruction lines, and the path acquires a call. These observations
do not measure runtime or establish a gain. The outline was deferred; production
frontier source is unchanged. The fresh original-binary replication did not
reproduce a stable long independent churn penalty.

## Constructor draft refinement

The first checked-metadata constructor draft passed 42 selected release tests
and nine debug tests, with source and logs saved in `constructor_v1/`. Its
balanced public concat bypasses recursive join and duplicate overflow checks.
However, the compiler left its shared constructor tail outlined, so existing
recursive joins acquired an extra call from `make_internal` into that tail.
The selected refinement explicitly inlines the shared tail. This is a
code-generation requirement to preserve the existing call path; it does not
claim the tail's static body is small (121 instructions in the first draft).

`build_concat_probe.py` and `constructor_abi_build_manifest.json` preserve the
separate five-argument and three-argument helper builds. With the inline tail,
`make_internal` returns to 151 instructions and no longer calls another
constructor. The five-argument unbalanced helper still prevents a tail call
from `join` on Windows x64: its fifth argument needs a stack slot. Passing only
the arena and two node IDs, then reloading weights in the unbalanced helper,
changes that edge into a tail jump. The wrapper has 82 rather than 92 static
instructions; the unbalanced body has 247 rather than 241. The selected version
avoids the added call/frame, but retains separate code blocks and weight loads;
these assembly observations do not establish the runtime tradeoff.

The final production candidate and its full raw dump, normalized relevant
functions, source/library hashes and test logs are under `constructor_final/`.
It passed 42 selected release tests and nine debug tests. These cover full node
sequences across balanced/unbalanced trees and repeats, three hash profiles,
eager/lazy mode, length/weight rejection before allocation, the existing
published fingerprint/edit-layout goldens, and all targeted exact-comparison
regressions. No arithmetic, fingerprint defaults, metadata layout or exact
comparator was changed in this round. `concat-final.patch` isolates the source
change from the frozen round-five candidate. The later measured source has one
test-only Clippy cleanup: the reference node ID is stored directly instead of
wrapping it in `Some` and immediately unwrapping it. The runtime prefix is
unchanged. The [full validation](../VALIDATION.md) passes 221 release tests and
all tooling/platform checks; final timing is recorded in the parent report.

## Actual final-comparison churn builds

After the nine-pair final measurement, `inspect_final_churn.py` checked the
retained final-comparison churn executable hashes and disassembled their actual
rlibs. `final_churn_comparison/manifest.json` records the inputs and results.
The normalized `exact_lcp`, recursive comparator and frontier expansion bodies
are identical between the round-five baseline and the constructor candidate:
507, 546 and 69 instruction lines respectively. They also match the retained
round-five candidate bodies. Thus this constructor change did not alter those
comparator instruction bodies in these builds. This differs from the earlier
released-0.3.3 versus round-five comparison, where `exact_lcp` did change.

No benchmark was executed or rebuilt for this check. Instruction identity does
not establish identical linked placement, data addresses, cache state or whole
edit-and-comparison cost. In particular, it neither explains nor dismisses the
100K independent churn timing increase observed in the final nine-pair cohort.
