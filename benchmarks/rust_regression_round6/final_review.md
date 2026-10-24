# Independent final result and documentation review

The primary cohort, separate confirmation and current documentation pass review.
No material numerical, provenance, correctness or interpretive issue remains.
This review performed no builds, tests or timings.

[final_review.py](final_review.py) and [final_review.json](final_review.json)
independently verify the initial 36 processes, 6,588 samples, 66 paired statistic
vectors and CSV rows, four archived executables, frozen inputs and current source
identity. [confirmation_review.py](confirmation_review.py) and
[confirmation_review.json](confirmation_review.json) independently verify the
separate 24 processes, 4,392 samples and 66 vectors/CSV rows, including the exact
same four executable hashes. Both recompute medians, paired ratios, ranges,
direction counts, order/seed estimates and exploratory log-t intervals from raw
outputs. No phase is pooled or substituted for another.

All 66 rendered confirmation table rows match their summary. The study's
primary/confirmation tables, six candidate-cohort Hashrope/Ropey rows, root
README, Rust README, PARITY and CHANGELOG agree with the measured values. The
documentation distinguishes ratios of paired times from ratios of displayed
medians, and preserves the mixed results:

- Concat ratios are 0.8550 and 0.8777, with 8/9 and 5/6 pairs faster. The latter
  interval includes one; this supports a targeted measured improvement without
  establishing a universal bound.
- Independent 100K exact edit/comparison changes from 1.0438 (8/9 slower) to
  0.9834 (4/6 slower). Confirmation process medians nevertheless increase from
  18.113 to 18.406 microseconds. These different estimators are both retained;
  the initial negative observation is not dismissed as noise.
- Independent 1M edit/comparison ratios are 1.0138 and 0.9604. Independent 1 MiB
  LCP reverses from 0.9801 to 1.0680 (5/6 slower in confirmation), and the 64 KiB
  and late-mismatch LCP estimates also reverse. No confirmation slowdown interval
  lies wholly above one, which does not establish equivalence or zero cost.

The report correctly treats unchanged comparator instruction bodies as weaker
than identical linked placement, memory state or effects of changed edits.
Controlled slot/affinity/A-A evidence and the retained-original replication are
kept separate from the final optimization experiment. No timing cause or broad
performance-parity claim is established from these controls.

All seven final validation log hashes and exit codes were checked. Parsed test
totals are 221 passed, zero failed or ignored; all recorded current source/test
hashes match, including the six production files in the measured candidate.
The first validation failure remains preserved and is a test-only Clippy issue.
The production constructor refinement separately passed source review: checked
length then weight precede allocation, the profile guard and hash operands are
unchanged, and both unbalanced-helper callers establish the same balance
condition before entering the retained branch. Local Windows validation and
wasm compilation are properly distinguished from fresh cross-platform CI or
performance evidence. Versions remain 0.3.3 and the changes remain unreleased.
