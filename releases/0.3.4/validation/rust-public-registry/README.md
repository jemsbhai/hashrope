# Public Rust 0.3.4 verification

Both crates were published after all 14 coordinated CI jobs passed. Public
archives match the clean reviewed packages byte for byte:

| Package | SHA-256 |
|---|---|
| hashrope 0.3.4 | `2c9b71b750b9f9bdf9ec291bb5a5d614e7a6df1d1f39d13f10f9f171864053a3` |
| hashrope-bio 0.3.4 | `beeae45eafa54b6e2856bd66839c4521159a1ac0e88f890290c3a2317a406e11` |

The consumer was copied into a fresh temporary directory with an empty Cargo
home and target directory. Every dependency resolved from crates.io; no local
paths or patches were used. The downloaded archive checksums match the table.
All five default-mode and seven all-feature tests passed, including deliberate
fingerprint collisions, profile guards, cross-profile equality, range bounds,
compressed repeats, biology exon/codon comparisons and huge bigint repeats.

See [evidence.json](evidence.json), [default tests](default-tests.log),
[bigint tests](bigint-tests.log), [resolved dependencies](dependency-tree.log)
and [Cargo.lock](Cargo.lock). The publication records preserve package and dry-run
checks, upload commands and independent public archive downloads for
[core](../rust-core-publication/results.json) and
[biology](../rust-bio-publication/results.json).

`verify.py` reproduces the registry-only verification using the reviewed archive
hashes as arguments. These are correctness checks, not new timing measurements.
