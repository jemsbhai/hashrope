# Downstream compatibility checks

Run the portable source compatibility gate from any directory with Python 3.10+
and Cargo on `PATH`. Git is required only when using a local source checkout.
The gate writes logs and `report.json` into the requested output directory.

## Public CI

```console
python packages/rust/scripts/verify_downstream.py --consumer hashrope-bio --output-dir downstream-results
```

The public gate downloads `jemsbhai/hashrope-bio` at commit
`206d1af783ca4c34626cfe1f9ce9a7bd7d4275cb` without credentials. It extracts only
selected regular source files into a fresh temporary directory and checks the
unchanged Rust library and all six binary targets through an isolated Cargo
manifest. Original consumer manifests and source checksums are retained in the
result artifacts. The temporary directory is removed after the run.

This command works on Linux, macOS, and Windows; it does not depend on a user's
checkout layout or shell. The crate's normal test matrix also exercises the
public node-layout and lazy-hash contracts through its synthetic downstream
contract tests. Those tests can run publicly without private consumer code.

## Authorized local cdh-sort gate

`jemsbhai/cdh-sort` was not retrievable through GitHub's unauthenticated repository,
commit, or archive endpoints when this gate was added. Public CI does not receive
its source or credentials. Use an already authorized local Git checkout:

```console
python packages/rust/scripts/verify_downstream.py --cdh-source /path/to/cdh-sort --output-dir downstream-results
```

On Windows, quote the normal Windows path. The candidate crate is located relative
to this script by default; `--candidate /path/to/hashrope/packages/rust` overrides
it. `--cargo /path/to/cargo` can select a Cargo executable explicitly.

The script exports exactly commit
`6d40877ba7675e2a3b41dda5ed524b7644ce2ddb` using `git archive`. It ignores dirty and
untracked files and never checks out or edits the original repository. In the
temporary copy it adds a crates.io patch for the candidate and updates only that
copy's lockfile to the candidate's current version. Cargo metadata then confirms
that the local candidate was resolved. Gates include:

- All unit, integration, and property tests, including compressed LZ77 construction,
  byte access, invariants, comparison, and sorting.
- Compilation of all targets, including downstream benchmarks and examples.
- Bounded lazy-construction benchmark smoke cases at compression ratios 10, 50,
  200, and 1000. These are correctness smoke runs, not performance measurements.

The earlier exploratory audit included preexisting dirty cdh-sort files. Its test
counts therefore differ from this clean pinned-commit gate. The reproducible gate
is the release evidence; it does not depend on those local edits.

`--consumer cdh-sort` runs only this consumer. The default `--consumer all` runs
both consumers and fails if a requested source is unavailable; it never silently
skips cdh-sort. `--bio-source /path/to/hashrope-bio` exports the pinned biology
commit locally instead of downloading it. `--offline` disables Cargo network
access but does not disable the public snapshot download; use both local source
options when source downloads must also be avoided.

## Existing biology limitations

The original hashrope-bio manifest requires `hashrope = "0.2.1"` with a stale local
path. That requirement excludes every 0.3.x version, including the candidate. The
harness selects the candidate explicitly to test unchanged consumer source; a
passing gate does not upgrade the real dependency or prove that Cargo will select
the new crate under the old requirement. Downstream adoption needs a separate
manifest change in hashrope-bio.

The pinned `gene_diff.rs` documentation example passes `&mut arena` twice in one
`diff_genes` call and fails Rust error E0499. The same failure was reproduced with
the published Hashrope 0.3.1 crate. The gate runs `cargo test --all-targets`, which
excludes doctests, and records this existing documentation issue in `report.json`.
It does not modify the example or treat that existing issue as a new regression.

## Evidence and failure behavior

`report.json` records source commits, archive checksums, individual source-file
checksums, the candidate version and source hashes, toolchain versions, consumer
results, and known limitations. Resolved Cargo lockfiles are retained alongside
the report. Sources are pinned; when an upstream lockfile is absent, Cargo resolves
the declared dependency ranges for that run. Cargo output is streamed to the console and saved in separate logs.
Every Cargo command must succeed. Changes to the candidate manifest or Rust source
during a run invalidate the result and require a rerun. Failed downloads, missing
pinned commits, wrong dependency resolution, and test failures return a nonzero
exit code. No source files are published, no consumer repository is modified, and
no registry credentials are needed for the public gate.

