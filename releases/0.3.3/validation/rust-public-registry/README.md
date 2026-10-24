# Public Rust registry consumer validation

Both published crates passed consumer tests from a fresh external project and
empty Cargo home. The consumer pins `hashrope = "=0.3.3"` and
`hashrope-bio = "=0.3.3"`; it has no path dependencies or patches. Cargo metadata
confirms that every dependency came from the public crates.io registry. The
resolved lockfile, full command logs, test source and verifier are retained here.

- Default features: **5 tests passed**.
- Optional bigint feature: **7 tests passed**, including the same five tests.
- Core public archive SHA-256:
  `f830cd18595989dadcb9b8d39c10a3807c54b3f89bfd1f8355b63badfc9a88f1`.
- Bio public archive SHA-256:
  `45a1c594d2a13ff3fdab1efe896821da6f29687e474e316362a5fdfea9ef2e45`.

The tests exercise the default byte and encoded-token collisions, exact whole and
range comparisons, cross-profile comparisons, range bounds, frozen-profile
guards, and gene/exon/codon change detection. Bigint checks include independent
1,300-bit repeat lengths and shifted ranges without materialization or added
rope nodes. These are installation and contract smoke tests, not performance
measurements or a replacement for the full package suites.

## Archive provenance

The first attempt deliberately stopped before executing tests because the core
download did not match the precommit `cargo package` archive hash supplied to the
verifier. The initial metadata and evidence are preserved under
`initial-expected-hash-check/`.

Inspection found that Cargo publish had retained its actual upload under
`target/package/tmp-crate/`. The downloaded core and bio archives match those
actual local upload files byte for byte. Comparing every member of the reviewed
precommit core archive with the downloaded one found exactly one changed file:
`.cargo_vcs_info.json`. It changed from commit `bd8f9ca` with `dirty: true` to the
clean release commit `3dda4b5`. Every source, test, manifest and documentation
member is byte-identical. The precommit core archive hash was
`b9159ef8c3fcb668203b8e8561aa9202e9acf757a074dbdb82e1e62e28ebc136`.
`artifact-comparison.json` records this comparison and every published member
hash; the second fresh registry-only run used the verified actual-upload hashes.

The public dependency resolution selected num-bigint 0.4.8, num-integer 0.1.47,
num-traits 0.2.19 and autocfg 1.5.1. No local dependency source was substituted.

To repeat in a fresh temporary directory:

```text
python verify.py --core-sha256 f830cd18595989dadcb9b8d39c10a3807c54b3f89bfd1f8355b63badfc9a88f1 --bio-sha256 45a1c594d2a13ff3fdab1efe896821da6f29687e474e316362a5fdfea9ef2e45
```

`SHA256SUMS.json` covers the retained machine-readable evidence, commands and
consumer sources; this explanatory README is separate.
