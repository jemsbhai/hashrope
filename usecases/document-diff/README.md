# Use Case 1: Document Diff-Location

**Working title**: Exact Document Comparison with Persistent Hash Ropes

## Summary

Compare two document versions using exact byte equality and longest-common-prefix queries. `rope_equal` and `rope_lcp` can skip shared immutable ranges in versions produced through persistent edits, and can prove equality of repeated ranges without expanding every copy. Independently loaded documents may require scanning their common content. Locating every changed region also requires a diff or alignment algorithm; one LCP query only finds the first divergence.

## Core hashrope features exploited

- `rope_equal` / `rope_lcp` — exact comparisons that reuse shared structure and verify unshared bytes
- `split` / `concat` — persistent edits with incremental fingerprint maintenance; logarithmic tree traversal assumes bounded leaves, while boundary copying and arithmetic add costs
- `RepeatNode` — boilerplate/template clauses represented by a shared child and one repeat node
- Root and substring fingerprints — fast candidate filters, followed by exact verification before declaring content unchanged; matching fingerprints alone do not prove identity

## Key benchmark numbers (from FINDINGS.md)

These historical measurements cover fingerprint and storage primitives, not exact document differencing. See [consumer exactness reruns](../../benchmarks/consumer_exactness/) for the initial correction and [exact-path optimization measurements](../../benchmarks/exact_optimization/) for the current comparison paths. Exact costs depend on shared structure, common-prefix length, and leaf sizes.

- substr_hash: 124 ns warm at 10K leaves in the measured configuration
- concat (rejoin): 55–80 ns, 5,405× faster than rehash at 100K
- RepeatNode: 228 bytes for the measured child-plus-repeat representation; object sizes depend on implementation and version

## Subsumes

- Section 3 "Recursive Financial Statement Auditing"
- Section 9 "Template Verification"
- Section 9 "Financial OCR Deduplication" (partially — the dedup aspect may warrant its own experiments)

## Target venues

- ACM SIGMOD / VLDB (data management)
- EDBT (extending database technology)
- ACM DocEng (document engineering)
- USENIX ATC (systems)

## Proof-of-concept outline

1. Two versions of a PDF contract → extract text
2. Test both independently built chunked ropes and revisions derived through persistent edits
3. Use `rope_equal` and `rope_lcp` to verify equality and find divergence, then use a diff/alignment method to identify changed paragraphs or clauses
4. Report correctness against an independent byte/text oracle and compare end-to-end timings with a conventional exact diff; report construction and fingerprint maintenance separately

## Status

Planning phase.
