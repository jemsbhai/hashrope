# Use Case 1: Document Diff-Location

**Working title**: Sub-logarithmic Document Differencing via Polynomial Hash Ropes

## Summary

Given two versions of a document (financial statement, legal contract, regulatory filing), locate the exact changed regions in O(log w) time using `substr_hash` binary search, vs O(n) full-text comparison.

## Core hashrope features exploited

- `substr_hash` at O(log w) — binary search to pinpoint changed regions
- `concat` at O(log w) — incremental hash maintenance as pages are ingested
- `RepeatNode` — boilerplate/template clauses stored in O(1) space regardless of repetition count
- Root hash comparison at O(1) — instant identity check before diving into diffs

## Key benchmark numbers (from FINDINGS.md)

- substr_hash: 124 ns warm at 10K leaves, O(log w) confirmed
- concat (rejoin): 55–80 ns, 5,405× faster than rehash at 100K
- RepeatNode: 228 bytes constant regardless of repetition count

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
2. Build hashrope from each version (chunked by paragraph or page)
3. Compare root hashes → if different, binary search via substr_hash
4. Report: exact paragraph/clause that changed, with timing comparison against naive diff

## Status

Planning phase.
