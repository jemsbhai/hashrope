# Use Case 3: RAG Index Incremental Updates

**Working title**: Incremental Corpus Fingerprints and Exact Reuse for Retrieval-Augmented Generation

## Summary

Maintain a persistent hash rope over a chunked RAG corpus. Edits update its noncryptographic fingerprint without rehashing unchanged subtrees. Before reusing indexed content, use exact byte comparison: equal root or substring fingerprints are only candidate matches. Revisions derived through persistent edits can skip shared immutable ranges; independently ingested equal content may need a linear scan.

## Core hashrope features exploited

- `split` / `concat` — replace a chunk while sharing unchanged nodes; logarithmic tree traversal assumes bounded leaves, with additional work for new bytes, leaf boundaries, and arithmetic
- `substr_hash` — compute range fingerprints from maintained metadata, then verify candidate matches before reuse
- `rope_equal` / `rope_lcp` — exact content comparison that can skip shared immutable ranges and verified periodic spans
- Root hash at O(1) — a candidate fingerprint, not a corpus identity guarantee
- Sliding-window `append_copy` — consume decompressor back-references while maintaining fingerprints; decoding the compressed format remains a separate task

## Key benchmark numbers (from FINDINGS.md)

These historical primitive timings do not measure exact changed-chunk detection or a complete RAG indexing pipeline. See [consumer exactness reruns](../../benchmarks/consumer_exactness/) for the initial correction and [exact-path optimization measurements](../../benchmarks/exact_optimization/) for the current comparison paths, including the cost of inspecting unshared content.

- concat (rejoin): 55–80 ns, 5,405× faster than rehash at 100K
- substr_hash: 124 ns warm at 10K leaves
- Sliding window post-arena: 91–151 MB/s

## Relationship to other use cases

- Could potentially be subsumed into the document-diff paper if the RAG-specific experiments don't warrant a standalone contribution
- If standalone: the compressed-document ingestion via sliding window is the differentiating angle from document-diff
- The "KV Cache Optimization" and "Federated Learning" use cases from the master document are weaker fits and should be discussed as speculative future work only

## Target venues

- ACL / EMNLP (NLP — if framed around RAG pipeline efficiency)
- NeurIPS / ICML (ML systems — if framed around training/serving infrastructure)
- MLSys (ML systems conference)
- SIGIR (information retrieval)

## Proof-of-concept outline

1. Build a corpus of ~10K–100K text chunks (e.g., Wikipedia articles chunked at 512 tokens)
2. Build hashrope over the corpus
3. Simulate document updates: modify 1%, 5%, 10% of chunks
4. Benchmark exact changed-chunk identification against a byte-comparison oracle, separately for persistent revisions and independently rebuilt inputs
5. Benchmark incremental fingerprint updates after chunk replacement (`split` + `concat`) against recomputation; keep this separate from exact reuse decisions
6. Demonstrate compressed-document ingestion through a decoder integrated with the sliding window

## Status

Planning phase. Note: if experiments show this is primarily an engineering contribution rather than a research contribution, consider merging into the document-diff paper as an application section.
