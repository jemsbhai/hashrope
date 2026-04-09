# Use Case 3: RAG Index Incremental Updates

**Working title**: Incremental Knowledge Base Integrity for Retrieval-Augmented Generation via Hash Ropes

## Summary

Maintain a hash rope over a chunked RAG corpus. When documents are added, modified, or deleted, update the integrity hash in O(log w) instead of rehashing the entire corpus. Use `substr_hash` to answer "what changed since last indexing?" queries without touching unchanged chunks.

## Core hashrope features exploited

- `split`/`concat` at O(log w) — replace a single chunk without rebuilding the index
- `substr_hash` at O(log w) — "has chunk range X–Y changed?" without reading those chunks
- Root hash at O(1) — instant corpus-level identity check
- Sliding window `append_copy` — process compressed documents without decompression (unique angle)

## Key benchmark numbers (from FINDINGS.md)

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
4. Benchmark: time to identify changed chunks via substr_hash binary search vs naive full-corpus rehash
5. Benchmark: time to update the integrity hash after chunk replacement (split + concat) vs rebuild
6. Demonstrate compressed-document ingestion via sliding window on gzip'd documents

## Status

Planning phase. Note: if experiments show this is primarily an engineering contribution rather than a research contribution, consider merging into the document-diff paper as an application section.
