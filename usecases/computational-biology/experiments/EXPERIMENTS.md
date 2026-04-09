# Computational Biology — Experiment Protocol

**Version**: 0.1 (draft)
**Date**: 2026-04-08
**Status**: Design phase — no results collected yet
**Depends on**: hashrope Python package (PyPI) and/or Rust crate

---

## 1. Claims

Each claim maps to one or more experiments. Claims are scoped to what
hashrope's exact-identity primitives can actually deliver — no overclaiming.

| ID | Claim | Experiments |
|----|-------|-------------|
| CB-C1 | substr_hash on a genome-scale rope (~3.2 Gbp) answers arbitrary region identity queries in O(log w) time, substantially faster than reading + hashing the region bytes | E-G1 |
| CB-C2 | RepeatNode compresses known tandem repeat loci (STRs/VNTRs) to O(1) space and compares repeat counts in O(1) time | E-G2 |
| CB-C3 | Binary search via substr_hash localizes a single-nucleotide mutation in an N-bp region in O(log N) hash comparisons | E-G3 |
| CB-C4 | Rope construction from chunked FASTA is a practical one-time cost that amortizes over many queries | E-G4 |
| CB-C5 | substr_hash on serialized protein structures enables O(log w) frame-to-frame identity testing on MD trajectories | E-P1 |
| CB-C6 | RepeatNode detects and compresses periodic conformations in MD trajectories | E-P2 |
| CB-C7 | Hash comparison on serialized SMILES/InChI strings provides O(1) exact compound lookup in a chemical database | E-D1 |
| CB-C8 | Hashrope-based drug resistance panel checks N known mutation sites in O(N · log w) total time | E-D2 |

---

## 2. Baselines

| Baseline | Tool/Method | Why |
|----------|-------------|-----|
| samtools faidx | samtools 1.x | Standard random-access on indexed FASTA; extracts region then hashes |
| Direct memcmp | hand-rolled | Read both regions into memory, compare byte-by-byte |
| SHA-256 region hash | sha2 crate / hashlib | Cryptographic hash of extracted region bytes — apples-to-apples on hash quality |
| BLAST+ | NCBI BLAST 2.x | Standard sequence similarity search — NOT a direct competitor (approximate vs exact) but included to contextualize speed |
| MDAnalysis | Python MDAnalysis | Standard MD trajectory analysis; frame extraction + comparison |
| RDKit | Python RDKit | Cheminformatics toolkit; canonical SMILES generation and comparison |
| Naive Python dict | hand-rolled | Hash each SMILES with Python's built-in hash, store in dict for O(1) lookup |

**Important framing note**: BLAST/HMMer are *approximate* similarity tools. We do NOT claim hashrope replaces them. Baselines are chosen to compare against tools that answer the *same question* (exact identity), not different questions (approximate similarity). Any paper must make this distinction explicitly.

---

## 3. Datasets

### 3.1 Genomics

| Dataset | Source | Size | Purpose |
|---------|--------|------|---------|
| GRCh38 (human reference genome) | NCBI/Ensembl | ~3.1 GB (FASTA) | Primary genome-scale benchmark |
| GRCh38 chr22 | NCBI | ~51 MB | Fast iteration during development |
| GRCh37 (hg19) | NCBI | ~3.1 GB | Version comparison (GRCh37 vs GRCh38 diff localization) |
| Pan troglodytes (chimp) reference | Ensembl | ~3.0 GB | Cross-species divergence localization |
| ClinVar VCF | NCBI ClinVar | ~100 MB | Known pathogenic variants — ground truth for mutation localization |
| UCSC STR catalog | UCSC Genome Browser | ~50 MB | Known tandem repeat loci — ground truth for RepeatNode compression |
| HIV-1 HXB2 reference | NCBI | 9.7 Kbp | Drug resistance mutation panel target |

### 3.2 Protein Structure / MD

| Dataset | Source | Size | Purpose |
|---------|--------|------|---------|
| PDB subset (top 1000 by resolution) | RCSB PDB | ~2 GB (PDB format) | Structural comparison benchmarks |
| AlphaFold human proteome | AlphaFold DB | ~12 GB (mmCIF) | Large-scale structural database integrity |
| MD trajectory — lysozyme in water | GROMACS tutorial / MDAnalysis test data | ~500 MB (.xtc) | Trajectory indexing benchmark |
| MD trajectory — protein folding (Trp-cage) | D.E. Shaw Research / MDDB | ~5 GB | Periodic conformation detection |
| Alanine dipeptide trajectory | Standard benchmark | ~50 MB | Small-scale validation of frame comparison |

### 3.3 Drug Discovery / Cheminformatics

| Dataset | Source | Size | Purpose |
|---------|--------|------|---------|
| PubChem canonical SMILES (subset) | PubChem FTP | ~1M compounds, ~200 MB | Compound lookup benchmark |
| ChEMBL approved drugs | ChEMBL | ~2,500 compounds | Small focused set for drug resistance panel |
| HIV Drug Resistance Database | Stanford HIVDB | ~5,000 sequences | Resistance mutation ground truth |
| WHO catalog of TB resistance mutations | WHO | ~1,500 mutations | Second resistance panel target |

### 3.4 Reproducibility

- All random data generated with seeded RNGs (seed=42)
- Genome files identified by MD5/SHA-256 checksum logged in results
- PDB/trajectory files identified by PDB ID or DOI
- Exact software versions recorded (samtools, BLAST, MDAnalysis, RDKit)

---

## 4. Experiments — Phase 1: Genomics

### E-G1: Genome-Scale substr_hash Query Performance

**Claim**: CB-C1
**Question**: How fast can hashrope answer "is region [start, start+L) identical between two ropes?" compared to reading + hashing the bytes?

**Method**:
1. Load GRCh38 chromosome 22 (~51 MB) as chunked FASTA
2. Build hashrope with chunk sizes c ∈ {256, 1024, 4096, 16384} bytes
3. For query region sizes L ∈ {100, 500, 1000, 5000, 10000, 50000, 100000} bp:
   a. Generate 1000 random (start, L) pairs within chromosome bounds
   b. Time `substr_hash(rope, start, L)` for each query — measure median, p5, p95
   c. Baseline 1: `samtools faidx chr22.fa chr22:{start}-{start+L}` → pipe to SHA-256
   d. Baseline 2: read region bytes from mmap'd file → `PolynomialHash.hash(bytes)`
   e. Baseline 3: read region bytes → memcmp against second copy
4. Repeat on full GRCh38 (~3.1 GB) with chunk size 4096 only

**Output**:
- Table: query time (ns) vs region size L, for each method
- Log-log plot: time vs L, showing O(log w) for hashrope vs O(L) for baselines
- Crossover point: at what L does hashrope become faster than direct read+hash?

**Expected outcome**: 
- hashrope: ~200–400 ns per query regardless of L (O(log w), not O(L))
- samtools: dominated by I/O, likely 10–100 µs per query
- Direct read+hash: O(L), crossing over with hashrope around L ≈ 500–2000 bp
- The advantage grows with L. At L=100Kbp, expect 100–1000× speedup.

**One-time cost to report**: rope construction time for the full genome. At ~215 MB/s (from FINDINGS.md E1), expect ~14 seconds for 3.1 GB. This amortizes over thousands of queries.

---

### E-G2: Tandem Repeat Compression via RepeatNode

**Claim**: CB-C2
**Question**: How effectively does RepeatNode compress known tandem repeat loci, and how fast is repeat-count comparison?

**Method**:
1. Download UCSC STR catalog for GRCh38 — extract loci with known motif and repeat count
2. For each locus (motif M of length d, repeat count q):
   a. Build `RepeatNode(Leaf(M), q)` — measure construction time and memory (node count × 64 bytes)
   b. Build naive rope: materialize M×q into a single Leaf — measure memory
   c. Compute compression ratio: naive_bytes / RepeatNode_bytes
3. Repeat-count comparison benchmark:
   a. For q1 ∈ {10, 20, 30, 40, 50} and q2 = q1 ± {0, 1, 5, 10}:
   b. Build RepeatNode(M, q1) and RepeatNode(M, q2)
   c. Time: hash comparison (O(1))
   d. Baseline: materialize both, compare byte-by-byte
4. Clinically relevant loci: benchmark on HTT (Huntington's CAG repeat, normal: 10–35, pathogenic: 36+),
   FMR1 (Fragile X CGG repeat, normal: 5–44, pathogenic: 200+), ATXN1, ATXN3, DMPK

**Output**:
- Table: locus, motif, q, RepeatNode bytes, naive bytes, compression ratio
- Table: repeat-count comparison time — hash vs materialize+compare
- Box plot: compression ratios across all UCSC STR catalog loci

**Expected outcome**:
- RepeatNode: 228 bytes constant regardless of q (from FINDINGS.md E6)
- Pathogenic HTT (q=70, motif=CAG): naive = 210 bytes, RepeatNode = 228 bytes — comparable. But FMR1 (q=800, motif=CGG): naive = 2,400 bytes, RepeatNode = 228 bytes — 10.5× savings.
- The real win is at large q. FMR1 full mutation (q=2000+): naive = 6KB+, RepeatNode = 228 bytes.
- Hash comparison: O(1) at ~1 ns vs O(q·d) for byte comparison.

---

### E-G3: Mutation Localization via Binary Search

**Claim**: CB-C3
**Question**: Given that a sample differs from the reference somewhere in a region of N bp, how many substr_hash comparisons does it take to pinpoint the mutation, and how long does this take?

**Method**:
1. Build reference rope from chr22 (4096-byte chunks)
2. Create mutant copy: introduce a single synthetic SNP at a known position
3. Build mutant rope from modified sequence
4. Binary search algorithm:
   ```
   lo, hi = 0, N
   while hi - lo > 1:
       mid = (lo + hi) // 2
       if substr_hash(ref, lo, mid-lo) == substr_hash(mut, lo, mid-lo):
           lo = mid
       else:
           hi = mid
   return lo  # mutation position
   ```
5. For region sizes N ∈ {1000, 10000, 100000, 1000000, 10000000, 51000000 (full chr22)}:
   a. Place SNP at random position within region
   b. Run binary search — count comparisons, measure total wall-clock time
   c. Verify: reported position matches ground truth
   d. Repeat 100 times per N (different random positions)
   e. Baseline: linear scan with per-base comparison
6. Extend to multiple mutations: introduce k ∈ {1, 2, 5, 10} SNPs, measure time to locate all

**Output**:
- Table: N, comparisons (should be ⌈log₂ N⌉), total time, time per comparison
- Comparison vs linear scan time
- Plot: localization time vs region size (log-log), expect slope ≈ 0 (log(log w))
- Multi-mutation: time vs k, expect O(k · log N · log w)

**Expected outcome**:
- For N = 51M (full chr22): ⌈log₂(51M)⌉ = 26 comparisons × ~300 ns = ~8 µs total
- Linear scan baseline: O(N) = ~51M byte comparisons = ~50 ms
- Speedup: ~6,000× for single mutation localization on chr22

---

### E-G4: Rope Construction Cost and Amortization

**Claim**: CB-C4
**Question**: What is the one-time cost of building a genome-scale rope, and how many queries amortize it?

**Method**:
1. Load GRCh38 FASTA, parse into chunks of size c ∈ {256, 1024, 4096, 16384}
2. Time: `Arena::new()` + sequential `from_bytes` + `concat` for all chunks
3. Measure: wall-clock time, peak memory, final node count, tree height
4. Compute: amortization point = construction_time / per_query_savings
   (where per_query_savings = baseline_query_time - hashrope_query_time)
5. Also time: loading the pre-built rope from a serialized format (if we implement serialization)

**Output**:
- Table: chunk size, construction time, peak memory, node count, height
- Amortization curve: total time (construction + N queries) vs N, for hashrope vs baseline
- Crossover point: N* where hashrope total time < baseline total time

**Expected outcome**:
- At 4096-byte chunks: ~800K leaves, construction ~14s (at ~215 MB/s)
- Per-query savings: ~10–100 µs (depending on region size)
- Amortization at ~10 µs savings per query: N* ≈ 1.4M queries
- At ~100 µs savings per query (large regions): N* ≈ 140K queries
- For a genomics pipeline processing thousands of samples against one reference, this amortizes easily.

---

## 5. Experiments — Phase 2: Protein Structure

### E-P1: MD Trajectory Frame Comparison

**Claim**: CB-C5
**Question**: Can substr_hash detect identical frames in an MD trajectory faster than reading + comparing the frame data?

**Method**:
1. Load alanine dipeptide trajectory (~50 MB, ~10K frames) via MDAnalysis
2. Serialize each frame: flatten atomic coordinates to bytes (float32 × 3 × N_atoms)
3. Build hashrope: each frame is a leaf, concatenated sequentially
4. For frame pairs (i, j) where i ∈ {0, 100, 500, 1000, 5000} and j sweeps all frames:
   a. Time: `substr_hash(rope, i*frame_size, frame_size) == substr_hash(rope, j*frame_size, frame_size)`
   b. Baseline: extract frames via MDAnalysis, compute np.array_equal or RMSD
5. Scale up: lysozyme trajectory (~500 MB, ~50K frames)
6. Scale up: Trp-cage folding trajectory (~5 GB, ~500K frames) — if available

**Output**:
- Table: trajectory size, frames, per-comparison time (hashrope vs MDAnalysis extract+compare)
- Heatmap: frame-vs-frame identity matrix (computed via substr_hash) — visualize periodic returns
- Time to compute full N×N identity matrix: hashrope vs naive

**Expected outcome**:
- hashrope per-comparison: ~200–400 ns (O(log w))
- MDAnalysis extract + compare: ~10–100 µs (I/O dominated)
- Full 10K×10K matrix: hashrope ~20 ms vs MDAnalysis ~1000s (many hours)
- The early-exit optimization (2.1 ns when query covers full node) means most "are these two frames identical?" checks resolve in nanoseconds when frames differ at the first subtree

**Important caveat**: Floating-point serialization must be deterministic. If the same conformation produces slightly different float representations (due to integration noise), all frames will appear unique. This experiment measures *bitwise* frame identity, which is most meaningful for:
- Detecting exact checkpoint duplicates
- Verifying trajectory integrity (was this frame corrupted?)
- Finding periodic orbits only when the simulation revisits *exactly* the same coordinates

For approximate structural similarity (RMSD < threshold), hashrope is not the right tool. State this honestly.

---

### E-P2: Periodic Conformation Detection via RepeatNode

**Claim**: CB-C6
**Question**: Can RepeatNode detect and compress periodic trajectories?

**Method**:
1. Construct a *synthetic* periodic trajectory: frame A → B → A → B → ... for 10K cycles
   (Use two real frames from alanine dipeptide as A and B)
2. Build rope naively (20K leaves) vs build rope with RepeatNode detection:
   - Scan for repeat patterns: if frames i..i+k equal frames i+k..i+2k, create RepeatNode
3. Measure: memory (node count × 64 bytes) for naive vs RepeatNode representation
4. Measure: hash computation time for the full trajectory
5. Extend: detect partial periodicity (A→B→A→B→C→A→B→A→B→C→...) — period of 5 frames

**Output**:
- Table: trajectory length, naive nodes, RepeatNode nodes, memory savings
- Timing: hash of full trajectory (RepeatNode vs materialized)

**Expected outcome**:
- Perfect periodicity (A→B × 10K): RepeatNode uses ~5 nodes vs 20K naive nodes. Memory: ~320 bytes vs ~1.3 MB.
- Near-zero hash computation time for the RepeatNode version (O(log q) = O(log 10K) ≈ 14 steps)

**Honest limitation**: Real MD trajectories are *almost never* exactly periodic. This experiment demonstrates the *capability* on synthetic data. The practical value is in detecting exact duplicates (e.g., simulation restart artifacts) and in representing *intentionally constructed* repeat structures (e.g., homopolymer simulations).

---

## 6. Experiments — Phase 3: Drug Discovery

### E-D1: Chemical Compound Exact Lookup

**Claim**: CB-C7
**Question**: How fast is hash-based exact compound lookup vs RDKit canonical SMILES comparison?

**Method**:
1. Download PubChem canonical SMILES for 1M compounds
2. Build hashrope: sort SMILES lexicographically, concatenate as length-prefixed byte strings, build rope
3. Also build: Python dict mapping hash(SMILES) → compound_id using PolynomialHash
4. Query benchmark — for 10,000 random query compounds (5,000 in database, 5,000 not):
   a. Hashrope: compute hash of query SMILES, compare against stored hashes — O(1) per lookup
   b. Python dict with built-in hash: same approach but with Python's hash()
   c. RDKit: canonical SMILES generation + string comparison against sorted list
   d. Naive: linear scan of all 1M SMILES strings
5. Measure: query time (median, p5, p95), false positive rate (hash collision check)

**Output**:
- Table: method, query time, false positive rate
- Collision analysis: among 1M compounds, how many hash collisions with Mersenne-61?

**Expected outcome**:
- Hash lookup: ~10–50 ns per query (hash computation + dict lookup)
- RDKit canonical comparison: ~1–10 µs per query (canonicalization overhead)
- Linear scan: ~100 ms per query (1M string comparisons)
- Collision probability: with Mersenne-61 (2^61 - 1 ≈ 2.3 × 10^18), for 1M compounds, 
  expected collisions ≈ C(10^6, 2) / 2^61 ≈ 2.2 × 10^-7. Effectively zero.

**Honest note**: The Python built-in dict will likely be competitive with hashrope for simple O(1) lookup. Hashrope's advantage emerges when you need *structural queries* — "is this compound a substring of that polymer?", "do compounds 500–600 in the sorted list match this batch?" — where substr_hash at O(log w) provides value a flat dict cannot.

---

### E-D2: Drug Resistance Mutation Panel

**Claim**: CB-C8
**Question**: How fast can hashrope check N known resistance mutation sites on a gene sequence?

**Method**:
1. Download HIV-1 HXB2 reference sequence (9,719 bp) — focus on reverse transcriptase gene (RT, ~1.7 Kbp)
2. Build hashrope from reference RT sequence (single leaf or chunked at 64-byte boundaries)
3. Define resistance panel: top 50 NNRTI + NRTI resistance positions from Stanford HIVDB
   (e.g., K103N, M184V, Y181C — each is a single amino acid = 3 nucleotide positions)
4. For each panel position p_i (i = 1..50):
   a. Time: `substr_hash(ref_rope, p_i * 3, 3)` — hash of the reference codon at position p_i
5. For a "patient sample":
   a. Introduce known resistance mutations at k ∈ {0, 1, 5, 10, 20} positions
   b. Build patient rope
   c. Time: check all 50 panel positions — `substr_hash(patient_rope, p_i*3, 3) == substr_hash(ref_rope, p_i*3, 3)` for each i
   d. Any mismatch = potential resistance mutation
   e. Baseline: extract codons from both sequences, compare strings
6. Scale: extend to full-length HIV genome (9.7 Kbp) with expanded panel (200 positions)
7. Scale: TB resistance panel — rpoB gene (~3.5 Kbp), 50 known rifampicin resistance positions

**Output**:
- Table: panel size, total check time (hashrope vs string extraction + comparison)
- Per-position time breakdown
- Accuracy: false negative rate (should be 0 — exact match is deterministic)
- Scaling plot: panel size vs total time, expect linear in panel size but with very small constant

**Expected outcome**:
- 50 positions × ~100–200 ns per substr_hash = ~5–10 µs total for full panel
- String extraction baseline: 50 positions × ~1–5 µs (seek + copy + compare) = ~50–250 µs
- Speedup: ~10–50× for the panel check itself
- The speedup is modest because the gene is small (1.7 Kbp). The real value is:
  (a) The rope can be pre-built once per reference and reused across thousands of patient samples
  (b) When extended to whole-genome resistance checking (e.g., TB with ~4.4 Mbp genome), the O(log w) advantage compounds

**Honest limitation**: In clinical practice, the bottleneck is *sequencing* (hours), not *analysis* (milliseconds). The hashrope advantage matters for high-throughput surveillance pipelines processing thousands of samples per day, where aggregate analysis time becomes significant.

---

## 7. Experiments — Phase 4: Comparative Genomics

### E-CG1: Cross-Species Divergence Localization

**Claim**: Extends CB-C3 to cross-species comparison
**Question**: How fast can hashrope localize divergent regions between human and chimp genomes?

**Method**:
1. Load GRCh38 chr22 and Pan_tro chr22 (syntenic region)
2. Build ropes for both (4096-byte chunks)
3. Top-level: compare root hashes — expect different (human/chimp ~1.3% divergent)
4. Binary search to enumerate *all* divergent blocks:
   ```
   def find_divergent_regions(ref, sample, start, length, min_block):
       if length <= min_block:
           yield (start, length)
           return
       mid = length // 2
       h_ref_left = substr_hash(ref, start, mid)
       h_sam_left = substr_hash(sample, start, mid)
       if h_ref_left != h_sam_left:
           yield from find_divergent_regions(ref, sample, start, mid, min_block)
       h_ref_right = substr_hash(ref, start+mid, length-mid)
       h_sam_right = substr_hash(sample, start+mid, length-mid)
       if h_ref_right != h_sam_right:
           yield from find_divergent_regions(ref, sample, start+mid, length-mid, min_block)
   ```
5. min_block ∈ {100, 1000, 10000} bp
6. Baseline: nucmer (MUMmer4) — standard whole-genome alignment tool
7. Measure: time to enumerate all divergent regions, number of regions found

**Output**:
- Table: min_block, divergent regions found, total time (hashrope vs nucmer)
- Visualization: chromosome ideogram with divergent regions marked
- Precision check: do hashrope-detected divergent regions overlap with nucmer-detected variants?

**Expected outcome**:
- At ~1.3% divergence, ~98.7% of substr_hash comparisons will match (early exit at O(1) or O(log w))
- The ~1.3% divergent regions will be recursively subdivided to min_block resolution
- Total comparisons: O(D · log(N/min_block)) where D = number of divergent blocks
- For chr22 (~51 Mbp), expect ~66K divergent 100bp blocks × log₂(510K) ≈ 19 levels = ~1.25M hash comparisons × ~300 ns = ~375 ms
- nucmer baseline: ~10–60 seconds for full chr22 alignment
- Hashrope advantage: ~30–160× faster, but *only answers identity/divergence*, not alignment

**Critical distinction for the paper**: nucmer produces *alignments* (insertions, deletions, rearrangements). Hashrope produces *divergence locations* (this block differs, that block matches). These are different outputs. Hashrope is not a replacement for alignment tools — it's a fast pre-filter that says "look here, not there."

---

## 8. Environment

All benchmarks on the same machine (from hashrope FINDINGS.md):
- CPU: [to be filled — run `wmic cpu get name`]
- RAM: 64 GB
- GPU: NVIDIA RTX 4090 (not used — CPU only)
- OS: Windows [version to be filled]
- Python: 3.10+ with hashrope from PyPI
- Rust: rustc 1.94.0 (if using Rust implementation)
- Compiler flags: `--release` (opt-level 3) for Rust
- Key dependencies: samtools, BLAST+, MDAnalysis, RDKit, MUMmer4

**Memory note**: Full GRCh38 rope at 4096-byte chunks:
- ~800K leaves × 64 bytes/node ≈ ~100 MB node metadata + ~3.1 GB leaf data = ~3.2 GB total
- Fits comfortably in 64 GB RAM. Includes room for a second rope (mutant/sample) simultaneously.

**Disk note**: GRCh38 FASTA + index: ~3.5 GB. PDB subset: ~2 GB. PubChem SMILES: ~200 MB. Total: ~6 GB.

---

## 9. Reporting Standards

- All times: median ± MAD from ≥100 iterations (or ≥10 for slow operations)
- Throughput: MB/s or queries/sec as appropriate
- Memory: peak RSS or analytical estimate (node_count × 64 + payload bytes)
- Statistical tests: Wilcoxon signed-rank for paired comparisons (hashrope vs baseline on same queries)
- Effect sizes reported alongside p-values
- All query positions/parameters logged for reproducibility
- Raw timing data archived in `results/` as CSV or JSON
- No cherry-picking: all experiments reported including negative results
- Negative results (e.g., "hashrope slower than baseline for small regions") reported honestly with analysis

---

## 10. Known Limitations (to state in any publication)

1. **Exact identity only**: hashrope detects whether two regions are *bitwise identical*. It cannot detect similarity, homology, or structural equivalence. It complements, not replaces, BLAST/HMMer/RMSD-based tools.

2. **Non-cryptographic hash**: Mersenne-61 polynomial hash is *not* collision-resistant against adversarial manipulation. For data integrity in untrusted environments (e.g., shared genomic databases), a cryptographic hash must be layered on top. For corruption detection and internal consistency checks, the polynomial hash is sufficient.

3. **Sequencing error sensitivity**: A single base error causes a hash mismatch. The binary search will localize the error, not the mutation. This approach works best on:
   - Reference-to-reference comparisons (no sequencing errors)
   - High-accuracy sequencing (PacBio HiFi Q30+, Illumina with error correction)
   - Consensus sequences (assembled, not raw reads)

4. **Floating-point determinism**: MD trajectory frame comparison requires bitwise-identical serialization. Different MD engines, compilers, or even run orders may produce slightly different float representations for the same physical conformation. The experiment must control for this.

5. **One-time construction cost**: Building a genome-scale rope takes ~14 seconds. This is practical for reference genomes (built once, queried many times) but impractical for single-use comparisons where the construction time exceeds the analysis time.

6. **Chunk size sensitivity**: The choice of chunk size affects tree depth, construction time, and query performance. The experiments must characterize this tradeoff rather than reporting only the best chunk size.

7. **No SIMD optimization (yet)**: Current polynomial hash is byte-at-a-time. SIMD unrolling (planned Phase 2 of hashrope) would improve construction throughput from ~215 MB/s to potentially ~1–2 GB/s, substantially reducing the amortization point.

---

## 11. Execution Order

1. **E-G4** (rope construction) — establishes baseline costs, validates setup
2. **E-G1** (region queries) — core claim, highest priority
3. **E-G3** (mutation localization) — builds on E-G1, demonstrates application
4. **E-G2** (tandem repeats) — independent experiment, uses RepeatNode
5. **E-D2** (resistance panel) — small-scale, fast to run, clinically motivating
6. **E-D1** (compound lookup) — independent, requires PubChem download
7. **E-P1** (MD frames) — requires MDAnalysis setup, larger data
8. **E-P2** (periodic detection) — synthetic data, builds on E-P1
9. **E-CG1** (cross-species) — requires chimp genome download, largest experiment

**Rationale**: Genomics experiments first (smallest setup cost, strongest claims). Drug discovery next (fast experiments, good narrative). Protein structure last (largest data, most setup). Cross-species comparison as the capstone experiment.

---

## 12. Paper Structure (Preliminary)

If targeting a single comprehensive paper:

1. **Introduction**: The O(log w) substr_hash gap — what no existing bioinformatics tool provides
2. **Background**: Polynomial hash algebra, rope data structure, RepeatNode
3. **Methods**: Hashrope construction from biological sequences
4. **Experiments**:
   - 4.1 Genome-scale region queries (E-G1, E-G4)
   - 4.2 Mutation localization (E-G3)
   - 4.3 Tandem repeat compression (E-G2)
   - 4.4 Drug resistance panels (E-D2)
   - 4.5 Cross-species divergence (E-CG1)
5. **Discussion**: When to use hashrope vs BLAST/samtools/alignment (the complementarity argument)
6. **Limitations**: Section 10 above
7. **Future work**: MD trajectory indexing (E-P1/E-P2), compound databases (E-D1), SIMD optimization

If results warrant two papers:
- Paper A (Bioinformatics): Sections 4.1–4.3, 4.5 — pure genomics
- Paper B (JCIM or Briefings): Sections 4.4, E-D1, E-P1 — applications in drug discovery and structural biology
