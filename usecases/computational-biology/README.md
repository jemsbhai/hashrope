# Use Case 2: Computational Biology & Drug Discovery

**Working title**: Compositional Polynomial Hashing for Sequence Identity, Structural Comparison, and Simulation Integrity in Computational Biology

---

## Scope

This use case covers hashrope applications across multiple sub-domains of computational biology where the core primitives — O(log w) `substr_hash`, O(1) `hash_concat`, O(log q) `RepeatNode`, and bounded-memory sliding window — provide advantages over existing tools.

The common thread: biological data is overwhelmingly **sequential and string-like** (DNA, RNA, amino acid sequences, simulation trajectories, serialized atomic coordinates). Hashrope is a string data structure. The fit is natural.

**Important constraint**: Hashrope provides **exact identity testing**, not similarity/approximate matching. Every sub-domain analysis below must be evaluated against this constraint. Where the domain primarily needs approximate matching (e.g., BLAST-style homology search), hashrope is complementary at best, not a replacement.

---

## Sub-domain 1: Genomics & Sequence Analysis

### 1a. Reference Genome Region Comparison

**Fit: STRONG**

A reference genome (GRCh38) is ~3.2 billion base pairs — a byte string. Build a hash rope with 4KB chunks (~800K leaves). Any arbitrary region query — "is the patient's chr7:1,200,000–1,205,000 identical to the reference?" — is a single `substr_hash` call at O(log w), roughly 200–300 ns extrapolated from FINDINGS.md benchmarks.

**What hashrope adds over existing tools:**
- samtools faidx: O(n) — must read and hash the region bytes. hashrope: O(log w) — never touches the bytes.
- Existing tools answer "give me the bytes at this region." Hashrope answers "are these two regions identical?" without reading either.

**Where it fits in the pipeline:**
- After alignment (minimap2/BWA-MEM finds candidate positions), hashrope provides fast exact-match verification.
- Variant calling: binary search via substr_hash to locate where a sample diverges from the reference. At O(log w) per comparison step, localizing a mutation in a 100Kbp region takes ~17 steps = ~5 µs total.

### 1b. Tandem Repeat Analysis

**Fit: STRONG (RepeatNode is purpose-built for this)**

Tandem repeats (STRs, VNTRs) are regions where a short motif repeats many times. In the human genome, ~3% is tandem repeats. The RepeatNode represents "CAGCAGCAG...CAG" (q copies of "CAG") in exactly 228 bytes regardless of q. This is not just memory savings — the hash is computed in O(log q), so comparing two individuals' repeat counts at a locus is:
1. Build RepeatNode(CAG, q1) and RepeatNode(CAG, q2)
2. Compare hashes in O(1)
3. If different, q1 ≠ q2 (no false negatives for exact match)

**Clinical relevance:** Huntington's disease, Fragile X syndrome, and ~50 other disorders are caused by tandem repeat expansions. Fast repeat-length comparison is clinically meaningful.

### 1c. Metagenomics — Taxonomic Classification

**Fit: MODERATE**

Environmental DNA samples contain fragments from thousands of species. Current approach: k-mer matching against reference databases (Kraken2, Centrifuge). Hashrope alternative: build ropes for each reference genome, use `substr_hash` to test if an environmental read exactly matches any reference region. This is faster per-query than full alignment, but:

- Limitation: environmental reads have sequencing errors (~0.1–1% error rate). Exact matching fails on reads with even one error. This limits the approach to high-accuracy sequencing (PacBio HiFi, ~Q30+).
- Advantage: for *confirmed exact matches*, hashrope provides the match location via binary search in O(log w) — faster than re-running alignment.

### 1d. Epigenomics — Methylation Pattern Comparison

**Fit: MODERATE**

CpG methylation data can be represented as a binary overlay on the genome (methylated/unmethylated at each CpG site). Building a rope of methylation states and using substr_hash to compare methylation patterns between samples at specific genomic regions is technically sound. RepeatNode could represent uniformly methylated/unmethylated regions (CpG islands, CpG deserts) efficiently.

- Limitation: methylation data is typically fractional (0–100% methylation per CpG from bisulfite sequencing), not binary. Discretization to binary loses information.
- Advantage: for comparing fully-methylated vs fully-unmethylated calls at single-CpG resolution, exact matching is appropriate.

---

## Sub-domain 2: Proteomics & Protein Structure

### 2a. Protein Sequence Motif Search

**Fit: STRONG**

Protein sequences are strings over a 20-letter amino acid alphabet. The polynomial hash works on any byte sequence. Build a rope from a protein database (UniProt: ~250M sequences, ~80GB total sequence data). Use substr_hash to answer: "does this motif appear at this position in this protein?" in O(log w).

**Where it fits:**
- PROSITE pattern matching: after a regex/HMM identifies candidate motifs, hashrope verifies exact identity of the matched region against a reference motif in O(log w).
- Domain boundary verification: given a predicted domain at residues 50–200, compare against the canonical domain sequence via substr_hash.

**Key limitation:** Most protein searches require *similarity* (BLAST, HMMer), not exact identity. Proteins diverge evolutionarily but retain function. Hashrope is useful for the subset of tasks requiring exact match (database integrity, clone detection, contamination checking).

### 2b. Protein-Protein Interaction Interface Comparison

**Fit: MODERATE-TO-STRONG (depends on representation)**

Protein-protein interactions occur at binding interfaces. These interfaces can be characterized by:
- **Contact residue sequences**: the ordered list of residues at the interface → a string. Two interfaces with identical contact residue sequences share the same structural hash.
- **Structural fingerprints**: serialize the 3D coordinates of interface atoms into a byte array. Exact hash match = identical geometry.

**Where hashrope adds value:**
- Scanning a database of known interfaces (PDB: ~200K structures): build ropes for each interface fingerprint, compare a query interface via root hash in O(1).
- Identifying conserved interaction motifs: substr_hash on interface sequences to find shared sub-motifs in O(log w).

**Key limitation:** Biologically equivalent interfaces are rarely *exactly* identical — they're similar within RMSD thresholds. Hashrope catches exact matches only. However, discretized/binned representations (e.g., contact maps at 1-Å bins) can turn approximate similarity into exact identity for practical purposes.

### 2c. Molecular Dynamics Trajectory Management

**Fit: STRONG (this is the RepeatNode + substr_hash sweet spot)**

Protein folding / MD simulations generate trajectories: ordered sequences of frames, each frame being a serialized array of atomic coordinates (~10K–100K atoms × 3 floats × 4 bytes = 120KB–1.2MB per frame). Simulations run for millions of frames.

**What hashrope enables:**
1. **Checkpoint indexing**: Build a rope of serialized frames. `substr_hash` at frame boundaries allows "is frame 500,000 identical to frame 300,000?" in O(log w) without reading either frame. This detects periodic orbits (protein returning to a previous conformation) in sub-microsecond time.
2. **RepeatNode for periodic trajectories**: If a protein oscillates between two conformations (as many do), the trajectory has repeat structure. RepeatNode represents "conformation A → B → A → B → ... × 10,000 cycles" in O(1) space instead of materializing 20,000 frames.
3. **Trajectory integrity**: Root hash verifies that a multi-terabyte trajectory file hasn't been corrupted or truncated. substr_hash localizes corruption to a specific frame range.
4. **Simulation comparison**: Two simulations of the same protein with different force fields — binary search via substr_hash to find where the trajectories diverge.

**Data scale**: A 1 µs all-atom MD simulation of a 50K-atom protein at 1 ps resolution = 1M frames × ~600KB = ~600GB. This is well within the range where O(log w) vs O(n) matters enormously.

### 2d. AlphaFold Structure Database Integrity

**Fit: STRONG**

The AlphaFold Protein Structure Database contains predicted structures for ~200M proteins. Each structure is a set of 3D coordinates + confidence scores. Building a hash rope over the database enables:
- **Version control**: when AlphaFold releases updated predictions, binary search via substr_hash identifies which structures changed.
- **Deduplication**: identical predicted structures (e.g., same protein from different species with identical sequences) detected via O(1) hash comparison.
- **Integrity**: researchers downloading subsets can verify against the root hash that their copy is unaltered.

---

## Sub-domain 3: Drug Discovery

### 3a. Virtual Screening — Molecular Fingerprint Matching

**Fit: WEAK for primary screening, MODERATE for verification**

Virtual screening searches for drug candidates by comparing molecular fingerprints (e.g., ECFP4, MACCS keys) against a target. This is fundamentally a *similarity* search (Tanimoto coefficient), which hashrope cannot accelerate.

**Where hashrope fits:**
- **Hit verification**: after a similarity search returns candidate molecules, verify exact fingerprint identity against a reference compound library. "Is this hit the same molecule we already have in our compound library?" → O(1) hash comparison.
- **Compound library deduplication**: pharmaceutical compound libraries (millions of molecules) often contain duplicates from different vendors. Build a rope of sorted fingerprints, detect duplicates via hash comparison.
- **NOT** the primary screening step. Be honest about this.

### 3b. Molecular Docking — Binding Pose Comparison

**Fit: MODERATE**

Molecular docking generates many candidate binding poses (orientations of a drug molecule in a protein binding site). Each pose is a set of atomic coordinates.

**Where hashrope fits:**
- **Pose deduplication**: docking software often generates near-identical poses. Discretize coordinates (e.g., round to 0.1 Å), serialize, hash → identical hashes = duplicate poses. Eliminates redundant downstream scoring.
- **Cross-docking comparison**: "does molecule A adopt the same pose in protein X as molecule B?" → hash comparison of discretized pose coordinates.

**Key limitation:** Pose comparison typically uses RMSD thresholds, not exact identity. The discretization trick helps but introduces binning artifacts.

### 3c. Drug Resistance Mutation Tracking

**Fit: STRONG**

This connects directly to Sub-domain 1a. Drug resistance (HIV, cancer, bacteria) arises from specific mutations in target genes. Monitoring resistance:
1. Build a rope of the wild-type gene sequence.
2. Patient sample → substr_hash at known resistance positions → O(log w) per position.
3. Hash mismatch at position X = potential resistance mutation.
4. Faster than re-sequencing or alignment for routine monitoring at *known* positions.

**Clinical workflow**: For a panel of 50 known resistance sites in HIV reverse transcriptase (~1.7Kbp gene), checking all 50 sites: 50 × ~100 ns = ~5 µs total. Current PCR-based assays take hours.

Note: the hashrope check operates on already-sequenced data. The bottleneck in clinical workflows is the sequencing itself, not the analysis. But for high-throughput surveillance (thousands of samples), the analysis time adds up.

### 3d. SMILES/InChI String Comparison for Chemical Databases

**Fit: STRONG**

Chemical compounds are represented as SMILES or InChI strings. PubChem has ~110M compounds, each with a canonical SMILES string. Build a rope over the database:
- **Exact compound lookup**: hash a query SMILES, compare against database hashes in O(1).
- **Substructure queries** (limited): if the substructure query maps to a substring of the SMILES, substr_hash can test identity. But SMILES substrings don't always correspond to chemical substructures (branching notation breaks this). InChI layers are more amenable — the connectivity layer can be compared independently via substr_hash.
- **Database integrity**: root hash over the entire PubChem database verifies no records were altered.

---

## Sub-domain 4: Phylogenetics & Evolutionary Biology

### 4a. Whole-Genome Alignment Verification

**Fit: STRONG**

Comparative genomics aligns whole genomes of different species. After alignment (via progressiveMauve, mugsy, etc.), verification: "are these two aligned regions truly identical?" → substr_hash in O(log w). For closely related species (human vs chimp, ~98.7% identity), the binary search localizes the ~1.3% divergent regions efficiently.

### 4b. Horizontal Gene Transfer Detection

**Fit: MODERATE**

Genes acquired from other species via horizontal transfer should have sequence identity with the donor. Build ropes of candidate genomes, use substr_hash to test whether a region of genome A is identical to any region of genome B. This is an O(log w) exact-match test, complementing the approximate phylogenetic methods (GC content anomaly, codon usage bias) typically used.

---

## Unified Experiment Plan

### Phase 1: Genomic Foundation (Proof-of-concept)
1. Download GRCh38 reference genome (FASTA, ~3.1 GB)
2. Build hashrope with 4KB chunks (~800K leaves)
3. Benchmark `substr_hash` queries at various region sizes (100bp, 1Kbp, 10Kbp, 100Kbp)
4. Compare against: samtools faidx + full hash, and direct memcmp
5. Tandem repeat compression: measure RepeatNode savings on known STR loci
6. Mutation localization: introduce synthetic SNPs, time the binary search

### Phase 2: Protein Structure
7. Load PDB structures, serialize atomic coordinates
8. Build rope per-protein, benchmark frame-to-frame comparison via substr_hash
9. MD trajectory indexing: load a real trajectory (GROMACS .xtc), build rope of frames
10. Measure: time to detect periodic conformations vs naive frame-by-frame comparison

### Phase 3: Drug Discovery Integration
11. Build rope over PubChem SMILES subset (~1M compounds)
12. Benchmark exact compound lookup via hash comparison
13. Compound library deduplication: measure dedup rate and time
14. Drug resistance panel: build rope of HIV RT gene, benchmark 50-site resistance check

### Phase 4: Comparative Genomics
15. Build ropes for human + chimp chromosomes
16. Binary search for divergent regions via substr_hash
17. Compare localization time against nucmer/MUMmer

---

## Key benchmark numbers (from hashrope FINDINGS.md)

| Operation | Performance | Scaling |
|-----------|------------|---------|
| substr_hash (warm, 10K leaves) | 124 ns | O(log w) |
| substr_hash (warm, 1M leaves) | 284 ns | O(log w) |
| Full-node early exit | 2.1 ns constant | O(1) |
| concat (rejoin) | 55–80 ns | O(log w) |
| RepeatNode | 228 bytes, 262–363 ns | O(1) space, O(log q) time |
| Memory overhead (4KB chunks) | 1.03× | Near-zero |

---

## Target venues

**Primary:**
- Bioinformatics (Oxford) — methodological paper with genomics + proteomics experiments
- PLoS Computational Biology — broader scope, drug discovery integration
- RECOMB — conference, algorithmic focus

**Secondary (sub-domain specific):**
- BMC Bioinformatics — genomics focus
- WABI (Workshop on Algorithms in Bioinformatics) — algorithmic focus
- Journal of Chemical Information and Modeling (JCIM) — drug discovery / cheminformatics angle
- Journal of Molecular Biology — structural comparison angle
- Briefings in Bioinformatics — review/methods paper covering all sub-domains

**Strategy:** One comprehensive paper covering Sub-domains 1–3 with Sub-domain 4 as future work, targeting Bioinformatics or PLoS Comp Bio. Alternatively, two papers: (1) genomics + sequence analysis for Bioinformatics, (2) structural biology + drug discovery for JCIM. The choice depends on how strong the experimental results are in each sub-domain.

---

## Honest limitations to acknowledge in any paper

1. **Not a similarity tool**: hashrope tests exact identity, not homology. It complements BLAST/HMMer, does not replace them.
2. **Non-cryptographic hash**: polynomial hash over Mersenne-61 is not collision-resistant against adversarial manipulation. For data integrity in untrusted environments, a cryptographic hash layer is needed on top.
3. **Sequencing error sensitivity**: any base-level error in sequencing data causes a hash mismatch. The approach works best on high-accuracy data (PacBio HiFi, Illumina at Q30+) or on reference-to-reference comparisons.
4. **Discretization artifacts**: when applying to continuous data (3D coordinates, methylation fractions), discretization is required, introducing binning artifacts and potential false negatives near bin boundaries.
5. **Sequential throughput**: at ~215 MB/s, initial rope construction from a 3.1 GB genome takes ~14 seconds. This is a one-time cost amortized over many queries, but should be reported honestly.

---

## Status

Planning phase. Folder renamed from `genomic-region-comparison` to `computational-biology` to reflect expanded scope.
