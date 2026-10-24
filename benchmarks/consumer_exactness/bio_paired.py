"""Nine alternating paired runs on stored HIV RT data, without registry installs.

Old consumer source comes from a pinned bio commit; both arms use current
Hashrope, the same prebuilt ropes, warm caches and data. This isolates the cost
of verifying consumer decisions, not unrelated core version changes. Run from
any directory. Results go to a new explicitly named JSON; historical data stays
untouched. No clinical treatment conclusions are drawn from this timing test.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import statistics
import subprocess
import sys
import time
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CORE = Path(os.environ.get("HASHROPE_BENCH_CORE", ROOT))
BIO = ROOT / "usecases/computational-biology/code"
BASELINE_REF = "206d1af783ca4c34626cfe1f9ce9a7bd7d4275cb"
sys.path[:0] = [str(CORE / "packages/python"), str(BIO / "src")]
from hashrope import PolynomialHash, Leaf, rope_concat
import hashrope
assert Path(hashrope.__file__).resolve() == (CORE / "packages/python/hashrope/__init__.py").resolve()
assert Path(sys.modules["hashrope.rope"].__file__).resolve() == (CORE / "packages/python/hashrope/rope.py").resolve()
from hashrope_bio.genomics import gene_diff as current_genes
from hashrope_bio.genomics import region_query as current_regions
from hashrope_bio.cheminformatics import resistance as current_panels


def old_module(relative, name, baseline_ref):
    source = subprocess.check_output(["git", "show", f"{baseline_ref}:{relative}"], cwd=BIO)
    module = types.ModuleType(name)
    sys.modules[name] = module
    exec(compile(source, f"{baseline_ref}:{relative}", "exec"), module.__dict__)
    return module, hashlib.sha256(source).hexdigest()


def fasta(path):
    rows, header, parts = [], None, []
    for line in path.read_bytes().splitlines():
        if line.startswith(b">"):
            if header is not None:
                rows.append((header, b"".join(parts).upper()))
            header, parts = line[1:].decode("ascii"), []
        elif line:
            parts.append(line.strip())
    if header is not None:
        rows.append((header, b"".join(parts).upper()))
    return rows


def build(data, h):
    root = None
    for offset in range(0, len(data), 64):
        root = rope_concat(root, Leaf(data[offset:offset + 64], h), h)
    return root


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=128)
    parser.add_argument("--runs", type=int, default=9)
    parser.add_argument("--baseline-ref", default=BASELINE_REF)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    baseline_ref = subprocess.check_output(["git", "rev-parse", "--verify", f"{args.baseline_ref}^{{commit}}"], cwd=BIO, text=True).strip()
    old_genes, gh = old_module("src/hashrope_bio/genomics/gene_diff.py", "bio_old_genes", baseline_ref)
    old_regions, rh = old_module("src/hashrope_bio/genomics/region_query.py", "bio_old_regions", baseline_ref)
    old_panels, ph = old_module("src/hashrope_bio/cheminformatics/resistance.py", "bio_old_panels", baseline_ref)
    files = [BIO / "data/hiv1_hxb2.fa", BIO / "data/hivdb_rt_expanded_aligned.fa"]
    reference = fasta(files[0])[0][1][2549:4229]
    all_patients = fasta(files[1])
    patients = all_patients[:args.samples]
    h = PolynomialHash()
    panel = current_panels.HIV_RT_NNRTI_PANEL + current_panels.HIV_RT_NRTI_PANEL
    cases = []
    build_start = time.perf_counter_ns()
    for header, sample in patients:
        match = re.search(r"RT:(\d+)-(\d+)", header)
        if match is None:
            raise ValueError("Missing RT coverage in stored FASTA header")
        first, last = map(int, match.groups())
        ref = reference[(first - 1) * 3:min(last * 3, len(reference))]
        n = min(len(ref), len(sample))
        ref, sample = ref[:n], sample[:n]
        if not n:
            continue
        relative_panel = [current_panels.ResistanceSite(
            site.gene, site.position - first + 1, site.wildtype_aa,
            site.mutant_aa, site.drug_class, site.annotation)
            for site in panel if first <= site.position <= last and (site.position - first + 1) * 3 <= n]
        gene = current_genes.GeneRegion("RT covered region", "HXB2 RT", "HIV-1", "+", 0, n)
        cases.append((build(ref, h), build(sample, h), ref, sample, relative_panel, gene))
    # Independent equal controls exercise verification beyond early mismatches.
    for _, _, ref, _, sites, gene in list(cases[:16]):
        cases.append((build(ref, h), build(ref, h), ref, ref, sites, gene))
    construction_ms = (time.perf_counter_ns() - build_start) / 1e6

    def run(kind, arm):
        output = []
        for a, b, raw_a, raw_b, sites, gene in cases:
            if kind == "regions":
                module = old_regions if arm == "old_hash_only" else current_regions
                output.append(module.regions_identical(a, b, 0, len(raw_a), h))
            elif kind == "gene_reports":
                module = old_genes if arm == "old_hash_only" else current_genes
                output.append(module.diff_genes(a, b, [gene], h, include_exons=False).gene_diffs[0].body_changed)
            else:
                module = old_panels if arm == "old_hash_only" else current_panels
                output.extend(result.is_mutant for result in module.check_resistance_panel(a, b, sites, h))
        return output

    expected = {
        "regions": [a == b for _, _, a, b, _, _ in cases],
        "gene_reports": [a != b for _, _, a, b, _, _ in cases],
        "panels": [a[(site.position - 1) * 3:site.position * 3] != b[(site.position - 1) * 3:site.position * 3]
                   for _, _, a, b, sites, _ in cases for site in sites],
    }
    results = {}
    for kind in expected:
        timings = {"old_hash_only": [], "exact": []}
        warm_times = {}
        for arm in timings:
            start = time.perf_counter_ns()
            assert run(kind, arm) == expected[kind]
            warm_times[arm] = time.perf_counter_ns() - start
        repetitions = max(1, min(1000, 5_000_000 // max(1, warm_times["old_hash_only"])))
        for iteration in range(args.runs):
            order = ("old_hash_only", "exact") if iteration % 2 == 0 else ("exact", "old_hash_only")
            for arm in order:
                start = time.perf_counter_ns()
                for _ in range(repetitions):
                    result = run(kind, arm)
                timings[arm].append((time.perf_counter_ns() - start) / 1e6 / repetitions)
                assert result == expected[kind], (kind, arm, iteration)
        medians = {arm: statistics.median(values) for arm, values in timings.items()}
        results[kind] = {"operations_per_run": len(expected[kind]), "all_byte_oracles_match": True,
                         "workload_repetitions_per_timed_sample": repetitions,
                         "milliseconds": timings, "median_ms": medians,
                         "exact_over_old_ratio": medians["exact"] / medians["old_hash_only"]}
        print(kind, results[kind]["median_ms"], flush=True)

    output = {
        "purpose": "Exact consumer migration on stored real HIV data; no clinical efficacy claim",
        "comparison_modes": ["old_hash_only", "exact"], "paired_runs": args.runs,
        "sampling": "nine within-process paired runs, not nine independent workers",
        "order": "alternating; one warmup each; construction excluded; same current core in both arms",
        "python": sys.version, "platform": platform.platform(), "hashrope_file": hashrope.__file__,
        "bio_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=BIO, text=True).strip(),
        "baseline_ref": baseline_ref,
        "old_source_sha256": {"gene_reports": gh, "regions": rh, "panels": ph},
        "exact_source_sha256": {kind: hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()
                                for kind, module in [("gene_reports", current_genes),
                                                     ("regions", current_regions), ("panels", current_panels)]},
        "core_source_sha256": {name: hashlib.sha256((CORE / "packages/python/hashrope" / name).read_bytes()).hexdigest()
                               for name in ["rope.py", "polynomial_hash.py"]},
        "corpora": [{"file": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
                     "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in files],
        "corpus_patient_records": len(all_patients), "patient_records_used": len(patients),
        "selection": "first FASTA records in stored order plus 16 independently built equal controls",
        "cases": len(cases), "construction_ms": construction_ms, "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
