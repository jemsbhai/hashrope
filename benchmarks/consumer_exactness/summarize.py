"""Aggregate fresh-process medians; preserve raw data and historical outputs."""
import json
import math
from pathlib import Path
import statistics

HERE = Path(__file__).resolve().parent
RUNS = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((HERE / "paper_results").glob("run_s*_i*.json"))]
assert len(RUNS) == 9, f"Expected nine complete independent processes, got {len(RUNS)}"


def stats(values):
    return {"mean": statistics.mean(values), "std": statistics.stdev(values) if len(values)>1 else 0,
            "n": len(values), "values": values}


def slope(xs, ys):
    xx, yy = list(map(math.log, xs)), list(map(math.log, ys))
    xm, ym = statistics.mean(xx), statistics.mean(yy)
    return sum((x-xm)*(y-ym) for x,y in zip(xx,yy)) / sum((x-xm)**2 for x in xx)


def aggregate(section, key, arms):
    rows = {}
    for run in RUNS:
        for row in run[section]:
            k = key(row)
            rows.setdefault(k, []).append(row)
    return {k: {arm: stats([r[arm]["median_ms"] for r in group]) for arm in arms} for k,group in rows.items()}


summary = {"comparison_modes": ["old_hash_only", "exact_structural"], "fresh_processes": len(RUNS),
           "exp005": aggregate("exp005", lambda r:f"{r['bytes']}/{r['fraction']}", ["old_hash_only","exact","brute"]),
           "exp017": aggregate("exp017_controlled", lambda r:str(r["expected_tokens"]), ["old_hash_only","exact","radix","numpy"])}
# R2 has three queries per process. First take their mean within each process;
# uncertainty remains across nine independent processes, not 27 queries.
summary["exp015_r2_cpu"] = {}
for length in [1000,16000,64000,256000,1000000]:
    row = {}
    for arm in ["old_hash_only","exact","radix","numpy"]:
        row[arm] = stats([statistics.mean(c[arm]["median_ms"] for c in r["exp015_r2_cpu"]
                         if c["expected_tokens"] == length) for r in RUNS])
    row["query_rope_build"] = stats([statistics.mean(c["query_rope_build_ms"] for c in r["exp015_r2_cpu"]
                                    if c["expected_tokens"] == length) for r in RUNS])
    summary["exp015_r2_cpu"][str(length)] = row
sizes = [64000,256000,1000000,2000000,4000000,8000000]
summary["exp005_slope_f05"] = {arm: slope(sizes,[summary["exp005"][f"{n}/0.5"][arm]["mean"] for n in sizes])
                              for arm in ["old_hash_only","exact","brute"]}
summary["exp005_slope_f05_per_run"] = {arm:stats([slope(sizes,[next(c[arm]["median_ms"] for c in r["exp005"]
                                  if c["bytes"]==n and c["fraction"]==.5) for n in sizes]) for r in RUNS])
                                  for arm in ["old_hash_only","exact","brute"]}
grid = sorted(map(int,summary["exp017"]))
stable = [n for n in grid if all(summary["exp017"][str(m)]["exact"]["mean"] + summary["exp017"][str(m)]["exact"]["std"]
                < summary["exp017"][str(m)]["radix"]["mean"] - summary["exp017"][str(m)]["radix"]["std"] for m in grid if m>=n)]
summary["exp017_stable_separation_from_tokens"] = min(stable) if stable else None
cross = summary["exp017_stable_separation_from_tokens"]
summary["exp017_original_crossover_gate_pass"] = bool(cross) and all(
    summary["exp017"][str(n)]["radix"]["mean"] <= summary["exp017"][str(n)]["exact"]["mean"] * 1.001
    for n in grid if n < cross)
summary["exp017_paired_wins_2M"] = sum(next(c for c in r["exp017_controlled"] if c["expected_tokens"]==2000000)["exact"]["median_ms"]
    < next(c for c in r["exp017_controlled"] if c["expected_tokens"]==2000000)["radix"]["median_ms"] for r in RUNS)
sources = set.intersection(*(set(r.get("source_sha256",{})) for r in RUNS))
summary["timed_worker_sources_consistent"] = {p: len({r["source_sha256"][p] for r in RUNS}) == 1 for p in sorted(sources)}
summary["correctness"] = {"all_timing_results_asserted_against_oracle": True,
    "exp005_cells": sum(len(r["exp005"]) for r in RUNS),
    "exp017_cells": sum(len(r["exp017_controlled"]) for r in RUNS),
    "r2_queries": sum(len(r["exp015_r2_cpu"]) for r in RUNS)}

EXTENDED = [json.loads(p.read_text(encoding="utf-8"))
            for p in sorted((HERE / "competitive_extended_results").glob("run_s*_i*.json"))]
if EXTENDED:
    assert len(EXTENDED) == 9, "Wait until all extended competitive workers finish"
    summary["exp017_real_pairs"] = {}
    for dataset in ["sharegpt_sample.jsonl", "lmsys_sample.jsonl"]:
        summary["exp017_real_pairs"][dataset] = {
            arm: stats([r["real_pairs"][dataset][f"{arm}_ms_stats"]["mean"] for r in EXTENDED])
            for arm in ["hashrope", "radix", "flat"]}
        summary["exp017_real_pairs"][dataset]["mean_lcp_tokens"] = statistics.mean(
            r["real_pairs"][dataset]["lcp_stats"]["mean"] for r in EXTENDED)
        summary["exp017_real_pairs"][dataset]["exact_faster_pairs"] = sum(
            p["hashrope_median_ms"] < p["radix_median_ms"]
            for r in EXTENDED for p in r["real_pairs"][dataset]["pair_results"])
    summary["exp017_k_sweep"] = {
        key: {arm: stats([r["k_sweep"][key][f"{arm}_median_ms"] for r in EXTENDED])
              for arm in ["hashrope", "radix", "flat"]}
        for key in EXTENDED[0]["k_sweep"]}
    summary["extended_correctness"] = {
        "real_pairs": sum(r["real_pairs"][d]["n_pairs"] for r in EXTENDED for d in r["real_pairs"]),
        "all_real_pairs_correct": all(r["real_pairs"][d]["all_correct"] for r in EXTENDED for d in r["real_pairs"]),
        "all_k_cells_correct": all(c["correctness_ok"] for r in EXTENDED for c in r["k_sweep"].values()),
        "source_hashes_consistent": all(r["source_sha256"] == EXTENDED[0]["source_sha256"] for r in EXTENDED)}


def cell(stats):
    return f"{stats['mean']:.4f} ± {stats['std']:.4f}"


lines = ["# Exact consumer reruns: 2026-10-23", "",
         "This correction preserves fingerprints and evaluates exact consumer decisions. All historical paper result files remain untouched.", "",
         "## Protocol and scope", "",
         "Paper CPU measurements use three corpus seeds × three fresh interpreter processes, five warmed query repetitions per cell, "
         "and mean ± sample standard deviation across the nine process medians. Old and exact arms run on the same independently built ropes; "
         "the old arm is loaded from paper commit `c9d8902e0f443f1a5ba12da78f087d5dbd2cad68`. "
         "Both use current Hashrope, isolating consumer semantics. No shared-root shortcut is supplied. "
         "The original 4096-byte leaf builder, controlled divergent-tail construction, vendored radix implementation, "
         "and NumPy int64 comparator are retained. Old/exact order reverses on alternating invocations. "
         "This differs from the earlier 512-byte-leaf microbenchmarks; leaf size materially changes costs.", "",
         "All corpora were re-tokenized offline using the cached GPT-2 tokenizer. Their stored canonical token prefixes match exactly. "
         "Full corpus hashes, source hashes, runtime versions, per-query samples and oracle expectations are in `paper_results/`. "
         "Each observed result was asserted against its oracle before timing data were accepted. "
         "The pre-tokenization manifest predates a final core range-wrapper refactor; each timed worker independently records its actual sources. "
         f"Timed worker source consistency: **{all(summary['timed_worker_sources_consistent'].values())}**. "
         "The workstation is Windows 11, Intel Core i9-14900HX, Python 3.12.2; no CPU affinity/frequency locking was used. "
         "Team CPU-heavy tests were held during timed workers after the startup barrier; a sub-second preliminary biology probe may overlap the earliest work. "
         "These measurements are local CPU results, not a reproduction of the historical GPU server.", "",
         "Construction is excluded from query latency, as in the original methodology. "
         "EXP015 query-rope construction is reported separately. Its three queries per prefix are averaged within each process before computing uncertainty. "
         + ("The full original ShareGPT/LMSYS real-pair and bounded K-candidate sweeps were rerun separately in nine fresh processes. "
          "They retain the original fixed timed arm order: radix, exact, NumPy."
          if EXTENDED else "The ShareGPT/LMSYS real-pair and K-candidate sweeps have not yet completed."), "",
         "## EXP005: exact LCP correctness survives; logarithmic exact-query claim does not", "",
         "Times below are milliseconds, at a mismatch halfway through independent corpus-derived buffers.", "",
         "| Bytes | Old hash-only | Exact | Python byte oracle |", "|---:|---:|---:|---:|"]
for n in sizes:
    r=summary["exp005"][f"{n}/0.5"]
    lines.append(f"| {n:,} | {cell(r['old_hash_only'])} | {cell(r['exact'])} | {cell(r['brute'])} |")
sl=summary["exp005_slope_f05"]
lines += ["", f"Log-log slopes over 64K–8M bytes: old hash-only **{sl['old_hash_only']:.3f}**, exact **{sl['exact']:.3f}**, "
          f"byte oracle **{sl['brute']:.3f}**. The historical ≤0.3 exact-scaling gate "
          f"{'fails' if sl['exact']>.3 else 'passes this finite grid only'}; there is no general logarithmic exactness guarantee. "
          "Zero hash calls now measures fingerprint independence, not work. Fractions 0 and 0.99 are also present in the raw and summary data.", "",
          "## EXP017: controlled competitive win remains, with different crossover", "",
          "| Shared tokens | Old hash-only ms | Exact ms | Vendored radix ms | NumPy ms |", "|---:|---:|---:|---:|---:|"]
for n in grid:
    r=summary["exp017"][str(n)]
    lines.append(f"| {n:,} | {cell(r['old_hash_only'])} | {cell(r['exact'])} | {cell(r['radix'])} | {cell(r['numpy'])} |")
r=summary["exp017"]["2000000"]
lines += ["", f"At 2M tokens exact beats radix in **{summary['exp017_paired_wins_2M']}/9** paired runs. "
          f"The ratio of mean radix time to mean exact time is **{r['radix']['mean']/r['exact']['mean']:.2f}×**. "
          f"Mean ± 1 SD bands remain separated from **{summary['exp017_stable_separation_from_tokens']:,} tokens** upward on this grid. "
          "This is a measured CPU implementation result, not evidence for sublinear exact comparison. NumPy remains the flat-array reference. "
          f"The original strict crossover gate **{'passes' if summary['exp017_original_crossover_gate_pass'] else 'fails'}**: "
          "it also demands that radix's mean win at every smaller point; at 4K the means favor exact but the bands overlap. "
          "Thus this is a qualified long-context result, not a renewed full historical B1 promotion.", ""]
if EXTENDED:
    lines += ["### Real conversational pairs and K candidates", "",
              "The original 200 pairs per dataset per seed were rerun with exact, radix, and NumPy helpers. "
              "The table averages the per-pair median across each process's 200 pairs, then reports mean ± sample SD across nine fresh processes. "
              "Cached/query ropes are independently constructed. The same pairs recur across invocations of a seed; the 3,600 comparisons are not 3,600 independent samples.", "",
              "| Dataset | Mean shared tokens | Exact ms/pair | Radix ms/pair | NumPy ms/pair |",
              "|---|---:|---:|---:|---:|"]
    for dataset, row in summary["exp017_real_pairs"].items():
        lines.append(f"| {dataset} | {row['mean_lcp_tokens']:.1f} | {cell(row['hashrope'])} | {cell(row['radix'])} | {cell(row['flat'])} |")
    lines += ["", "For K candidates, radix performs one tree lookup; Hashrope and NumPy perform K pairwise comparisons. "
              "The original workload shares token values but builds separate Hashrope roots; it does not add a Hashrope candidate index or share prefix nodes. "
              "Times below cover the whole lookup across K candidates, with construction excluded.", "",
              "| Shared tokens / K | Exact ms | Radix ms | NumPy ms |", "|---|---:|---:|---:|"]
    for key, row in summary["exp017_k_sweep"].items():
        length, k = key.split("/")
        lines.append(f"| {int(length):,} / {k} | {cell(row['hashrope'])} | {cell(row['radix'])} | {cell(row['flat'])} |")
    lines += ["", "All real-pair and K-cell oracle gates passed. These cells test the corrected implementation directly; "
              "the historical hash-only penalties should not be carried forward as exact performance estimates.", ""]
lines += ["## EXP015: CPU identification only", "",
          "| Cached tokens | Old hash-only ms | Exact ms | Radix ms | NumPy ms | Query-rope build ms |",
          "|---:|---:|---:|---:|---:|---:|"]
for n,r in summary["exp015_r2_cpu"].items():
    lines.append(f"| {int(n):,} | {cell(r['old_hash_only'])} | {cell(r['exact'])} | {cell(r['radix'])} | {cell(r['numpy'])} | {cell(r['query_rope_build'])} |")
r2_large = summary["exp015_r2_cpu"]["1000000"]
lines += ["", f"At 1M cached tokens, query-rope construction alone costs **{cell(r2_large['query_rope_build'])} ms**, "
          f"excluded from the **{cell(r2_large['exact'])} ms** exact query. "
          "This eager construction cost is part of any end-to-end integration that builds each incoming query from flat tokens. "
          "Prebuilt exact identification can remain competitive while eager query-rope construction is much more expensive. "
          "The 4-A100 model-serving, TTFT and NVML energy measurements were not rerun; their historical engine-cache observations "
          "must not be presented as measurements of the corrected exact identification path. No new GPU/energy claim is supported here.", ""]

if (HERE / "edit_results.json").exists():
    edit=json.loads((HERE/"edit_results.json").read_text())
    summary["edit_comparison"]={}
    for n in [10000,100000,1000000]:
        for shared in [False,True]:
            key=f"{n}/{shared}"
            summary["edit_comparison"][key]={arm:stats([c['median_ms_per_cycle'] for run in edit['runs'] for c in run['cells']
              if c['bytes']==n and c['shared_history']==shared and c['arm']==arm]) for arm in ['hashrope_exact','hashrope_candidate','ropey_exact','ropey_rehash_candidate']}
    lines += ["## New edit-plus-target comparison workload", "",
              "This is separate from EXP020's maintained-fingerprint read. Each batch performs 50 seeded ASCII churn edits and "
              "compares each edited value with a prebuilt target: 50% equal, 50% different at the last byte. "
              "Independent targets and targets sharing edit ancestry are separate regimes. Target construction/destruction and warmup are untimed. "
              "Seeds 1/2/3 × three fresh processes, three batches per process; per-process medians are aggregated below. "
              "Every 50-bit equality-result sequence is checked exactly. Arm order rotates across seed/invocation. "
              "Ropey 1.6.1 uses its own exact equality; rehash uses Hashrope's polynomial hash over materialized Ropey bytes.", "",
              "| Bytes / target history | Hashrope exact ms/cycle | Ropey exact | Hashrope candidate | Ropey rehashed candidate |",
              "|---|---:|---:|---:|---:|"]
    for key,r in summary["edit_comparison"].items():
        n, shared = key.split("/")
        history = "shared ancestry" if shared == "True" else "independent"
        lines.append(f"| {int(n):,} / {history} | {cell(r['hashrope_exact'])} | {cell(r['ropey_exact'])} | {cell(r['hashrope_candidate'])} | {cell(r['ropey_rehash_candidate'])} |")
    independent = summary["edit_comparison"]["1000000/False"]
    shared = summary["edit_comparison"]["1000000/True"]
    lines += ["", "At 1M bytes, the ratio of mean Hashrope exact time to mean candidate time is "
              f"**{independent['hashrope_exact']['mean']/independent['hashrope_candidate']['mean']:.2f}×** for independent targets and "
              f"**{shared['hashrope_exact']['mean']/shared['hashrope_candidate']['mean']:.2f}×** for targets with shared ancestry. "
              "These ratios include the same edits in both arms; they measure this explicit workload, not a general equality bound. "
              "Ropey exact equality is faster in both 1M-byte regimes: Hashrope exact takes "
              f"**{independent['hashrope_exact']['mean']/independent['ropey_exact']['mean']:.2f}×** its time for independent targets and "
              f"**{shared['hashrope_exact']['mean']/shared['ropey_exact']['mean']:.2f}×** with shared ancestry."]
    lines += ["", "A maintained fingerprint read is still cheap, but it is not an exact identity decision. "
              "The historical EXP020 speedup applies to its original maintained-versus-recomputed fingerprint workload; "
              "exact equality needs the separate target-dependent measurements above.", ""]
if (HERE/"bio_results.json").exists():
    bio=json.loads((HERE/"bio_results.json").read_text())
    summary["bio"] = bio
    lines += ["## Biology consumers", "", "`bio_results.json` and `bio_paired.py` compare old/exact region, gene-report and resistance-panel "
              f"decisions on {bio['patient_records_used']:,} stored HIV patient records plus 16 independent equal controls against byte oracles. "
              "These are nine paired repetitions in one process, not nine independent process estimates. "
              "Construction is excluded; fixed 64-byte leaves and the data selection are recorded. "
              "The results do not reproduce the old chr22 gene-scale speedup claims or establish clinical efficacy.", "",
              "| Whole corpus batch | Old hash-only median ms | Exact median ms | Exact / old |",
              "|---|---:|---:|---:|"]
    for name, row in bio["results"].items():
        lines.append(f"| {name.replace('_', ' ')} | {row['median_ms']['old_hash_only']:.4f} | {row['median_ms']['exact']:.4f} | {row['exact_over_old_ratio']:.2f}× |")
    lines += ["", f"Constructing the {bio['cases']:,} cases took {bio['construction_ms']:.2f} ms and is excluded from these batch timings. "
              "The longer exact decisions are a real correction cost on these independently built inputs; tiny-codon panel queries amortize that cost differently.", ""]
lines += ["## Unchanged or unavailable claims", "",
          "EXP019 branch sharing/COW, raw edits, and compression do not depend on treating matching fingerprints as content equality. "
          "Their structural claims are not invalidated by this collision; their full historical benchmarks were not rerun here. "
          "Changes to default hash parameters, all large-scale/real-serving effects, and an end-to-end production reuse integration remain outside these measurements.", "",
          "## Reproduce", "", "The Python environment needs NumPy, Torch, Transformers and a locally cached GPT-2 tokenizer; "
          "tokenization deliberately runs offline. Rust requires Cargo and rustc. From the Hashrope checkout, "
          "capture the pinned paper baseline into a new temporary folder:", "",
          "```powershell", "$paperRepo = 'E:/data/code/claudecode/hashropeictaipaper'",
          "$paperBaseline = Join-Path $env:TEMP ('hashrope-paper-baseline-' + [guid]::NewGuid().ToString('N'))",
          "New-Item -ItemType Directory -Path $paperBaseline | Out-Null",
          "$paperArchive = Join-Path $paperBaseline 'baseline.zip'",
          "git -C $paperRepo archive --format=zip --output=$paperArchive c9d8902e0f443f1a5ba12da78f087d5dbd2cad68",
          "Expand-Archive -LiteralPath $paperArchive -DestinationPath $paperBaseline",
          "python benchmarks/consumer_exactness/paper_cpu.py --paper $paperRepo --baseline $paperBaseline",
          "python benchmarks/consumer_exactness/competitive_extended.py --paper $paperRepo",
          "python benchmarks/consumer_exactness/edit_compare.py",
          "python benchmarks/consumer_exactness/bio_paired.py --samples 2000 --output benchmarks/consumer_exactness/bio_results.json",
          "python benchmarks/consumer_exactness/summarize.py", "```", "",
          "The Rust runner uses the checked-in dependency lock, fetches pinned Ropey 1.6.1 when necessary, and records Cargo.lock, binary hash and compile-time source hashes. "
          "A reused --binary must match its build manifest and current sources. "
          "Token arrays are ignored, regenerable preparation artifacts; their SHA-256 values remain in provenance. "
          "No historical latest JSON or figure is overwritten by these commands.", ""]
(HERE/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
(HERE/"README.md").write_text("\n".join(lines),encoding="utf-8")
print(HERE/"README.md")
