"""Aggregate paired frozen-core measurements without modifying prior reports."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
from provenance import verify_worker_sources

HERE=Path(__file__).resolve().parent


def stats(values):
    return {"mean":statistics.mean(values),"std":statistics.stdev(values) if len(values)>1 else 0,
            "median":statistics.median(values),"min":min(values),"max":max(values),"n":len(values),"values":values}


def fmt(value):
    return f"{value['mean']:.4f} ± {value['std']:.4f}"


def slope(xs,ys):
    x,y=list(map(math.log,xs)),list(map(math.log,ys))
    xm,ym=statistics.mean(x),statistics.mean(y)
    return sum((a-xm)*(b-ym) for a,b in zip(x,y))/sum((a-xm)**2 for a in x)


def validate_context(kind, workers, results, allow_unbound=False):
    """Reject mixed workloads or context files whose content does not match its ID."""
    fingerprints={worker.get("context_fingerprint") for worker in workers}
    if fingerprints=={None} and allow_unbound:
        return "pending final provenance verification"
    assert None not in fingerprints and len(fingerprints)==1,(kind,fingerprints)
    fingerprint=fingerprints.pop()
    document=json.loads((results/"contexts"/(fingerprint+".json")).read_text())
    assert document.pop("fingerprint")==fingerprint
    assert hashlib.sha256(json.dumps(document,sort_keys=True).encode()).hexdigest()==fingerprint
    assert document["kind"]==kind
    for label,digest in document["core_archives"].items():
        assert hashlib.sha256((HERE/(label+".zip")).read_bytes()).hexdigest()==digest
    for worker in workers:
        verify_worker_sources(document,worker)
    return fingerprint


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--partial",action="store_true")
    parser.add_argument("--results-root",type=Path,default=HERE/"results")
    parser.add_argument("--output-dir",type=Path,default=HERE)
    args=parser.parse_args()
    results=args.results_root.resolve()
    output=args.output_dir.resolve()
    output.mkdir(parents=True,exist_ok=True)
    baseline_manifest=json.loads((HERE/"baseline_manifest.json").read_text())
    optimized_manifest=json.loads((HERE/"optimized_manifest.json").read_text())
    core_expected={"optimized":optimized_manifest["core_sha256"]}
    core_expected["baseline"]={name:baseline_manifest["files"]["hashrope/"+name]["sha256"] for name in core_expected["optimized"]}
    summary={"baseline":"completed exact implementation before optimization","optimized":"frozen reviewed exact implementation",
             "uncertainty":"mean and sample SD across fresh-process medians; bio uses within-process paired medians",
             "results_root":str(results),"context_fingerprints":{}}
    runs={}
    for kind in ["rust","paper","extended","bio"]:
        cohort={label:[json.loads(p.read_text()) for p in sorted((results/label/kind).glob("run_*.json"))]
                for label in ["baseline","optimized"]}
        expected=1 if kind=="bio" else 9
        if all(len(cohort[label])==expected for label in cohort):
            runs[kind]=cohort
        elif not args.partial:
            raise AssertionError((kind,{label:len(v) for label,v in cohort.items()}))
    for kind,cohort in runs.items():
        summary["context_fingerprints"][kind]=validate_context(kind,[w for group in cohort.values() for w in group],results,args.partial)
        summary[kind]={}
        for label,workers in cohort.items():
            assert all(w["core_sha256"]==core_expected[label] and w["optimization_core"]==label for w in workers)
            summary[kind][label]={}
            if kind=="rust":
                for n in [10000,100000,1000000]:
                    for shared in [False,True]:
                        key=f"{n}/{shared}"
                        summary[kind][label][key]={arm:stats([c["median_ms_per_cycle"] for w in workers for c in w["cells"]
                            if c["bytes"]==n and c["shared_history"]==shared and c["arm"]==arm])
                            for arm in ["hashrope_exact","ropey_exact","hashrope_candidate","ropey_rehash_candidate"]}
            elif kind=="paper":
                for section,arms,key in [
                    ("exp005",["old_hash_only","exact","brute"],lambda c:f"{c['bytes']}/{c['fraction']}"),
                    ("exp017_controlled",["old_hash_only","exact","radix","numpy"],lambda c:str(c["expected_tokens"]))]:
                    summary[kind][label][section]={}
                    for row in workers[0][section]:
                        name=key(row)
                        summary[kind][label][section][name]={arm:stats([next(c for c in w[section] if key(c)==name)[arm]["median_ms"] for w in workers]) for arm in arms}
                r2={}
                for n in [1000,16000,64000,256000,1000000]:
                    r2[str(n)]={arm:stats([statistics.mean(c[arm]["median_ms"] for c in w["exp015_r2_cpu"] if c["expected_tokens"]==n) for w in workers]) for arm in ["old_hash_only","exact","radix","numpy"]}
                    r2[str(n)]["query_rope_build"]=stats([statistics.mean(c["query_rope_build_ms"] for c in w["exp015_r2_cpu"] if c["expected_tokens"]==n) for w in workers])
                summary[kind][label]["r2"]=r2
                sizes=[64000,256000,1000000,2000000,4000000,8000000]
                summary[kind][label]["exp005_f05_slope"]=slope(sizes,[summary[kind][label]["exp005"][f"{n}/0.5"]["exact"]["mean"] for n in sizes])
                summary[kind][label]["exp005_historical_slope_gate_passes"]=summary[kind][label]["exp005_f05_slope"]<=0.3
                controlled=summary[kind][label]["exp017_controlled"]
                grid=sorted(map(int,controlled))
                stable=[n for n in grid if all(controlled[str(m)]["exact"]["mean"]+controlled[str(m)]["exact"]["std"]
                    < controlled[str(m)]["radix"]["mean"]-controlled[str(m)]["radix"]["std"] for m in grid if m>=n)]
                crossover=min(stable) if stable else None
                summary[kind][label]["strict_crossover"]={"separation_from_tokens":crossover,"passes":bool(crossover) and all(
                    controlled[str(n)]["radix"]["mean"] <= controlled[str(n)]["exact"]["mean"]*1.001 for n in grid if n<crossover)}
                summary[kind][label]["strict_crossover"]["at_grid_lower_boundary"]=crossover==min(grid)
                summary[kind][label]["paired_2M_wins"]=sum(next(c for c in w["exp017_controlled"] if c["expected_tokens"]==2000000)["exact"]["median_ms"]
                    < next(c for c in w["exp017_controlled"] if c["expected_tokens"]==2000000)["radix"]["median_ms"] for w in workers)
            elif kind=="extended":
                summary[kind][label]["real_pairs"]={d:{arm:stats([w["real_pairs"][d][f"{arm}_ms_stats"]["mean"] for w in workers]) for arm in ["hashrope","radix","flat"]}
                    for d in workers[0]["real_pairs"]}
                summary[kind][label]["k_sweep"]={k:{arm:stats([w["k_sweep"][k][f"{arm}_median_ms"] for w in workers]) for arm in ["hashrope","radix","flat"]}
                    for k in workers[0]["k_sweep"]}
                assert all(c["all_correct"] for w in workers for c in w["real_pairs"].values())
                assert all(c["correctness_ok"] for w in workers for c in w["k_sweep"].values())
            else:
                summary[kind][label]=workers[0]
    lines=["# Exact comparison optimization: measured correction costs", "",
           "The frozen baseline is the already-correct exact implementation, including the prior collision fixes. "
           "Fingerprints and defaults are identical. These results measure optimization, separately from the earlier hash-only-to-exact migration.", "",
           "## Evidence and method", "",
           "The archives and manifests capture both complete measured core versions. "
           "`unchanged_code.json` verifies six hash/sliding/export modules, 33 Python definitions outside `_leaf_lcp`, "
           "and 31 Rust construction/hash/edit functions stayed unchanged. "
           "[VALIDATION.md](VALIDATION.md) records 639 passing tests, the no_std release check, and explicit experimental exclusions.", "",
           "Rust profiling of independent equal 1M-byte trees found 483 query allocations requesting 184,176 bytes, "
           "192,304 ordered-memo key comparisons and zero memo hits. The optimized ordinary-tree path performs zero heap allocations on this fixture, "
           "uses no memo entries, and reduces normalization visits from 19,534 to 15,626 while still checking 1M bytes. "
           "The compressed-repeat/DAG engine remains available through the bounded fallback. "
           "`rust_profile*.json` records diagnostic counters and ablations; their single-process timings are explanatory, not the final performance estimates.", "",
           "Python profiling isolated the unequal leaf's byte-by-byte scan. Native comparisons now narrow that unequal interval. "
           "The tree/cache/periodic logic is unchanged. [python_profile.md](python_profile.md) discloses the local tradeoff: "
           "a 64 KiB leaf mismatching at byte 64 increased from 2.161 to 4.091 µs, while its last-byte mismatch fell from 1,979.246 to 7.180 µs. "
           "Those are diagnostic offset samples; not every mismatch position improves.", "",
           "Final measurements alternate frozen baseline and optimized cores by seed/invocation, in fresh processes. "
           "The same data, builders, sizes, match/mismatch mix, inner arm order, warmups, oracles and repetitions are retained. "
           "Tables report mean ± sample SD over nine process medians, in milliseconds unless labeled otherwise; speedups are ratios of means. "
           "Real-pair and R2 queries are first averaged within each process. Construction is excluded from query latency. "
           "Biology retains nine paired repetitions within one process per core and its original adaptive batching. "
           "Microbenchmarks retain three alternating rounds and seven samples per round. "
           "Team load-generating tests were paused; this Windows 11 / Intel Core i9-14900HX workstation had no CPU affinity or frequency locking.", ""]
    if "rust" in summary:
        lines += ["## Rust: edit plus exact target comparison", "",
                  "Each batch performs 50 seeded churn edits and compares each edited value with a prebuilt target, "
                  "50% equal and 50% different at the last byte. All 50 boolean outcomes are checked. "
                  "Fixtures are ASCII, so Hashrope byte offsets and Ropey character offsets coincide. "
                  "Target preparation and destruction remain outside timing. Ropey is pinned to 1.6.1.", "",
                  "| Bytes / target history | Baseline exact | Optimized exact | Speedup | Ropey exact, optimized cohort |",
                  "|---|---:|---:|---:|---:|"]
        for key,b in summary["rust"]["baseline"].items():
            o=summary["rust"]["optimized"][key]
            n,shared=key.split("/")
            lines.append(f"| {int(n):,} / {'shared ancestry' if shared=='True' else 'independent'} | {fmt(b['hashrope_exact'])} | {fmt(o['hashrope_exact'])} | {b['hashrope_exact']['mean']/o['hashrope_exact']['mean']:.2f}× | {fmt(o['ropey_exact'])} |")
        b,o=summary["rust"]["baseline"],summary["rust"]["optimized"]
        independent,shared=o["1000000/False"],o["1000000/True"]
        lines += ["", f"At 1M bytes, optimized Hashrope remains **{independent['hashrope_exact']['mean']/independent['ropey_exact']['mean']:.2f}× slower than Ropey** on independent targets, "
                  f"and is **{shared['ropey_exact']['mean']/shared['hashrope_exact']['mean']:.2f}× faster** with shared ancestry. "
                  "This is an improvement with explicit limits, not a universal equality win. "
                  f"Independent candidate-only control time is {b['1000000/False']['hashrope_candidate']['mean']:.5f}→{o['1000000/False']['hashrope_candidate']['mean']:.5f} ms; "
                  f"shared candidate control is {b['1000000/True']['hashrope_candidate']['mean']:.5f}→{o['1000000/True']['hashrope_candidate']['mean']:.5f} ms.", ""]
    if "paper" in summary:
        b,o=summary["paper"]["baseline"],summary["paper"]["optimized"]
        lines += ["## Paper CPU consumers", "", "EXP-005, mismatching halfway through independently built corpus-derived ropes:", "",
                  "| Bytes | Baseline exact | Optimized exact |", "|---:|---:|---:|"]
        for n in [64000,256000,1000000,2000000,4000000,8000000]:
            lines.append(f"| {n:,} | {fmt(b['exp005'][f'{n}/0.5']['exact'])} | {fmt(o['exp005'][f'{n}/0.5']['exact'])} |")
        lines += ["", f"The exact f=0.5 log-log slope fitted to the six cohort means is {b['exp005_f05_slope']:.3f}→{o['exp005_f05_slope']:.3f}. "
                  f"The former slope ≤ 0.3 gate {'passes' if o['exp005_historical_slope_gate_passes'] else 'fails'} for the optimized core. "
                  "Optimization does not establish a logarithmic exact-comparison bound. Fractions 0 and 0.99 are retained in the raw/summary data.", "",
                  "Controlled EXP-017 prefix queries:", "", "| Shared tokens | Baseline exact | Optimized exact | Radix, optimized cohort | NumPy |", "|---:|---:|---:|---:|---:|"]
        for n,r in o["exp017_controlled"].items():
            lines.append(f"| {int(n):,} | {fmt(b['exp017_controlled'][n]['exact'])} | {fmt(r['exact'])} | {fmt(r['radix'])} | {fmt(r['numpy'])} |")
        gate=o["strict_crossover"]
        lines += ["",f"Optimized exact wins {o['paired_2M_wins']}/9 paired comparisons at 2M tokens. "
                  f"The original strict crossover criterion {'passes' if gate['passes'] else 'fails'} on this grid; "
                  f"separated mean ± 1 SD bands start at {gate['separation_from_tokens']} tokens. "
                  "This empirical grid result is separate from the retired logarithmic hash-call guard."]
        if gate["at_grid_lower_boundary"]:
            lines += ["The separation threshold equals the smallest measured point. The criterion's below-threshold condition is therefore vacuous; "
                      "this shows superiority over the sampled grid, not a measured interior crossover."]
        lines += ["", "EXP-015 R2 CPU identification:", "", "| Cached tokens | Baseline exact | Optimized exact | Query build, optimized cohort |", "|---:|---:|---:|---:|"]
        for n,r in o["r2"].items():
            lines.append(f"| {int(n):,} | {fmt(b['r2'][n]['exact'])} | {fmt(r['exact'])} | {fmt(r['query_rope_build'])} |")
        lines += ["", "Eager query construction is part of an integration that starts with flat tokens; the small prebuilt-query latency does not include it. "
                  "No 4-A100 serving, TTFT or NVML energy measurement was rerun, and no new GPU claim follows.", ""]
    if "extended" in summary:
        b,o=summary["extended"]["baseline"],summary["extended"]["optimized"]
        lines += ["## Real conversational pairs and candidate sets", "", "Original 200 pairs per dataset per seed, with independently built ropes:", "",
                  "| Dataset | Baseline exact | Optimized exact | Radix, optimized cohort | NumPy |", "|---|---:|---:|---:|---:|"]
        for name,r in o["real_pairs"].items():
            lines.append(f"| {name} | {fmt(b['real_pairs'][name]['hashrope'])} | {fmt(r['hashrope'])} | {fmt(r['radix'])} | {fmt(r['flat'])} |")
        lines += ["", "The optimized aggregate means reverse the previous real-pair radix advantage on both datasets; NumPy remains faster. "
                  "These are query-only results for the selected conversational pairs, not a universal indexing advantage."]
        lines += ["", "K candidates retain one radix lookup versus K pairwise Hashrope/NumPy queries; prefix values are shared but the Hashrope roots are independently constructed.", "",
                  "| Shared tokens / K | Baseline exact | Optimized exact | Radix, optimized cohort | NumPy |", "|---|---:|---:|---:|---:|"]
        for name,r in o["k_sweep"].items():
            lines.append(f"| {name} | {fmt(b['k_sweep'][name]['hashrope'])} | {fmt(r['hashrope'])} | {fmt(r['radix'])} | {fmt(r['flat'])} |")
        k100=o['k_sweep']['64000/100']
        lines += ["", f"At 64K/K=100, exact remains {k100['hashrope']['mean']/k100['radix']['mean']:.2f}× slower than the indexed radix lookup. "
                  "At 512K/K=10 the optimized exact mean increases, with overlapping sample variability; this cell does not establish an optimization gain."]
    if "bio" in summary:
        b,o=summary["bio"]["baseline"],summary["bio"]["optimized"]
        lines += ["", "## Biology consumers", "", "All 2,000 stored HIV records plus 16 independent equal controls, checked against byte oracles. "
                  "These are whole-corpus median batch times over nine within-process paired repetitions per core, not nine fresh-process estimates.", "",
                  "| Workload | Baseline exact ms | Optimized exact ms | Speedup | Inner repetitions baseline / optimized |", "|---|---:|---:|---:|---:|"]
        for name,r in o["results"].items():
            before=b["results"][name]
            lines.append(f"| {name} | {before['median_ms']['exact']:.4f} | {r['median_ms']['exact']:.4f} | {before['median_ms']['exact']/r['median_ms']['exact']:.2f}× | {before['workload_repetitions_per_timed_sample']} / {r['workload_repetitions_per_timed_sample']} |")
        lines += ["", "Adaptive batching is unchanged from the earlier protocol; values normalize the recorded inner repetitions. "
                  "The gene-report median increases, while the old-hash control also varies substantially. These samples do not support a uniform biology-consumer speedup. "
                  "Construction is excluded. These timings do not validate clinical efficacy or reproduce chromosome-scale claims.", ""]
    micro_path=results/"microbench.json"
    if micro_path.exists():
        micro=json.loads(micro_path.read_text())
        assert micro["core_sha256"]==core_expected
        summary["context_fingerprints"]["micro"]=validate_context("micro",[micro],results,args.partial)
        summary["micro"]={}
        summary["micro_primary_aggregation"]="median of three round medians [minimum round median, maximum round median]"
        summary["micro_memory"]={}
        lines += ["## Existing structural and operation microbenchmarks", "", "Same fixtures, three alternating rounds × seven samples. "
                  "Values below are the median [minimum, maximum] of the three within-process medians, matching the original microbenchmark aggregation. "
                  "These round ranges are not confidence intervals. Full memory/scratch samples remain in `results/microbench.json`.", ""]
        for language in ["python","rust"]:
            data=micro[language+"_runs"]
            metrics={label:{m['name']:stats([next(c for c in r['measurements'] if c['name']==m['name'])['median_ns'] for r in data if r['variant']==label])
                     for m in next(r for r in data if r['variant']==label)['measurements']} for label in ['baseline','optimized']}
            summary["micro"][language]=metrics
            summary["micro_memory"][language]={label:[{"memory":r["memory"],"scratch_memory":r.get("scratch_memory",[])}
                for r in data if r["variant"]==label] for label in ["baseline","optimized"]}
            lines += [f"{language.capitalize()} (µs per operation):", "", "| Fixture | Baseline | Optimized |", "|---|---:|---:|"]
            for name,b in metrics['baseline'].items():
                o=metrics['optimized'][name]
                before=f"{b['median']/1000:.4f} [{b['min']/1000:.4f}, {b['max']/1000:.4f}]"
                after=f"{o['median']/1000:.4f} [{o['min']/1000:.4f}, {o['max']/1000:.4f}]"
                lines.append(f"| {name} | {before} | {after} |")
        rust_memory=summary["micro_memory"]["rust"]
        assert all(r["memory"]==rust_memory["baseline"][0]["memory"] for group in rust_memory.values() for r in group)
        memory=rust_memory["baseline"][0]["memory"]
        lines += ["",f"Rust layout controls are identical: NodeInner {memory['node_inner_bytes']} bytes, "
                  f"Arena inline {memory['arena_inline_bytes']} bytes, {memory['nodes_for_64KiB']} nodes for 64 KiB. "
                  "Python retained/root-object measurements and comparison scratch samples are retained in the summary and raw data. "
                  "Hash/edit source is unchanged; control timings still vary with measurement conditions and compiled layout. "
                  "Source equivalence does not establish identical timing. "
                  "The 2 TiB paths retain their compact behavior, but Rust DAG LCP, split and substring controls show higher optimized medians in these rounds. "
                  "Python early/late/shared and same-root measurements also move in both directions. This is not a universal no-regression result.", ""]
    elif not args.partial:
        raise AssertionError("Microbenchmarks not complete")
    lines += ["## Provenance and reproduction", ""]
    if (HERE/"final_source_verification.json").exists():
        lines += ["[final_source_verification.json](final_source_verification.json) records the post-timing source check against the measured snapshot, including any documentation-only edits.", ""]
    paper_fingerprint=summary["context_fingerprints"].get("paper")
    variants=[]
    if paper_fingerprint and paper_fingerprint!="pending final provenance verification":
        variants=json.loads((results/"contexts"/(paper_fingerprint+".json")).read_text()).get("untimed_source_variants",[])
    harness_note=("The first baseline paper worker ran before three untimed startup import-path-check statements were added. "
                  "Both exact harness versions are preserved in `paper_cpu_pre_import_guards.txt` and `paper_cpu_import_guards.txt`, matching the recorded source hashes. "
                  "The paper context explicitly binds both versions, verifies the two timed nested-function ASTs are identical and every remaining source byte matches, "
                  "and verifies that worker's recorded import path selects the baseline snapshot. All later workers assert the package/module paths at startup. "
                  if variants else "Worker package/module paths are checked against their selected source snapshot. ")
    lines += ["The original in-memory driver is preserved in `executed_matrix_source.txt`; `run_matrix.py` adds hardened guards for future runs. "
              "`results/contexts/*.json` binds each result's executed driver, timed harness/helpers, data hashes, workloads, runtimes and core archives. "
              "The initial `execution_context_*.json` files are retained as captured; the final paper context additionally records the verified startup-only source variant. "
              "Each final result references its context fingerprint; Rust additionally records compiled-binary hashes and Cargo.lock. "
              +harness_note+
              "Future resume attempts reject driver, harness, data, workload or core mismatches. "
              "Use a new output directory for a fresh reproduction with the hardened driver:", "",
              "```powershell", "python benchmarks/exact_optimization/run_matrix.py --output-root benchmarks/exact_optimization/reproduction_results",
              "python benchmarks/exact_optimization/summarize.py --results-root benchmarks/exact_optimization/reproduction_results --output-dir benchmarks/exact_optimization/reproduction_report", "```", "",
              "The runner restores saved source archives and verifies/reuses the prior token caches, regenerating missing caches with the offline GPT-2 tokenizer. "
              "Python requires NumPy, Torch and Transformers; Rust requires Cargo/rustc and the locked Ropey dependencies. "
              "`preserved_results.json` verifies 25 original consumer_exactness and comparison_safety report/result files remain byte-for-byte preserved. "
              "Post-timing documentation-only source edits, if any, are recorded separately with AST/body equivalence to the measured snapshot.", ""]
    if output!=HERE:
        artifact_prefix=Path(os.path.relpath(HERE,output)).as_posix()
        lines=[line.replace("](VALIDATION.md)",f"]({artifact_prefix}/VALIDATION.md)")
                   .replace("](python_profile.md)",f"]({artifact_prefix}/python_profile.md)")
                   .replace("](final_source_verification.json)",f"]({artifact_prefix}/final_source_verification.json)") for line in lines]
    (output/"summary.json").write_text(json.dumps(summary,indent=2))
    (output/"README.md").write_text("\n".join(lines),encoding="utf-8")
    print(output/"README.md")


if __name__=="__main__":
    main()
