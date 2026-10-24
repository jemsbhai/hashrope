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
    assert {n:h for n,h in core_expected['baseline'].items() if n.startswith('packages/python/')} == {
        n:h for n,h in core_expected['optimized'].items() if n.startswith('packages/python/')}, 'Python control implementation changed'
    summary={"baseline":"completed round-two exact implementation","optimized":"frozen reviewed round-three Rust implementation",
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
            if kind!='bio':
                seeds=[1,2,3] if kind=='rust' else [42,43,44]
                assert {(w['seed'],w['invocation']) for w in workers}=={(seed,invocation) for seed in seeds for invocation in range(3)},(kind,label)
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
    micro_path=results/"microbench.json"
    if micro_path.exists():
        micro=json.loads(micro_path.read_text())
        assert micro["core_sha256"]==core_expected
        summary["context_fingerprints"]["micro"]=validate_context("micro",[micro],results,args.partial)
        summary["micro"]={}
        summary["micro_memory"]={}
        summary["micro_primary_aggregation"]="median of three round medians [minimum round median, maximum round median]"
        for language in ["python","rust"]:
            data=micro[language+"_runs"]
            assert all(sum(r["variant"]==label for r in data)==3 for label in ["baseline","optimized"])
            summary["micro"][language]={label:{m["name"]:stats([next(c for c in r["measurements"] if c["name"]==m["name"])["median_ns"]
                for r in data if r["variant"]==label]) for m in next(r for r in data if r["variant"]==label)["measurements"]}
                for label in ["baseline","optimized"]}
            summary["micro_memory"][language]={label:[{"memory":r["memory"],"scratch_memory":r.get("scratch_memory",[])}
                for r in data if r["variant"]==label] for label in ["baseline","optimized"]}
    elif not args.partial:
        raise AssertionError("Microbenchmarks not complete")

    summary['complete_matrix']=all(kind in summary for kind in ['rust','paper','extended','bio','micro'])
    lines=["# Edit and metadata optimization: round three", "",
           "The baseline is the completed second optimization round. This investigation isolates edit and metadata costs behind the remaining Rust-versus-Ropey gap, especially on small inputs. "
           "The [round-one report](../exact_optimization/README.md) and [round-two report](../exact_optimization_round2/README.md), their raw results and all original source archives remain preserved.", "",
           "## Method and scope", "",
           "Both measured core versions are immutable archives with file manifests. The same consumer runners, corpus inputs, builders, sizes, "
           "seed/invocation order, match/mismatch mix, warmups and fixed repetition counts are retained. "
           "Nine-process tables report mean ± sample SD across process medians, in milliseconds. Real-pair and R2 queries are first averaged within each process. "
           "Speedup ratios use the reported means. Biology retains nine paired repetitions within one process per core, with recorded adaptive batching; its inner counts may differ across cores. "
           "Microbenchmarks report the median of three round medians and their minimum/maximum; these ranges are not confidence intervals.", "",
           "Python is unchanged between the two core snapshots and is rerun as a control. Differences in Python consumer timings therefore do not measure the Rust optimization. "
           "Preparation remains outside query timing. CPU-heavy team work is paused for final timing, on the same Windows 11 / Intel Core i9-14900HX workstation, "
           "without CPU affinity or frequency locking. No GPU serving, TTFT, NVML energy or chromosome-scale experiment is rerun.", ""]
    if not summary['complete_matrix']:
        lines += ['This is a partial aggregate while the remaining final benchmark cohorts run. Only complete cohorts are included below.', '']
    if 'rust' in summary:
        rb,ro=summary['rust']['baseline'],summary['rust']['optimized']
        independent,shared='1000000/False','1000000/True'
        independent_gain=1-ro[independent]['hashrope_exact']['mean']/rb[independent]['hashrope_exact']['mean']
        shared_gain=1-ro[shared]['hashrope_exact']['mean']/rb[shared]['hashrope_exact']['mean']
        change=lambda gain:f"falls {gain:.1%}" if gain>=0 else f"rises {-gain:.1%}"
        small_b,small_o=rb['10000/False']['hashrope_exact']['mean'],ro['10000/False']['hashrope_exact']['mean']
        small_shared_b,small_shared_o=rb['10000/True']['hashrope_exact']['mean'],ro['10000/True']['hashrope_exact']['mean']
        small_wins=[sum(x>y for x,y in zip(rb[k]['hashrope_exact']['values'],ro[k]['hashrope_exact']['values'])) for k in ['10000/False','10000/True']]
        outcome=(f"At 10K bytes, independent-target edit-and-compare changes from {small_b*1000:.3f} to {small_o*1000:.3f} µs "
                 f"({small_b/small_o:.2f}× speedup), and shared history from {small_shared_b*1000:.3f} to {small_shared_o*1000:.3f} µs "
                 f"({small_shared_b/small_shared_o:.2f}×), with {small_wins[0]}/9 and {small_wins[1]}/9 paired improvements respectively. "
                 f"Exact comparisons plus edits still take {small_o/ro['10000/False']['ropey_exact']['mean']:.2f}× and "
                 f"{small_shared_o/ro['10000/True']['ropey_exact']['mean']:.2f}× Ropey's time at 10K bytes. "
                 f"In the repeated 1M-byte workload, independent-target exact time {change(independent_gain)} "
                 f"and is {ro[independent]['hashrope_exact']['mean']/ro[independent]['ropey_exact']['mean']:.2f}× Ropey's time. "
                 f"Shared-history time {change(shared_gain)} and takes {ro[shared]['hashrope_exact']['mean']/ro[shared]['ropey_exact']['mean']:.1%} of Ropey's time in this fixture.")
        if 'micro' in summary:
            key='compare/self_concat_DAG_2TiB/exact_lcp'
            mb,mo=summary['micro']['rust']['baseline'][key],summary['micro']['rust']['optimized'][key]
            outcome+=f" The 2 TiB self-concatenation DAG's LCP changes from {mb['median']/1000:.3f} to {mo['median']/1000:.3f} µs (median of three round medians)."
        lines[4:4]=[outcome, '']
    if (HERE/"VALIDATION.md").exists():
        lines += ["[VALIDATION.md](VALIDATION.md) records the correctness checks and explicit experimental exclusions.", ""]
    if (HERE/"unchanged_code.json").exists():
        checks=json.loads((HERE/"unchanged_code.json").read_text())
        lines += [f"[unchanged_code.json](unchanged_code.json) verifies {len(checks['module_checksums'])} complete modules, "
                  f"{len(checks['python_ast'])} Python definitions and storage/default declarations against the baseline. "
                  f"It audits {len(checks['rust_functions'])} existing Rust production functions, recording {len(checks['intentional_changes'])} selected changes and "
                  f"{len(checks['new_functions'])} addition{'s' if len(checks['new_functions'])!=1 else ''} explicitly; all other production function bodies and existing signatures remain unchanged. "
                  "Source checks establish change scope; differential baseline tests establish the preserved fingerprint semantics.", ""]
    if (HERE/'edit_profile.md').exists():
        lines += ['## Profiling evidence', '',
                  '[edit_profile.md](edit_profile.md) separates eager/lazy edits, exact comparisons, combined operations and fingerprint primitives. '
                  'Its three-seed process diagnostics retain the original sizes and add a clearly labeled exploratory 1K cell. '
                  'Separately prepared phases have different locality and retained-state behavior; their timings cannot be added or subtracted for exact attribution.', '',
                  'At 10K independent targets, diagnostic eager editing falls from 5.065 to 1.983 µs while comparison-only changes from 2.238 to 2.116 µs. '
                  'The exact comparator is unchanged. Leaf hashing still processes about 709 bytes in 6.92 calls per edit; eliminating redundant source-leaf copies removes '
                  '2.96 fresh allocation requests and 692.40 requested bytes per edit. Node creation and arena growth do not decrease. '
                  'Those allocation figures describe requested traffic, not peak or retained memory. Warmed power-cache misses are zero, so cold exponentiation is not the measured bottleneck.', '',
                  'Default-profile eight-byte polynomial grouping reduces the diagnostic 512-byte hash from 2,168 to 169 ns, but the one-byte hash increases '
                  'from 1.830 to 2.253 ns (0.423 ns). The generic/custom-profile path remains scalar. '
                  '[REVIEW.md](REVIEW.md) gives the arithmetic bounds and independent compatibility review. '
                  '[hash_kernel_diagnostic.json](hash_kernel_diagnostic.json) and [edit_variants.json](edit_variants.json) preserve primitive/combined ablations, '
                  'including the rejected extra modular-arithmetic specialization. Their diagnostic repetitions and source provenance remain distinct from the final matrix below.', '']
    if (HERE/'DAG_PROFILE.md').exists():
        lines += ['## Profiling evidence', '',
                  '[DAG_PROFILE.md](DAG_PROFILE.md) records the compressed-DAG investigation and rejected per-recursion shortcut. '
                  'The dispatcher-only ablation does not explain the previous 15% regression, so its cause remains unresolved. '
                  'The selected query-level structural-period proof reduces the equal-DAG diagnostic from 404 normalization visits to 84 normalization plus 82 period-proof visits, '
                  '81 memo lookups to one and seven query allocations to one, while verifying the same two-byte witness. '
                  'These work counters explain a new optimization; single-process diagnostic timings are separate from the final cohorts below.', '']
    if (HERE/'fast_path_findings.md').exists():
        lines += ['[fast_path_findings.md](fast_path_findings.md) records the independent-range investigation: direct-child traversal, portable eight-byte mismatch localization, '
                  'and a fixed-stack frontier for cross-arena spans of at least 8 KiB. Same-arena comparisons retain recursive shared-subtree skips. '
                  'The frontier uses 1,600 bytes of frontier storage on the stack on this 64-bit target, excluding other locals, and no heap allocation; '
                  'its diagnostic runs are separate from the final matrix. A 512-byte leaf with a first-byte mismatch measured about 0.7 ns slower in the combined diagnostic variant, '
                  'so the leaf-localization improvement is not universal.', '']
    lines += ['Exact comparisons remain independent of fingerprints. Existing accepted custom profiles can produce different stored fingerprints for identical bytes in different tree shapes: '
              'with modulus 17 and base 2, bytes [0, 13] hash to 16 as one leaf and 0 as a concatenation under the legacy reduction. '
              'Rejecting equality on a fingerprint mismatch would therefore break existing exact semantics. This round preserves those profiles and fingerprint values.', '']
    if "rust" in summary:
        b,o=summary["rust"]["baseline"],summary["rust"]["optimized"]
        lines += ["## Rust edit plus exact target comparison", "",
                  "Each timed batch performs 50 seeded churn edits and compares the edited rope with a prebuilt target: 50% equal and 50% different at the last byte. "
                  "Every boolean result is checked. ASCII fixtures make byte and Ropey character offsets coincide. "
                  "Target preparation/destruction is excluded. Ropey is pinned to 1.6.1. Independent targets and shared ancestry remain separate regimes. "
                  "Both fingerprint-candidate arms use Hashrope from their measured core version, including Ropey's rehash-candidate arm, so they can benefit from this round's hash optimization. "
                  "Ropey exact uses its unchanged equality implementation.", "",
                  "| Bytes / history | Baseline exact | Optimized exact | Baseline / optimized | Paired improvements | Ropey, optimized cohort |",
                  "|---|---:|---:|---:|---:|---:|"]
        for key,before in b.items():
            after=o[key]; n,shared=key.split("/")
            wins=sum(x>y for x,y in zip(before['hashrope_exact']['values'],after['hashrope_exact']['values']))
            after['paired_improvements']=wins
            lines.append(f"| {int(n):,} / {'shared' if shared=='True' else 'independent'} | {fmt(before['hashrope_exact'])} | "
                         f"{fmt(after['hashrope_exact'])} | {before['hashrope_exact']['mean']/after['hashrope_exact']['mean']:.2f}× | {wins}/9 | {fmt(after['ropey_exact'])} |")
        for key in ["1000000/False","1000000/True"]:
            case=o[key]
            ratio=case['hashrope_exact']['mean']/case['ropey_exact']['mean']
            lines += ["",f"At 1M bytes ({'shared ancestry' if key.endswith('True') else 'independent targets'}), optimized exact Hashrope/Ropey is {ratio:.3f}. "
                      f"The separate fingerprint-candidate control changes from {b[key]['hashrope_candidate']['mean']:.5f} to {case['hashrope_candidate']['mean']:.5f} ms/cycle."]
        small_before,small_after=b['10000/False'],o['10000/False']
        lines += ['', f"For 10K independent targets, the edit-plus-fingerprint-candidate arm is {fmt(small_before['hashrope_candidate'])}→{fmt(small_after['hashrope_candidate'])} ms, "
                  f"versus {fmt(small_after['ropey_exact'])} ms for Ropey edit plus exact comparison. The decomposition profiles separately measure edit and metadata costs.", '',
                  "These fixture-specific ratios and paired win counts do not establish statistical significance, a universal equality advantage or identical timings on every path.", ""]
        lines += ['Fingerprint-candidate arms from the same runs (milliseconds per cycle); these are distinct from exact equality:', '',
                  '| Bytes / history | Hashrope candidate, baseline | Hashrope candidate, optimized | Ropey rehash candidate, baseline | Ropey rehash candidate, optimized |',
                  '|---|---:|---:|---:|---:|']
        for key,before in b.items():
            after=o[key]; n,shared=key.split('/')
            lines.append(f"| {int(n):,} / {'shared' if shared=='True' else 'independent'} | {fmt(before['hashrope_candidate'])} | "
                         f"{fmt(after['hashrope_candidate'])} | {fmt(before['ropey_rehash_candidate'])} | {fmt(after['ropey_rehash_candidate'])} |")
        lines.append('')
    if "paper" in summary:
        b,o=summary["paper"]["baseline"],summary["paper"]["optimized"]
        lines += ["## Unchanged Python paper controls", "",
                  "EXP-005, mismatch halfway through independently built corpus-derived ropes; all other fractions remain in the raw results:", "",
                  "| Bytes | Baseline core | Optimized core |", "|---:|---:|---:|"]
        for n in [64000,256000,1000000,2000000,4000000,8000000]:
            lines.append(f"| {n:,} | {fmt(b['exp005'][f'{n}/0.5']['exact'])} | {fmt(o['exp005'][f'{n}/0.5']['exact'])} |")
        lines += ["", f"Log-log slopes fitted to the six cohort means are {b['exp005_f05_slope']:.3f} and {o['exp005_f05_slope']:.3f}. "
                  f"The historical ≤0.3 gate {'passes' if b['exp005_historical_slope_gate_passes'] else 'fails'} for the baseline and "
                  f"{'passes' if o['exp005_historical_slope_gate_passes'] else 'fails'} for the optimized cohort. "
                  "This gate remains descriptive; no logarithmic exact-LCP bound follows.", "",
                  "Controlled EXP-017:", "", "| Shared tokens | Baseline exact | Optimized-core exact | Radix, optimized cohort | NumPy, optimized cohort |",
                  "|---:|---:|---:|---:|---:|"]
        for n,row in o['exp017_controlled'].items():
            lines.append(f"| {int(n):,} | {fmt(b['exp017_controlled'][n]['exact'])} | {fmt(row['exact'])} | {fmt(row['radix'])} | {fmt(row['numpy'])} |")
        gate=o['strict_crossover']
        lines += ["",f"Optimized-cohort exact wins {o['paired_2M_wins']}/9 paired comparisons against Radix at 2M tokens. "
                  f"The literal original crossover-grid condition {'passes' if gate['passes'] else 'fails'}; stable mean ± one sample SD band separation starts at {gate['separation_from_tokens']} tokens. "
                  "These are dispersion bands, not confidence intervals or a significance test."]
        if gate['at_grid_lower_boundary']:
            lines += ["That is the smallest measured point, so the below-threshold condition is vacuous; an interior crossover is not observed."]
        lines += ["", "EXP-015 R2 CPU identification:", "", "| Cached tokens | Baseline exact | Optimized-core exact | Radix, optimized cohort | Excluded query build, optimized cohort |",
                  "|---:|---:|---:|---:|---:|"]
        for n,row in o['r2'].items():
            lines.append(f"| {int(n):,} | {fmt(b['r2'][n]['exact'])} | {fmt(row['exact'])} | {fmt(row['radix'])} | {fmt(row['query_rope_build'])} |")
        lines += ["", "An integration starting with flat tokens must account for query construction; the prebuilt-query figures exclude it.", ""]
    if "extended" in summary:
        b,o=summary['extended']['baseline'],summary['extended']['optimized']
        lines += ["## Unchanged Python real-pair and candidate-set controls", "",
                  "Original 200 ShareGPT/LMSYS pairs per dataset/seed. Inner timing order remains fixed: radix, exact Hashrope, NumPy. "
                  "K candidates retain one indexed radix lookup versus K pairwise comparisons with independently built Hashrope roots.", ""]
        for section,title in [('real_pairs','Conversational pairs'),('k_sweep','Shared tokens / K')]:
            lines += [title+':', '', '| Cell | Baseline exact | Optimized-core exact | Radix, optimized cohort | NumPy, optimized cohort |', '|---|---:|---:|---:|---:|']
            for name,row in o[section].items():
                lines.append(f"| {name} | {fmt(b[section][name]['hashrope'])} | {fmt(row['hashrope'])} | {fmt(row['radix'])} | {fmt(row['flat'])} |")
            lines.append('')
    if "bio" in summary:
        b,o=summary['bio']['baseline'],summary['bio']['optimized']
        lines += ['## Unchanged Python biology controls', '',
                  'All 2,000 stored HIV records plus 16 independent equal controls, checked against byte oracles. Whole-corpus median batch times below are milliseconds, '
                  'from nine within-process paired repetitions; they are not nine fresh-process estimates. Adaptive inner counts are normalized in the timings.', '',
                  '| Workload | Baseline exact | Optimized-core exact | Old-hash control, baseline / optimized | Inner repetitions, baseline / optimized |',
                  '|---|---:|---:|---:|---:|']
        for name,row in o['results'].items():
            before=b['results'][name]
            lines.append(f"| {name} | {before['median_ms']['exact']:.4f} | {row['median_ms']['exact']:.4f} | "
                         f"{before['median_ms']['old_hash_only']:.4f} / {row['median_ms']['old_hash_only']:.4f} | "
                         f"{before['workload_repetitions_per_timed_sample']} / {row['workload_repetitions_per_timed_sample']} |")
        lines += ['', 'Construction is excluded. These are timing controls, not clinical-efficacy or chromosome-scale measurements.', '']
    if 'micro' in summary:
        lines += ['## Structural and operation microbenchmarks', '',
                  'Same fixtures, three alternating rounds × seven samples. Values are median [minimum, maximum] of three round medians, in microseconds. '
                  'These round ranges are not confidence intervals; all individual samples and scratch-memory measurements remain in results/microbench.json.', '',
                  'Rust microbench independent roots are built separately within one Arena and use the recursive same-arena exact path. '
                  'Their 1 MiB input repeats the byte sequence 0 through 255; the independent-target churn arm uses 1M edited ASCII bytes and separately allocated target arenas, '
                  'with the frontier path for long cross-arena ranges. These comparison timings are not interchangeable or additive. '
                  'The construction microbench includes Arena::new and arena destruction inside timing, whereas churn excludes initial setup.', '']
        b,o=summary['micro']['rust']['baseline'],summary['micro']['rust']['optimized']
        show=lambda k,label:f"{label} {b[k]['median']/1000:.3f}→{o[k]['median']/1000:.3f} µs"
        lines += ['The changed Rust construction/hash/edit paths improve in these fixtures: '+', '.join([
            show('core/from_bytes_64KiB','64 KiB full arena construction/destruction'),show('core/split_chunked_64KiB','split'),
            show('core/substr_warm_64KiB','warm substring fingerprint')])+'. '
            'The exact engines are unchanged, and their timing movements are mixed: '+', '.join([
            show('compare/1048576/independent_equal/exact_lcp','independent 1 MiB LCP'),
            show('compare/1048576/independent_equal/exact_equal','equality'),
            show('compare/1048576/late_mismatch/exact_lcp','late-mismatch LCP')])+'. '
            'Overlapping round ranges and unchanged engine source do not establish either a new exact-comparison optimization or timing equivalence. '
            'Compressed DAG/Repeat cases remain symbolic. Hash-only LCP rows are retained performance controls and remain collision-unsafe.', '']
        for language,cohorts in summary['micro'].items():
            lines += [language.capitalize()+':', '', '| Fixture | Baseline | Optimized |', '|---|---:|---:|']
            for name,before in cohorts['baseline'].items():
                after=cohorts['optimized'][name]
                form=lambda v:f"{v['median']/1000:.4f} [{v['min']/1000:.4f}, {v['max']/1000:.4f}]"
                lines.append(f"| {name} | {form(before)} | {form(after)} |")
            lines.append('')
        memories=summary['micro_memory']['rust']
        assert all(row['memory']==memories['baseline'][0]['memory'] for rows in memories.values() for row in rows)
        memory=memories['baseline'][0]['memory']
        lines += [f"Rust layout remains NodeInner {memory['node_inner_bytes']} bytes / Arena inline {memory['arena_inline_bytes']} bytes, "
                  f"excluding heap allocations, with {memory['nodes_for_64KiB']} nodes for 64 KiB. Source compatibility does not imply identical operation timing. "
                  "Compressed-path, construction/hash/edit and unchanged Python results must be assessed with their recorded round variability.", '']
    lines += ['## Provenance and reproduction', '',
              'Every result binds its immutable core, driver, consumer helpers, input hashes, workload and runtime context. '
              'Resume and aggregation reject mismatched identities. Rust results also reference compiled-binary hashes and locked dependency manifests. '
              'Python worker imports select and check the intended frozen package. Previous reports/results are retained unchanged.', '']
    if (HERE/'preserved_results.json').exists():
        preserved=json.loads((HERE/'preserved_results.json').read_text())
        assert all(row['unchanged'] for row in preserved.values())
        lines += [f"[preserved_results.json](preserved_results.json) verifies {len(preserved)} prior benchmark source/artifact files remain byte-identical.", '']
    if (HERE/'final_source_verification.json').exists():
        lines += ['[final_source_verification.json](final_source_verification.json) verifies the final working source against the measured optimized archive, including any documentation-only edits.', '']
    if (HERE/'final_audit.json').exists():
        lines += ['[final_audit.json](final_audit.json) records the complete matrix, context identities, verified resume and alternate-output aggregation checks. '
                  'The final revised paper compiled successfully; the completion record is in VALIDATION.md and the paper LOGBOOK.', '']
    lines += ['Use a fresh output directory for a full reproduction:', '', '```powershell',
              'python benchmarks/exact_optimization_round3/run_matrix.py --output-root benchmarks/exact_optimization_round3/reproduction_results',
              'python benchmarks/exact_optimization_round3/provenance.py --results-root benchmarks/exact_optimization_round3/reproduction_results',
              'python benchmarks/exact_optimization_round3/summarize.py --results-root benchmarks/exact_optimization_round3/reproduction_results --output-dir benchmarks/exact_optimization_round3/reproduction_report',
              '```', '',
              'The runner restores saved core archives, verifies corpus-token caches and can regenerate missing caches with the offline GPT-2 tokenizer. '
              'The original local paper/biology datasets remain required. Python requires NumPy, Torch and Transformers; Rust requires Cargo/rustc and the locked Ropey dependencies.', '']
    if output!=HERE:
        relative=Path(os.path.relpath(HERE,output)).as_posix()
        for name in ['VALIDATION.md','unchanged_code.json','final_source_verification.json','final_audit.json','DAG_PROFILE.md','fast_path_findings.md','preserved_results.json',
                     'edit_profile.md','REVIEW.md','hash_kernel_diagnostic.json','edit_variants.json']:
            lines=[line.replace(f']({name})',f']({relative}/{name})') for line in lines]
        lines=[line.replace('](../exact_optimization/README.md)',f']({relative}/../exact_optimization/README.md)') for line in lines]
        lines=[line.replace('](../exact_optimization_round2/README.md)',f']({relative}/../exact_optimization_round2/README.md)') for line in lines]
    (output/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    (output/'README.md').write_text('\n'.join(lines),encoding='utf-8')
    print(output/'README.md')


if __name__=='__main__':
    main()
