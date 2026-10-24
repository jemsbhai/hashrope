"""Render results.json as a concise Markdown report with no external dependencies."""
import json
from pathlib import Path
import statistics


directory = Path(__file__).resolve().parent
results = json.loads((directory / "results.json").read_text(encoding="utf-8"))
lines = ["# Exact comparison and profile-check benchmarks", "",
         "Measured on 2026-10-23. These are local microbenchmarks, not universal performance guarantees.", "",
         f"Baseline: `{results['baseline_ref']}`. Current source SHA-256 values and all raw samples are in `results.json`.",
         f"The measured sources remained unchanged during the run: **{results['source_unchanged_during_run']}**.", "",
         f"Environment: {results.get('cpu', {}).get('Name', 'CPU not captured')}; "
         f"{results['platform']}; Python {results['python']}; {results['rustc']}. "
         "Rust uses Cargo release defaults. No CPU pinning, frequency locking, or process isolation was applied.", "",
         f"{results['rounds']} rounds alternate baseline/current and current/baseline. "
         "Each round has seven samples after warm-up. Tables give the median of the round medians; "
         "brackets show the minimum and maximum round median. Times are microseconds per operation. "
         "Fixture creation is excluded unless the operation is construction. GC is disabled only during Python timing. "
         "All comparison fixtures are checked against known exact answers before timing.", "",
         "## Existing operations", "",
         "Python default construction remains one leaf; its other nonrepeat fixtures are manually built "
         "from 512-byte leaves. Rust construction already chunks to 512-byte leaves. "
         "Python construction reuses a hasher; Rust construction includes a fresh arena and its disposal. "
         "Do not compare absolute timings across languages. "
         "Concat/split retain their newly allocated nodes in Rust's arena; Python discards each returned result. "
         "Substring and exponent caches are warmed. These results do not cover lazy first-query materialization.", ""]


def times(runtime, variant, name):
    return [record["median_ns"] / 1000 for run in results[runtime + "_runs"]
            if run["variant"] == variant for record in run["measurements"] if record["name"] == name]


def formatted(values):
    return f"{statistics.median(values):.3f} [{min(values):.3f}–{max(values):.3f}]"


for runtime in ("python", "rust"):
    lines += [f"### {runtime.title()}", "", "| Operation | Baseline µs | Current µs | Current / baseline |",
              "|---|---:|---:|---:|"]
    names = [x["name"] for x in results[runtime + "_runs"][0]["measurements"] if x["name"].startswith("core/")]
    for name in names:
        before, after = times(runtime, "baseline", name), times(runtime, "current", name)
        lines.append(f"| {name.removeprefix('core/')} | {formatted(before)} | {formatted(after)} | "
                     f"{statistics.median(after)/statistics.median(before):.3f}× |")
    lines.append("")

lines += ["## Exact comparisons", "",
          "Hash-only LCP is the existing prefix-hash binary-search strategy using current code; "
          "its speed does not make it an exact correctness substitute. The equality API is also timed "
          "in the raw results. Byte inputs deterministically cycle through 0–255; results need validation "
          "on application traces before deployment decisions. "
          "Shared late mismatch reuses the original prefix after replacing its last byte. "
          "The Repeat case compares independent `ab × 2^40` with `abab × 2^39`; "
          "the DAG case independently self-concatenates two seeds 40 times. Both represent 2 TiB without expansion.", ""]
for runtime in ("python", "rust"):
    lines += [f"### {runtime.title()}", "", "| Input | Exact LCP µs | Hash-only LCP µs |", "|---|---:|---:|"]
    for label in ["1048576/same_root", "1048576/independent_equal", "1048576/early_mismatch",
                  "1048576/late_mismatch", "1048576/shared_late_mismatch", "repeat_2TiB", "self_concat_DAG_2TiB"]:
        lines.append(f"| {label} | {formatted(times(runtime, 'current', 'compare/' + label + '/exact_lcp'))} | "
                     f"{formatted(times(runtime, 'current', 'compare/' + label + '/hash_only_lcp'))} |")
    lines.append("")

lines += ["## Memory", "",
          "Python `tracemalloc` starts after creating the input bytes. Construction figures include "
          "new nodes, chunk slices, hasher/cache, and profiles; the single-leaf input buffer is reused "
          "and therefore excluded. These are live traced Python allocations, not process RSS.", "",
          "| Python construction | Baseline retained bytes | Current retained bytes | Change |", "|---|---:|---:|---:|"]
for name in ("single_leaf_64KiB", "chunked_64KiB"):
    values = {}
    for variant in ("baseline", "current"):
        values[variant] = statistics.median([x["retained_bytes"] for run in results["python_runs"]
            if run["variant"] == variant for item in run["memory"] if item["name"] == name for x in item["runs"]])
    lines.append(f"| {name} | {values['baseline']:.0f} | {values['current']:.0f} | "
                 f"{values['current'] - values['baseline']:+.0f} ({100*(values['current']/values['baseline']-1):+.2f}%) |")
lines += ["", "The profile adds one pointer per Python node; the shared profile/hasher has a small fixed cost. "
          "Scratch peaks below include temporary memo/stack entries and slice buffers. "
          "They exclude existing ropes and warm caches.", "", "| Python query | Current peak scratch bytes |", "|---|---:|"]
for record in next(x for x in results["python_runs"] if x["variant"] == "current")["scratch_memory"]:
    values = [item["peak_bytes"] for run in results["python_runs"] if run["variant"] == "current"
              for item in run["scratch_memory"] if item["name"] == record["name"]]
    lines.append(f"| {record['name'].removeprefix('compare/')} | {statistics.median(values):.0f} |")
lines += ["", "Rust reports layout/node counts only, not allocator totals or exact-query scratch peaks:", "",
          "| Rust metric | Baseline | Current |", "|---|---:|---:|"]
for metric in ("node_inner_bytes", "arena_inline_bytes", "nodes_for_64KiB"):
    before = next(x for x in results["rust_runs"] if x["variant"] == "baseline")["memory"][metric]
    after = next(x for x in results["rust_runs"] if x["variant"] == "current")["memory"][metric]
    lines.append(f"| {metric} | {before} | {after} |")
lines += ["", "## Reproduce", "", "From the repository root:", "", "```powershell",
          f"python benchmarks/comparison_safety/run_benchmarks.py --baseline-ref {results['baseline_ref']}",
          "python benchmarks/comparison_safety/summarize.py", "```", "",
          "The runner extracts the baseline into a fresh temporary directory and builds separate Rust binaries. "
          "It does not change either checkout or inspect the computational-biology submodule. "
          "The individual `bench_python.py` and `bench_rust.py` launchers also accept `--package PATH` "
          "and `--exact`, allowing focused reruns. Remove temporary benchmark directories manually when no longer needed.", ""]
(directory / "README.md").write_text("\n".join(lines), encoding="utf-8")
print(directory / "README.md")
