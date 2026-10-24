"""Build isolated timing/counter diagnostics; never modify the production core.

Example: python edit_profile.py --core-root .cache/baseline/hashrope --label baseline
The default output is beside this script. Existing results are never replaced.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import statistics
import subprocess
import tempfile
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
COUNTERS = ["hash_calls", "hashed_bytes", "concat_hash_calls", "power_calls",
            "power_inline_hits", "power_tree_hits", "power_misses",
            "split_leaf_clone_bytes", "split_leaf_calls"]
ALLOCATIONS = ["allocations", "allocation_bytes", "reallocations", "reallocation_requested_bytes", "deallocations"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_once(text, old, new):
    assert text.count(old) == 1, (old, text.count(old))
    return text.replace(old, new, 1)


def instrument(core):
    source = core / "src/polynomial_hash.rs"
    code = source.read_text(encoding="utf-8")
    insertions = [
        ("pub fn hash(&self, data: &[u8]) -> u64 {", "crate::diagnostic::count(0, 1); crate::diagnostic::count(1, data.len() as u64);"),
        ("pub fn hash_concat(&mut self, h_a: u64, len_b: u64, h_b: u64) -> u64 {", "crate::diagnostic::count(2, 1);"),
        ("pub fn power(&mut self, n: u64) -> u64 {", "crate::diagnostic::count(3, 1);"),
        ("if (n as usize) < INLINE_CACHE_SIZE {", "crate::diagnostic::count(4, 1);"),
        ("if let Some(&v) = self.overflow_cache.get(&n) {", "crate::diagnostic::count(5, 1);"),
    ]
    for marker, counter in insertions:
        code = replace_once(code, marker, marker + "\n        #[cfg(feature = \"diagnostic-counters\")] { " + counter + " }")
    code = replace_once(code, "// Binary exponentiation", "#[cfg(feature = \"diagnostic-counters\")] crate::diagnostic::count(6, 1);\n        // Binary exponentiation")
    source.write_text(code, encoding="utf-8")
    source = core / "src/rope.rs"
    code = source.read_text(encoding="utf-8")
    start = code.index("    fn split_inner(")
    end = code.index("    fn split_repeat(", start)
    body = code[start:end]
    # Detect the baseline clone explicitly. Optimized borrowed-leaf variants
    # have no clone charge; their leaf call count is still collected below.
    cloned = "match self.nodes[id as usize].clone()" in body
    marker = "NodeInner::Leaf { data, .. } => {"
    if marker in body:
        counters = "crate::diagnostic::count(8, 1);"
        if cloned:
            counters += " crate::diagnostic::count(7, data.len() as u64);"
        body = replace_once(body, marker, marker + "\n                #[cfg(feature = \"diagnostic-counters\")] { " + counters + " }")
    elif "Kind::Leaf => {" in body:
        body = replace_once(body, "Kind::Leaf => {", "Kind::Leaf => {\n                #[cfg(feature = \"diagnostic-counters\")] crate::diagnostic::count(8, 1);")
    else:
        raise AssertionError("Review new split representation before collecting counters")
    source.write_text(code[:start] + body + code[end:], encoding="utf-8")
    source = core / "src/lib.rs"
    source.write_text(source.read_text(encoding="utf-8") + '\n#[cfg(feature = "diagnostic-counters")]\npub mod diagnostic;\n', encoding="utf-8")
    (core / "src/diagnostic.rs").write_text('''use core::sync::atomic::{AtomicBool, AtomicU64, Ordering::Relaxed};
static ACTIVE: AtomicBool = AtomicBool::new(false);
static VALUES: [AtomicU64; 9] = [const { AtomicU64::new(0) }; 9];
pub fn count(index: usize, amount: u64) {
    if ACTIVE.load(Relaxed) { VALUES[index].fetch_add(amount, Relaxed); }
}
pub fn start() {
    for value in &VALUES { value.store(0, Relaxed); }
    ACTIVE.store(true, Relaxed);
}
pub fn end() -> [u64; 9] {
    ACTIVE.store(false, Relaxed);
    core::array::from_fn(|i| VALUES[i].load(Relaxed))
}
''', encoding="utf-8")
    manifest = core / "Cargo.toml"
    manifest.write_text(manifest.read_text(encoding="utf-8") + '\n[features]\ndiagnostic-counters=[]\n', encoding="utf-8")
    return {"split_leaf_clone_present": cloned}


def build(core_root, label, output):
    core_root = core_root.resolve()
    sources = [HERE / "edit_profile.rs", Path(__file__), core_root / "packages/rust/Cargo.toml",
               *sorted((core_root / "packages/rust/src").glob("*.rs"))]
    source_hashes = {str(p): sha(p) for p in sources}
    project = Path(tempfile.mkdtemp(prefix="hashrope-edit-profile-"))
    copied = project / "core"
    shutil.copytree(core_root / "packages/rust/src", copied / "src")
    shutil.copyfile(core_root / "packages/rust/Cargo.toml", copied / "Cargo.toml")
    features = instrument(copied)
    (project / "src").mkdir()
    shutil.copyfile(HERE / "edit_profile.rs", project / "src/main.rs")
    (project / "Cargo.toml").write_text(
        '[package]\nname="edit-profile"\nversion="0.0.0"\nedition="2021"\n'
        '[dependencies]\nhashrope={path="core"}\nropey="=1.6.1"\n'
        '[features]\ncounters=["hashrope/diagnostic-counters"]\n', encoding="utf-8")
    original_lock = ROOT / "benchmarks/consumer_exactness/edit_compare.Cargo.lock"
    (project / "Cargo.lock").write_text(original_lock.read_text().replace('name = "exact-edit-comparison"', 'name = "edit-profile"'), encoding="utf-8")
    binaries = {}
    for mode in ["timing", "counters"]:
        command = ["cargo", "build", "--release", "--locked", "--manifest-path", str(project / "Cargo.toml")]
        if mode == "counters": command += ["--features", "counters"]
        subprocess.run(command, check=True)
        destination = project / f"edit_profile_{mode}.exe"
        shutil.copyfile(project / "target/release/edit-profile.exe", destination)
        binaries[mode] = str(destination)
    assert source_hashes == {str(p): sha(p) for p in sources}, "Input source changed during diagnostic build"
    metadata = {"core_root": str(core_root), "project": str(project), "binaries": binaries,
                "binary_sha256": {mode: sha(Path(path)) for mode, path in binaries.items()},
                "source_sha256": source_hashes, "instrumentation": features,
                "instrumented_core_sha256": {str(p.relative_to(copied)): sha(p) for p in copied.rglob("*") if p.is_file()},
                "cargo_lock": (project / "Cargo.lock").read_text(), "cargo_lock_sha256": sha(project / "Cargo.lock"),
                "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip()}
    (output / f"edit_profile.{label}.build.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    with zipfile.ZipFile(output / f"edit_profile.{label}.sources.zip", "x", zipfile.ZIP_DEFLATED) as archive:
        for name in ["edit_profile.py", "edit_profile.rs"]:
            archive.write(HERE / name, name)
        archive.write(output / f"edit_profile.{label}.build.json", "build.json")
        for file in [core_root / "packages/rust/Cargo.toml", *sorted((core_root / "packages/rust/src").glob("*.rs"))]:
            archive.write(file, "core/" + file.relative_to(core_root / "packages/rust").as_posix())
    return metadata


def parse(text):
    rows = []
    for line in text.splitlines():
        values = line.split("\t")
        assert len(values) == 22, values
        row = dict(zip(["bytes", "shared", "arm", "phase", "sample", "operations", "nanoseconds", "new_nodes",
                        *ALLOCATIONS, *COUNTERS], values))
        row = {key: value if key in ["arm", "phase"] else value == "true" if key == "shared" else int(value)
               for key, value in row.items()}
        row["ns_per_operation"] = row["nanoseconds"] / row["operations"]
        rows.append(row)
    return rows


def summarize(result):
    cells = {}
    for run in result["runs"]:
        for row in run["rows"]:
            key = f"{row['bytes']}/{row['shared']}/{row['arm']}/{row['phase']}"
            cell = cells.setdefault(key, {"timing": {}, "counters": {}})
            cell[run["mode"]].setdefault(str(run["seed"]), []).append(row)
    output = {}
    for key, cell in cells.items():
        medians = [statistics.median(row["ns_per_operation"] for row in rows) for rows in cell["timing"].values()]
        counts = [rows[0] for rows in cell["counters"].values()]
        output[key] = {"mean_process_median_ns": statistics.mean(medians),
                       "std_process_median_ns": statistics.stdev(medians) if len(medians) > 1 else 0,
                       "process_medians_ns": medians,
                       "mean_counters_per_operation": {name: statistics.mean(row[name] / row["operations"] for row in counts)
                                                       for name in ["new_nodes", *ALLOCATIONS, *COUNTERS]}}
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--core-root", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument("--reuse-build", action="store_true")
    parser.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--sizes", type=int, nargs="+", default=[10000, 100000, 1000000])
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    destination = output / f"edit_profile.{args.label}.json"
    if destination.exists(): raise SystemExit(f"Preserve existing diagnostic: {destination}")
    if args.reuse_build:
        metadata = json.loads((output / f"edit_profile.{args.label}.build.json").read_text())
        assert Path(metadata["core_root"]) == args.core_root.resolve()
        assert all(sha(Path(p)) == digest for p, digest in metadata["source_sha256"].items())
        assert all(sha(Path(metadata["binaries"][mode])) == digest for mode, digest in metadata["binary_sha256"].items())
    else:
        metadata = build(args.core_root, args.label, output)
    if args.build_only: return
    result = {"label": args.label, "build": metadata, "seeds": args.seeds, "sizes": args.sizes, "samples_per_process": args.samples,
              "method": "Diagnostic decomposition; 50 seeded edits, sizes 10K/100K/1M, independent and shared targets; timing/counter binaries separate.",
              "limitations": [
                  "Comparison-only retains all 50 edited roots and repeats their 50 half-equal comparisons three times; this changes locality and Ropey COW during untimed preparation.",
                  "Phase times are separate experiments, not additive components or a subtraction-based attribution of combined time.",
                  "Lazy Hashrope is a diagnostic that suppresses hash materialization, not the eager public workload; differences also include cache/allocation effects.",
                  "Allocation bytes are requested sizes, not peak live bytes, RSS or retained arena capacity; reallocations include the entire requested new block size.",
                  "Counter-build elapsed times are instrumented and excluded from timing summaries; one counter sample per seed.",
                  "Three process medians across three workload seeds are diagnostic dispersion, not confirmatory confidence intervals.",
              ], "runs": []}
    for mode in ["timing", "counters"]:
        for seed in args.seeds:
            print(f"edit diagnostic {args.label}: {mode} seed={seed}", flush=True)
            raw = subprocess.check_output([metadata["binaries"][mode], str(seed), str(args.samples if mode == "timing" else 1), ",".join(map(str, args.sizes))], text=True)
            raw_path = output / f"edit_profile.{args.label}.{mode}.s{seed}.tsv"
            raw_path.write_text(raw, encoding="utf-8")
            result["runs"].append({"mode": mode, "seed": seed, "raw_sha256": sha(raw_path), "rows": parse(raw)})
    result["summary"] = summarize(result)
    destination.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(destination)


if __name__ == "__main__":
    main()
