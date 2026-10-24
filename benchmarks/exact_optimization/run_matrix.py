"""Alternate baseline/optimized cores on unchanged consumer workloads.

Final timings start only after --freeze has captured the reviewed optimized
core. Historical measurements are never overwritten. Every worker is a fresh
process. --resume accepts only results for the same immutable core snapshots.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CONSUMER = ROOT / "benchmarks/consumer_exactness"
PAPER = ROOT.parent / "hashropeictaipaper"
RESULTS = HERE / "results"
BASELINE = HERE / ".cache/baseline/hashrope"
OPTIMIZED = HERE / ".cache/optimized/hashrope"
HISTORICAL_REF = "c9d8902e0f443f1a5ba12da78f087d5dbd2cad68"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def core_hashes(root):
    files = sorted((root / "packages/rust/src").glob("*.rs")) + sorted((root / "packages/python/hashrope").glob("*.py"))
    return {p.relative_to(root).as_posix(): sha(p) for p in files}


def freeze():
    archive = HERE / "optimized.zip"
    if archive.exists():
        raise SystemExit("An optimized snapshot already exists; preserve it instead of silently replacing")
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as output:
        for relative in ["packages/python/hashrope", "packages/rust/src", "packages/rust/Cargo.toml", "packages/rust/README.md"]:
            source = ROOT / relative
            for file in sorted(source.rglob("*")) if source.is_dir() else [source]:
                if file.is_file() and "__pycache__" not in file.parts:
                    output.write(file, "hashrope/" + file.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(archive) as source:
        source.extractall(OPTIMIZED.parent)
    manifest = {"archive_sha256": sha(archive), "core_sha256": core_hashes(OPTIMIZED),
                "snapshot_matches_working_core": core_hashes(OPTIMIZED) == core_hashes(ROOT)}
    assert manifest["snapshot_matches_working_core"]
    (HERE / "optimized_manifest.json").write_text(json.dumps(manifest, indent=2))
    print("Frozen optimized core:", manifest["archive_sha256"])


def restore():
    if not BASELINE.exists():
        with zipfile.ZipFile(HERE / "baseline.zip") as archive:
            archive.extractall(BASELINE.parents[0])
    if not OPTIMIZED.exists():
        with zipfile.ZipFile(HERE / "optimized.zip") as archive:
            archive.extractall(OPTIMIZED.parent)
    manifest = json.loads((HERE / "baseline_manifest.json").read_text())
    assert sha(HERE / "baseline.zip") == manifest["archive_sha256"]
    for name, info in manifest["files"].items():
        if name.startswith("hashrope/packages/"):
            assert sha(BASELINE.parent / name) == info["sha256"], name
    optimized = json.loads((HERE / "optimized_manifest.json").read_text())
    assert sha(HERE / "optimized.zip") == optimized["archive_sha256"]
    assert core_hashes(OPTIMIZED) == optimized["core_sha256"]


def prepare():
    historical = HERE / ".cache/hash_only_paper/src"
    historical.mkdir(parents=True, exist_ok=True)
    (historical / "lcp.py").write_bytes(subprocess.check_output(["git", "-C", str(PAPER), "show", f"{HISTORICAL_REF}:src/lcp.py"]))
    provenance = json.loads((CONSUMER / "paper_results/provenance.json").read_text())
    for seed in [42,43,44]:
        path = CONSUMER / f"paper_results/tokens_s{seed}.npy"
        if not path.exists():
            env = os.environ.copy()
            env.update(HF_HUB_OFFLINE="1", TRANSFORMERS_OFFLINE="1")
            command = "import sys,numpy as np;sys.path.insert(0,sys.argv[1]);from scripts.exp017_bench import tokenize_corpus;from pathlib import Path;np.save(sys.argv[3],np.asarray(tokenize_corpus(Path(sys.argv[1]),int(sys.argv[2])),dtype=np.int64))"
            subprocess.run([sys.executable, "-c", command, str(PAPER), str(seed), str(path)], env=env, check=True)
        assert sha(path) == provenance["data"][str(seed)]["tokens_sha256"]


def rust_build(label, core):
    output = subprocess.check_output([sys.executable, str(CONSUMER / "edit_compare.py"), "--build-only", "--core-root", str(core)], text=True)
    binary = Path(output.strip().splitlines()[-1])
    build = json.loads(binary.with_suffix(".build.json").read_text())
    directory = RESULTS / label
    directory.mkdir(parents=True, exist_ok=True)
    builds = directory / "rust_builds"
    builds.mkdir(exist_ok=True)
    (builds / (build["binary_sha256"] + ".json")).write_text(json.dumps(build, indent=2))
    return binary, build["binary_sha256"]


def microbench(cores, resume, context):
    target = RESULTS / "microbench.json"
    manifests = {label: core_hashes(core) for label,core in cores.items()}
    if target.exists():
        if not resume:
            raise SystemExit(f"Refusing to overwrite {target}")
        previous = json.loads(target.read_text())
        assert previous["core_sha256"] == manifests
        assert previous.get("context_fingerprint") == context["fingerprint"], "Resume harness/workload/data changed"
        return
    scripts = ROOT / "benchmarks/comparison_safety"
    sys.path.insert(0,str(scripts))
    from bench_rust import parse_output
    binaries = {label: Path(subprocess.check_output([sys.executable,str(scripts / "bench_rust.py"),
        "--package",str(core / "packages/rust"),"--exact","--build-only"],text=True).strip())
        for label,core in cores.items()}
    result = {"rounds":3,"samples_per_round":7,"core_sha256":manifests,"context_fingerprint":context["fingerprint"],
              "harness_sha256":{p.name:sha(p) for p in scripts.glob("bench_*") if p.is_file()},
              "binary_sha256":{label:sha(path) for label,path in binaries.items()},
              "python_runs":[],"rust_runs":[]}
    for round_id in range(3):
        order = ["baseline","optimized"] if round_id % 2 == 0 else ["optimized","baseline"]
        for label in order:
            print(f"Microbench Python: {label} round={round_id}",flush=True)
            run=json.loads(subprocess.check_output([sys.executable,str(scripts / "bench_python.py"),
                "--package",str(cores[label] / "packages/python"),"--exact"],text=True,timeout=180))
            run.update(variant=label,round=round_id)
            result["python_runs"].append(run)
        for label in order:
            print(f"Microbench Rust: {label} round={round_id}",flush=True)
            run=parse_output(subprocess.check_output([str(binaries[label])],text=True,timeout=180))
            run.update(variant=label,round=round_id)
            result["rust_runs"].append(run)
    target.write_text(json.dumps(result,indent=2))


def main():
    global RESULTS
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--kind", choices=["all","paper","extended","rust","bio","micro"], default="all")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--output-root", type=Path, default=RESULTS)
    args = parser.parse_args()
    RESULTS = args.output_root
    if args.freeze:
        freeze()
        return
    restore()
    prepare()
    if args.prepare_only:
        return
    cores = {"baseline": BASELINE, "optimized": OPTIMIZED}
    from provenance import make,save
    for kind in (["rust","paper","extended","bio","micro"] if args.kind == "all" else [args.kind]):
        context = make(kind,Path(__file__))
        save(context,RESULTS)
        if kind == "micro":
            microbench(cores,args.resume,context)
            continue
        binaries = {label: rust_build(label, core) for label,core in cores.items()} if kind == "rust" else {}
        seeds = [1,2,3] if kind == "rust" else [42,43,44]
        pairs = [(seed, invocation) for seed in seeds for invocation in range(3)] if kind != "bio" else [(42,0)]
        for seed, invocation in pairs:
            order = ["baseline", "optimized"] if (seed + invocation) % 2 == 0 else ["optimized", "baseline"]
            for label in order:
                directory = RESULTS / label / kind
                directory.mkdir(parents=True, exist_ok=True)
                target = directory / f"run_s{seed}_i{invocation}.json"
                if target.exists():
                    if not args.resume:
                        raise SystemExit(f"Refusing to overwrite {target}")
                    previous = json.loads(target.read_text())
                    assert previous.get("core_sha256") == core_hashes(cores[label]), "Resume result uses a different or unidentified core"
                    assert previous.get("context_fingerprint") == context["fingerprint"], "Resume harness/workload/data changed"
                    continue
                env = os.environ.copy()
                env["HASHROPE_BENCH_CORE"] = str(cores[label])
                print(f"{kind}: {label} seed={seed} invocation={invocation}", flush=True)
                if kind == "paper":
                    command = [sys.executable, str(CONSUMER / "paper_cpu.py"), "--paper", str(PAPER),
                        "--baseline", str(HERE / ".cache/hash_only_paper"), "--output", str(directory),
                        "--tokens", str(CONSUMER / "paper_results"), "--worker", "--seed", str(seed), "--invocation", str(invocation)]
                    subprocess.run(command, env=env, check=True)
                elif kind == "extended":
                    subprocess.run([sys.executable, str(CONSUMER / "competitive_extended.py"), "--paper", str(PAPER),
                        "--tokens", str(CONSUMER / "paper_results"), "--output", str(directory), "--worker",
                        "--seed", str(seed), "--invocation", str(invocation)], env=env, check=True)
                elif kind == "bio":
                    subprocess.run([sys.executable, str(CONSUMER / "bio_paired.py"), "--samples", "2000", "--runs", "9", "--output", str(target)], env=env, check=True)
                else:
                    binary, binary_sha = binaries[label]
                    output = subprocess.check_output([str(binary),str(seed),str(invocation)], text=True)
                    cells = {}
                    for line in output.splitlines():
                        _, n, shared, arm, _, ms = line.split("\t")
                        cell = cells.setdefault(f"{n}/{shared}/{arm}", {"bytes":int(n),"shared_history":shared=="true", "arm":arm,"samples_ms_per_cycle":[]})
                        cell["samples_ms_per_cycle"].append(float(ms))
                    for cell in cells.values():
                        cell["median_ms_per_cycle"] = statistics.median(cell["samples_ms_per_cycle"])
                    target.write_text(json.dumps({"seed":seed,"invocation":invocation,"binary_sha256":binary_sha,"cells":list(cells.values())},indent=2))
                data = json.loads(target.read_text())
                data["optimization_core"] = label
                data["core_sha256"] = core_hashes(cores[label])
                data["context_fingerprint"] = context["fingerprint"]
                for name,digest in data.get("source_sha256",{}).items():
                    if name in context["source_sha256"]:
                        assert digest==context["source_sha256"][name],name
                data["process_order"] = order
                target.write_text(json.dumps(data, indent=2))


if __name__ == "__main__":
    main()
