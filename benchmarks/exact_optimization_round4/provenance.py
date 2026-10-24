"""Bind measurements/resume to actual driver, harness, inputs and workloads."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PAPER=ROOT.parent/"hashropeictaipaper"
CONSUMER=ROOT/"benchmarks/consumer_exactness"
BIO=ROOT/"usecases/computational-biology/code"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make(kind,driver):
    files=[driver,Path(__file__)]
    data=[]
    common=[PAPER/"src/lcp.py",PAPER/"src/competitive.py",PAPER/"src/flatten.py",
            PAPER/"scripts/exp017_bench.py",PAPER/"third_party/sglang_radix_cache/radix_cache.py"]
    if kind=="rust":
        files += [CONSUMER/name for name in ["edit_compare.py","edit_compare.rs","edit_compare.Cargo.lock"]]
        config={"seeds":[1,2,3],"invocations":3,"samples":3,"cycles":50,"equal_fraction":0.5,"sizes":[10000,100000,1000000],"history":[False,True]}
    elif kind=="paper":
        files += [CONSUMER/"paper_cpu.py",PAPER/"src/exp015_identify.py",*common]
        data += [PAPER/f"data/raw/corpus_s{s}.txt" for s in [42,43,44]]
        data += [CONSUMER/f"paper_results/tokens_s{s}.npy" for s in [42,43,44]]
        config={"seeds":[42,43,44],"invocations":3,"reps":5,"byte_sizes":[64000,256000,1000000,2000000,4000000,8000000],"fractions":[0,.5,.99],"controlled_tokens":[1000,4000,16000,64000,128000,256000,512000,1000000,2000000],"r2_tokens":[1000,16000,64000,256000,1000000],"queries_per_r2_cell":3,"historical_helper_commit":"c9d8902e0f443f1a5ba12da78f087d5dbd2cad68"}
    elif kind=="extended":
        files += [CONSUMER/"competitive_extended.py",PAPER/"src/realpairs.py",*common]
        data += [PAPER/"data/canonical"/name for name in ["sharegpt_sample.jsonl","lmsys_sample.jsonl"]]
        data += [CONSUMER/f"paper_results/tokens_s{s}.npy" for s in [42,43,44]]
        config={"seeds":[42,43,44],"invocations":3,"reps":5,"real_pairs_per_dataset":200,"k_grid":{"64000":[1,10,100],"512000":[1,10]},"timed_arm_order":["radix","exact","numpy"]}
    elif kind=="bio":
        files += [CONSUMER/"bio_paired.py",*sorted((BIO/"src").rglob("*.py"))]
        data += [BIO/"data"/name for name in ["hiv1_hxb2.fa","hivdb_rt_expanded_aligned.fa"]]
        config={"patients":2000,"equal_controls":16,"within_process_paired_runs":9,"inner_repetitions":"original adaptive warmup method, recorded separately per workload/core","historical_consumer_commit":"206d1af783ca4c34626cfe1f9ce9a7bd7d4275cb"}
    elif kind=="micro":
        files += [ROOT/"benchmarks/comparison_safety"/name for name in ["bench_python.py","bench_rust.py","bench_rust.rs"]]
        config={"rounds":3,"samples_per_round":7,"exact_enabled_for_both":True}
    else:
        raise ValueError(kind)
    sources={str(p.resolve()):sha(p) for p in files}
    # Driver filename is irrelevant; actual captured contents determine identity.
    sources["executed_driver"]=sources.pop(str(driver.resolve()))
    environment_names=["PYTHONHASHSEED","OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS"]
    if kind in {"rust","micro"}:
        environment_names += ["RUSTFLAGS","CARGO_ENCODED_RUSTFLAGS","CARGO_BUILD_TARGET","CARGO_PROFILE_RELEASE_OPT_LEVEL","CARGO_PROFILE_RELEASE_LTO","CARGO_PROFILE_RELEASE_CODEGEN_UNITS"]
    runtime={"python_executable":sys.executable,"python_build":list(platform.python_build()),"python_implementation":platform.python_implementation(),
             "platform":platform.platform(),"machine":platform.machine(),"processor_identifier":os.environ.get("PROCESSOR_IDENTIFIER"),
             "environment":{name:os.environ.get(name) for name in environment_names}}
    if kind in {"rust","micro"}:
        runtime["rustc_verbose"]=subprocess.check_output(["rustc","--version","--verbose"],text=True).strip()
        runtime["cargo_version"]=subprocess.check_output(["cargo","--version"],text=True).strip()
    plan=json.loads((HERE/"measurement_plan.json").read_text())
    document={"kind":kind,"source_sha256":sources,"data_sha256":{str(p.resolve()):sha(p) for p in data},
              "workload":config,"python":platform.python_version(),
              "runtime":runtime,
              "runtime_packages":{name:importlib.metadata.version(name) for name in ["numpy","torch","transformers"]},
              "study_mode":plan["study_mode"],"measurement_plan_sha256":sha(HERE/"measurement_plan.json"),
              "manifest_sha256":{name:sha(HERE/name) for name in ["baseline_manifest.json","reference_manifest.json","optimized_manifest.json"]},
              "core_archive_files":plan["core_archive_files"],
              "core_archives":{label:sha(HERE/name) for label,name in plan["core_archive_files"].items()},
              "preserved_working_baseline_sha256":sha(HERE/"baseline.zip")}
    assert document["core_archives"]==plan["core_archives"]
    document["fingerprint"]=hashlib.sha256(json.dumps(document,sort_keys=True).encode()).hexdigest()
    return document


def save(document,directory):
    destination=directory/"contexts"/(document["fingerprint"]+".json")
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(document,indent=2))


def verify_worker_sources(document,worker):
    for name,digest in worker.get("source_sha256",{}).items():
        if name in document["source_sha256"]:
            assert digest==document["source_sha256"][name],name


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--results-root",type=Path,default=HERE/"results")
    args=parser.parse_args()
    for kind in ["rust","paper","extended","bio","micro"]:
        document=make(kind,HERE/"run_matrix.py")
        captured=json.loads((args.results_root/"contexts"/(document["fingerprint"]+".json")).read_text())
        assert captured==document,"Harness/input changed during timing: "+kind
        paths=[args.results_root/"microbench.json"] if kind=="micro" else sorted(args.results_root.glob("*/"+kind+"/run_*.json"))
        assert len(paths)==(1 if kind=="micro" else 2 if kind=="bio" else 18),(kind,len(paths))
        for path in paths:
            worker=json.loads(path.read_text())
            assert worker["context_fingerprint"]==document["fingerprint"]
            verify_worker_sources(document,worker)
            if kind=="rust":
                build=json.loads((path.parent.parent/"rust_builds"/(worker["binary_sha256"]+".json")).read_text())
                assert build["binary_sha256"]==worker["binary_sha256"]
                verify_worker_sources(document,build)
                assert hashlib.sha256(build["cargo_lock"].encode()).hexdigest()==build["cargo_lock_sha256"]
                assert build["cargo_lock_sha256"]==document["source_sha256"][str((CONSUMER/"edit_compare.Cargo.lock").resolve())]
                for relative,digest in worker["core_sha256"].items():
                    if relative.startswith("packages/rust/"):
                        recorded=[value for name,value in build["source_sha256"].items() if name.replace("\\","/").endswith("/"+relative)]
                        assert recorded==[digest],relative
        print("Verified",kind)


if __name__=="__main__":
    main()
