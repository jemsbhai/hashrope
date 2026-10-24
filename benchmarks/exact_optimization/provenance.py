"""Bind measurements/resume to actual driver, harness, inputs and workloads."""
import argparse
import ast
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PAPER=ROOT.parent/"hashropeictaipaper"
CONSUMER=ROOT/"benchmarks/consumer_exactness"
BIO=ROOT/"usecases/computational-biology/code"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make(kind,driver):
    files=[driver]
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
    document={"kind":kind,"source_sha256":sources,"data_sha256":{str(p.resolve()):sha(p) for p in data},
              "workload":config,"python":platform.python_version(),
              "runtime_packages":{name:importlib.metadata.version(name) for name in ["numpy","torch","transformers"]},
              "core_archives":{"baseline":sha(HERE/"baseline.zip"),"optimized":sha(HERE/"optimized.zip")}}
    document["fingerprint"]=hashlib.sha256(json.dumps(document,sort_keys=True).encode()).hexdigest()
    return document


def save(document,directory):
    destination=directory/"contexts"/(document["fingerprint"]+".json")
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(document,indent=2))


def source_variant(document):
    """Recover the one observed pre-guard harness exactly, without changing timings."""
    source=CONSUMER/"paper_cpu.py"
    canonical=source.read_bytes()
    block=('    selected_core = Path(os.environ.get("HASHROPE_BENCH_CORE", Path(__file__).resolve().parents[2]))\n'
           '    assert Path(hashrope.__file__).resolve() == (selected_core / "packages/python/hashrope/__init__.py").resolve()\n'
           '    assert Path(sys.modules["hashrope.rope"].__file__).resolve() == (selected_core / "packages/python/hashrope/rope.py").resolve()\n')
    assert canonical.count(block.encode())==1
    earlier=canonical.replace(block.encode(),b"",1)
    earlier_sha=hashlib.sha256(earlier).hexdigest()
    assert earlier_sha=="456ba33f4a2050d7828e0207935f454ba693f337ba780766420dd2d9c7ee4001"
    assert hashlib.sha256(canonical).hexdigest()=="807bed81945cd9f7fff369eaec42d9b39773a1373c1a561519a5e8cdd96f4e72"
    artifact=HERE/"paper_cpu_pre_import_guards.txt"
    canonical_artifact=HERE/"paper_cpu_import_guards.txt"
    if artifact.exists():
        assert artifact.read_bytes()==earlier
    else:
        artifact.write_bytes(earlier)
    if canonical_artifact.exists():
        assert canonical_artifact.read_bytes()==canonical
    else:
        canonical_artifact.write_bytes(canonical)
    def timed_functions(raw):
        return {n.name:ast.dump(n,include_attributes=False) for n in ast.walk(ast.parse(raw))
                if isinstance(n,ast.FunctionDef) and n.name in {"measure","token_cell"}}
    assert timed_functions(earlier)==timed_functions(canonical)
    assert len(timed_functions(earlier))==2
    record={"source_path":str(source.resolve()),"canonical_sha256":sha(source),
            "variant_sha256":earlier_sha,"source_artifact":artifact.name,"canonical_source_artifact":canonical_artifact.name,
            "difference":"Canonical adds only one selected-core assignment and two import-path assertions before fixture preparation; every remaining byte is identical.",
            "removed_startup_block":block,"timed_nested_function_asts_identical":["measure","token_cell"],
            "allowed_worker":{"optimization_core":"baseline","seed":42,"invocation":0},
            "expected_import_path":str((HERE/".cache/baseline/hashrope/packages/python/hashrope/__init__.py").resolve())}
    result=dict(document)
    result["untimed_source_variants"]=[record]
    result.pop("fingerprint")
    result["fingerprint"]=hashlib.sha256(json.dumps(result,sort_keys=True).encode()).hexdigest()
    return result


def verify_worker_sources(document,worker):
    for name,digest in worker.get("source_sha256",{}).items():
        if name not in document["source_sha256"] or digest==document["source_sha256"][name]:
            continue
        variants=[v for v in document.get("untimed_source_variants",[]) if v["source_path"]==name and v["variant_sha256"]==digest]
        assert len(variants)==1,name
        variant=variants[0]
        assert all(worker.get(k)==v for k,v in variant["allowed_worker"].items())
        assert Path(worker["hashrope_file"]).resolve()==Path(variant["expected_import_path"]).resolve()
        recovered=(HERE/variant["source_artifact"]).read_bytes()
        assert hashlib.sha256(recovered).hexdigest()==digest
        canonical=(HERE/variant["canonical_source_artifact"]).read_bytes()
        assert hashlib.sha256(canonical).hexdigest()==variant["canonical_sha256"]
        assert canonical.replace(variant["removed_startup_block"].encode(),b"",1)==recovered


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--capture",action="store_true")
    parser.add_argument("--stamp",action="store_true")
    args=parser.parse_args()
    directory=HERE/"results"
    driver=HERE/"executed_matrix_source.txt"
    for kind in ["rust","paper","extended","bio","micro"]:
        document=make(kind,driver)
        record=HERE/("execution_context_"+kind+".json")
        if args.capture:
            assert not record.exists(),record
            record.write_text(json.dumps(document,indent=2))
            save(document,directory)
        if args.stamp:
            captured=json.loads(record.read_text())
            assert captured==document,"Harness/input changed during timing: "+kind
            if kind=="paper":
                document=source_variant(document)
                save(document,directory)
            paths=[directory/"microbench.json"] if kind=="micro" else sorted(directory.glob("*/"+kind+"/run_*.json"))
            assert len(paths)==(1 if kind=="micro" else 2 if kind=="bio" else 18),(kind,len(paths))
            for path in paths:
                value=json.loads(path.read_text())
                # The worker records independently observed source hashes.
                verify_worker_sources(document,value)
                value["context_fingerprint"]=document["fingerprint"]
                path.write_text(json.dumps(value,indent=2))
            print("Verified and bound",kind)


if __name__=="__main__":
    main()
