"""Controlled DAG diagnostics, never instrumentation in production.

The bypass ablation removes only the ordinary-tree dispatch. Counter builds
explain work; their instrumented times are not performance estimates.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NAMES = ["lcp_calls", "normalize_visits", "memo_lookups", "memo_hits", "memo_inserts", "leaf_pairs", "leaf_bytes", "periodic_witnesses", "structural_period_visits"]


def build(label, core, bypass=False, counters=False):
    project = HERE / ".cache" / ("dag_profile_" + label + ("_counters" if counters else ""))
    library = project / "hashrope"
    shutil.copytree(core / "packages/rust/src", library / "src", dirs_exist_ok=True)
    for name in ["Cargo.toml", "README.md"]:
        shutil.copyfile(core / "packages/rust" / name, library / name)
    path = library / "src/rope.rs"
    source = path.read_text()
    if bypass:
        start = source.index("fn exact_lcp(")
        end = source.index("\nstruct FastComparison", start)
        source = source[:start] + """fn exact_lcp(left: &Arena, right: &Arena, a: RangeStart, b: RangeStart, length: u64) -> u64 {
    ExactComparison::new(left, right).lcp(a, b, length)
}
""" + source[end:]
    main = (HERE / "dag_profile.rs").read_text()
    if counters:
        source += """
use core::sync::atomic::{AtomicU64, Ordering};
static DAG_PROFILE: [AtomicU64; 9] = [const { AtomicU64::new(0) }; 9];
fn dag_tick(n: usize, v: u64) { DAG_PROFILE[n].fetch_add(v, Ordering::Relaxed); }
pub fn dag_reset() { for x in &DAG_PROFILE { x.store(0, Ordering::Relaxed); } }
pub fn dag_read() -> [u64;9] { core::array::from_fn(|i| DAG_PROFILE[i].load(Ordering::Relaxed)) }
"""
        method = "lcp_inner" if "fn lcp_inner(" in source else "lcp"
        signature = f"fn {method}(&mut self, a: RangeStart, b: RangeStart, length: u64) -> u64 {{"
        source=source.replace(signature, signature + "\n        dag_tick(0,1);")
        source=source.replace("fn normalize(arena: &Arena, mut range: RangeStart, length: u64) -> RangeStart {\n        loop {", "fn normalize(arena: &Arena, mut range: RangeStart, length: u64) -> RangeStart {\n        loop { dag_tick(1,1);")
        source=source.replace("if self.equal_ranges.contains(&key) {", "dag_tick(2,1); if self.equal_ranges.contains(&key) { dag_tick(3,1);")
        source=source.replace("self.equal_ranges.insert(key);", "dag_tick(4,1); self.equal_ranges.insert(key);")
        source=source.replace("self.equal_ranges.insert(original);", "dag_tick(4,1); self.equal_ranges.insert(original);")
        source=source.replace("if x != y {", "dag_tick(5,1); dag_tick(6,remaining); if x != y {")
        source=source.replace("if witness < remaining {", "if witness < remaining { dag_tick(7,1);")
        source=source.replace("if witness < length {", "if witness < length { dag_tick(7,1);")
        source=source.replace("let mut repeated = false;\n        loop {", "let mut repeated = false;\n        loop { dag_tick(8,1);")
        main=main.replace("// COUNTER_RESET", "hashrope::rope::dag_reset();").replace("COUNTER_READ", "hashrope::rope::dag_read()")
    else:
        main=main.replace("COUNTER_READ", "[0u64;9]")
    path.write_text(source)
    (project / "src").mkdir(exist_ok=True)
    (project / "src/main.rs").write_text(main)
    (project / "Cargo.toml").write_text('[package]\nname="dag-profile"\nversion="0.0.0"\nedition="2021"\n[dependencies]\nhashrope={path="hashrope"}\n')
    subprocess.run(["cargo","build","--release","--manifest-path",str(project / "Cargo.toml")],check=True,capture_output=True)
    output=subprocess.check_output([str(project / "target/release/dag-profile.exe")],text=True)
    rows=[json.loads(line) for line in output.splitlines()]
    for row in rows: row["counters"]=dict(zip(NAMES,row["counters"]))
    return {"label":label,"instrumented":counters,"source_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"rows":rows}


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--core",type=Path,default=HERE / ".cache/baseline/hashrope")
    p.add_argument("--label",default="before")
    p.add_argument("--bypass",action="store_true")
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    result={"method":__doc__,"variants":[build(args.label,args.core,args.bypass,c) for c in [False,True]]}
    args.output.write_text(json.dumps(result,indent=2))
    print(args.output)
