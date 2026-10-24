"""Instrument only an extracted baseline copy; no production instrumentation.

Counters quantify BTree memo, normalization and leaf work. A global allocator
counts query-only allocations. Timing comes from a separate uninstrumented
binary, with diagnostic no-memo / no-byte-check ablations on ordinary equal
trees only; those ablations are not correctness-preserving implementations.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
BASE = HERE / ".cache/baseline/hashrope"
NAMES = ["lcp_calls", "normalize_calls", "normalize_visits", "memo_lookups", "memo_hits",
         "memo_insert_calls", "memo_unique_inserts", "memo_key_comparisons",
         "leaf_pairs", "leaf_bytes", "shared_skips", "periodic_witnesses"]


def build(variant, optimized_core=None):
    project = HERE / ".cache" / ("profile_" + variant)
    project.mkdir(parents=True, exist_ok=True)
    library = project / "hashrope"
    core = optimized_core if variant.startswith("optimized") else BASE
    shutil.copytree(core / "packages/rust/src", library / "src", dirs_exist_ok=True)
    shutil.copyfile(core / "packages/rust/Cargo.toml", library / "Cargo.toml")
    shutil.copyfile(core / "packages/rust/README.md", library / "README.md")
    path = library / "src/rope.rs"
    source = path.read_text()
    instrumented = variant.endswith("counters")
    if instrumented:
        prelude = "\nuse core::sync::atomic::{AtomicU64, Ordering};\nstatic PROFILE: [AtomicU64; 12] = [const { AtomicU64::new(0) }; 12];\nfn tick(n: usize, value: u64) { PROFILE[n].fetch_add(value, Ordering::Relaxed); }\npub fn profile_reset() { for x in &PROFILE { x.store(0, Ordering::Relaxed); } }\npub fn profile_read() -> [u64;12] { core::array::from_fn(|i| PROFILE[i].load(Ordering::Relaxed)) }\n"
        source += prelude
        source = source.replace("type EqualRange = (NodeId, u64, NodeId, u64, u64);", """
#[derive(Clone,Copy,PartialEq,Eq)]
struct EqualRange((NodeId,u64,NodeId,u64,u64));
impl Ord for EqualRange { fn cmp(&self, other: &Self) -> core::cmp::Ordering {
    tick(7,1); self.0.cmp(&other.0)
} }
impl PartialOrd for EqualRange { fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> { Some(self.cmp(other)) } }
""")
        source = source.replace("(a.node, a.offset, b.node, b.offset, length)\n", "EqualRange((a.node, a.offset, b.node, b.offset, length))\n")
        source = source.replace("fn normalize(arena: &Arena, mut range: RangeStart, length: u64) -> RangeStart {\n        loop {", "fn normalize(arena: &Arena, mut range: RangeStart, length: u64) -> RangeStart {\n        tick(1,1);\n        loop {\n            tick(2,1);")
        source = source.replace("fn lcp(&mut self, a: RangeStart, b: RangeStart, length: u64) -> u64 {", "fn lcp(&mut self, a: RangeStart, b: RangeStart, length: u64) -> u64 {\n        tick(0,1);")
        source = source.replace("fn lcp(&mut self, mut a: RangeStart, mut b: RangeStart, length: u64) -> Result<u64, u64> {", "fn lcp(&mut self, mut a: RangeStart, mut b: RangeStart, length: u64) -> Result<u64, u64> {\n        tick(0,1);")
        source = source.replace("fn normalize(arena: &Arena, mut range: RangeStart, length: u64, visits_left: &mut usize) -> Option<RangeStart> {\n        loop {", "fn normalize(arena: &Arena, mut range: RangeStart, length: u64, visits_left: &mut usize) -> Option<RangeStart> {\n        tick(1,1);\n        loop {\n            tick(2,1);")
        source = source.replace("return Ok(if x == y {", "tick(8,1); tick(9,remaining);\n                return Ok(if x == y {")
        source = source.replace("if self.same_arena && a.node == b.node && a.offset == b.offset {", "if self.same_arena && a.node == b.node && a.offset == b.offset {\n                tick(10,1);")
        source = source.replace("if self.equal_ranges.contains(&key) {", "tick(3,1);\n            if self.equal_ranges.contains(&key) {\n                tick(4,1);")
        source = source.replace("if x != y {", "tick(8,1); tick(9,remaining);\n                    if x != y {")
        source = source.replace("if witness < remaining {", "if witness < remaining {\n                        tick(11,1);")
        source = source.replace("self.equal_ranges.insert(key);", "tick(5,1); if self.equal_ranges.insert(key) { tick(6,1); }")
        source = source.replace("self.equal_ranges.insert(original);", "tick(5,1); if self.equal_ranges.insert(original) { tick(6,1); }")
    elif variant in ["no_memo", "no_memo_no_leaf_bytes"]:
        source = source.replace("self.equal_ranges.contains(&key)", "false")
        source = source.replace("self.equal_ranges.insert(key);", "")
        source = source.replace("self.equal_ranges.insert(original);", "")
        if variant == "no_memo_no_leaf_bytes":
            source = source.replace("if x != y {", "if false {")
    elif variant == "no_leaf_bytes":
        source = source.replace("if x != y {", "if false {")
    path.write_text(source)
    main = (HERE / "rust_profile.rs").read_text().replace("PROFILE_COUNTERS", "true" if instrumented else "false")
    if not instrumented:
        main = main.replace("hashrope::rope::profile_reset();", "")
        main = main.replace("hashrope::rope::profile_read()", "[0u64;12]")
    (project / "src").mkdir(exist_ok=True)
    (project / "src/main.rs").write_text(main)
    (project / "Cargo.toml").write_text('[package]\nname="exact-profile"\nversion="0.0.0"\nedition="2021"\n[dependencies]\nhashrope={path="hashrope"}\n')
    subprocess.run(["cargo", "build", "--release", "--manifest-path", str(project / "Cargo.toml")], check=True, capture_output=True)
    output = subprocess.check_output([str(project / "target/release/exact-profile.exe")], text=True)
    rows = []
    for line in output.splitlines():
        row = json.loads(line)
        row["counters"] = dict(zip(NAMES, row["counters"]))
        rows.append(row)
    return {"variant": variant, "core_rope_sha256": hashlib.sha256((core / "src/rope.rs" if core.name == "rust" else core / "packages/rust/src/rope.rs").read_bytes()).hexdigest(),
            "instrumented_source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "harness_sha256": hashlib.sha256(main.encode()).hexdigest(), "rows": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--optimized-core", type=Path)
    parser.add_argument("--variants", default="baseline,counters,no_memo,no_leaf_bytes,no_memo_no_leaf_bytes")
    parser.add_argument("--output", type=Path, default=HERE / "rust_profile.json")
    args = parser.parse_args()
    result = {"method": __doc__, "variants": [build(v,args.optimized_core) for v in args.variants.split(",")]}
    args.output.write_text(json.dumps(result, indent=2))
    print(args.output)
