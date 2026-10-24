"""Ablate traversal and byte-localization changes on immutable diagnostic copies.

All variants remain exact. Fixture construction is excluded. Timings are
diagnostics, not the final unchanged comparison matrix. --build-only prepares
all binaries without occupying a timing window. --run executes existing builds.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
BASELINE = HERE / ".cache/baseline/hashrope"
VARIANTS = ["baseline", "children", "words", "children_words"]

OLD_PREFIX = """            let piece = ExactComparison::boundary(self.left, a)
                .min(ExactComparison::boundary(self.right, b)).min(remaining);
            debug_assert!(piece > 0 && piece < remaining);
            match self.lcp(a, b, piece) {"""
NEW_PREFIX = """            let a_boundary = ExactComparison::boundary(self.left, a);
            let b_boundary = ExactComparison::boundary(self.right, b);
            let piece = a_boundary.min(b_boundary).min(remaining);
            debug_assert!(piece > 0 && piece < remaining);
            let prefix_a = match &self.left.nodes[a.node as usize] {
                NodeInner::Internal { left, .. } => RangeStart { node: *left, offset: a.offset },
                _ => a,
            };
            let prefix_b = match &self.right.nodes[b.node as usize] {
                NodeInner::Internal { left, .. } => RangeStart { node: *left, offset: b.offset },
                _ => b,
            };
            match self.lcp(prefix_a, prefix_b, piece) {"""
OLD_ADVANCE = """            a.offset += piece;
            b.offset += piece;
            matched += piece;"""
NEW_ADVANCE = """            a = match &self.left.nodes[a.node as usize] {
                NodeInner::Internal { right, .. } if piece == a_boundary => RangeStart { node: *right, offset: 0 },
                _ => RangeStart { offset: a.offset + piece, ..a },
            };
            b = match &self.right.nodes[b.node as usize] {
                NodeInner::Internal { right, .. } if piece == b_boundary => RangeStart { node: *right, offset: 0 },
                _ => RangeStart { offset: b.offset + piece, ..b },
            };
            matched += piece;"""
WORD_HELPER = """
// The caller has already proved equal lengths and unequal bytes.
fn first_different_byte(left: &[u8], right: &[u8]) -> usize {
    let mut offset = 0;
    while left.len() - offset >= 8 {
        let a = u64::from_le_bytes(left[offset..offset + 8].try_into().unwrap());
        let b = u64::from_le_bytes(right[offset..offset + 8].try_into().unwrap());
        let different = a ^ b;
        if different != 0 { return offset + (different.trailing_zeros() / 8) as usize; }
        offset += 8;
    }
    offset + left[offset..].iter().zip(&right[offset..]).position(|(a,b)| a != b).unwrap()
}
"""


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(variant):
    project = HERE / ".cache" / ("fast_path_" + variant)
    library = project / "hashrope"
    shutil.copytree(BASELINE / "packages/rust/src", library / "src", dirs_exist_ok=True)
    shutil.copyfile(BASELINE / "packages/rust/Cargo.toml", library / "Cargo.toml")
    shutil.copyfile(BASELINE / "packages/rust/README.md", library / "README.md")
    path = library / "src/rope.rs"
    source = path.read_text()
    start = source.index("impl FastComparison<'_> {")
    end = source.index("struct ExactComparison<'a> {", start)
    fast = source[start:end]
    if "children" in variant:
        assert fast.count(OLD_PREFIX) == 1 and fast.count(OLD_ADVANCE) == 1
        fast = fast.replace(OLD_PREFIX, NEW_PREFIX).replace(OLD_ADVANCE, NEW_ADVANCE)
    if "words" in variant:
        old = "x.iter().zip(y).position(|(x, y)| x != y).unwrap()"
        assert fast.count(old) == 1
        fast = fast.replace(old, "first_different_byte(x, y)") + WORD_HELPER + "\n"
    if "frontier" in variant:
        signature = "fn lcp(&mut self, mut a: RangeStart, mut b: RangeStart, length: u64) -> Result<u64, u64> {"
        assert fast.count(signature) == 1
        condition = "length >= 8192" if "cross" not in variant else "!self.same_arena && length >= 8192"
        fast = fast.replace(signature, signature + f"\n        if {condition} {{ return self.frontier_lcp(a, b, length); }}")
        fast += (HERE / "fast_path_frontier.rs.txt").read_text() + "\n"
    source = source[:start] + fast + source[end:]
    path.write_text(source)
    (project / "src").mkdir(exist_ok=True)
    shutil.copyfile(HERE / "fast_path_profile.rs", project / "src/main.rs")
    (project / "Cargo.toml").write_text('[package]\nname="fast-path-profile"\nversion="0.0.0"\nedition="2021"\n[dependencies]\nhashrope={path="hashrope"}\n')
    subprocess.run(["cargo", "build", "--release", "--offline", "--manifest-path", str(project / "Cargo.toml")], capture_output=True, check=True)
    binary = project / "target/release/fast-path-profile.exe"
    metadata = {"variant": variant, "baseline_rope_sha256": sha(BASELINE / "packages/rust/src/rope.rs"),
                "variant_rope_sha256": sha(path), "harness_sha256": sha(HERE / "fast_path_profile.rs"),
                "binary_sha256": sha(binary), "binary": str(binary)}
    binary.with_suffix(".json").write_text(json.dumps(metadata, indent=2))
    return metadata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--output", type=Path, default=HERE / "fast_path_diagnostics.json")
    parser.add_argument("--variants", default=",".join(VARIANTS))
    args = parser.parse_args()
    builds = [json.loads((HERE / ".cache" / ("fast_path_" + variant) / "target/release/fast-path-profile.json").read_text())
              if args.run else build(variant) for variant in args.variants.split(",")]
    if args.build_only:
        print(json.dumps(builds, indent=2))
        return
    if args.output.exists():
        raise SystemExit(f"Refusing to overwrite {args.output}")
    result = {"purpose": __doc__, "order": [], "runs": [], "builds": builds}
    for round_id in range(3):
        order = builds if round_id % 2 == 0 else list(reversed(builds))
        result["order"].append([b["variant"] for b in order])
        for metadata in order:
            binary = Path(metadata["binary"])
            assert sha(binary) == metadata["binary_sha256"]
            rows = [json.loads(line) for line in subprocess.check_output([str(binary)], text=True).splitlines()]
            result["runs"].append({"variant": metadata["variant"], "round": round_id, "rows": rows})
            print(metadata["variant"], round_id, flush=True)
    args.output.write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
