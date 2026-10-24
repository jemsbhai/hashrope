"""Build two constructor drafts to inspect ABI boundaries; never time them."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--dumpbin", type=Path, required=True)
    args = p.parse_args()
    root = args.root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    codegen = root / "benchmarks/rust_regression_round6/codegen"
    snapshot = codegen / "constructor_v2/rope.rs.txt"
    baseline = root / "benchmarks/rust_regression_round6/.cache/baseline/hashrope/packages/rust"
    workload = root / "benchmarks/rust_lcp_round5/workloads/edit_compare.rs"
    manifest = {"script_sha256": sha(Path(__file__)), "workload_sha256": sha(workload),
                "dumpbin_sha256": sha(args.dumpbin), "rope_snapshot_sha256": sha(snapshot),
                "arms": {}}
    for arm in ("five_args", "three_args"):
        dest = output / arm
        core = dest / "hashrope"
        shutil.copytree(baseline, core, ignore=shutil.ignore_patterns("target"))
        rope = core / "src/rope.rs"
        text = snapshot.read_bytes()
        if arm == "three_args":
            newline = b"\r\n" if b"\r\n" in text else b"\n"
            substitutions = [
                (b"self.join_unbalanced(left, right, wl, wr)", b"self.join_unbalanced(left, right)"),
                (b"self.join_unbalanced(l, r, wl, wr)", b"self.join_unbalanced(l, r)"),
                (b"fn join_unbalanced(&mut self, left: NodeId, right: NodeId, wl: u64, wr: u64) -> NodeId {",
                 b"fn join_unbalanced(&mut self, left: NodeId, right: NodeId) -> NodeId {" + newline +
                 b"        let wl = self.node_weight(left);" + newline +
                 b"        let wr = self.node_weight(right);")]
            for old, new in substitutions:
                assert text.count(old) == 1
                text = text.replace(old, new)
        rope.write_bytes(text)
        (dest / "src").mkdir()
        shutil.copy2(workload, dest / "src/main.rs")
        cargo = dest / "Cargo.toml"
        cargo.write_text('[package]\nname="exact-edit-comparison"\nversion="0.0.0"\nedition="2021"\n'
                         '[dependencies]\nhashrope={path="hashrope"}\nropey="=1.6.1"\n')
        lock = root / "benchmarks/rust_lcp_round5/.cache/build/candidate_churn/Cargo.lock"
        shutil.copy2(lock, dest / "Cargo.lock")
        command = ["cargo", "build", "--release", "--offline", "--locked", "--manifest-path", str(cargo)]
        result = subprocess.run(command, capture_output=True, text=True)
        (dest / "build.log").write_text(result.stdout + result.stderr)
        assert result.returncode == 0, result.stderr
        binary = dest / "target/release/exact-edit-comparison.exe"
        rlibs = list((binary.parent / "deps").glob("libhashrope-*.rlib"))
        assert len(rlibs) == 1
        dump = subprocess.run([str(args.dumpbin), "/DISASM:NOBYTES", str(rlibs[0])],
                              capture_output=True, check=True)
        (dest / "core.dump.txt").write_bytes(dump.stdout)
        manifest["arms"][arm] = {"rope_sha256": sha(rope),
            "cargo_manifest_sha256": sha(cargo), "cargo_lock_sha256": sha(dest / "Cargo.lock"),
            "command": command, "binary_sha256": sha(binary), "rlib_sha256": sha(rlibs[0]),
            "dump_sha256": sha(dest / "core.dump.txt")}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
