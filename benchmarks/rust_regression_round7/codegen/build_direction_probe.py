"""Build frozen-current and direction-only join variants; never time them."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def direction_variant(data):
    newline = b"\r\n" if b"\r\n" in data else b"\n"
    changes = [
        (b"self.join_unbalanced(left, right)", b"self.join_unbalanced(left, right, wl > wr)"),
        (b"self.join_unbalanced(l, r)", b"self.join_unbalanced(l, r, wl > wr)"),
        (newline.join([b"    // Pass only node IDs so join can tail-call this branch on the Windows x64",
                       b"    // ABI. Passing both weights would add a stack argument and an extra call."]),
         newline.join([b"    // Pass the known direction in the fourth Windows x64 register argument.",
                       b"    // This retains tail calls and avoids looking up either weight again."])),
        (newline.join([b"    fn join_unbalanced(&mut self, left: NodeId, right: NodeId) -> NodeId {",
                       b"        let wl = self.node_weight(left);", b"        let wr = self.node_weight(right);",
                       b"        if wl > wr {"]),
         newline.join([b"    fn join_unbalanced(&mut self, left: NodeId, right: NodeId, left_heavy: bool) -> NodeId {",
                       b"        if left_heavy {"]))]
    for old, new in changes:
        assert data.count(old) == 1, old
        data = data.replace(old, new)
    return data


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--dumpbin", type=Path, required=True)
    args = p.parse_args()
    root = args.root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    study = root / "benchmarks/rust_regression_round7"
    core_source = study / ".cache/baseline/hashrope/packages/rust"
    workload = root / "benchmarks/rust_regression_round6/final_comparison/workloads/edit_compare.rs"
    record = {"script_sha256": sha(Path(__file__)), "baseline_manifest_sha256": sha(study / "baseline_manifest.json"),
              "workload_sha256": sha(workload), "dumpbin_sha256": sha(args.dumpbin), "arms": {}}
    for arm in ("current", "direction"):
        dest = output / arm
        core = dest / "hashrope"
        shutil.copytree(core_source, core, ignore=shutil.ignore_patterns("target"))
        rope = core / "src/rope.rs"
        before = sha(rope)
        if arm == "direction":
            rope.write_bytes(direction_variant(rope.read_bytes()))
        (dest / "src").mkdir()
        shutil.copy2(workload, dest / "src/main.rs")
        cargo = dest / "Cargo.toml"
        cargo.write_text('[package]\nname="exact-edit-comparison"\nversion="0.0.0"\nedition="2021"\n'
                         '[dependencies]\nhashrope={path="hashrope"}\nropey="=1.6.1"\n')
        lock = root / "benchmarks/rust_regression_round6/final_comparison/.cache/build/candidate_churn/Cargo.lock"
        shutil.copy2(lock, dest / "Cargo.lock")
        command = ["cargo", "build", "--release", "--offline", "--locked", "--manifest-path", str(cargo)]
        result = subprocess.run(command, capture_output=True, text=True)
        (dest / "build.log").write_text(result.stdout + result.stderr)
        assert result.returncode == 0, result.stderr
        binary = dest / "target/release/exact-edit-comparison.exe"
        rlibs = list((binary.parent / "deps").glob("libhashrope-*.rlib"))
        assert len(rlibs) == 1
        dump = subprocess.run([str(args.dumpbin), "/DISASM:NOBYTES", str(rlibs[0])], capture_output=True, check=True)
        (dest / "core.dump.txt").write_bytes(dump.stdout)
        record["arms"][arm] = {"rope_before_sha256": before, "rope_sha256": sha(rope),
            "cargo_manifest_sha256": sha(cargo), "cargo_lock_sha256": sha(dest / "Cargo.lock"),
            "command": command, "binary_sha256": sha(binary), "rlib_sha256": sha(rlibs[0]),
            "dump_sha256": sha(dest / "core.dump.txt")}
    (output / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
