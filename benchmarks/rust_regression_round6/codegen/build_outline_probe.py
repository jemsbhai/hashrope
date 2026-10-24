"""Build (never time) unchanged round-five and frontier-outline-only cores."""
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
    baseline = root / "benchmarks/rust_regression_round6/.cache/baseline/hashrope/packages/rust"
    workload = root / "benchmarks/rust_lcp_round5/workloads/edit_compare.rs"
    manifest = {"script_sha256": sha(Path(__file__)), "workload_sha256": sha(workload),
                "dumpbin_sha256": sha(args.dumpbin), "arms": {}}
    for arm in ("unchanged", "outlined"):
        dest = output / arm
        core = dest / "hashrope"
        shutil.copytree(baseline, core, ignore=shutil.ignore_patterns("target"))
        rope = core / "src/rope.rs"
        before = sha(rope)
        if arm == "outlined":
            text = rope.read_bytes()
            needle = b"    fn frontier_lcp(&mut self, a: RangeStart, b: RangeStart, length: u64) -> Result<u64, u64> {"
            assert text.count(needle) == 1
            # Preserve every other source byte, including existing CRLFs.
            newline = b"\r\n" if b"\r\n" in text else b"\n"
            rope.write_bytes(text.replace(needle, b"    #[inline(never)]" + newline + needle))
        (dest / "src").mkdir()
        shutil.copy2(workload, dest / "src/main.rs")
        cargo = dest / "Cargo.toml"
        cargo.write_text('[package]\nname="exact-edit-comparison"\nversion="0.0.0"\nedition="2021"\n'
                         '[dependencies]\nhashrope={path="hashrope"}\nropey="=1.6.1"\n')
        lock_source = root / "benchmarks/rust_lcp_round5/.cache/build/candidate_churn/Cargo.lock"
        shutil.copy2(lock_source, dest / "Cargo.lock")
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
        manifest["arms"][arm] = {"rope_before_sha256": before, "rope_sha256": sha(rope),
            "cargo_manifest_sha256": sha(cargo), "cargo_lock_sha256": sha(dest / "Cargo.lock"),
            "command": command, "binary_sha256": sha(binary), "rlib_sha256": sha(rlibs[0]),
            "dump_sha256": sha(dest / "core.dump.txt")}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
