"""Inspect actual final-comparison churn rlibs without rebuilding/executing."""
import argparse
import difflib
import hashlib
import json
from pathlib import Path
import subprocess

from inspect_retained import extract


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dumpbin", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    build_path = args.root / "benchmarks/rust_regression_round6/final_comparison/build.json"
    build = json.loads(build_path.read_text())
    record = {"script_sha256": sha(Path(__file__).read_bytes()),
              "extractor_sha256": sha(Path(__file__).with_name("inspect_retained.py").read_bytes()),
              "build_manifest_sha256": sha(build_path.read_bytes()),
              "dumpbin_sha256": sha(args.dumpbin.read_bytes()), "inputs": {}, "comparison": {},
              "limits": "COFF instruction comparison, not linked code/data placement. Rust/LLVM hashes and anonymous data symbols normalized; branch offsets/registers retained. No timing inference."}
    functions = ["rope9exact_lcp", "FastComparison13lcp_recursive", "ExactFrontier6expand"]
    bodies = {}
    for arm in ("baseline", "candidate"):
        binary = Path(build[arm]["churn"]["binary"])
        assert sha(binary.read_bytes()) == build[arm]["churn"]["binary_sha256"]
        rlibs = list((binary.parent / "deps").glob("libhashrope-*.rlib"))
        assert len(rlibs) == 1, rlibs
        command = [str(args.dumpbin), "/DISASM:NOBYTES", str(rlibs[0])]
        result = subprocess.run(command, capture_output=True, check=True)
        dump = result.stdout.decode().replace("\r\n", "\n")
        (args.output / f"{arm}.dump.txt").write_text(dump)
        record["inputs"][arm] = {"binary": str(binary), "binary_sha256": sha(binary.read_bytes()),
            "rlib": str(rlibs[0]), "rlib_sha256": sha(rlibs[0].read_bytes()), "command": command,
            "dump_sha256": sha(dump.encode())}
        bodies[arm] = {}
        for function in functions:
            body = extract(dump, function)
            bodies[arm][function] = body
            (args.output / f"{arm}.{function}.asm.txt").write_text(body)
    for function in functions:
        old, new = (bodies[arm][function] for arm in ("baseline", "candidate"))
        record["comparison"][function] = {"same_normalized_instructions": old == new,
            "baseline_lines": len(old.splitlines()), "candidate_lines": len(new.splitlines()),
            "baseline_sha256": sha(old.encode()), "candidate_sha256": sha(new.encode())}
        (args.output / f"{function}.diff.txt").write_text("".join(difflib.unified_diff(
            old.splitlines(True), new.splitlines(True), fromfile="baseline", tofile="candidate")))
    (args.output / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record["comparison"], indent=2))


if __name__ == "__main__":
    main()
