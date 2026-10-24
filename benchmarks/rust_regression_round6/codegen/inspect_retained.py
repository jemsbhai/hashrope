"""Disassemble retained round-five rlibs; no builds or benchmark execution.

Normalization keeps instruction offsets, branch targets and registers. It only
removes Rust/LLVM symbol hash suffixes and anonymous data symbol identifiers.
This establishes instruction equality, not data-target or linked-layout equality.
"""
import argparse
import difflib
import hashlib
import json
from pathlib import Path
import re
import subprocess

FUNCTIONS = ["Arena6concat", "Arena4join", "Arena13make_internal",
             "ExactFrontier6expand", "rope9exact_lcp",
             "FastComparison13lcp_recursive", "Arena11split_inner"]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def extract(text, name):
    sections = re.findall(r"(?m)^([^\s].*):\n((?:  [^\n]+\n)+)", text)
    found = [body for symbol, body in sections
             if symbol.startswith("_ZN8hashrope") and name in symbol]
    assert len(found) == 1, (name, len(found))
    body = "\n".join(line.strip() for line in found[0].splitlines()) + "\n"
    body = re.sub(r"17h[0-9a-f]{16}", "17hHASH", body)
    body = re.sub(r"\.llvm\.\d+", ".llvm.HASH", body)
    body = re.sub(r"anon\.[0-9a-f]{32}\.\d+", "anon.HASH.N", body)
    return body


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--dumpbin", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    build_path = args.root / "benchmarks/rust_lcp_round5/build.json"
    build = json.loads(build_path.read_text())
    manifest = {"script_sha256": sha(Path(__file__).read_bytes()),
                "build_manifest_sha256": sha(build_path.read_bytes()),
                "dumpbin": str(args.dumpbin),
                "dumpbin_sha256": sha(args.dumpbin.read_bytes()),
                "inputs": {}, "comparisons": {},
                "limits": __doc__}
    all_functions = {}
    for cohort in ("micro", "churn"):
        for arm in ("baseline", "candidate"):
            record = build[arm][cohort]
            binary = Path(record["binary"])
            binary_sha = sha(binary.read_bytes())
            assert binary_sha == record["binary_sha256"], (arm, cohort)
            rlibs = list((binary.parent / "deps").glob("libhashrope-*.rlib"))
            assert len(rlibs) == 1, rlibs
            rlib = rlibs[0]
            command = [str(args.dumpbin), "/DISASM:NOBYTES", str(rlib)]
            proc = subprocess.run(command, capture_output=True, check=True)
            text = proc.stdout.decode("utf-8", errors="replace").replace("\r\n", "\n")
            key = f"{cohort}.{arm}"
            (args.output / f"{key}.dump.txt").write_text(text, encoding="utf-8")
            manifest["inputs"][key] = {"binary": str(binary), "binary_sha256": binary_sha,
                "rlib": str(rlib), "rlib_sha256": sha(rlib.read_bytes()), "command": command,
                "dump_sha256": sha(text.encode()), "stderr": proc.stderr.decode(errors="replace")}
            all_functions[key] = {}
            for name in FUNCTIONS:
                body = extract(text, name)
                all_functions[key][name] = body
                (args.output / f"{cohort}.{name}.{arm}.asm.txt").write_text(body, encoding="utf-8")
        manifest["comparisons"][cohort] = {}
        for name in FUNCTIONS:
            old = all_functions[f"{cohort}.baseline"][name]
            new = all_functions[f"{cohort}.candidate"][name]
            manifest["comparisons"][cohort][name] = {
                "same_normalized_instructions": old == new,
                "baseline_lines": len(old.splitlines()), "candidate_lines": len(new.splitlines()),
                "baseline_sha256": sha(old.encode()), "candidate_sha256": sha(new.encode())}
            delta = difflib.unified_diff(old.splitlines(True), new.splitlines(True),
                                        fromfile="baseline", tofile="candidate")
            (args.output / f"{cohort}.{name}.diff.txt").write_text("".join(delta), encoding="utf-8")
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest["comparisons"], indent=2))


if __name__ == "__main__":
    main()
