"""Run unchanged Python consumer suites and biology Rust before final timing."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PAPER = ROOT.parent / "hashropeictaipaper"
BIO = ROOT / "usecases/computational-biology/code"
PYTHON = ROOT / "packages/python"


def core_hashes():
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for directory, pattern in [(ROOT/"packages/rust/src", "*.rs"), (PYTHON/"hashrope", "*.py")]
            for p in sorted(directory.glob(pattern))}


def run_check(name, cwd, command):
    before = core_hashes()
    env = os.environ.copy()
    env.update(PYTHONPATH=os.pathsep.join([str(PYTHON), str(BIO/"src")]),
               HF_HUB_OFFLINE="1", TRANSFORMERS_OFFLINE="1")
    started = time.time()
    result = subprocess.run(command, cwd=cwd, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    (HERE/(name+".log")).write_text(result.stdout, encoding="utf-8")
    after = core_hashes()
    record = {"command": command, "cwd": str(cwd), "seconds": time.time()-started,
              "returncode": result.returncode, "core_before": before, "core_after": after,
              "core_unchanged_during_check": before == after, "log": name+".log"}
    print(name, result.returncode, result.stdout.strip().splitlines()[-1], flush=True)
    return name, record


def main():
    expected = str((PYTHON/"hashrope/__init__.py").resolve())
    entry = ("from pathlib import Path;import hashrope;"
             f"assert Path(hashrope.__file__).resolve()==Path({expected!r}),hashrope.__file__;"
             "import pytest;raise SystemExit(pytest.main(['-q','tests']))")
    py = [sys.executable, "-c", entry]
    commands = [("python_core", PYTHON, py), ("paper", PAPER, py), ("bio_python", BIO, py),
                ("bio_rust", BIO/"rust", ["cargo", "test", "--release"])]
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = dict(pool.map(lambda args: run_check(*args), commands))
    (HERE/"python_and_bio.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    assert all(row["returncode"] == 0 and row["core_unchanged_during_check"] for row in results.values()), results


if __name__ == "__main__":
    main()
