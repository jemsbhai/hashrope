"""Alternate already-built diagnostic cores; keep confirmatory matrix separate."""
import argparse
import json
from pathlib import Path
import subprocess

from edit_profile import HERE, parse, sha, summarize


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--labels", nargs=2, default=["baseline_final", "selected"])
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument("--sizes", type=int, nargs="+", default=[1000, 10000, 100000, 1000000])
    args = parser.parse_args()
    output = args.output_dir.resolve()
    results = {}
    template = json.loads((HERE / "edit_profile.baseline.json").read_text())
    for label in args.labels:
        destination = output / f"edit_profile.{label}.json"
        if destination.exists(): raise SystemExit(f"Preserve existing result: {destination}")
        metadata = json.loads((output / f"edit_profile.{label}.build.json").read_text())
        assert all(sha(Path(path)) == digest for path, digest in metadata["source_sha256"].items())
        assert all(sha(Path(metadata["binaries"][mode])) == digest for mode, digest in metadata["binary_sha256"].items())
        result = {key: value for key, value in template.items() if key not in ["runs", "summary", "build", "label"]}
        result.update(label=label, build=metadata, sizes=args.sizes, runs=[],
                      method="Same diagnostic decomposition, with alternating core order by seed; 1K is an exploratory extra cell.",
                      pair_driver_sha256=sha(Path(__file__)),
                      source_archive_sha256=sha(output / f"edit_profile.{label}.sources.zip"))
        results[label] = result
    for mode in ["timing", "counters"]:
        for seed in [1, 2, 3]:
            order = args.labels if seed % 2 else list(reversed(args.labels))
            for label in order:
                print(f"Paired diagnostic {mode}: {label}, seed={seed}", flush=True)
                binary = results[label]["build"]["binaries"][mode]
                path = output / f"edit_profile.{label}.{mode}.s{seed}.tsv"
                if path.exists(): raise SystemExit(f"Preserve existing raw diagnostic: {path}")
                raw = subprocess.check_output([binary, str(seed), "3" if mode == "timing" else "1",
                                               ",".join(map(str, args.sizes))], text=True)
                path.write_text(raw, encoding="utf-8")
                results[label]["runs"].append({"mode": mode, "seed": seed, "core_order": order,
                                               "raw_sha256": sha(path), "rows": parse(raw)})
    for label, result in results.items():
        result["summary"] = summarize(result)
        (output / f"edit_profile.{label}.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print("Paired decomposition complete.")


if __name__ == "__main__":
    main()
