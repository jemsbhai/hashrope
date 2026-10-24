"""Inspect separately emitted assembly; never rebuilds or replaces measured binaries."""
import difflib
import hashlib
import json
from pathlib import Path
import re
import zipfile

HERE = Path(__file__).resolve().parent
DIRECTORY = HERE / ".cache/assembly/release/deps"
PARTS = ["14FastComparison13lcp_recursive", "5Arena10lcp_across", "5Arena12equal_across",
         "9exact_lcp", "15ExactComparison3lcp", "15ExactComparison9lcp_inner", "15ExactComparison9normalize"]


def normalized_function(source, part):
    match = re.search(r"^_ZN26hashrope_diagnostic_slot_[ab]4rope" + part + r"17h[0-9a-f]+E:$", source, re.M)
    assert match, part
    end = source.index(".seh_endproc", match.start())
    body = source[match.start():end]
    body = re.sub(r"hashrope_diagnostic_slot_[ab]", "hashrope_diagnostic_SLOT", body)
    body = re.sub(r"17h[0-9a-f]{16}E", "17hSYMBOL_HASH_E", body)
    body = re.sub(r"\.LBB\d+_", ".LBBFUNCTION_", body)
    body = re.sub(r"\.Lfunc_begin\d+", ".Lfunc_beginFUNCTION", body)
    # Anonymous globals remain distinct within a function. Renumber by first
    # appearance; this normalizes relocations, not their target contents.
    anonymous = {}
    def rename(match):
        return anonymous.setdefault(match.group(), f"ANONYMOUS_RELOCATION_{len(anonymous)}")
    body = re.sub(r"__unnamed_\d+", rename, body)
    body = re.sub(r"anon\.[0-9a-f]{32}\.\d+", rename, body)
    body = re.sub(r"\.Ltmp\d+", rename, body)
    return body


def main():
    paths = {slot: next(DIRECTORY.glob(f"hashrope_diagnostic_slot_{slot}-*.s")) for slot in ["a", "b"]}
    source = {slot: path.read_text() for slot, path in paths.items()}
    output = {"scope": "Separately emitted release library assembly, old_new slot build. Not a final linked-binary identity proof.",
              "normalization": "Crate-slot names, symbol hashes, local labels, and anonymous-relocation labels; relocation target contents and final placement not compared.",
              "assembly_sha256": {slot: hashlib.sha256(path.read_bytes()).hexdigest() for slot, path in paths.items()},
              "functions": []}
    diffs = []
    for part in PARTS:
        bodies = {slot: normalized_function(text, part) for slot, text in source.items()}
        output["functions"].append({"name": part, "normalized_equal": bodies["a"] == bodies["b"],
            "line_counts": {slot: len(body.splitlines()) for slot, body in bodies.items()},
            "normalized_sha256": {slot: hashlib.sha256(body.encode()).hexdigest() for slot, body in bodies.items()}})
        if bodies["a"] != bodies["b"]:
            diffs.extend(difflib.unified_diff(bodies["a"].splitlines(True), bodies["b"].splitlines(True), fromfile=f"round2/{part}", tofile=f"round3/{part}"))
    (HERE / "assembly_audit.json").write_text(json.dumps(output, indent=2), encoding="utf-8")
    (HERE / "assembly_differences.txt").write_text("".join(diffs), encoding="utf-8")
    with zipfile.ZipFile(HERE / "comparison_assembly.zip", "x", zipfile.ZIP_DEFLATED) as archive:
        for slot, path in paths.items():
            archive.write(path, f"slot_{slot}.s")
        archive.write(Path(__file__), Path(__file__).name)
        archive.write(HERE / "assembly_audit.json", "assembly_audit.json")
        archive.write(HERE / "assembly_differences.txt", "assembly_differences.txt")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
