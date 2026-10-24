"""Read retained round-six rlibs/PEs: no compilation or benchmark execution.

COFF instruction comparison removes symbol hashes and anonymous data symbols;
branch offsets and registers remain. PE placement uses an exact 64-byte prologue
signature, requires a unique .text match and an x64 runtime-function table entry.
It records linked placement, not ASLR-resolved runtime addresses or causation.
"""
import argparse
import difflib
import hashlib
import json
from pathlib import Path
import re
import struct
import subprocess

FUNCTIONS = ["rope9exact_lcp", "FastComparison13lcp_recursive", "ExactFrontier6expand",
             "Arena6concat", "Arena4join", "Arena15join_unbalanced", "Arena11split_inner"]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def sections(text):
    result = {}
    for match in re.finditer(r"(?m)^(_ZN8hashrope[^\n]*):\n", text):
        start = match.end()
        lines = []
        for line in text[start:].splitlines():
            # Unwind funclets can immediately follow without a blank line.
            # Stop at their unindented labels instead of counting them in the
            # primary function's instruction body.
            if not line or not line.startswith(" "):
                break
            lines.append(line)
        result[match.group(1)] = "\n".join(lines)
    return result


def select(text, name):
    found = [body for symbol, body in sections(text).items() if name in symbol]
    assert len(found) <= 1, (name, len(found))
    return found[0] if found else None


def normalized(body):
    body = "\n".join(line.strip() for line in body.splitlines()) + "\n"
    body = re.sub(r"17h[0-9a-f]{16}", "17hHASH", body)
    body = re.sub(r"\.llvm\.\d+", ".llvm.HASH", body)
    return re.sub(r"anon\.[0-9a-f]{32}\.\d+", "anon.HASH.N", body)


def prologue(body):
    data = bytearray()
    for line in body.splitlines():
        instruction = re.match(r"  ([0-9A-F]{16}): ((?:[0-9A-F]{2} )*[0-9A-F]{2})", line)
        if instruction:
            assert int(instruction[1], 16) == len(data), line
            # A relative call or named relocation in a signature would not
            # preserve its object bytes at final link. These selected prologues
            # should reach 64 bytes before any such instruction.
            if "call " in line or "_ZN" in line or "anon." in line:
                break
            data.extend(bytes.fromhex(instruction[2]))
        else:
            continuation = re.fullmatch(r"\s+((?:[0-9A-F]{2}\s*)+)", line)
            assert continuation, line
            data.extend(bytes.fromhex(continuation[1]))
        if len(data) >= 64:
            return bytes(data[:64])
    return bytes(data)


def pe_map(data):
    pe = struct.unpack_from("<I", data, 0x3c)[0]
    assert data[pe:pe + 4] == b"PE\0\0"
    machine, count = struct.unpack_from("<HH", data, pe + 4)
    assert machine == 0x8664
    opt_size = struct.unpack_from("<H", data, pe + 20)[0]
    opt = pe + 24
    assert struct.unpack_from("<H", data, opt)[0] == 0x20b
    image_base = struct.unpack_from("<Q", data, opt + 24)[0]
    entries = []
    for i in range(count):
        start = opt + opt_size + i * 40
        name = data[start:start + 8].rstrip(b"\0").decode()
        size, rva, raw_size, raw = struct.unpack_from("<IIII", data, start + 8)
        entries.append({"name": name, "size": size, "rva": rva, "raw_size": raw_size, "raw": raw})
    def offset(rva):
        for entry in entries:
            if entry["rva"] <= rva < entry["rva"] + max(entry["size"], entry["raw_size"]):
                return entry["raw"] + rva - entry["rva"]
        raise AssertionError(rva)
    exception_rva, exception_size = struct.unpack_from("<II", data, opt + 112 + 3 * 8)
    runtime = {}
    if exception_rva:
        start = offset(exception_rva)
        assert exception_size % 12 == 0
        for i in range(exception_size // 12):
            begin, end, unwind = struct.unpack_from("<III", data, start + i * 12)
            runtime[begin] = {"end_rva": end, "unwind_rva": unwind}
    return image_base, next(entry for entry in entries if entry["name"] == ".text"), runtime


def placement(data, signature):
    base, section, runtime = pe_map(data)
    code = data[section["raw"]:section["raw"] + section["raw_size"]]
    hits = []
    position = code.find(signature)
    while position >= 0:
        hits.append(section["rva"] + position)
        position = code.find(signature, position + 1)
    assert len(signature) >= 24, len(signature)
    assert len(hits) == 1, hits
    rva = hits[0]
    assert rva in runtime, (hex(rva), "not a runtime-function entry")
    return {"signature_bytes": len(signature), "signature_hex": signature.hex(),
            "rva": rva, "rva_hex": hex(rva), "preferred_address_hex": hex(base + rva),
            "mod64": rva % 64, "mod4096": rva % 4096,
            "runtime_function_size": runtime[rva]["end_rva"] - rva,
            "unwind_rva": runtime[rva]["unwind_rva"]}


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
              "build_manifest_sha256": sha(build_path.read_bytes()),
              "dumpbin_sha256": sha(args.dumpbin.read_bytes()), "inputs": {},
              "comparison": {}, "limits": __doc__}
    bodies = {}
    for cohort in ("micro", "churn"):
        for arm in ("baseline", "candidate"):
            key = f"{cohort}.{arm}"
            binary = Path(build[arm][cohort]["binary"])
            exe = binary.read_bytes()
            assert sha(exe) == build[arm][cohort]["binary_sha256"]
            rlibs = list((binary.parent / "deps").glob("libhashrope-*.rlib"))
            assert len(rlibs) == 1, rlibs
            dumps = {}
            for mode, suffix in (("/DISASM:NOBYTES", "normalized_source"), ("/DISASM", "bytes")):
                command = [str(args.dumpbin), mode, str(rlibs[0])]
                result = subprocess.run(command, capture_output=True, check=True)
                dumps[suffix] = result.stdout.decode().replace("\r\n", "\n")
                (args.output / f"{key}.{suffix}.dump.txt").write_text(dumps[suffix])
            record["inputs"][key] = {"binary": str(binary), "binary_sha256": sha(exe),
                "rlib": str(rlibs[0]), "rlib_sha256": sha(rlibs[0].read_bytes()),
                "dump_sha256": {k: sha(v.encode()) for k, v in dumps.items()}, "functions": {}}
            bodies[key] = {}
            for function in FUNCTIONS:
                source = select(dumps["normalized_source"], function)
                if source is None:
                    record["inputs"][key]["functions"][function] = None
                    continue
                body = normalized(source)
                bodies[key][function] = body
                (args.output / f"{key}.{function}.asm.txt").write_text(body)
                signature = prologue(select(dumps["bytes"], function))
                try:
                    location = placement(exe, signature)
                except AssertionError as error:
                    # Short leaf helpers without .pdata entries are not forced
                    # into a possibly incorrect linked-symbol identification.
                    location = {"unresolved": str(error), "signature_hex": signature.hex()}
                record["inputs"][key]["functions"][function] = {
                    "instruction_lines": len(body.splitlines()), "normalized_sha256": sha(body.encode()),
                    "placement": location}
        record["comparison"][cohort] = {}
        for function in FUNCTIONS:
            old = bodies[f"{cohort}.baseline"].get(function)
            new = bodies[f"{cohort}.candidate"].get(function)
            record["comparison"][cohort][function] = None if old is None or new is None else old == new
            if old is not None and new is not None:
                (args.output / f"{cohort}.{function}.diff.txt").write_text("".join(difflib.unified_diff(
                    old.splitlines(True), new.splitlines(True), fromfile="baseline", tofile="candidate")))
    (args.output / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record["comparison"], indent=2))


if __name__ == "__main__":
    main()
