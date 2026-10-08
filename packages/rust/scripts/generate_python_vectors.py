"""Regenerate valid-input Rust conformance fixtures with published hashrope 0.2.2.

Run with the extracted, checksum-verified PyPI wheel on PYTHONPATH. This script
does not install packages or alter Python implementation files. Integer fields
are decimal strings and bytes are hex so every language can read them exactly.
"""

import argparse
import hashlib
import json
import zipfile
from pathlib import Path

import hashrope as hr


def decimal(value):
    return str(value)


def snapshot(node):
    return {"length": decimal(hr.rope_len(node)), "hash": decimal(hr.rope_hash(node))}


def profile(prime, base):
    h = hr.PolynomialHash(prime, base)
    samples = [b"", b"\0", b"\0\0", b"\xff", b"hello world", bytes(range(256)),
               "Aé🙂e\u0301\0漢字".encode("utf-8"), b"\xff\xfe\x80\xc0\xaf\x00",
               bytes((i * 73 + 19) % 256 for i in range(513)),
               bytes((i * 137 + 31) % 256 for i in range(4097))]
    rows = []
    for data in samples:
        node = hr.rope_from_bytes(data, h)
        assert hr.rope_hash(node) == h.hash(data)
        cuts = sorted({0, len(data), len(data) // 2, min(1, len(data)), min(3, len(data)), min(512, len(data))})
        splits = []
        for cut in cuts:
            left, right = hr.rope_split(node, cut, h)
            assert hr.rope_to_bytes(left) == data[:cut]
            assert hr.rope_to_bytes(right) == data[cut:]
            splits.append({"at": decimal(cut), "left": snapshot(left), "right": snapshot(right)})
        ranges = sorted({(0, 0), (0, len(data)), (len(data), 0),
                         (len(data) // 3, len(data) // 2), (min(1, len(data)), min(4, max(0, len(data) - 1)))})
        substrings = [{"start": decimal(start), "length": decimal(length),
                       "hash": decimal(hr.rope_substr_hash(node, start, length, h))}
                      for start, length in ranges]
        for sub in substrings:
            start, length = int(sub["start"]), int(sub["length"])
            assert int(sub["hash"]) == h.hash(data[start:start + length])
        rows.append({"bytes_hex": data.hex(), **snapshot(node), "splits": splits, "substrings": substrings})

    repeats = []
    pattern = b"a\0\xff\xc3\xa9"
    child = hr.rope_from_bytes(pattern, h)
    for count in [0, 1, 2, 17, 1000000, (1 << 61), (1 << 80) + 17]:
        node = hr.rope_repeat(child, count, h)
        length = hr.rope_len(node)
        ranges = [(0, 0), (0, length)]
        if count:
            ranges += [(1, min(13, length - 1)), (length - 1, 1)]
        if count > 2:
            ranges += [(3, length - 7)]
        queries = [{"start": decimal(start), "length": decimal(size),
                    "hash": decimal(hr.rope_substr_hash(node, start, size, h))}
                   for start, size in ranges]
        repeats.append({"pattern_hex": pattern.hex(), "count": decimal(count),
                        **snapshot(node), "substrings": queries})

    powers = [{"exponent": decimal(n), "value": decimal(h.power(n))}
              for n in [0, 1, 255, 256, 257, 1 << 32, (1 << 32) + 1, (1 << 64) - 1, (1 << 80) + 17]]
    phis = [{"count": decimal(q), "alpha": decimal(a), "value": decimal(hr.phi(q, a, prime))}
            for q in [0, 1, 2, 63, 1 << 32, (1 << 80) + 17]
            for a in [0, 1, base, prime - 1]]
    moduli = [{"input": decimal(a), "value": decimal(hr.mersenne_mod(a, prime))}
              for a in [0, 1, 256, (1 << 128) - 1, (1 << 1200) + 991]]
    overlap = []
    for length in [1, 5, 6, 27, 1000000, (1 << 80) + 17]:
        rem = length % len(pattern)
        overlap.append({"pattern_hex": pattern.hex(), "length": decimal(length),
                        "hash": decimal(h.hash_overlap(h.hash(pattern), len(pattern), length, h.hash(pattern[:rem])))})

    # Persistent edit composition: insert, delete, replace, append. No mutation API
    # is invented: these are the same split/concat primitives exposed in Python.
    data = "original é🙂 bytes".encode("utf-8")
    node = hr.rope_from_bytes(data, h)
    edits = [{"kind": "from", "bytes_hex": data.hex(), **snapshot(node)}]
    for start, removed, inserted in [(0, 0, b"prefix:"), (11, 3, b"\xff\x00"), (4, 7, b""), (2, 0, b"INSERT")]:
        left, rest = hr.rope_split(node, start, h)
        _, right = hr.rope_split(rest, removed, h)
        middle = hr.rope_from_bytes(inserted, h)
        node = hr.rope_concat(hr.rope_concat(left, middle, h), right, h)
        data = data[:start] + inserted + data[start + removed:]
        assert hr.rope_to_bytes(node) == data
        edits.append({"kind": "splice", "start": decimal(start), "removed": decimal(removed),
                      "insert_hex": inserted.hex(), "bytes_hex": data.hex(), **snapshot(node)})

    sw = hr.SlidingWindow(8, 4, prime, base)
    decoded = bytearray()
    stream = []
    operations = [("bytes", b"abcdefgh"), ("copy", (3, 3)), ("copy", (4, 13)),
                  ("bytes", bytes(range(32))), ("bytes", b""), ("copy", (8, 9)),
                  ("copy", (1, 17)), ("bytes", "é🙂".encode("utf-8"))]
    for kind, value in operations:
        if kind == "bytes":
            sw.append_bytes(value)
            decoded.extend(value)
            op = {"kind": kind, "bytes_hex": value.hex()}
        else:
            offset, length = value
            sw.append_copy(offset, length)
            for _ in range(length):
                decoded.append(decoded[-offset])
            op = {"kind": kind, "offset": decimal(offset), "length": decimal(length)}
        assert sw.current_hash() == h.hash(bytes(decoded))
        stream.append({**op, "hash": decimal(sw.current_hash()), "position": decimal(sw.pos),
                       "window_length": decimal(sw.window_len)})
    return {"prime": decimal(prime), "base": decimal(base), "bytes": rows,
            "repeats": repeats, "powers": powers, "phis": phis, "moduli": moduli,
            "overlaps": overlap, "edits": edits, "stream": stream}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel", type=Path, required=True, help="Downloaded hashrope-0.2.2 wheel for provenance validation")
    parser.add_argument("--check", action="store_true", help="Verify fixture is reproducible without writing")
    args = parser.parse_args()
    assert hr.__version__ == "0.2.2", hr.__version__
    digest = hashlib.sha256(args.wheel.read_bytes()).hexdigest()
    assert digest == "710112e3aa14a991918a15331c1f418ce9d1c9d93c34a9da09e22be33c68a2d2", digest
    with zipfile.ZipFile(args.wheel) as wheel:
        for name in ["__init__.py", "polynomial_hash.py", "rope.py", "sliding.py"]:
            loaded = Path(hr.__file__).parent / name
            assert loaded.read_bytes() == wheel.read(f"hashrope/{name}"), f"Import did not use the published wheel: {loaded}"
    result = {
        "schema": 1,
        "reference": {"package": "hashrope", "version": hr.__version__, "wheel_sha256": digest,
                      "source": "https://pypi.org/project/hashrope/0.2.2/",
                      "contract": "Valid byte-oriented operations; decimal integer strings and hex bytes. This is a test fixture, not a rope serialization format."},
        "profiles": [profile(p, b) for p, b in [(3, 2), (7, 3), (31, 3),
                     (hr.MERSENNE_61, 131), (hr.MERSENNE_61, 257),
                     (hr.MERSENNE_127, 131), (hr.MERSENNE_127, hr.MERSENNE_127 - 2), ((1 << 521) - 1, 131)]]}
    output = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "python-0.2.2.json"
    text = json.dumps(result, indent=2, ensure_ascii=True) + "\n"
    if args.check:
        assert output.read_text(encoding="utf-8") == text, "Fixture differs from published Python reference"
        print("Published Python vectors match checked-in fixture")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text, encoding="utf-8", newline="\n")
        print(f"Wrote {len(result['profiles'])} profiles to {output}")


if __name__ == "__main__":
    main()
