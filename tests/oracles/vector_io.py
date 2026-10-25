#!/usr/bin/env python3
"""Physics tests of VectorAscii2Binary and VectorBinary2Ascii against the documented binary vector format.

    python3 tests/oracles/vector_io.py ascii2binary PROGRAM WORK_DIR
    python3 tests/oracles/vector_io.py binary2ascii PROGRAM WORK_DIR

The format (docs/reference/data-formats.md, RealVector::WriteVector): a 4-byte little-endian int, the
dimension, then the components as little-endian 8-byte IEEE doubles. This script reads and writes that layout
itself with Python's struct module; it shares no code with DiagHam.

ascii2binary: writes a text vector (one number per line, chosen to need every bit: 1/3, 1e-300, -0.1, a
17-digit value, an integer above 2^53), runs the program, and requires the binary file to hold exactly the
doubles Python reads from the same text (correctly rounded, so bit for bit).
binary2ascii: writes the same doubles in the binary layout, runs the program, and requires each printed value
to equal the double to 14 significant digits, upstream's output precision (relative 1e-13).
Exit status 0 on pass, 1 on a failed check.
"""
import struct, subprocess, sys
from pathlib import Path

TEXT = ["0.5", "-1.25", "3e-05", "123456.789", "0.33333333333333331", "1e-300", "-0.1",
        "9007199254740993", "-2.7182818284590451", "6.02214076e+23"]
VALUES = [float(t) for t in TEXT]


def write_vec(path, values):
    path.write_bytes(struct.pack("<i", len(values)) + struct.pack(f"<{len(values)}d", *values))


def read_vec(path):
    data = path.read_bytes()
    (n,) = struct.unpack_from("<i", data, 0)
    if n < 0 or len(data) != 4 + 8 * n:
        raise ValueError(f"{path}: header says {n} components, file has {len(data)} bytes")
    return list(struct.unpack_from(f"<{n}d", data, 4))


def main():
    mode, program, work = sys.argv[1], sys.argv[2], Path(sys.argv[3])
    work.mkdir(parents=True, exist_ok=True)
    if mode == "ascii2binary":
        (work / "vector.txt").write_text("\n".join(TEXT) + "\n")
        r = subprocess.run([program, "-i", "vector.txt", "-o", "vector.vec"], cwd=work, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout, r.stderr); print("FAIL: program exited with", r.returncode); return 1
        got = read_vec(work / "vector.vec")
        bad = [(t, g) for t, g, v in zip(TEXT, got, VALUES) if struct.pack("<d", g) != struct.pack("<d", v)]
        print(f"{len(got)} components read from the binary file; expected {len(VALUES)}; bit mismatches: {len(bad)}")
        for t, g in bad:
            print(f"  text {t} -> binary {g!r}, Python reads {float(t)!r}")
        ok = len(got) == len(VALUES) and not bad
    elif mode == "binary2ascii":
        write_vec(work / "vector.vec", VALUES)
        r = subprocess.run([program, "-i", "vector.vec", "-s"], cwd=work, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout, r.stderr); print("FAIL: program exited with", r.returncode); return 1
        printed = [float(x) for x in r.stdout.split()]
        worst = max((abs(p - v) / max(abs(v), 1e-300) for p, v in zip(printed, VALUES)), default=0.0)
        print(f"{len(printed)} values printed; expected {len(VALUES)}; max relative diff {worst:.3e} (tolerance 1e-13)")
        ok = len(printed) == len(VALUES) and worst <= 1e-13
    else:
        print("usage: vector_io.py ascii2binary|binary2ascii PROGRAM WORK_DIR"); return 2
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
