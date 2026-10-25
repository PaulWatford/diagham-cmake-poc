#!/usr/bin/env python3
"""Spectrum of the total angular momentum L^2 on the sphere, from exact state counting (no DiagHam, no numpy).

    python3 tests/oracles/sphere_l2.py --write tests/data/sphere_l2
    python3 tests/oracles/sphere_l2.py --check tests/data/sphere_l2   (selftest.sphere_l2_oracle)

In the sector of total momentum Lz = M (the program's -z is 2M), every multiplet with L >= |M| contributes one
state, with L^2 eigenvalue L(L + 1). The number of multiplets of total L is d(L) - d(L + 1), where d(M) is the
dimension of the Lz = M sector (counted by tests/oracles/dimensions.py). FQHESphereL2Diagonalize diagonalises
L^2 in the sector, so its spectrum is this list.
"""
import argparse, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dimensions import sphere_lz_dimension  # noqa: E402


def l2_spectrum(n, two_s, two_lz, fermion):
    """Sorted L(L+1) values, one per multiplet with 2L >= |2Lz|, 2L of the same parity as 2Lz."""
    rows = []
    two_l = abs(two_lz)
    while True:
        d = sphere_lz_dimension(n, two_s, two_l, fermion)
        if d == 0:
            break
        mult = d - sphere_lz_dimension(n, two_s, two_l + 2, fermion)
        rows += [(two_l / 2) * (two_l / 2 + 1)] * mult
        two_l += 2
    return sorted(rows)


CASES = {
    "fermions_n4_2s9_lz0": (4, 9, 0, True),
    "fermions_n3_2s6_2lz2": (3, 6, 2, True),
    "bosons_n4_2s6_lz0": (4, 6, 0, False),
    "bosons_n3_2s5_2lz3": (3, 5, 3, False),
}


def render(name, rows):
    return f"# {name}; from tests/oracles/sphere_l2.py, do not edit\n# L(L+1)\n" + "".join(f"{e:.15g}\n" for e in rows)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--write"); ap.add_argument("--check")
    a = ap.parse_args()
    files = {n + ".dat": render(n, l2_spectrum(*c)) for n, c in CASES.items()}
    if a.write:
        Path(a.write).mkdir(parents=True, exist_ok=True)
        for n, t in files.items():
            (Path(a.write) / n).write_text(t)
        print(f"wrote {len(files)} files")
    if a.check:
        bad = [n for n, t in files.items() if not (Path(a.check) / n).exists() or (Path(a.check) / n).read_text() != t]
        for n in bad:
            print("differs", n)
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
