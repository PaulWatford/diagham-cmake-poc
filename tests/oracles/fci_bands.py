#!/usr/bin/env python3
"""Band structures of the FCI lattice models from their published Bloch Hamiltonians (no DiagHam).

    python3 tests/oracles/fci_bands.py --write tests/data/fci_bands
    python3 tests/oracles/fci_bands.py --check tests/data/fci_bands      (selftest.fci_bands_oracle)

Checkerboard lattice (Sun, Gu, Katsura & Das Sarma, PRL 106, 236803 (2011)), two bands:
    h(k) = d1(k) 1 + Re B1(k) sigma_x + Im B1(k) sigma_y + d3(k) sigma_z
    B1 = 4 t1 [cos(kx/2) cos(ky/2) cos(phi) + i sin(kx/2) sin(ky/2) sin(phi)],  phi = pi/4
    d1 = 4 tpp cos kx cos ky,   d3 = mu_s + 2 t2 (cos kx - cos ky)
    E_{0,1}(k) = d1 -/+ sqrt(|B1|^2 + d3^2),  k = 2 pi (i/Nx, j/Ny).
The program's default hoppings t1 = 1, t2 = 1/(2+sqrt 2), tpp = 1/(2+2 sqrt 2) are the paper's
flat-band point. The expected file lists kx ky E_0 E_1 per momentum, in the order the program's
--export-onebodytext writes them; the tests compare columns E_0 and E_1 (check_spectrum spectrum).
"""
import argparse, math, sys
from pathlib import Path


def checkerboard(nx, ny, t1=1.0, t2=1.0 / (2 + math.sqrt(2)), tpp=1.0 / (2 + 2 * math.sqrt(2)), mu_s=0.0):
    rows = []
    for i in range(nx):
        for j in range(ny):
            kx, ky = 2 * math.pi * i / nx, 2 * math.pi * j / ny
            b1 = complex(4 * t1 * math.cos(kx / 2) * math.cos(ky / 2) * math.cos(math.pi / 4),
                         4 * t1 * math.sin(kx / 2) * math.sin(ky / 2) * math.sin(math.pi / 4))
            d1 = 4 * tpp * math.cos(kx) * math.cos(ky)
            d3 = mu_s + 2 * t2 * (math.cos(kx) - math.cos(ky))
            g = math.sqrt(abs(b1) ** 2 + d3 * d3)
            rows.append((i, j, d1 - g, d1 + g))
    return rows


CASES = [
    ("checkerboard_3x3", lambda: checkerboard(3, 3), "default hoppings (flat-band point), 3x3 unit cells"),
    ("checkerboard_4x4", lambda: checkerboard(4, 4), "default hoppings, 4x4"),
    ("checkerboard_4x3_t2_0.2_tpp_0.1", lambda: checkerboard(4, 3, t2=0.2, tpp=0.1), "away from the flat-band point"),
]


def render_all():
    return {f"{name}.dat": f"# kx ky E_0 E_1 -- {note}; from the Bloch Hamiltonian in tests/oracles/fci_bands.py, do not edit\n"
            + "".join(f"{i} {j} {e0:.15g} {e1:.15g}\n" for i, j, e0, e1 in f()) for name, f, note in CASES}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", metavar="DIR")
    ap.add_argument("--check", metavar="DIR")
    a = ap.parse_args()
    files = render_all()
    if a.write:
        Path(a.write).mkdir(parents=True, exist_ok=True)
        for n, t in files.items():
            (Path(a.write) / n).write_text(t)
        print(f"wrote {len(files)} files to {a.write}")
    if a.check:
        bad = 0
        for n, t in files.items():
            p = Path(a.check) / n
            if not p.exists():
                bad += 1; print(f"MISSING {p}"); continue
            av = [[float(x) for x in l.split()[2:]] for l in p.read_text().splitlines() if l and not l.startswith("#")]
            bv = [[float(x) for x in l.split()[2:]] for l in t.splitlines() if l and not l.startswith("#")]
            if len(av) != len(bv) or max(abs(x - y) for ra, rb in zip(av, bv) for x, y in zip(ra, rb)) > 1e-12:
                bad += 1; print(f"MISMATCH {p}")
        print(f"{len(files)} files, {bad} mismatches")
        sys.exit(1 if bad else 0)
    if not a.write and not a.check:
        for n, t in files.items():
            print(f"== {n}\n" + "\n".join(t.splitlines()[:4]))


if __name__ == "__main__":
    main()
