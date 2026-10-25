#!/usr/bin/env python3
"""Independent exact diagonalisation of the Hofstadter model on a torus (numpy, no DiagHam).

    python3 tests/oracles/hofstadter.py --write tests/data/hofstadter
    python3 tests/oracles/hofstadter.py --check tests/data/hofstadter   (selftest.hofstadter_oracle)

Model (DiagHam FCIHofstadterModel --real-space --full-momentum, square lattice):
- Lx = x * X by Ly = y * Y sites, periodic, flux phi = q / (X * Y) quanta per plaquette;
- H = -t sum (e^{i A} c+_r' c_r + h.c.) over nearest neighbours, t = 1, in the Landau gauge: the hop
  (i, j) -> (i + 1, j) carries e^{2 pi i phi j}, the hop (i, j) -> (i, j + 1) carries 1. This gauge is
  periodic when phi * Ly is an integer, and it reproduces the program's single-particle spectrum
  (checked: --singleparticle-spectrum, 4 x 6 sites, phi = 1/3; the other Landau gauge gives a different
  torus spectrum, so the holonomies matter);
- hard-core bosons: no other term (the on-site U cannot act); fermions: U sum over nearest-neighbour bonds
  of n_r n_r' (the program's --u-potential for fermions, default 1).
Rows written: E for every eigenvalue of the full Fock space with N particles.

Uses: the one-particle spectra are the physics goldens of the program's lowest band (its default band
projection; the lowest band is the lowest x * y levels, separated by a gap at these fluxes). The hard-core
boson spectra are what --real-space should give, and do not match it (open finding U45): the program's
ground state lies below the rigorous bound 2 e_0 set by its own one-particle spectrum.
"""
import argparse, cmath, itertools, math, sys
from pathlib import Path


def hopping(Lx, Ly, phi):
    """[(a, b, amplitude)]: c+_a c_b terms of the hopping Hamiltonian, both directions."""
    idx = lambda i, j: (i % Lx) * Ly + (j % Ly)
    terms = []
    for i in range(Lx):
        for j in range(Ly):
            a = idx(i, j)
            for b, phase in ((idx(i + 1, j), cmath.exp(2j * math.pi * phi * j)), (idx(i, j + 1), 1.0)):
                terms.append((b, a, -phase))
                terms.append((a, b, -phase.conjugate() if isinstance(phase, complex) else -phase))
    return terms


def bonds(Lx, Ly):
    idx = lambda i, j: (i % Lx) * Ly + (j % Ly)
    return [(idx(i, j), idx(i + 1, j)) for i in range(Lx) for j in range(Ly)] + \
           [(idx(i, j), idx(i, j + 1)) for i in range(Lx) for j in range(Ly)]


def spectrum(n, Lx, Ly, phi, fermion, U=0.0):
    import numpy as np
    N = Lx * Ly
    states = [sum(1 << s for s in c) for c in itertools.combinations(range(N), n)]
    index = {s: k for k, s in enumerate(states)}
    H = np.zeros((len(states), len(states)), dtype=complex)
    terms = hopping(Lx, Ly, phi); bl = bonds(Lx, Ly)
    for k, s in enumerate(states):
        for a, b, amp in terms:                       # c+_a c_b
            if not (s >> b) & 1 or (a != b and (s >> a) & 1):
                continue
            t = s & ~(1 << b)
            sign = 1
            if fermion:
                sign *= (-1) ** bin(t & ((1 << b) - 1)).count("1")
                sign *= (-1) ** bin(t & ((1 << a) - 1)).count("1")
            t |= 1 << a
            H[index[t], k] += sign * amp
        if U:
            H[k, k] += U * sum(1 for p, q in bl if (s >> p) & 1 and (s >> q) & 1)
    return sorted(float(e) for e in np.linalg.eigvalsh(H))


CASES = {
    # 4 x 2 cells of 1 x 3 sites, one flux quantum per cell: 24 sites, phi = 1/3
    "hardcore_bosons_n2_x4_y2_X1_Y3_q1": lambda: spectrum(2, 4, 6, 1 / 3, fermion=False),
    "fermions_n1_x4_y2_X1_Y3_q1": lambda: spectrum(1, 4, 6, 1 / 3, fermion=True),
    "fermions_n1_x3_y2_X1_Y4_q1": lambda: spectrum(1, 3, 8, 1 / 4, fermion=True),
    "fermions_n1_x2_y3_X1_Y4_q1": lambda: spectrum(1, 2, 12, 1 / 4, fermion=True),
    # 3 x 2 cells of 1 x 4 sites, phi = 1/4: 24 sites; three hard-core bosons
    "hardcore_bosons_n3_x3_y2_X1_Y4_q1": lambda: spectrum(3, 3, 8, 1 / 4, fermion=False),
}


def render(name, rows):
    return f"# {name}; from tests/oracles/hofstadter.py, do not edit\n# E\n" + "".join(f"{e:.15g}\n" for e in rows)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--write"); ap.add_argument("--check")
    a = ap.parse_args()
    files = {n + ".dat": render(n, f()) for n, f in CASES.items()}
    if a.write:
        Path(a.write).mkdir(parents=True, exist_ok=True)
        for n, t in files.items():
            (Path(a.write) / n).write_text(t)
        print(f"wrote {len(files)} files")
    if a.check:
        bad = 0
        for n, t in files.items():
            p = Path(a.check) / n
            old = [float(l) for l in p.read_text().splitlines() if l and not l.startswith("#")] if p.exists() else None
            new = [float(l) for l in t.splitlines() if l and not l.startswith("#")]
            if old is None or len(old) != len(new) or any(abs(x - y) > 1e-9 for x, y in zip(old, new)):
                print("differs", n); bad += 1
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
