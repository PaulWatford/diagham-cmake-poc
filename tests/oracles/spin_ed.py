#!/usr/bin/env python3
"""Independent dense exact diagonalisation of spin chains (numpy, no symmetry).

    python3 tests/oracles/spin_ed.py --write tests/data/spin_ed
    python3 tests/oracles/spin_ed.py --check tests/data/spin_ed    (selftest.spin_ed_oracle, needs numpy)

H = sum_i [ J (Sx_i Sx_i+1 + Sy_i Sy_i+1) + Jz Sz_i Sz_i+1 ] + J2 sum_i Sz_i Sz_i+2 + hz sum_i Sz_i
on L spins of length s, open or periodic; the full (2s+1)^L matrix is built
with Kronecker products and diagonalised with numpy.linalg.eigvalsh. The
ground energy of each case is written with a tolerance; the tests compare
the DiagHam program's lowest eigenvalue (Sz = 0 sector) with it.
DiagHam options used: GenericOpenSpinChain -j J -z (Jz - J) --hz-value hz;
GenericPeriodicSpinChain --nn-coupling (Jz - J) --nnn-coupling J2 (its
next-nearest term is Sz Sz only, established by matching the three
possibilities against the program).
"""
import argparse, math, sys
from pathlib import Path


def ground_energy(L, two_s, J=1.0, Jz=None, J2=0.0, hz=0.0, periodic=False):
    import numpy as np
    Jz = J if Jz is None else Jz
    d = two_s + 1
    m = np.arange(two_s / 2, -two_s / 2 - 1, -1.0)
    sz = np.diag(m)
    sp = np.zeros((d, d))
    for i in range(1, d):
        sp[i - 1, i] = math.sqrt((two_s / 2) * (two_s / 2 + 1) - m[i] * (m[i] + 1))
    sx = (sp + sp.T) / 2
    sy = (sp - sp.T) / 2j
    I = np.eye(d)

    def op(o, i):
        r = np.array([[1.0 + 0j]])
        for j in range(L):
            r = np.kron(r, o if j == i else I)
        return r
    Sx = [op(sx, i) for i in range(L)]
    Sy = [op(sy, i) for i in range(L)]
    Sz = [op(sz, i) for i in range(L)]
    H = np.zeros((d ** L, d ** L), dtype=complex)
    bonds = [(i, i + 1) for i in range(L - 1)] + ([(L - 1, 0)] if periodic else [])
    for i, j in bonds:
        H += J * (Sx[i] @ Sx[j] + Sy[i] @ Sy[j]) + Jz * Sz[i] @ Sz[j]
    if J2:
        nnn = [(i, i + 2) for i in range(L - 2)] + ([(L - 2, 0), (L - 1, 1)] if periodic else [])
        for i, j in nnn:
            H += J2 * Sz[i] @ Sz[j]
    for i in range(L):
        H += hz * Sz[i]
    return float(np.linalg.eigvalsh(H)[0])


CASES = [
    # name, kwargs, tolerance, note
    ("open_heisenberg_L8", dict(L=8, two_s=1), 1e-12, "spin-1/2 open chain, 256 states"),
    ("open_spin1_heisenberg_L6", dict(L=6, two_s=2), 1e-11, "spin-1 open chain, 729 states"),
    ("open_xxz_L8_jz1.5_hz0.3", dict(L=8, two_s=1, Jz=1.5, hz=0.3), 1e-12, "open XXZ chain with a field (Sz=0 sector ground state)"),
    ("ring_xxz_L8_jz1.5", dict(L=8, two_s=1, Jz=1.5, periodic=True), 1e-12, "periodic XXZ ring"),
    ("ring_j1j2_L8_j2_0.4", dict(L=8, two_s=1, J2=0.4, periodic=True), 1e-12, "periodic ring with a next-nearest Sz Sz coupling"),
    ("ring_xxz_L10_jz0.5", dict(L=10, two_s=1, Jz=0.5, periodic=True), 1e-12, "periodic XXZ ring, 1024 states"),
]


def render():
    lines = ["# name  value  tolerance  # numpy dense ED, tests/oracles/spin_ed.py; do not edit"]
    for name, kw, tol, note in CASES:
        lines.append(f"{name} {ground_energy(**kw):.16g} {tol:g}  # {note}")
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", metavar="DIR")
    ap.add_argument("--check", metavar="DIR")
    a = ap.parse_args()
    text = render()
    if a.write:
        Path(a.write).mkdir(parents=True, exist_ok=True)
        (Path(a.write) / "values.txt").write_text(text)
        print(f"wrote {len(CASES)} values to {a.write}/values.txt")
    if a.check:
        p = Path(a.check) / "values.txt"
        old = p.read_text().splitlines() if p.exists() else []
        new = text.splitlines()
        ok = len(old) == len(new)
        if ok:
            for o, n in zip(old[1:], new[1:]):
                oname, oval = o.split()[0], float(o.split()[1])
                nname, nval, ntol = n.split()[0], float(n.split()[1]), float(n.split()[2])
                if oname != nname or abs(oval - nval) > ntol:
                    ok = False
        print("values.txt", "matches" if ok else "MISMATCH")
        sys.exit(0 if ok else 1)
    if not a.write and not a.check:
        print(text)


if __name__ == "__main__":
    main()
