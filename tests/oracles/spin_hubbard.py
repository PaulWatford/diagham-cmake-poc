#!/usr/bin/env python3
"""Independent answers for the spin-chain and Hubbard goldens.

    python3 tests/oracles/spin_hubbard.py --write tests/data/spin_hubbard
    python3 tests/oracles/spin_hubbard.py --check tests/data/spin_hubbard   (selftest.spin_hubbard_oracle)

Closed forms (no numerics beyond a sum):
- XX ring, J sum_i (Sx Sx + Sy Sy): Jordan-Wigner free fermions with hopping
  J/2; N = L/2 fermions with periodic (N odd) or antiperiodic (N even)
  momenta; E0 = -J sum over the N lowest cos k.
- AKLT ring, sum_i [S_i.S_j + (1/3)(S_i.S_j)^2] for spin 1: E0 = -2L/3
  exactly (each bond projector term is -2/3 on the valence-bond state).
- Haldane-Shastry ring, sum_{i<j} S_i.S_j / d(i,j)^2 with the chord distance
  d = (L/pi) sin(pi|i-j|/L): E0 = -(pi^2/24) (L + 5/L) (Haldane 1988).
- Square-lattice tight binding, t = 1, periodic: eps(k) = -2(cos kx + cos ky);
  at U = 0 the N-electron ground energy is twice the sum of the N/2 lowest
  eps (both spins).
Independent numerics (numpy, dense diagonalisation of the full 2^L or 3^L
matrix, no symmetry, nothing shared with DiagHam):
- spin-1/2 Heisenberg ring L=8, spin-1 Heisenberg ring L=6,
- transverse-field Ising ring L=8: SpinChainXYZ uses Pauli matrices for the
  couplings and S = sigma/2 for the field, H = Jx sum sigma^x sigma^x + h sum S^z
  (established by matching, then checked at 1e-12; the generic XYZ case
  with three couplings did not match any simple convention and is left out).
"""
import argparse, itertools, math, sys
from pathlib import Path


def xx_ring(L, J=1.0):
    N = L // 2
    shift = 0.5 if N % 2 == 0 else 0.0
    e = sorted(-math.cos(2 * math.pi * (n + shift) / L) for n in range(L))
    return J * sum(e[:N])


def aklt_ring(L):
    return -2.0 * L / 3.0


def haldane_shastry(L):
    return -(math.pi ** 2 / 24.0) * (L + 5.0 / L)


def tight_binding_square(Lx, Ly, N, t=1.0):
    eps = sorted(-2 * t * (math.cos(2 * math.pi * i / Lx) + math.cos(2 * math.pi * j / Ly)) for i in range(Lx) for j in range(Ly))
    return 2 * sum(eps[:N // 2])


def tight_binding_square_sectors(Lx, Ly, N, t=1.0):
    """U = 0 ground energy of every total-momentum sector (kx, ky), N/2 up + N/2 down electrons:
    the lowest sum of single-particle energies over the fillings whose momenta add up to (kx, ky)."""
    orbitals = [(i, j, -2 * t * (math.cos(2 * math.pi * i / Lx) + math.cos(2 * math.pi * j / Ly)))
                for i in range(Lx) for j in range(Ly)]
    best = {}
    for up in itertools.combinations(orbitals, N // 2):
        for down in itertools.combinations(orbitals, N - N // 2):
            k = (sum(o[0] for o in up + down) % Lx, sum(o[1] for o in up + down) % Ly)
            e = sum(o[2] for o in up + down)
            if k not in best or e < best[k]:
                best[k] = e
    return [(kx, ky, best[(kx, ky)]) for kx in range(Lx) for ky in range(Ly)]


def _numpy_spin_chain(L, s2, terms, field=0.0, coupling_scale=1.0):
    import numpy as np
    d = s2 + 1
    m = np.arange(s2 / 2, -s2 / 2 - 1, -1.0)
    sz = np.diag(m)
    sp = np.zeros((d, d))
    for i in range(1, d):
        sp[i - 1, i] = math.sqrt((s2 / 2) * (s2 / 2 + 1) - m[i] * (m[i] + 1))
    sx = (sp + sp.T) / 2
    sy = (sp - sp.T) / 2j
    I = np.eye(d)

    def op(o, i):
        r = np.array([[1.0 + 0j]])
        for j in range(L):
            r = np.kron(r, o if j == i else I)
        return r
    S = {"x": [op(sx, i) for i in range(L)], "y": [op(sy, i) for i in range(L)], "z": [op(sz, i) for i in range(L)]}
    H = np.zeros((d ** L, d ** L), dtype=complex)
    for axis, J in terms.items():
        for i in range(L):
            H += coupling_scale * J * S[axis][i] @ S[axis][(i + 1) % L]
    for i in range(L):
        H += field * S["z"][i]
    return float(np.linalg.eigvalsh(H)[0])


CASES = [
    ("heisenberg_ring_L8", lambda: _numpy_spin_chain(8, 1, {"x": 1, "y": 1, "z": 1}), 1e-12, "numpy dense ED, 256 states"),
    ("spin1_heisenberg_ring_L6", lambda: _numpy_spin_chain(6, 2, {"x": 1, "y": 1, "z": 1}), 1e-11, "numpy dense ED, 729 states"),
    ("tfim_ring_L8_h0.5", lambda: _numpy_spin_chain(8, 1, {"x": 1}, field=0.5, coupling_scale=4), 1e-11, "numpy dense ED; Jx sigma^x sigma^x + h S^z"),
    ("xx_ring_L6", lambda: xx_ring(6), 1e-12, "free fermions: -sum of 3 lowest cos k, periodic k"),
    ("xx_ring_L8", lambda: xx_ring(8), 1e-12, "free fermions: -sum of 4 lowest cos k, antiperiodic k"),
    ("xx_ring_L10", lambda: xx_ring(10), 1e-12, "free fermions: -sum of 5 lowest cos k, periodic k"),
    ("aklt_ring_L6", lambda: aklt_ring(6), 1e-12, "-2L/3 exactly"),
    ("aklt_ring_L8", lambda: aklt_ring(8), 1e-12, "-2L/3 exactly"),
    ("haldane_shastry_L6", lambda: haldane_shastry(6), 1e-12, "-(pi^2/24)(L + 5/L)"),
    ("haldane_shastry_L8", lambda: haldane_shastry(8), 1e-12, "-(pi^2/24)(L + 5/L)"),
    ("haldane_shastry_L10", lambda: haldane_shastry(10), 1e-12, "-(pi^2/24)(L + 5/L)"),
    ("hubbard_2x4_N8_U0", lambda: tight_binding_square(2, 4, 8), 1e-10, "tight binding: 2 x (-4-2-2+0)"),
    ("hubbard_4x2_N8_U0", lambda: tight_binding_square(4, 2, 8), 1e-10, "tight binding, same spectrum"),
    ("hubbard_3x2_N6_U0", lambda: tight_binding_square(3, 2, 6), 1e-10, "tight binding: 2 x (-4-1-1)"),
]


def render():
    lines = ["# name  value  tolerance  # how the value is known; generated by tests/oracles/spin_hubbard.py, do not edit"]
    for name, f, tol, note in CASES:
        lines.append(f"{name} {f():.16g} {tol:g}  # {note}")
    return "\n".join(lines) + "\n"


SECTOR_FILE = "hubbard_2x4_N8_U0_sectors.dat"


def render_sectors():
    return ("# kx ky E: tight-binding ground energy of each total-momentum sector of the 2x4 Hubbard model at U = 0,\n"
            "# 4 up + 4 down electrons (tests/oracles/spin_hubbard.py tight_binding_square_sectors; do not edit)\n"
            + "".join(f"{kx} {ky} {e:.16g}\n" for kx, ky, e in tight_binding_square_sectors(2, 4, 8)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", metavar="DIR")
    ap.add_argument("--check", metavar="DIR")
    a = ap.parse_args()
    text = render()
    if a.write:
        Path(a.write).mkdir(parents=True, exist_ok=True)
        (Path(a.write) / "values.txt").write_text(text)
        (Path(a.write) / SECTOR_FILE).write_text(render_sectors())
        print(f"wrote {len(CASES)} values to {a.write}/values.txt and {SECTOR_FILE}")
    if a.check:
        # numeric comparison within each case's tolerance: the last digits of the numpy
        # results depend on the BLAS build and thread count, the physics does not
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
                    print(f"MISMATCH {oname}: committed {oval!r}, recomputed {nval!r}, tolerance {ntol:g}")
        print("values.txt", "matches within tolerance" if ok else "MISMATCH")
        ps = Path(a.check) / SECTOR_FILE
        sectors_ok = ps.exists() and ps.read_text() == render_sectors()
        print(SECTOR_FILE, "matches" if sectors_ok else "MISMATCH")
        sys.exit(0 if (ok and sectors_ok) else 1)
    if not a.write and not a.check:
        print(text)


if __name__ == "__main__":
    main()
