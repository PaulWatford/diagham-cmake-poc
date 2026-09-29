#!/usr/bin/env python3
"""Independent band-projected many-body ED of interacting lattice models (numpy, no DiagHam).

    python3 tests/oracles/fci_manybody.py --write tests/data/fci_manybody
    python3 tests/oracles/fci_manybody.py --check tests/data/fci_manybody   (selftest.fci_manybody_oracle)

Model (checkerboard lattice, DiagHam FCICheckerboardLatticeModel --single-band):
- sites A at integer cell positions R, B at R + (1/2, 1/2), nx x ny cells, periodic;
- one-body Bloch matrix [[d1 + d3, b], [b*, d1 - d3]] of fci_bands.checkerboard (its eigenvalues are the
  program's band energies; the many-body spectrum does not depend on the phase/sublattice conventions of
  the eigenvectors: b -> b* is complex conjugation and A <-> B leaves the interaction invariant);
- Bloch functions in the position gauge: amplitude a e^{ik.R} / sqrt(Nc) on A(R), b e^{ik.(R + (1/2,1/2))} / sqrt(Nc)
  on B, k = 2 pi (i/nx, j/ny) with i, j in 0..n-1;
- interaction U * sum over the four nearest-neighbour A-B bonds (B at R + (1/2, 1/2) + (+-1/2, +-1/2)) of
  n_A n_B, projected onto the lowest band: H = sum_k eps_k n_k + U sum_bonds P n_A n_B P, in the Fock space
  of the Nc band states (dense, the C(Nc, N) states).
Also the exactly solvable atomic limit of FCIAtomicLimitLatticeModel (atomic_limit).
Rows written: kx ky E for every eigenvalue of every total-momentum sector (kx, ky) = sum of the k indices mod n.
"""
import argparse, itertools, math, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lattice_fermions import _add, _remove  # noqa: E402


def checkerboard_band(nx, ny, t1, t2, tpp, mu_s=0.0):
    import numpy as np
    out = []
    for i in range(nx):
        for j in range(ny):
            kx, ky = 2 * math.pi * i / nx, 2 * math.pi * j / ny
            b = complex(4 * t1 * math.cos(kx / 2) * math.cos(ky / 2) * math.cos(math.pi / 4),
                        4 * t1 * math.sin(kx / 2) * math.sin(ky / 2) * math.sin(math.pi / 4))
            d1 = 4 * tpp * math.cos(kx) * math.cos(ky)
            d3 = mu_s + 2 * t2 * (math.cos(kx) - math.cos(ky))
            w, v = np.linalg.eigh(np.array([[d1 + d3, b], [b.conjugate(), d1 - d3]]))
            out.append(((i, j), (kx, ky), w[0], v[:, 0]))
    return out


def checkerboard_manybody(nx, ny, n_part, U, t1=1.0, t2=1.0 / (2 + math.sqrt(2)), tpp=1.0 / (2 + 2 * math.sqrt(2)), mu_s=0.0, flat=False):
    """[(kx, ky, E)] for every eigenstate of the N-particle band-projected problem."""
    import numpy as np
    bands = checkerboard_band(nx, ny, t1, t2, tpp, mu_s)
    nc = nx * ny
    # site amplitudes: phiA[a][(x, y)], phiB[a][(x, y)] with B position (x + 1/2, y + 1/2)
    phiA, phiB = [], []
    for (_, (kx, ky), _, v) in bands:
        pa, pb = {}, {}
        for x in range(nx):
            for y in range(ny):
                pa[x, y] = v[0] * np.exp(1j * (kx * x + ky * y)) / math.sqrt(nc)
                pb[x, y] = v[1] * np.exp(1j * (kx * (x + .5) + ky * (y + .5))) / math.sqrt(nc)
        phiA.append(pa)
        phiB.append(pb)
    # bonds: A(x, y) with the B sites at position (x + 1/2 + sx/2, y + 1/2 + sy/2): cell (x + (1 + sx)/2, ...) ->
    # B cell index c means position c + 1/2, so the four B cells are x-1, x and y-1, y
    W = np.zeros((nc,) * 4, dtype=complex)         # coefficient of a+_a a+_b a_d a_c
    for x in range(nx):
        for y in range(ny):
            for bx in (x - 1, x):
                for by in (y - 1, y):
                    Bc = (bx % nx, by % ny)
                    fa = np.array([phiA[a][x, y] for a in range(nc)])
                    fb = np.array([phiB[a][Bc] for a in range(nc)])
                    # position phases must use the un-wrapped B position: phi at cell bx (not wrapped) equals the wrapped one
                    # because k*N is a multiple of 2 pi
                    W += U * np.einsum("a,b,d,c->abcd", fa.conj(), fb.conj(), fb, fa)
    states = list(itertools.combinations(range(nc), n_part))
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)), dtype=complex)
    nzW = [(a, b, c, d, W[a, b, c, d]) for a in range(nc) for b in range(nc) for c in range(nc) for d in range(nc)
           if abs(W[a, b, c, d]) > 1e-14]
    eps = [0.0 if flat else b[2] for b in bands]
    for s in states:
        col = idx[s]
        H[col, col] += sum(eps[a] for a in s)
        for a, b, c, d, w in nzW:
            r = _remove(s, c)
            if r is None:
                continue
            r2 = _remove(r[1], d)
            if r2 is None:
                continue
            r3 = _add(r2[1], b)
            if r3 is None:
                continue
            r4 = _add(r3[1], a)
            if r4 is None:
                continue
            H[idx[r4[1]], col] += w * r[0] * r2[0] * r3[0] * r4[0]
    assert np.abs(H - H.conj().T).max() < 1e-10, np.abs(H - H.conj().T).max()
    mom = [tuple(sum(bands[a][0][i] for a in s) % n for i, n in enumerate((nx, ny))) for s in states]
    rows = []
    for m in sorted(set(mom)):
        sel = [i for i, mm in enumerate(mom) if mm == m]
        Hs = H[np.ix_(sel, sel)]
        assert np.abs(H[np.ix_(sel, [i for i in range(len(states)) if i not in sel])]).max(initial=0) < 1e-9
        for e in np.linalg.eigvalsh(Hs):
            rows.append((m[0], m[1], float(e)))
    return sorted(rows, key=lambda r: r[2])


def atomic_limit(ncell, n_part, U, boson=True):
    """[(0, 0, E)] for FCIAtomicLimitLatticeModel --single-band: every cell hosts one localised orbital and only
    same-cell pairs of bosons interact; the convention of the program is 4 U per pair
    (E = 4 U sum_c n_c (n_c - 1) / 2). Fermions are inert (a spinless orbital cannot interact with itself): E = 0.
    Momentum labels are not recorded (a count per sector needs translation eigenstates): kx = ky = 0."""
    occs = itertools.combinations_with_replacement(range(ncell), n_part) if boson else itertools.combinations(range(ncell), n_part)
    rows = []
    for occ in occs:
        e = 4.0 * U * sum(occ.count(c) * (occ.count(c) - 1) / 2 for c in set(occ)) if boson else 0.0
        rows.append((0, 0, e))
    return sorted(rows, key=lambda r: r[2])


CASES = {
    "checkerboard_3x3_n2_u1": lambda: checkerboard_manybody(3, 3, 2, 1.0),
    "checkerboard_3x3_n2_u0": lambda: checkerboard_manybody(3, 3, 2, 0.0),
    "checkerboard_3x3_n2_flat": lambda: checkerboard_manybody(3, 3, 2, 1.0, flat=True),
    "atomic_boson_3x3_n2_u1.5": lambda: atomic_limit(9, 2, 1.5),
    "atomic_boson_3x3_n3_u1": lambda: atomic_limit(9, 3, 1.0),
    "atomic_boson_4x2_n2_u1": lambda: atomic_limit(8, 2, 1.0),
    "atomic_fermion_3x3_n2": lambda: atomic_limit(9, 2, 1.0, boson=False),
    "checkerboard_3x3_n3_flat": lambda: checkerboard_manybody(3, 3, 3, 1.0, flat=True),
    "checkerboard_4x3_n2_flat_t2_0.2_tpp_0.1": lambda: checkerboard_manybody(4, 3, 2, 1.0, t2=0.2, tpp=0.1, flat=True),
    "checkerboard_4x3_n2_u1_t2_0.2_tpp_0.1": lambda: checkerboard_manybody(4, 3, 2, 1.0, t2=0.2, tpp=0.1),
}


def render(rows, note):
    return f"# {note}; from tests/oracles/fci_manybody.py, do not edit\n# kx ky E\n" + "".join(f"{a} {b} {e:.15g}\n" for a, b, e in rows)


def render_all():
    return {name + ".dat": render(f(), name) for name, f in CASES.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write")
    ap.add_argument("--check")
    a = ap.parse_args()
    files = render_all()
    if a.write:
        Path(a.write).mkdir(parents=True, exist_ok=True)
        for n, t in files.items():
            (Path(a.write) / n).write_text(t)
        print(f"wrote {len(files)} files")
    if a.check:
        bad = 0
        for n, t in files.items():
            p = Path(a.check) / n
            if not p.exists():
                print("missing", n); bad += 1; continue
            old = [l.split() for l in p.read_text().splitlines() if l and not l.startswith("#")]
            new = [l.split() for l in t.splitlines() if l and not l.startswith("#")]
            if len(old) != len(new) or any(o[:2] != m[:2] or abs(float(o[2]) - float(m[2])) > 1e-9 for o, m in zip(old, new)):
                print("differs", n); bad += 1
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
