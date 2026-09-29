#!/usr/bin/env python3
"""Independent exact diagonalisation of lattice fermion models (numpy, no DiagHam).

    python3 tests/oracles/lattice_fermions.py --write tests/data/lattice_fermions
    python3 tests/oracles/lattice_fermions.py --check tests/data/lattice_fermions   (selftest.lattice_fermions_oracle)

Two builders:
- hubbard_spectrum(hop, U, n_up, n_down): spinful fermions on the sites of a
  Hermitian hopping matrix hop (complex allowed), H = sum_{ij,s} hop[i,j] c+_is c_js
  + U sum_i n_i,up n_i,down, on the explicit Fock basis of the (N_up, N_down) sector
  (occupied-mode tuples, fermion signs from the sorted order).
- quadratic_spectrum(L, h, d, mu): spinless fermions on L sites with hopping h[i,j]
  (Hermitian), pairing d[i,j] (antisymmetric, the term d[i,j] c_i c_j + h.c.) and
  on-site mu[i] n_i, in the full 2^L space through Jordan-Wigner matrices (the
  Kitaev chain).
Momentum is not resolved: the programs are run with their momentum quantum
numbers switched off, or the reference lists every momentum sector.
"""
import argparse, itertools, math, sys
from pathlib import Path


def _remove(block, j):
    if j not in block:
        return None
    i = block.index(j)
    return (-1) ** i, block[:i] + block[i + 1:]


def _add(block, j):
    if j in block:
        return None
    l = sorted(block + (j,))
    return (-1) ** l.index(j), tuple(l)


def hubbard_spectrum(hop, U, n_up, n_down, onsite=None):
    """Every eigenvalue of the (N_up, N_down) sector; modes are site + spin * nsites."""
    import numpy as np
    ns = len(hop)
    onsite = onsite if onsite is not None else [0.0] * ns
    ups = list(itertools.combinations(range(ns), n_up))
    downs = list(itertools.combinations(range(ns), n_down))
    states = [(u, d) for u in ups for d in downs]
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)), dtype=complex)
    nz = [(i, j, hop[i][j]) for i in range(ns) for j in range(ns) if abs(hop[i][j]) > 1e-15]
    for s in states:
        col = idx[s]
        u, d = s
        diag = U * len(set(u) & set(d)) + sum(onsite[i] for i in u) + sum(onsite[i] for i in d)
        H[col, col] += diag
        for which, block in ((0, u), (1, d)):
            for i, j, t in nz:
                r = _remove(block, j)
                if r is None:
                    continue
                a = _add(r[1], i)
                if a is None:
                    continue
                new = (a[1], d) if which == 0 else (u, a[1])
                H[idx[new], col] += t * r[0] * a[0]
    assert np.abs(H - H.conj().T).max() < 1e-10
    return np.sort(np.linalg.eigvalsh(H).real)


def jordan_wigner(L):
    """Annihilation operators c_i (2^L x 2^L) for L spinless modes."""
    import numpy as np
    sz = np.diag([1.0, -1.0])
    sm = np.array([[0.0, 1.0], [0.0, 0.0]])       # |1> -> |0> with basis (|0>, |1>) = (empty, occupied)
    I = np.eye(2)
    cs = []
    for i in range(L):
        mats = [sz] * i + [sm] + [I] * (L - i - 1)
        M = mats[0]
        for m in mats[1:]:
            M = np.kron(M, m)
        cs.append(M)
    return cs


def quadratic_spectrum(L, h, d, mu):
    """Every eigenvalue of H = sum h_ij c+_i c_j + sum_{i<j} (d_ij c_i c_j + h.c.) + sum mu_i n_i on 2^L states."""
    import numpy as np
    cs = jordan_wigner(L)
    H = np.zeros((2 ** L, 2 ** L), dtype=complex)
    for i in range(L):
        for j in range(L):
            if abs(h[i][j]) > 1e-15:
                H += h[i][j] * cs[i].conj().T @ cs[j]
            if i < j and abs(d[i][j]) > 1e-15:
                T = d[i][j] * cs[i] @ cs[j]
                H += T + T.conj().T
        H += mu[i] * cs[i].conj().T @ cs[i]
    assert np.abs(H - H.conj().T).max() < 1e-10
    return np.sort(np.linalg.eigvalsh(H).real)


# ------------------------------------------------------------------ lattices
def square_hopping(nx, ny, t=1.0, periodic=True):
    """Nearest-neighbour hopping -t on an nx x ny square lattice, site (x, y) -> x * ny + y (the bond set, each pair once)."""
    ns = nx * ny
    hop = [[0.0] * ns for _ in range(ns)]
    for x in range(nx):
        for y in range(ny):
            i = x * ny + y
            for dx, dy in ((1, 0), (0, 1)):
                X, Y = x + dx, y + dy
                if periodic:
                    X, Y = X % nx, Y % ny
                elif X >= nx or Y >= ny:
                    continue
                j = X * ny + Y
                if i == j:
                    continue
                hop[i][j] += -t
                hop[j][i] += -t
    return hop


def haldane_bloch(kx, ky, t1, t2, phi, t3=0.0, mu_s=0.0):
    """Haldane honeycomb Bloch Hamiltonian (2 x 2) at k, lattice vectors a1 = (1, 0), a2 = (1/2, sqrt3/2):
    nearest neighbours delta_1..3, next-nearest with the phase phi (in units of pi) in the Haldane sense,
    t3 on the third neighbours, mu_s the staggered potential. Used only through the many-body spectrum,
    which is invariant under the phase-sign and gauge conventions."""
    import numpy as np
    a1 = np.array([1.0, 0.0]); a2 = np.array([0.5, math.sqrt(3) / 2])
    k = np.array([kx, ky])
    d1 = [np.array([0.0, 1 / math.sqrt(3)]), np.array([0.5, -1 / (2 * math.sqrt(3))]), np.array([-0.5, -1 / (2 * math.sqrt(3))])]
    b = [a1, a2 - a1, -a2]               # next-nearest vectors (one orientation)
    f = sum(np.exp(1j * k @ v) for v in d1)
    hz = 2 * t2 * sum(math.sin(math.pi * phi) * 0 + np.sin(k @ v) for v in b) * 2 * math.sin(math.pi * phi)
    h0 = 2 * t2 * math.cos(math.pi * phi) * sum(np.cos(k @ v) for v in b)
    f3 = t3 * sum(np.exp(1j * k @ (-2 * v)) for v in d1)
    return np.array([[h0 + hz / 2 + mu_s, -t1 * f - f3], [-t1 * np.conj(f) - np.conj(f3), h0 - hz / 2 - mu_s]])


def bloch_to_real_space(bloch, nx, ny, norb, a1=(1.0, 0.0), a2=(0.5, math.sqrt(3) / 2)):
    """Real-space hopping matrix of a periodic nx x ny cluster from a Bloch Hamiltonian: the inverse Fourier
    transform over the cluster's momenta; site index = (x * ny + y) * norb + orbital."""
    import numpy as np
    a1 = np.array(a1); a2 = np.array(a2)
    # reciprocal vectors
    A = np.array([a1, a2])
    B = 2 * math.pi * np.linalg.inv(A).T
    ns = nx * ny * norb
    hop = np.zeros((ns, ns), dtype=complex)
    cells = [(x, y) for x in range(nx) for y in range(ny)]
    for mx in range(nx):
        for my in range(ny):
            k = B[0] * mx / nx + B[1] * my / ny
            Hk = bloch(k[0], k[1])
            for (x, y) in cells:
                for (X, Y) in cells:
                    phase = np.exp(1j * (k @ ((x - X) * a1 + (y - Y) * a2)))
                    i0 = (x * ny + y) * norb; j0 = (X * ny + Y) * norb
                    hop[i0:i0 + norb, j0:j0 + norb] += Hk * phase / (nx * ny)
    return hop


# ------------------------------------------------------------------ the programs' models
def haldane_bloch_diagham(x, y, t1, t2, phi, mu_s=0.0):
    """DiagHam's Haldane Bloch Hamiltonian (TightBindingModelHaldaneHoneycombLattice) at the reduced momenta
    x = k.a1, y = k.a2; phi in radians (the program's --phi is in radians unless --phase-in-pi is given).
    B1 = t1 (1 + e^{i(x+y)} + e^{iy}), d0 = 2 t2 cos(phi) (cos x + cos y + cos(x+y)),
    d3 = 2 t2 sin(phi) (sin x + sin y - sin(x+y)) + mu_s; H = [[d0 + d3, B1], [B1*, d0 - d3]]."""
    import numpy as np
    B1 = t1 * complex(1 + math.cos(x + y) + math.cos(y), math.sin(x + y) + math.sin(y))
    d0 = 2 * t2 * math.cos(phi) * (math.cos(x) + math.cos(y) + math.cos(x + y))
    d3 = 2 * t2 * math.sin(phi) * (math.sin(x) + math.sin(y) - math.sin(x + y)) + mu_s
    return np.array([[d0 + d3, B1], [np.conj(B1), d0 - d3]])


def cluster_from_bloch(bloch, nx, ny, norb):
    """Real-space one-body matrix of the periodic nx x ny cluster: the inverse Fourier transform of the Bloch
    Hamiltonian over the cluster momenta (x, y) = 2 pi (mx/nx, my/ny); site = (cx * ny + cy) * norb + orbital.
    This is the model the programs diagonalise (their real-space Hamiltonian comes from the same H(k));
    for clusters where periodic images of a bond coincide it differs from a bond list."""
    import numpy as np
    ns = nx * ny * norb
    hop = np.zeros((ns, ns), dtype=complex)
    for mx in range(nx):
        for my in range(ny):
            x, y = 2 * math.pi * mx / nx, 2 * math.pi * my / ny
            Hk = bloch(x, y)
            for cx in range(nx):
                for cy in range(ny):
                    for CX in range(nx):
                        for CY in range(ny):
                            ph = np.exp(1j * (x * (cx - CX) + y * (cy - CY)))
                            i0 = (cx * ny + cy) * norb
                            j0 = (CX * ny + CY) * norb
                            hop[i0:i0 + norb, j0:j0 + norb] += Hk * ph / (nx * ny)
    return hop


def ssh_hopping(cells, delta, t=1.0):
    """HubbardSSHModel: spinless fermions, cells x 2 sites, intra-cell hopping -t(1 - delta), inter-cell -t(1 + delta), periodic."""
    ns = 2 * cells
    hop = [[0.0] * ns for _ in range(ns)]
    for c in range(cells):
        a, b = 2 * c, 2 * c + 1
        hop[a][b] += -t * (1 - delta); hop[b][a] += -t * (1 - delta)
        nxt = (2 * c + 2) % ns
        hop[b][nxt] += -t * (1 + delta); hop[nxt][b] += -t * (1 + delta)
    return hop


# name, kind, parameters, note
CASES = [
    ("square_2x3_n6_u2_sz0", "hubbard", lambda: (square_hopping(2, 3, 1.0, True), 2.0, 3, 3), "HubbardSquareLatticeModel 2x3, N=6 (3 up 3 down), U=2, every momentum sector; the x bonds are doubled for nx=2 (400 states)"),
    ("haldane_3x3_n2_u0_sz0", "hubbard", lambda: (cluster_from_bloch(lambda x, y: haldane_bloch_diagham(x, y, 1.0, 1.0, 1.0 / 3), 3, 3, 2).tolist(), 0.0, 1, 1), "HubbardHaldaneLatticeModel 3x3, t1=t2=1, phi=1/3 rad, N=2 (1 up 1 down), U=0, every momentum sector (324 states)"),
    ("haldane_3x3_n2_u2_sz0", "hubbard", lambda: (cluster_from_bloch(lambda x, y: haldane_bloch_diagham(x, y, 1.0, 1.0, 1.0 / 3), 3, 3, 2).tolist(), 2.0, 1, 1), "the same with U=2"),
    ("haldane_3x3_n2_u2_mus0.1_sz0", "hubbard", lambda: ((cluster_from_bloch(lambda x, y: haldane_bloch_diagham(x, y, 1.0, 1.0, 1.0 / 3), 3, 3, 2) + __import__("numpy").diag([0.1, 0.0] * 9)).tolist(), 2.0, 1, 1), "the same with U=2 and --mu-s 0.1: the program puts mu_s on the A sublattice only"),
    ("ssh_4cells_n4_delta0.3", "hubbard", lambda: (ssh_hopping(4, 0.3), 0.0, 4, 0), "HubbardSSHModel 4 unit cells, N=4 spinless, delta=0.3, translations off (70 states)"),
]


def render_all():
    files = {}
    for name, kind, par, note in CASES:
        hop, U, n_up, n_down = par()
        e = hubbard_spectrum(hop, U, n_up, n_down)
        files[f"{name}_spectrum.dat"] = (f"# {note}; {len(e)} eigenvalues from tests/oracles/lattice_fermions.py\n"
                                         + "".join(f"{x:.15g}\n" for x in e))
    return files


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", metavar="DIR")
    ap.add_argument("--check", metavar="DIR")
    ap.add_argument("--tol", type=float, default=1e-9)
    a = ap.parse_args()
    files = render_all()
    if a.write:
        Path(a.write).mkdir(parents=True, exist_ok=True)
        for name, text in files.items():
            (Path(a.write) / name).write_text(text)
        print(f"wrote {len(files)} files to {a.write}")
    if a.check:
        bad = 0
        for name, text in files.items():
            p = Path(a.check) / name
            if not p.exists():
                bad += 1
                print(f"MISSING {p}")
                continue
            av = [float(l.split()[-1]) for l in p.read_text().splitlines() if l and not l.startswith("#")]
            bv = [float(l.split()[-1]) for l in text.splitlines() if l and not l.startswith("#")]
            if not (len(av) == len(bv) and max(abs(x - y) for x, y in zip(av, bv)) <= a.tol):
                bad += 1
                print(f"MISMATCH {p}")
        print(f"{len(files)} files, {bad} mismatches")
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
