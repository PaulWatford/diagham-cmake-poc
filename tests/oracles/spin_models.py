#!/usr/bin/env python3
"""Independent dense exact diagonalisation of spin models from explicit operator lists (numpy, no DiagHam).

    python3 tests/oracles/spin_models.py --write tests/data/spin_models
    python3 tests/oracles/spin_models.py --check tests/data/spin_models    (selftest.spin_models_oracle)

A model is a list of terms (coefficient, [(site, operator), ...]) on L spins of
length s (two_s = 2s), operators among x, y, z, +, -, and polynomials thereof
built by the helpers below; the full (2s+1)^L matrix is assembled with Kronecker
products and diagonalised with numpy. Every eigenvalue of the requested
sector is written: the Sz sector (total S^z), the momentum sector of a periodic
chain (translation by one site, eigenvalue exp(2 pi i k / L)), or the whole
spectrum. Conventions of the programs (which operator each option multiplies,
S or sigma, boundary terms) are fixed on two-spin cases and stated per case.
"""
import argparse, itertools, math, sys
from pathlib import Path


def spin_ops(two_s):
    import numpy as np
    d = two_s + 1
    m = np.arange(two_s / 2, -two_s / 2 - 1, -1.0)
    sz = np.diag(m)
    sp = np.zeros((d, d))
    for i in range(1, d):
        sp[i - 1, i] = math.sqrt((two_s / 2) * (two_s / 2 + 1) - m[i] * (m[i] + 1))
    sm = sp.T.copy()
    sx = (sp + sm) / 2
    sy = (sp - sm) / (2j)
    return {"x": sx, "y": sy, "z": sz, "+": sp, "-": sm, "1": np.eye(d)}


def build(L, two_s, terms):
    """Dense H from terms [(coeff, [(site, op or matrix), ...]), ...]; op is a key of spin_ops or a (d x d) matrix."""
    import numpy as np
    ops = spin_ops(two_s)
    d = two_s + 1
    H = np.zeros((d ** L, d ** L), dtype=complex)
    for coeff, factors in terms:
        mats = [np.eye(d, dtype=complex) for _ in range(L)]
        for site, op in factors:
            o = ops[op] if isinstance(op, str) else op
            mats[site % L] = mats[site % L] @ o
        M = mats[0]
        for k in range(1, L):
            M = np.kron(M, mats[k])
        H += coeff * M
    return H


def sz_total(L, two_s):
    import numpy as np
    return np.real(np.diag(build(L, two_s, [(1.0, [(i, "z")]) for i in range(L)])))


def translation(L, two_s):
    """Matrix of T: site i -> i + 1 on the product basis."""
    import numpy as np
    d = two_s + 1
    n = d ** L
    T = np.zeros((n, n))
    for idx in range(n):
        digits = []
        x = idx
        for _ in range(L):
            digits.append(x % d)
            x //= d
        digits = digits[::-1]              # digits[0] = site 0 (most significant in the Kronecker order)
        new = digits[-1:] + digits[:-1]    # site i -> i + 1
        j = 0
        for t in new:
            j = j * d + t
        T[j, idx] = 1.0
    return T


def spectrum(L, two_s, terms, sz=None, k=None, hermitian_tol=1e-10):
    """Every eigenvalue; sz = twice the total S^z of the sector (None = all), k = momentum sector (periodic chain)."""
    import numpy as np
    H = build(L, two_s, terms)
    assert np.abs(H - H.conj().T).max() < hermitian_tol, "not Hermitian"
    n = H.shape[0]
    B = np.eye(n, dtype=complex)
    if sz is not None:
        keep = np.isclose(2 * sz_total(L, two_s), sz)
        B = B[:, keep]
    if k is not None:
        T = translation(L, two_s)
        Tb = B.conj().T @ T @ B
        w, v = np.linalg.eig(Tb)
        sel = np.isclose(np.angle(w) * L / (2 * math.pi) % L, k % L, atol=1e-6) | np.isclose(np.angle(w) * L / (2 * math.pi) % L, (k % L) - L, atol=1e-6)
        if not sel.any():
            return np.zeros(0)
        q, _ = np.linalg.qr(v[:, sel])
        B = B @ q
    return np.sort(np.linalg.eigvalsh(B.conj().T @ H @ B).real)


# --------------------------------------------------------------- model builders (DiagHam's conventions per case)
def heisenberg_xyz(L, jx, jy, jz, periodic, h=0.0, hx=0.0, hy=0.0):
    terms = []
    bonds = range(L if periodic else L - 1)
    for i in bonds:
        terms += [(jx, [(i, "x"), (i + 1, "x")]), (jy, [(i, "y"), (i + 1, "y")]), (jz, [(i, "z"), (i + 1, "z")])]
    for i in range(L):
        if h:
            terms.append((h, [(i, "z")]))
        if hx:
            terms.append((hx, [(i, "x")]))
        if hy:
            terms.append((hy, [(i, "y")]))
    return terms


def square_lattice(nx, ny, jxy, jz, periodic=True, open_x=False):
    """S^x S^x + S^y S^y (jxy) and S^z S^z (jz) on the nearest-neighbour bonds of an nx x ny square lattice, site (x, y) -> x * ny + y."""
    terms = []
    bonds = set()
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
                bonds.add(tuple(sorted((i, j))))
    for i, j in sorted(bonds):
        terms += [(jxy, [(i, "x"), (j, "x")]), (jxy, [(i, "y"), (j, "y")]), (jz, [(i, "z"), (j, "z")])]
    return terms, sorted(bonds)


# --------------------------------------------------------------- the programs' conventions, each fixed on the cases below
def xyz_chain(L, jx, jy, jz, periodic, h=0.0):
    """SpinChainXYZ: H = -sum_bonds (Jx sigma^x sigma^x + Jy sigma^y sigma^y + Jz sigma^z sigma^z) + h sum S^z (Pauli
    couplings, the source negates them); with -b 1 the boundary xx and yy bonds are multiplied by the parity
    P = prod_i sigma^z_i (its "translation invariant boundary term"), the boundary zz bond is plain."""
    terms = []
    for i in range(L - 1):
        terms += [(-4 * jx, [(i, "x"), (i + 1, "x")]), (-4 * jy, [(i, "y"), (i + 1, "y")]), (-4 * jz, [(i, "z"), (i + 1, "z")])]
    if periodic:
        terms.append((-4 * jz, [(L - 1, "z"), (0, "z")]))
        parity = [(i, "z") for i in range(L)]
        for J, a in ((jx, "x"), (jy, "y")):
            terms.append((-4 * J * 2 ** L, [(L - 1, a), (0, a)] + parity))
    if h:
        terms += [(h, [(i, "z")]) for i in range(L)]
    return terms


def j1j2_chain(L, j1, j2):
    """PeriodicSpinChainJ1J2: J1 S_i.S_i+1 + J2 S_i.S_i+2, periodic, full Heisenberg for both."""
    terms = []
    for i in range(L):
        terms += [(j1, [(i, a), (i + 1, a)]) for a in "xyz"]
        terms += [(j2, [(i, a), (i + 2, a)]) for a in "xyz"]
    return terms


def double_triangle_chain(L, j1, j2, djz1, djz2):
    """DoubleTriangleSpinChain: periodic zigzag, J1 (nearest) and J2 (next-nearest) Heisenberg with
    the z couplings J1 + djz1 and J2 + djz2."""
    terms = []
    for i in range(L):
        terms += [(j1, [(i, "x"), (i + 1, "x")]), (j1, [(i, "y"), (i + 1, "y")]), (j1 + djz1, [(i, "z"), (i + 1, "z")])]
        terms += [(j2, [(i, "x"), (i + 2, "x")]), (j2, [(i, "y"), (i + 2, "y")]), (j2 + djz2, [(i, "z"), (i + 2, "z")])]
    return terms


def tfim_2d(nx, ny, jz, hx, hz, wrap):
    """TwoDimensionalTransverseFieldIsingModel: Jz sum_bonds S^z S^z + hx sum S^x + hz sum S^z, bonds along y then x,
    then the wrap-around bonds (0,k)-(nx-1,k) and (j,0)-(j,ny-1) (for nx = 2 or ny = 2 these double a bond).
    The source adds the wrap-around bonds whether or not --use-periodic is given (U35)."""
    idx = lambda j, k: j * ny + k
    bonds = []
    for j in range(nx):
        for k in range(1, ny):
            bonds.append((idx(j, k - 1), idx(j, k)))
    for k in range(ny):
        for j in range(1, nx):
            bonds.append((idx(j - 1, k), idx(j, k)))
    if wrap:
        for k in range(ny):
            bonds.append((idx(0, k), idx(nx - 1, k)))
        for j in range(nx):
            bonds.append((idx(j, 0), idx(j, ny - 1)))
    terms = [(jz, [(i, "z"), (j, "z")]) for i, j in bonds]
    terms += [(hx, [(i, "x")]) for i in range(nx * ny)] + [(hz, [(i, "z")]) for i in range(nx * ny)]
    return terms


def obrien_fendley_as_coded(L):
    """PeriodicSpinChainOBrienFendley with its default factors: sum_i [S^x S^x + S^y S^y]_{i,i+1}
    + 1/4 sum_i [(S^-_i S^+_i+1)^2 + h.c.] on spin-1 (SpinChainOBrienFendleyHamiltonian, diagonal part zero)."""
    ops = spin_ops(2)
    sp, sm = ops["+"], ops["-"]
    terms = []
    for i in range(L):
        terms += [(1.0, [(i, "x"), (i + 1, "x")]), (1.0, [(i, "y"), (i + 1, "y")])]
        terms += [(0.25, [(i, sm @ sm), (i + 1, sp @ sp)]), (0.25, [(i, sp @ sp), (i + 1, sm @ sm)])]
    return terms


def pair_projector_poly(J, two_s):
    """P_J for two spins s as a polynomial in S_i.S_j: prod_{J' != J} (S.S - x_J')/(x_J - x_J'), x_J = [J(J+1) - 2 s(s+1)]/2."""
    import numpy as np
    s = two_s / 2
    xs = {j: (j * (j + 1) - 2 * s * (s + 1)) / 2 for j in range(0, two_s + 1)}
    poly = np.poly1d([1.0])
    for jp, x in xs.items():
        if jp != J:
            poly = poly * np.poly1d([1.0, -x]) / (xs[J] - x)
    return poly


def projector_chain(L, two_s, projectors, periodic=True):
    """sum_bonds sum_J factor_J P_J(i, i+1) from the S.S polynomials (PeriodicSpinChainGeneralizedAKLT: P_3 + P_4 on spin 2)."""
    import numpy as np
    ops = spin_ops(two_s)
    terms = []
    for J, fac in projectors:
        poly = pair_projector_poly(J, two_s)
        for power, c in enumerate(poly.coeffs[::-1]):
            if abs(c) < 1e-14:
                continue
            for i in range(L if periodic else L - 1):
                if power == 0:
                    terms.append((fac * c, [(i, "1")]))
                    continue
                for combo in itertools.product("xyz", repeat=power):
                    A = np.eye(two_s + 1, dtype=complex)
                    B = np.eye(two_s + 1, dtype=complex)
                    for a in combo:
                        A = A @ ops[a]
                        B = B @ ops[a]
                    terms.append((fac * c, [(i, A), (i + 1, B)]))
    return terms


def disordered_chain(L, j, djz, fields):
    """GenericOpenSpinChainWithDisorder: J (S^x S^x + S^y S^y) + (J + djz) S^z S^z on an open chain, plus sum h_i S^z_i
    from the Disorder = ... line of --disorder-file."""
    return heisenberg_xyz(L, j, j, j + djz, False) + [(fields[i], [(i, "z")]) for i in range(L)]


def potts3_open(L, j, f):
    """Potts3ChainModel: H = -J sum (sigma_i sigma_i+1^+ + h.c.) - f sum (tau_i + tau_i^+), sigma = diag(1, w, w^2), tau the cyclic shift."""
    import numpy as np
    w = complex(math.cos(2 * math.pi / 3), math.sin(2 * math.pi / 3))
    sig = np.diag([1, w, w * w])
    tau = np.roll(np.eye(3), 1, axis=0)
    terms = []
    for i in range(L - 1):
        terms += [(-j, [(i, sig), (i + 1, sig.conj().T)]), (-j, [(i, sig.conj().T), (i + 1, sig)])]
    for i in range(L):
        terms += [(-f, [(i, tau)]), (-f, [(i, tau.T)])]
    return terms


# name, L, 2s, terms, sector (sz = 2 Sz or None, k or None), note
CASES = [
    ("xyz_open_L8_h0.2", 8, 1, lambda: xyz_chain(8, 1.0, 0.6, 0.3, False, h=0.2), None, None, "SpinChainXYZ open, Jx=1 Jy=0.6 Jz=0.3 (Pauli, negated), field 0.2 S^z; whole spectrum"),
    ("xyz_periodic_L8", 8, 1, lambda: xyz_chain(8, 1.0, 0.6, 0.3, True), None, None, "SpinChainXYZ -b 1: parity-twisted boundary bond; whole spectrum"),
    ("fullgeneric_open_L6_fields", 6, 1, lambda: heisenberg_xyz(6, 1.0, 0.6, 0.3, False, h=0.3, hx=0.2, hy=0.1), None, None, "FullGenericOpenSpinChain: S couplings and three fields; whole spectrum"),
    ("heisenberg2d_3x3_jz0.7_sz1", 9, 1, lambda: square_lattice(3, 3, 1.0, 0.7, periodic=True)[0], 1, None, "TwoDimensionalHeisenbergModel 3x3 periodic, J=1 (xx+yy), Jz=0.7, Sz=+1/2 sector (126 states)"),
    ("j1j2_L8_sz0_k0", 8, 1, lambda: j1j2_chain(8, 1.0, 0.4), 0, 0, "PeriodicSpinChainJ1J2 spin-1/2, J2=0.4, Sz=0, k=0"),
    ("j1j2_L8_sz0_k1", 8, 1, lambda: j1j2_chain(8, 1.0, 0.4), 0, 1, "k=1"),
    ("j1j2_L8_sz0_k3", 8, 1, lambda: j1j2_chain(8, 1.0, 0.4), 0, 3, "k=3"),
    ("doubletriangle_L8_sz0", 8, 1, lambda: double_triangle_chain(8, 1.0, 0.5, 0.2, 0.1), 0, None, "DoubleTriangleSpinChain J1=1 J2=0.5 djz1=0.2 djz2=0.1, Sz=0 (70 states)"),
    ("aklt_p3p4_spin2_L4_sz0_k0", 4, 4, lambda: projector_chain(4, 4, [(3, 1.0), (4, 1.0)]), 0, 0, "PeriodicSpinChainGeneralizedAKLT spin-2: P_3 + P_4 on every bond, Sz=0, k=0"),
    ("aklt_p3p4_spin2_L4_sz0_k2", 4, 4, lambda: projector_chain(4, 4, [(3, 1.0), (4, 1.0)]), 0, 2, "k=2"),
    ("tfim2d_2x3_periodic", 6, 1, lambda: tfim_2d(2, 3, 1.0, 0.7, 0.2, True), None, None, "TwoDimensionalTransverseFieldIsingModel 2x3 --use-periodic: 12 bonds (x bonds doubled for nx=2), hx=0.7, hz=0.2; whole spectrum"),
    ("tfim2d_3x3_periodic", 9, 1, lambda: tfim_2d(3, 3, 1.0, 0.7, 0.2, True), None, None, "3x3 periodic, 18 bonds; whole spectrum (512)"),
    ("tfim2d_2x3_open", 6, 1, lambda: tfim_2d(2, 3, 1.0, 0.7, 0.2, False), None, None, "2x3 with open boundaries (7 bonds): the program adds the wrap-around bonds anyway (U35)"),
    ("obrienfendley_spin1_L4_sz0", 4, 2, lambda: obrien_fendley_as_coded(4), 0, None, "PeriodicSpinChainOBrienFendley default factors, spin-1, Sz=0 (19 states)"),
    ("disorder_open_L6_sz0", 6, 1, lambda: disordered_chain(6, 1.0, 0.2, [0.3, -0.1, 0.2, 0.05, -0.25, 0.15]), 0, None, "GenericOpenSpinChainWithDisorder J=1 djz=0.2, fields from the disorder file, Sz=0 (20 states)"),
    ("potts3_open_L6", 6, 2, lambda: potts3_open(6, 1.0, 0.5), None, None, "Potts3ChainModel open, J=1 f=0.5; whole spectrum (729)"),
]


def render_all():
    files = {}
    for name, L, two_s, tf, sz, k, note in CASES:
        e = spectrum(L, two_s, tf(), sz=sz, k=k)
        files[f"{name}_spectrum.dat"] = (f"# L={L} 2s={two_s} sector 2Sz={sz} k={k}: {note}; {len(e)} eigenvalues from tests/oracles/spin_models.py\n"
                                         + "".join(f"{x:.15g}\n" for x in e))
    files["disorder_open_L6.dat"] = "Disorder = 0.3 -0.1 0.2 0.05 -0.25 0.15\n"
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
            if not name.endswith("_spectrum.dat"):
                ok = p.read_text() == text
            else:
                av = [float(l.split()[-1]) for l in p.read_text().splitlines() if l and not l.startswith("#")]
                bv = [float(l.split()[-1]) for l in text.splitlines() if l and not l.startswith("#")]
                ok = len(av) == len(bv) and max(abs(x - y) for x, y in zip(av, bv)) <= a.tol
            if not ok:
                bad += 1
                print(f"MISMATCH {p}")
        print(f"{len(files)} files, {bad} mismatches")
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
