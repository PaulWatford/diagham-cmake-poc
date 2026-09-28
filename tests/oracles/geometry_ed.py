#!/usr/bin/env python3
"""Independent exact diagonalisation of two-body pseudopotential Hamiltonians on the cylinder and the disk (numpy, no DiagHam).

    python3 tests/oracles/geometry_ed.py --write tests/data/geometry_ed
    python3 tests/oracles/geometry_ed.py --check tests/data/geometry_ed     (selftest.geometry_ed_oracle)

Cylinder (Landau gauge, circumference L along y, l_B = 1): the N_orb orbitals
have k_j = 2 pi j / L and guiding centre x_j = k_j; the torus matrix element
of tests/oracles/torus_ed.py with no periodic images and the q_x sum turned
into an integral:
    A(1234) = (1/(2 pi L)) int dq_x sum_{q_y = 2 pi t / L} V(q) e^{-q^2/2} cos(q_x (x_1 - x_3)),
    t = j_1 - j_4, momentum conservation j_1 + j_2 = j_3 + j_4,
    H = 1/2 sum A(1234) c+_1 c+_2 c_3 c_4 (c_4 acting first),
with V(q) e^{-q^2/2} -> 4 pi sum_m V_m L_m(q^2) e^{-q^2/2} for Haldane pseudopotentials
(the q = 0 term kept), as on the torus: a pair of relative angular momentum m costs 2 V_m.
DiagHam's cylinder has L = sqrt(2 pi r N_orb) for its option -r (checked at r = 1, 2, 1/2),
and it labels the momentum sectors by twice the momentum (its -y 2k is Ky = k here; -y also
needs --nbr-ky 1 to compute one sector only).

Disk (symmetric gauge): orbitals phi_m = z^m e^{-|z|^2/4} / sqrt(2 pi 2^m m!).
A pair of relative angular momentum m and centre-of-mass angular momentum M is
(z_1 - z_2)^m (z_1 + z_2)^M e^{-(|z_1|^2 + |z_2|^2)/4}, expanded with binomials
into phi_{m1} phi_{m2}; the pseudopotential projector is P_m = sum_M |m, M><m, M|
and H = 2 sum_m V_m sum_{i<j} P_m(ij) (a pair costs 2 V_m, as on the sphere).
The basis of N particles with total angular momentum L_z is every partition
of L_z into N orbital indices (no orbital cutoff, the program has none).
The disk Coulomb pseudopotentials are V_m = Gamma(m + 1/2) / (2 m!) in units of
e^2 / (eps l_B) (the relative wavefunction z^m e^{-|z|^2/8}, <1/r>).
"""
import argparse, itertools, math, sys
from collections import Counter
from pathlib import Path


def laguerre(m, x):
    if m == 0:
        return 1.0
    l0, l1 = 1.0, 1.0 - x
    for k in range(1, m):
        l0, l1 = l1, ((2 * k + 1 - x) * l1 - k * l0) / (k + 1)
    return l1


def _remove(block, j, fermion):
    if j not in block:
        return None
    l = list(block)
    i = l.index(j)
    amp = (-1) ** i if fermion else math.sqrt(block.count(j))
    del l[i]
    return amp, tuple(l)


def _add(block, j, fermion):
    if fermion and j in block:
        return None
    l = sorted(block + (j,))
    amp = (-1) ** l.index(j) if fermion else math.sqrt(l.count(j))
    return amp, tuple(l)


def _assemble(states, norb, V, fermion, conserve):
    """H = 1/2 sum V(1234) c+_1 c+_2 c_3 c_4 (c_4 first) on the given basis; conserve(j1, j2, j3, j4) selects the
    allowed j2 for given j1, j3, j4 (momentum / angular momentum conservation)."""
    import numpy as np
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for s in states:
        col = idx[s]
        for j4 in set(s):
            r4 = _remove(s, j4, fermion)
            for j3 in set(r4[1]):
                r3 = _remove(r4[1], j3, fermion)
                for j1 in range(norb):
                    j2 = conserve(j1, j3, j4)
                    if j2 is None or j2 < 0 or j2 >= norb:
                        continue
                    v = V(j1, j2, j3, j4)
                    if v == 0:
                        continue
                    a2 = _add(r3[1], j2, fermion)
                    if a2 is None:
                        continue
                    a1 = _add(a2[1], j1, fermion)
                    if a1 is None:
                        continue
                    H[idx[a1[1]], col] += 0.5 * r4[0] * r3[0] * a2[0] * a1[0] * v
    return np.linalg.eigvalsh(H)


# ------------------------------------------------------------------ cylinder
def cylinder_spectrum(n, norb, pp, fermion, ratio=1.0, ky=0, kernel=None, qmax=12.0, nq=4001):
    """Every eigenvalue of the total-momentum sector Ky = sum_j (j - (norb-1)/2) = ky (DiagHam's -y);
    pp = Haldane pseudopotentials, or kernel(q2) = V(q) e^{-q^2/2} directly (q = 0 kept)."""
    import numpy as np
    L = math.sqrt(2 * math.pi * norb * ratio)      # DiagHam: Length = sqrt(2 pi r N_orb), r its -r option
    if kernel is None:
        kernel = lambda q2: 4 * math.pi * math.exp(-q2 / 2) * sum(v * laguerre(m, q2) for m, v in enumerate(pp) if v)
    qx = np.linspace(-qmax, qmax, nq)
    dq = qx[1] - qx[0]
    cache = {}

    def A(j1, j2, j3, j4):
        key = (j1, j2, j3, j4)
        if key in cache:
            return cache[key]
        qy = 2 * math.pi * (j1 - j4) / L
        q2 = qx * qx + qy * qy
        f = np.array([kernel(x) for x in q2]) if not hasattr(kernel, "vectorized") else kernel(q2)
        v = np.sum(f * np.cos(qx * (2 * math.pi / L) * (j1 - j3))) * dq / (2 * math.pi * L)
        cache[key] = v
        return v

    gen = itertools.combinations if fermion else itertools.combinations_with_replacement
    shift = (norb - 1) / 2
    states = [c for c in gen(range(norb), n) if abs(sum(j - shift for j in c) - ky) < 1e-9]
    return _assemble(states, norb, A, fermion, lambda j1, j3, j4: j3 + j4 - j1)


# ------------------------------------------------------------------ disk
def disk_pair_state(m, M):
    """{(m1, m2): amplitude} of the normalised pair state (z1 - z2)^m (z1 + z2)^M in the orbital basis."""
    c = Counter()
    for a in range(m + 1):
        for b in range(M + 1):
            m1 = a + b
            m2 = m + M - m1
            c[(m1, m2)] += math.comb(m, a) * (-1) ** (m - a) * math.comb(M, b) * math.sqrt(2 ** m1 * math.factorial(m1) * 2 ** m2 * math.factorial(m2))
    norm = math.sqrt(sum(x * x for x in c.values()))
    return {k: x / norm for k, x in c.items() if x}


def disk_spectrum(n, lz, pp, fermion):
    """Every eigenvalue at total angular momentum lz for H = 2 sum_m V_m sum_{i<j} P_m(ij)."""
    norb = lz + 1
    pair = {}          # (m, total) -> amplitudes

    def V(j1, j2, j3, j4):
        tot = j1 + j2
        v = 0.0
        for m, Vm in enumerate(pp):
            if Vm == 0 or m > tot:
                continue
            key = (m, tot)
            if key not in pair:
                pair[key] = disk_pair_state(m, tot - m)
            p = pair[key]
            v += Vm * p.get((j1, j2), 0.0) * p.get((j3, j4), 0.0)
        return 2.0 * v

    gen = itertools.combinations if fermion else itertools.combinations_with_replacement
    states = [c for c in gen(range(norb), n) if sum(c) == lz]
    return _assemble(states, norb, V, fermion, lambda j1, j3, j4: j3 + j4 - j1)


def disk_coulomb_pseudopotentials(mmax):
    return [math.gamma(m + 0.5) / (2 * math.factorial(m)) for m in range(mmax + 1)]


# ------------------------------------------------------------------ exact first-quantised check (disk, three particles)
def disk_spectrum_exact3(lz, pp):
    """Every eigenvalue for three bosons at total angular momentum lz, computed in first quantisation
    with exact rational arithmetic: H = 2 sum_m V_m sum_{i<j} P_m(ij) where P_m(ij) keeps the
    (z_i - z_j)^m component of the pair (i, j) (the relative-angular-momentum decomposition is
    orthogonal, so no metric is needed); the basis is the monomial symmetric polynomials of degree
    lz. A second route to the same numbers as disk_spectrum, sharing nothing with it."""
    import numpy as np
    from fractions import Fraction
    from collections import defaultdict
    N = 3

    def add(p, q, c=1):
        out = defaultdict(Fraction, p)
        for k, v in q.items():
            out[k] += c * v
        return {k: v for k, v in out.items() if v != 0}

    def mul(p, q):
        out = defaultdict(Fraction)
        for k1, v1 in p.items():
            for k2, v2 in q.items():
                out[tuple(a + b for a, b in zip(k1, k2))] += v1 * v2
        return {k: v for k, v in out.items() if v != 0}

    def power(p, n, dim):
        r = {(0,) * dim: Fraction(1)}
        for _ in range(n):
            r = mul(r, p)
        return r

    def var(i):
        e = [0] * N
        e[i] = 1
        return {tuple(e): Fraction(1)}

    def project_pair(p, i, j, m):
        zi, zj = var(i), var(j)
        u, v = add(zi, zj), add(zi, zj, -1)
        half = Fraction(1, 2)
        U = {(0,) * N + (1, 0): Fraction(1)}
        Vv = {(0,) * N + (0, 1): Fraction(1)}
        zi_uv = {k: half * c for k, c in add(U, Vv).items()}
        zj_uv = {k: half * c for k, c in add(U, Vv, -1).items()}
        total = {}
        for k, c in p.items():
            rest = list(k)
            a, b = rest[i], rest[j]
            rest[i] = rest[j] = 0
            term = mul({tuple(rest) + (0, 0): c}, mul(power(zi_uv, a, N + 2), power(zj_uv, b, N + 2)))
            total = add(total, term)
        result = {}
        for k, c in total.items():
            if k[N + 1] != m:
                continue
            result = add(result, mul({k[:N]: c}, mul(power(u, k[N], N), power(v, m, N))))
        return result

    parts = [lam for lam in itertools.product(range(lz + 1), repeat=N) if sum(lam) == lz and list(lam) == sorted(lam, reverse=True)]
    basis = [{perm: Fraction(1) for perm in set(itertools.permutations(lam))} for lam in parts]
    M = np.zeros((len(basis), len(basis)))
    for col, b in enumerate(basis):
        hb = {}
        for m, Vm in enumerate(pp):
            if not Vm:
                continue
            for i, j in itertools.combinations(range(N), 2):
                hb = add(hb, project_pair(b, i, j, m), 2 * Fraction(Vm))
        for row, lead in enumerate(parts):
            M[row, col] = float(hb.get(tuple(lead), 0))
    return sorted(np.linalg.eigvals(M).real)


# ------------------------------------------------------------------ cases
# kind, name, parameters, note
CASES = [
    # cylinder: n, norb, pp, fermion, ratio, ky (this oracle's), note.  Program: -l norb-1, -r ratio, -y 2*ky --nbr-ky 1
    ("cyl", "fermions_n4_norb9_v1v3_r1_ky0", (4, 9, [0, 1, 0, 0.3], True, 1.0, 0), "V1 + V3, 9 orbitals, ratio 1, Ky = 0 (12 states)"),
    ("cyl", "fermions_n4_norb9_v1v3_r2_ky0", (4, 9, [0, 1, 0, 0.3], True, 2.0, 0), "ratio 2 (L = sqrt(2 pi N_orb r))"),
    ("cyl", "fermions_n4_norb9_v1v3_r0.5_ky0", (4, 9, [0, 1, 0, 0.3], True, 0.5, 0), "ratio 1/2"),
    ("cyl", "fermions_n4_norb9_v1v3_r1_ky2", (4, 9, [0, 1, 0, 0.3], True, 1.0, 2), "Ky = 2 (the program's -y 4: it labels sectors by twice the momentum)"),
    ("cyl", "bosons_n4_norb7_delta_r1_ky0", (4, 7, [1], False, 1.0, 0), "bosons, delta interaction = V0 (FQHECylinderBosonsDeltaInteraction), 18 states"),
    ("cyl", "bosons_n4_norb7_delta_r2_ky1", (4, 7, [1], False, 2.0, 1), "bosons delta, ratio 2, Ky = 1 (-y 2)"),
    ("cyl", "fermions_n4_norb9_laplacian_delta_r1_ky0", (4, 9, [0, 1], True, 1.0, 0), "fermions, Laplacian delta = V1 (FQHECylinderFermionsLaplacianDelta)"),
    # disk: n, lz, pp, fermion, note
    ("disk", "bosons_n3_lz6_v0", (3, 6, [1], False), "three bosons at Lz = 6 (the nu = 1/2 Laughlin point), V0 = 1: exact spectrum 0, 3/2, 15/8, 33/16, 9/4, 3, 6"),
    ("disk", "bosons_n4_lz12_v0_v2", (4, 12, [1, 0, 0.5], False), "four bosons at Lz = 12, V0 = 1, V2 = 1/2 (34 states)"),
    ("disk", "bosons_n3_lz9_v2", (3, 9, [0, 0, 1], False), "three bosons at Lz = 9, V2 = 1: four zero modes then 3/32, 3/4, 27/32, ..."),
]


def render_all():
    files = {}
    for kind, name, par, note in CASES:
        if kind == "cyl":
            n, norb, pp, fermion, ratio, ky = par
            e = cylinder_spectrum(n, norb, pp, fermion, ratio, ky)
            files[f"{name}_pp.dat"] = ("# generated by tests/oracles/geometry_ed.py; do not edit\nPseudopotentials = "
                                       + " ".join(f"{x:.16g}" for x in list(pp) + [0] * (norb - len(pp))) + "\n")
            head = f"# cylinder {'fermions' if fermion else 'bosons'} N={n} N_orb={norb} ratio={ratio} Ky={ky} (program -y {2 * ky}): {note}"
            files[f"{name}_spectrum.dat"] = head + f"; {len(e)} eigenvalues from tests/oracles/geometry_ed.py\n" + "".join(f"{2 * ky} {x:.15g}\n" for x in e)
        else:
            n, lz, pp, fermion = par
            e = disk_spectrum(n, lz, pp, fermion)
            files[f"{name}_pp.dat"] = ("# generated by tests/oracles/geometry_ed.py; do not edit\nPseudopotentials = "
                                       + " ".join(f"{x:.16g}" for x in list(pp) + [0] * (lz + 1 - len(pp))) + "\n")
            files[f"{name}_spectrum.dat"] = (f"# disk {'fermions' if fermion else 'bosons'} N={n} Lz={lz}: {note}; {len(e)} eigenvalues from tests/oracles/geometry_ed.py\n"
                                             + "".join(f"{lz} {x:.15g}\n" for x in e))
    v = disk_coulomb_pseudopotentials(8)
    files["disk_coulomb_2s8.txt"] = ("# Coulomb pseudopotentials on the disk, V_m = Gamma(m + 1/2) / (2 m!), units e^2/(eps l_B); tests/oracles/geometry_ed.py, do not edit\n"
                                     "Pseudopotentials = " + " ".join(f"{x:.15g}" for x in v) + "\n")
    return files


def selfcheck():
    """The two disk routes must agree (three-boson cases), a k-independent guard on the assembly."""
    worst = 0.0
    for kind, name, par, note in CASES:
        if kind == "disk" and par[0] == 3 and not par[3]:
            a = disk_spectrum(*par)
            b = disk_spectrum_exact3(par[1], par[2])
            worst = max(worst, max(abs(x - y) for x, y in zip(a, b)))
    return worst


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
        w = selfcheck()
        print(f"disk: second-quantised vs exact first-quantised (three bosons): max difference {w:.2g}")
        if w > 1e-9:
            bad += 1
        for name, text in files.items():
            p = Path(a.check) / name
            if not p.exists():
                bad += 1
                print(f"MISSING {p}")
                continue
            if name.endswith("_pp.dat") or name.endswith(".txt"):
                ok = p.read_text() == text
            else:
                av = [float(l.split()[-1]) for l in p.read_text().splitlines() if l and not l.startswith("#")]
                bv = [float(l.split()[-1]) for l in text.splitlines() if l and not l.startswith("#")]
                ok = len(av) == len(bv) and max(abs(x - y) for x, y in zip(av, bv)) <= a.tol
            if not ok:
                bad += 1
                print(f"MISMATCH {p}")
        print(f"{len(files)} files, {bad} problems")
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
