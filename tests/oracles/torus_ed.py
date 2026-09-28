#!/usr/bin/env python3
"""Independent exact diagonalisation in the lowest Landau level on a torus (numpy, no DiagHam).

    python3 tests/oracles/torus_ed.py --write tests/data/torus_ed
    python3 tests/oracles/torus_ed.py --check tests/data/torus_ed     (selftest.torus_ed_oracle, needs numpy)

Landau gauge, Nphi orbitals j = 0..Nphi-1 with k_j = 2 pi j / Ly, Lx Ly = 2 pi Nphi
(l_B = 1), aspect ratio r = Lx/Ly. Two-body Hamiltonian in the standard form
(Yoshioka 1984; Chakraborty & Pietilainen):

    H = 1/2 sum_{j1 j2 j3 j4} A(j1,j2,j3,j4) c+_{j1} c+_{j2} c_{j3} c_{j4},
    A = (1/(Lx Ly)) sum_{q != 0} V(q) exp(-q^2/2) exp(i q_x (X_{j1} - X_{j3})),
        q = 2 pi (s/Lx, t/Ly), t = j1 - j4 (mod Nphi), X_j = k_j.

Coulomb: V(q) = 2 pi / |q| (e^2/eps = 1), the q = 0 term dropped (neutralising
background). Haldane pseudopotentials V_m: V(q) exp(-q^2/2) is replaced by
4 pi sum_m V_m L_m(q^2) exp(-q^2/2) with the q = 0 term kept (it adds the constant
2 pi sum_m V_m N(N-1)/(2 pi Nphi), a "pair count" that the program keeps too).
The 4 pi is fixed by the two-particle problem: a pair of relative angular
momentum m then has energy 2 V_m, as in the sphere oracle (tests/oracles/sphere_ed.py).
Fermions and bosons; one Ky sector at a time; every eigenvalue is written.
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


def coulomb_kernel(q2):
    """V(q) exp(-q^2/2) for Coulomb"""
    return 2 * math.pi / math.sqrt(q2) * math.exp(-q2 / 2)


def pseudopotential_kernel(pp, scale=4 * math.pi):
    def kernel(q2):
        return scale * math.exp(-q2 / 2) * sum(v * laguerre(m, q2) for m, v in enumerate(pp) if v)
    return kernel


def spectrum(n, nphi, kernel, fermion=True, r=1.0, ky=0, keep_q0=False, smax=40, kmax=6):
    import numpy as np
    lxly = 2 * math.pi * nphi
    lx, ly = math.sqrt(lxly * r), math.sqrt(lxly / r)
    cache = {}

    def A(j1, j2, j3, j4):
        key = (j1, j2, j3, j4)
        if key in cache:
            return cache[key]
        tot = 0.0
        for k in range(-kmax, kmax + 1):
            t = (j1 - j4) + k * nphi
            for s in range(-smax, smax + 1):
                if s == 0 and t == 0 and not keep_q0:
                    continue
                qx, qy = 2 * math.pi * s / lx, 2 * math.pi * t / ly
                tot += kernel(qx * qx + qy * qy) * math.cos(qx * (2 * math.pi / ly) * (j1 - j3))
        v = tot / lxly
        cache[key] = v
        return v

    gen = itertools.combinations if fermion else itertools.combinations_with_replacement
    states = [c for c in gen(range(nphi), n) if sum(c) % nphi == ky]
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for s in states:
        c0 = Counter(s)
        for j3 in c0:
            for j4 in c0:
                c = Counter(c0)
                if fermion:
                    if j3 == j4:
                        continue
                    l = sorted(c.elements())
                    pre = (-1) ** l.index(j4); l.remove(j4)          # c_4 acts first
                    pre *= (-1) ** l.index(j3); l.remove(j3)
                else:
                    pre = math.sqrt(c[j4]); c[j4] -= 1
                    if c[j3] <= 0:
                        continue
                    pre *= math.sqrt(c[j3]); c[j3] -= 1
                    l = sorted(c.elements())
                for j1 in range(nphi):
                    j2 = (j3 + j4 - j1) % nphi
                    amp = pre * A(j1, j2, j3, j4)
                    if amp == 0:
                        continue
                    if fermion:
                        if j2 in l or j1 == j2:
                            continue
                        l2 = sorted(l + [j2]); amp *= (-1) ** l2.index(j2)
                        if j1 in l2:
                            continue
                        l1 = sorted(l2 + [j1]); amp *= (-1) ** l1.index(j1)
                        key = tuple(l1)
                    else:
                        cc = Counter(l)
                        amp *= math.sqrt(cc[j2] + 1); cc[j2] += 1
                        amp *= math.sqrt(cc[j1] + 1); cc[j1] += 1
                        key = tuple(sorted(cc.elements()))
                    H[idx[key], idx[s]] += 0.5 * amp
    return np.linalg.eigvalsh(H)


CASES = [
    # name, N, Nphi, fermion, ratio, ky, interaction ("coulomb" or pseudopotential list), note
    ("fermions_n3_nphi9_coulomb_ky0", 3, 9, True, 1.0, 0, "coulomb", "Coulomb, 84 states over Ky, the Ky=0 sector"),
    ("fermions_n4_nphi8_coulomb_ky0", 4, 8, True, 1.0, 0, "coulomb", "Coulomb, N=4"),
    ("fermions_n3_nphi9_coulomb_r2_ky0", 3, 9, True, 2.0, 0, "coulomb", "Coulomb, aspect ratio 2"),
    ("fermions_n4_nphi12_coulomb_ky0", 4, 12, True, 1.0, 0, "coulomb", "Coulomb at nu=1/3, N=4"),
    ("bosons_n3_nphi6_coulomb_ky0", 3, 6, False, 1.0, 0, "coulomb", "Coulomb, bosons"),
    ("fermions_n3_nphi9_v1_ky0", 3, 9, True, 1.0, 0, [0, 1], "V1 pseudopotential: Laughlin zero mode plus its spectrum"),
    ("fermions_n4_nphi12_v1v3_ky0", 4, 12, True, 1.0, 0, [0, 1, 0, 0.3], "V1 + V3"),
    ("bosons_n4_nphi8_v0v2_ky0", 4, 8, False, 1.0, 0, [1, 0, 0.5], "bosons, V0 + V2"),
]


def render_all(only=None):
    files = {}
    for name, n, nphi, fermion, r, ky, inter, note in CASES:
        if only and name != only:
            continue
        kernel = coulomb_kernel if inter == "coulomb" else pseudopotential_kernel(inter)
        e = spectrum(n, nphi, kernel, fermion, r, ky, keep_q0=(inter != "coulomb"))
        if inter != "coulomb":
            files[f"{name}_pp.dat"] = ("# generated by tests/oracles/torus_ed.py; do not edit\nPseudopotentials = " + " ".join(f"{x:.16g}" for x in inter) + "\n")
        files[f"{name}_spectrum.dat"] = (f"# {'fermions' if fermion else 'bosons'} N={n} Nphi={nphi} ratio={r} Ky={ky}: {note}; "
                                         f"{len(e)} eigenvalues from the independent ED in tests/oracles/torus_ed.py\n"
                                         + "".join(f"{ky} {x:.15g}\n" for x in e))
    return files


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", metavar="DIR")
    ap.add_argument("--check", metavar="DIR")
    ap.add_argument("--only")
    ap.add_argument("--tol", type=float, default=1e-9)
    a = ap.parse_args()
    files = render_all(a.only)
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
                bad += 1; print(f"MISSING {p}"); continue
            if name.endswith("_pp.dat"):
                ok = p.read_text() == text
            else:
                av = [float(l.split()[1]) for l in p.read_text().splitlines() if l and not l.startswith("#")]
                bv = [float(l.split()[1]) for l in text.splitlines() if l and not l.startswith("#")]
                ok = len(av) == len(bv) and max(abs(x - y) for x, y in zip(av, bv)) <= a.tol
            if not ok:
                bad += 1; print(f"MISMATCH {p}")
        print(f"{len(files)} files, {bad} mismatches")
        sys.exit(1 if bad else 0)
    if not a.write and not a.check:
        for name, text in files.items():
            print(f"== {name}\n" + "\n".join(text.splitlines()[:5]))


if __name__ == "__main__":
    main()
