#!/usr/bin/env python3
"""Independent exact diagonalisation of pseudopotential Hamiltonians in the lowest Landau level on the sphere.

    python3 tests/oracles/sphere_ed.py --write tests/data/sphere_ed
    python3 tests/oracles/sphere_ed.py --check tests/data/sphere_ed      (selftest.sphere_ed_oracle, needs numpy)

For N fermions or bosons in the 2S+1 orbitals of the lowest Landau level
with a two-body interaction given by its Haldane pseudopotentials V_m
(m = 2S - L, L the pair angular momentum), the Hamiltonian is built here
from scratch: pair projectors from Clebsch-Gordan coefficients (Wigner 3j
symbols computed exactly in tests/oracles/fqhe_sphere.py), second
quantisation on an explicit Fock basis of one Lz sector, dense
diagonalisation with numpy. Nothing is shared with DiagHam. DiagHam's
normalisation is H = 2 sum_{i<j} sum_m V_m P_m(ij) (a two-particle state of
relative angular momentum m has energy 2 V_m); that factor 2 is the only
convention taken from the program, established with the two-particle case.

Each case writes two files: <case>_pp.dat, the pseudopotential file the
program is given (so both sides use identical numbers -- for the Coulomb
cases those numbers come from this repository's own closed form, not from
CoulombPseudopotentials), and <case>_spectrum.dat, every eigenvalue of the
Lz sector as "Lz E" lines. The tests compare the program's full spectrum
of that sector with it (check_spectrum spectrum ..., 1e-9).
"""
import argparse, itertools, math, sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fqhe_sphere as fs  # noqa: E402


def cg(j1, m1, j2, m2, J, M):
    """Clebsch-Gordan <j1 m1 j2 m2 | J M>, doubled arguments."""
    if m1 + m2 != M:
        return 0.0
    return (-1) ** ((j1 - j2 + M) // 2) * math.sqrt(J + 1) * fs.wigner_3j(j1, j2, J, m1, m2, -M)


def spectrum(n, two_s, pp, fermion, two_lz=0, scale=2.0):
    import numpy as np
    norb = two_s + 1
    gen = itertools.combinations if fermion else itertools.combinations_with_replacement
    states = [c for c in gen(range(norb), n) if sum(2 * j - two_s for j in c) == two_lz]
    idx = {s: i for i, s in enumerate(states)}
    cache = {}

    def V(j1, j2, j3, j4):     # <j1 j2| V |j3 j4> = sum_m V_m <j1 j2|L M><j3 j4|L M>, L = 2S - m
        k = (j1, j2, j3, j4)
        if k in cache:
            return cache[k]
        v = 0.0
        M = (2 * j1 - two_s) + (2 * j2 - two_s)
        for m, Vm in enumerate(pp):
            J = 2 * two_s - 2 * m
            if abs(M) > J or Vm == 0:
                continue
            v += Vm * cg(two_s, 2 * j1 - two_s, two_s, 2 * j2 - two_s, J, M) * cg(two_s, 2 * j3 - two_s, two_s, 2 * j4 - two_s, J, M)
        cache[k] = v
        return v

    H = np.zeros((len(states), len(states)))
    for s in states:                      # H = 1/2 sum <12|V|34> a+_1 a+_2 a_4 a_3 (a_3 acts first)
        c0 = Counter(s)
        for j3 in c0:
            for j4 in c0:
                c = Counter(c0)
                if fermion:
                    if j3 == j4:
                        continue
                    l = sorted(c.elements())
                    pre = (-1) ** l.index(j3)
                    l.remove(j3)
                    pre *= (-1) ** l.index(j4)
                    l.remove(j4)
                else:
                    pre = math.sqrt(c[j3])
                    c[j3] -= 1
                    if c[j4] <= 0:
                        continue
                    pre *= math.sqrt(c[j4])
                    c[j4] -= 1
                    l = sorted(c.elements())
                M = (2 * j3 - two_s) + (2 * j4 - two_s)
                if (M + 2 * two_s) % 2:
                    continue
                for j1 in range(norb):
                    j2 = (M + 2 * two_s) // 2 - j1
                    if j2 < 0 or j2 >= norb:
                        continue
                    amp = pre * V(j1, j2, j3, j4)
                    if amp == 0:
                        continue
                    if fermion:
                        if j2 in l:
                            continue
                        l2 = sorted(l + [j2])
                        amp *= (-1) ** l2.index(j2)
                        if j1 in l2:
                            continue
                        l1 = sorted(l2 + [j1])
                        amp *= (-1) ** l1.index(j1)
                        key = tuple(l1)
                    else:
                        cc = Counter(l)
                        amp *= math.sqrt(cc[j2] + 1)
                        cc[j2] += 1
                        amp *= math.sqrt(cc[j1] + 1)
                        cc[j1] += 1
                        key = tuple(sorted(cc.elements()))
                    H[idx[key], idx[s]] += 0.5 * scale * amp
    return np.linalg.eigvalsh(H)


def coulomb(two_s):
    return fs.coulomb_pseudopotentials(two_s)


CASES = [
    # name, N, 2S, fermion, 2Lz, pseudopotentials, note
    ("fermions_n4_2s9_coulomb", 4, 9, True, 0, lambda: coulomb(9), "Coulomb (closed form from fqhe_sphere.py), Lz=0, 18 states"),
    ("fermions_n5_2s12_coulomb", 5, 12, True, 0, lambda: coulomb(12), "Coulomb, Lz=0"),
    ("fermions_n4_2s9_v1_v3", 4, 9, True, 0, lambda: [0, 1, 0, 0.3, 0, 0.1], "generic odd pseudopotentials V1=1, V3=0.3, V5=0.1, Lz=0"),
    ("fermions_n4_2s11_coulomb_lz2", 4, 11, True, 4, lambda: coulomb(11), "Coulomb, 2Lz=4 sector"),
    ("bosons_n4_2s6_coulomb", 4, 6, False, 0, lambda: coulomb(6), "Coulomb, bosons, Lz=0, 18 states"),
    ("bosons_n5_2s6_v0_v2", 5, 6, False, 0, lambda: [1, 0, 0.5, 0, 0.25], "generic even pseudopotentials V0=1, V2=0.5, V4=0.25"),
    ("bosons_n6_2s8_coulomb", 6, 8, False, 0, lambda: coulomb(8), "Coulomb, bosons N=6, Lz=0"),
]


def render_all():
    files = {}
    for name, n, two_s, fermion, two_lz, pp_f, note in CASES:
        pp = pp_f()
        files[f"{name}_pp.dat"] = (f"# {note}; generated by tests/oracles/sphere_ed.py, do not edit\n"
                                   "Pseudopotentials = " + " ".join(f"{x:.16g}" for x in pp) + "\n")
        e = spectrum(n, two_s, pp, fermion, two_lz)
        files[f"{name}_spectrum.dat"] = (f"# {'fermions' if fermion else 'bosons'} N={n} 2S={two_s} 2Lz={two_lz}: {note}; "
                                         f"{len(e)} eigenvalues from the independent ED in tests/oracles/sphere_ed.py\n"
                                         + "".join(f"{two_lz // 2 if two_lz % 2 == 0 else two_lz / 2} {x:.15g}\n" for x in e))
    return files


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", metavar="DIR")
    ap.add_argument("--check", metavar="DIR")
    ap.add_argument("--tol", type=float, default=1e-9, help="tolerance for --check on the spectra")
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
            if name.endswith("_pp.dat"):
                ok = p.read_text() == text
            else:   # spectra: numeric comparison (floating point may differ in the last digits between numpy builds)
                a_vals = [float(l.split()[1]) for l in p.read_text().splitlines() if l and not l.startswith("#")]
                b_vals = [float(l.split()[1]) for l in text.splitlines() if l and not l.startswith("#")]
                ok = len(a_vals) == len(b_vals) and max(abs(x - y) for x, y in zip(a_vals, b_vals)) <= a.tol
            if not ok:
                bad += 1
                print(f"MISMATCH {p}")
        print(f"{len(files)} files, {bad} mismatches")
        sys.exit(1 if bad else 0)
    if not a.write and not a.check:
        for name, text in files.items():
            print(f"== {name}\n" + "\n".join(text.splitlines()[:4]))


if __name__ == "__main__":
    main()
