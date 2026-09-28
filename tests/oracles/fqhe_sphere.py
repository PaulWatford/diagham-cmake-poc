#!/usr/bin/env python3
"""Independent answers for the FQHE-on-the-sphere goldens (no DiagHam, no numpy).

    python3 tests/oracles/fqhe_sphere.py --write tests/data/fqhe_sphere
    python3 tests/oracles/fqhe_sphere.py --check tests/data/fqhe_sphere   (selftest.fqhe_sphere_oracle)

1. Zero-mode counting. The zero-energy states of the (k+1)-body hard-core
   interaction (k=1: the V_0 / V_1 pseudopotential, Laughlin; k=2:
   Moore-Read; k=3: Read-Rezayi) on N particles in 2S+1 orbitals are in
   one-to-one correspondence with the (k, r)-admissible occupation
   configurations -- at most k particles in r consecutive orbitals, r = 2
   for bosons, r = 3 for the fermionic Laughlin state, r = 4 for the
   fermionic Moore-Read state (Read & Rezayi 1996; Bernevig & Haldane
   2008, the Jack construction). Counting admissible configurations per
   Lz is a finite enumeration; this file does it.
2. Coulomb pseudopotentials in the lowest Landau level on the sphere
   of radius R = sqrt(S) l_B: 1/|r1 - r2| = (1/R) sum_k P_k(cos theta12)
   and <(S S) L| P_k |(S S) L> = (-1)^(2S+L) (2S+1)^2 {S S L; S S k}
   (S k S; -S 0 S)^2 for monopole harmonics (Fano, Ortolani & Colombo
   1986), with V_m = V_(L = 2S - m). Wigner 3j and 6j symbols are computed
   here with the Racah formulas in exact integer arithmetic.
3. Particle entanglement spectrum of the Laughlin state: for N_A
   particles the rank of the reduced density matrix in each Lz_A sector
   equals the number of (1,3)-admissible configurations of N_A particles
   in the same 2S+1 orbitals (Sterdyniak, Regnault & Bernevig 2011).
"""
import argparse, sys
from fractions import Fraction
from math import factorial as fact, sqrt
from pathlib import Path


# ---- 1. admissible configurations -------------------------------------------------
def admissible_lz_counts(n, norb, k, r, fermion=False):
    """{2Lz: number of (k,r)-admissible occupations of n particles on norb orbitals}
    (fermions: at most one particle per orbital as well)."""
    out = {}
    max_occ = 1 if fermion else k

    def rec(j, left, occ):
        if j == norb:
            if left == 0:
                two_lz = sum(o * (2 * i - (norb - 1)) for i, o in enumerate(occ))
                out[two_lz] = out.get(two_lz, 0) + 1
            return
        for m in range(0, min(max_occ, left) + 1):
            if sum(occ[max(0, j - r + 1):j]) + m <= k:
                rec(j + 1, left - m, occ + [m])
    rec(0, n, [])
    return out


# ---- 2. Wigner symbols (doubled integer arguments) and Coulomb pseudopotentials ----
def _tri(a, b, c):
    return Fraction(fact((a + b - c) // 2) * fact((a - b + c) // 2) * fact((-a + b + c) // 2), fact((a + b + c) // 2 + 1))


def wigner_3j(j1, j2, j3, m1, m2, m3):
    if m1 + m2 + m3 != 0 or (j1 + j2 + j3) % 2 or j3 < abs(j1 - j2) or j3 > j1 + j2:
        return 0.0
    for j, m in ((j1, m1), (j2, m2), (j3, m3)):
        if abs(m) > j or (j - m) % 2:
            return 0.0
    pre = _tri(j1, j2, j3) * Fraction(fact((j1 + m1) // 2) * fact((j1 - m1) // 2) * fact((j2 + m2) // 2)
                                      * fact((j2 - m2) // 2) * fact((j3 + m3) // 2) * fact((j3 - m3) // 2))
    s = Fraction(0)
    kmin = max(0, (j2 - j3 - m1) // 2, (j1 - j3 + m2) // 2)
    kmax = min((j1 + j2 - j3) // 2, (j1 - m1) // 2, (j2 + m2) // 2)
    for k in range(kmin, kmax + 1):
        d = (fact(k) * fact((j1 + j2 - j3) // 2 - k) * fact((j1 - m1) // 2 - k) * fact((j2 + m2) // 2 - k)
             * fact((j3 - j2 + m1) // 2 + k) * fact((j3 - j1 - m2) // 2 + k))
        s += Fraction((-1) ** k, d)
    return (-1) ** ((j1 - j2 - m3) // 2) * sqrt(float(pre)) * float(s)


def wigner_6j(a, b, c, d, e, g):
    def ok(x, y, z):
        return (x + y + z) % 2 == 0 and abs(x - y) <= z <= x + y
    if not (ok(a, b, c) and ok(a, e, g) and ok(d, b, g) and ok(d, e, c)):
        return 0.0
    pre = sqrt(float(_tri(a, b, c) * _tri(a, e, g) * _tri(d, b, g) * _tri(d, e, c)))
    s = Fraction(0)
    zmin = max((a + b + c) // 2, (a + e + g) // 2, (d + b + g) // 2, (d + e + c) // 2)
    zmax = min((a + b + d + e) // 2, (a + c + d + g) // 2, (b + c + e + g) // 2)
    for z in range(zmin, zmax + 1):
        den = (fact(z - (a + b + c) // 2) * fact(z - (a + e + g) // 2) * fact(z - (d + b + g) // 2)
               * fact(z - (d + e + c) // 2) * fact((a + b + d + e) // 2 - z) * fact((a + c + d + g) // 2 - z)
               * fact((b + c + e + g) // 2 - z))
        s += Fraction((-1) ** z * fact(z + 1), den)
    return pre * float(s)


def coulomb_pseudopotentials(two_s, landau_level=0):
    """[V_0, V_1, ..., V_2l] in units of e^2/(epsilon l_B) for the Landau level n on a sphere of
    monopole strength Q = S: the orbitals are the monopole harmonics of angular momentum l = Q + n
    (2l + 1 of them), R = sqrt(Q) l_B, and V_m = <l l; L = 2l - m | 1/(R |r1 - r2|) | l l; L> from the
    multipole expansion of 1/|r1 - r2| with the monopole-harmonic Gaunt coefficients
    (Fano, Ortolani and Colombo, PRB 34, 2670 (1986)):
        V_L = (1/R) sum_k (-1)^(2l + L) (2l + 1)^2 {l l L; l l k} (l k l; -Q 0 Q)^2.
    n = 0 gives the lowest-Landau-level formula (l = Q)."""
    two_l = two_s + 2 * landau_level
    R = sqrt(two_s / 2)
    out = []
    for m in range(two_l + 1):
        L2 = 2 * (two_l - m)
        v = 0.0
        for k in range(two_l + 1):
            v += ((-1) ** ((2 * two_l + L2) // 2) * (two_l + 1) ** 2
                  * wigner_6j(two_l, two_l, L2, two_l, two_l, 2 * k) * wigner_3j(two_l, 2 * k, two_l, -two_s, 0, two_s) ** 2)
        out.append(v / R)
    return out


COULOMB_HIGHER_LL = [(8, 1), (6, 2), (10, 1)]   # (2S, Landau level): CoulombPseudopotentials --landau-level n


# ---- 3. the cases -------------------------------------------------------------------
ZERO_MODES = [
    # name, N, 2S, k, r, 2Lz, note
    ("bosons_laughlin_1_2_n6_2s10", 6, 10, 1, 2, 0, "bosonic Laughlin nu=1/2 at 2S=2(N-1): unique V_0 zero mode"),
    ("bosons_laughlin_1_2_n6_2s12", 6, 12, 1, 2, 0, "two extra flux = two quasiholes of angular momentum N/2: L=6,4,2,0"),
    ("bosons_laughlin_1_2_n7_2s13", 7, 13, 1, 2, 1, "one quasihole (one extra flux): a single L=N/2 multiplet, one state at 2Lz=1"),
    ("bosons_mooreread_n6_2s4", 6, 4, 2, 2, 0, "bosonic Moore-Read nu=1 at 2S=N-2: unique three-body hard-core zero mode"),
    ("bosons_mooreread_n6_2s5", 6, 5, 2, 2, 0, "one extra flux = two half-charge quasiholes: L=3,1 multiplets"),
    ("bosons_mooreread_n8_2s6", 8, 6, 2, 2, 0, "bosonic Moore-Read at N=8"),
    ("bosons_readrezayi3_n6_2s3", 6, 3, 3, 2, 0, "bosonic Read-Rezayi k=3 (four-body hard-core), one extra flux"),
    ("bosons_readrezayi3_n9_2s4", 9, 4, 3, 2, 0, "bosonic Read-Rezayi k=3 at 2S=2(N-3)/3: unique zero mode"),
    ("fermions_laughlin_1_3_n6_2s18", 6, 18, 1, 3, 0, "three extra flux = three quasiholes of the nu=1/3 Laughlin state"),
    ("fermions_mooreread_n6_2s9", 6, 9, 2, 4, 0, "fermionic Moore-Read nu=1/2 at 2S=2N-3: unique three-body hard-core zero mode"),
    ("fermions_mooreread_n6_2s10", 6, 10, 2, 4, 0, "fermionic Moore-Read, one extra flux: same counting as the bosonic case at 2S=5 (L=3,1)"),
]
COULOMB_2S = (6, 10, 15)
PES_LAUGHLIN = (6, 15)   # N, 2S of the Laughlin state whose PES is counted
PES_CASES = [
    # name, N_A, 2Lz_A (None = all sectors)
    ("pes_laughlin_n6_2s15_na2_all", 2, None),
    ("pes_laughlin_n6_2s15_na3_all", 3, None),
    ("pes_laughlin_n6_2s15_na3_2lz-1", 3, -1),
]


def render_all():
    files = {}
    lines = ["# zero-mode counts: name  2Lz  count  -- (k,r)-admissible configurations of N particles on 2S+1 orbitals",
             "# generated by tests/oracles/fqhe_sphere.py; do not edit"]
    for name, n, two_s, k, r, two_lz, note in ZERO_MODES:
        counts = admissible_lz_counts(n, two_s + 1, k, r, fermion=name.startswith("fermions"))
        lines.append(f"{name} {two_lz} {counts.get(two_lz, 0)}  # N={n} 2S={two_s} ({k},{r})-admissible, total {sum(counts.values())}: {note}")
    files["zero_modes.txt"] = "\n".join(lines) + "\n"
    for two_s in COULOMB_2S:
        v = coulomb_pseudopotentials(two_s)
        files[f"coulomb_2s{two_s}.txt"] = (f"# lowest-Landau-level Coulomb pseudopotentials on the sphere, 2S={two_s}, units e^2/(eps l_B), R=sqrt(S)\n"
                                           f"# from Wigner 3j/6j algebra (tests/oracles/fqhe_sphere.py); do not edit\n"
                                           "Pseudopotentials = " + " ".join(f"{x:.15g}" for x in v) + "\n")
    for two_s, ll in COULOMB_HIGHER_LL:
        v = coulomb_pseudopotentials(two_s, ll)
        files[f"coulomb_2s{two_s}_ll{ll}.txt"] = (f"# Coulomb pseudopotentials on the sphere in Landau level n={ll}, 2S={two_s} (l = S + n, {two_s + 2 * ll + 1} values), units e^2/(eps l_B), R=sqrt(S)\n"
                                                 f"# from Wigner 3j/6j algebra with the monopole-harmonic Gaunt coefficients (tests/oracles/fqhe_sphere.py); do not edit\n"
                                                 "Pseudopotentials = " + " ".join(f"{x:.15g}" for x in v) + "\n")
    n, two_s = PES_LAUGHLIN
    lines = [f"# particle entanglement spectrum of the nu=1/3 Laughlin state, N={n}, 2S={two_s}: name  N_A  2Lz_A|all  nonzero levels",
             "# = (1,3)-admissible configurations of N_A particles on 2S+1 orbitals; generated by tests/oracles/fqhe_sphere.py"]
    for name, na, two_lz in PES_CASES:
        counts = admissible_lz_counts(na, two_s + 1, 1, 3, fermion=True)
        val = sum(counts.values()) if two_lz is None else counts.get(two_lz, 0)
        lines.append(f"{name} {na} {'all' if two_lz is None else two_lz} {val}")
    files["pes_laughlin_n6_2s15.txt"] = "\n".join(lines) + "\n"
    return files


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", metavar="DIR")
    ap.add_argument("--check", metavar="DIR")
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
            if not p.exists() or p.read_text() != text:
                bad += 1
                print(f"MISMATCH {p}")
        print(f"{len(files)} files, {bad} mismatches")
        sys.exit(1 if bad else 0)
    if not a.write and not a.check:
        for name, text in files.items():
            print(f"== {name}\n{text}")


if __name__ == "__main__":
    main()
