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

spectrum_species generalises this to several species (SU(2), SU(3), SU(4) spin;
layers with inter-layer Coulomb 2 pi exp(-q d)/q), to Landau level n (form
factor [L_n(q^2/2)]^2 e^{-q^2/2}) and to the magnetic-translation sectors Kx
(eigenvalues of the translation of every particle by Nphi/gcd(N, Nphi)
orbitals inside a Ky sector). A pair of different species in relative
angular momentum m costs 2 V_m, like a same-species pair.
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


def landau_level_kernel(kernel, n):
    """The same interaction in Landau level n: the form factor exp(-q^2/2) becomes [L_n(q^2/2)]^2 exp(-q^2/2)."""
    if n == 0:
        return kernel
    return lambda q2: kernel(q2) * laguerre(n, q2 / 2) ** 2


def bilayer_coulomb_kernel(d):
    """Inter-layer Coulomb at separation d (magnetic lengths): V(q) = 2 pi exp(-q d) / q, with the form factor."""
    return lambda q2: 2 * math.pi / math.sqrt(q2) * math.exp(-math.sqrt(q2) * d) * math.exp(-q2 / 2)


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


def spectrum_species(numbers, nphi, kernels, fermion=True, r=1.0, ky=0, keep_q0=None, kx=None, smax=40, kmax=6):
    """k species (SU(k) spin, layers) with numbers[i] particles of species i in the same Nphi orbitals;
    kernels[(i, j)] (i <= j) is V(q) exp(-q^2/2) (times the Landau-level form factor) of an (i, j)
    pair, keep_q0[(i, j)] whether its q = 0 term is kept (pseudopotentials) or dropped (Coulomb).
    H = 1/2 sum_i A^{ii} c+c+cc + sum_{i<j} A^{ij} c+_i c+_j c_j c_i, the same matrix element as
    spectrum() for every pair. One Ky sector; with kx given, the sector of the magnetic translation
    T: j -> j + 1 (all particles) restricted to the Ky sector (T^d, d = Nphi / gcd(N, Nphi)), whose
    eigenvalues exp(2 pi i kx / g) with g = gcd(N, Nphi) label the Kx sectors; the labelling
    (which phase is called kx = 1) follows the program's, fixed once on the N = 4, Nphi = 8 case."""
    import numpy as np
    keep_q0 = keep_q0 or {}
    lxly = 2 * math.pi * nphi
    lx, ly = math.sqrt(lxly * r), math.sqrt(lxly / r)
    cache = {}

    def A(key, j1, j2, j3, j4):
        ck = (key, j1, j2, j3, j4)
        if ck in cache:
            return cache[ck]
        kernel = kernels[key]
        q0 = keep_q0.get(key, False)
        tot = 0.0
        for k in range(-kmax, kmax + 1):
            t = (j1 - j4) + k * nphi
            for s in range(-smax, smax + 1):
                if s == 0 and t == 0 and not q0:
                    continue
                qx, qy = 2 * math.pi * s / lx, 2 * math.pi * t / ly
                tot += kernel(qx * qx + qy * qy) * math.cos(qx * (2 * math.pi / ly) * (j1 - j3))
        v = tot / lxly
        cache[ck] = v
        return v

    k = len(numbers)
    gen = itertools.combinations if fermion else itertools.combinations_with_replacement
    blocks = [list(gen(range(nphi), n)) for n in numbers]
    states = [s for s in itertools.product(*blocks) if sum(sum(b) for b in s) % nphi == ky]
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for s in states:
        col = idx[s]
        for i in range(k):
            for j in range(i, k):
                if (i, j) not in kernels:
                    continue
                if i == j:
                    block = s[i]          # 1/2 sum A(1234) c+_1 c+_2 c_3 c_4, c_4 acting first
                    for j4 in set(block):
                        r4 = _remove(block, j4, fermion)
                        for j3 in set(r4[1]):
                            r3 = _remove(r4[1], j3, fermion)
                            for j1 in range(nphi):
                                j2 = (j3 + j4 - j1) % nphi
                                v = A((i, i), j1, j2, j3, j4)
                                if v == 0:
                                    continue
                                a2 = _add(r3[1], j2, fermion)
                                if a2 is None:
                                    continue
                                a1 = _add(a2[1], j1, fermion)
                                if a1 is None:
                                    continue
                                new = s[:i] + (a1[1],) + s[i + 1:]
                                H[idx[new], col] += 0.5 * r4[0] * r3[0] * a2[0] * a1[0] * v
                else:
                    # sum A(1234) c+_{1 i} c+_{2 j} c_{3 j} c_{4 i}: the momentum transfer t = j1 - j4 is
                    # within species i, the pair (2, 3) within species j; c_{4 i} acts first
                    for j4 in set(s[i]):
                        r4 = _remove(s[i], j4, fermion)
                        for j3 in set(s[j]):
                            r3 = _remove(s[j], j3, fermion)
                            for j1 in range(nphi):
                                j2 = (j3 + j4 - j1) % nphi
                                v = A((i, j), j1, j2, j3, j4)
                                if v == 0:
                                    continue
                                a2 = _add(r3[1], j2, fermion)
                                if a2 is None:
                                    continue
                                a1 = _add(r4[1], j1, fermion)
                                if a1 is None:
                                    continue
                                new = list(s)
                                new[i] = a1[1]
                                new[j] = a2[1]
                                H[idx[tuple(new)], col] += r4[0] * r3[0] * a2[0] * a1[0] * v
    if kx is None:
        return np.linalg.eigvalsh(H)
    # magnetic translation T^d within the Ky sector
    n = sum(numbers)
    g = math.gcd(n, nphi)
    d = nphi // g
    Td = np.zeros((len(states), len(states)), dtype=complex)
    for s, i in idx.items():
        new = []
        sgn = 1
        for b in s:
            shifted = [(j + d) % nphi for j in b]
            if fermion:
                # sign of the permutation that sorts the shifted block
                perm = sorted(range(len(shifted)), key=lambda t: shifted[t])
                inv = sum(1 for a in range(len(perm)) for c in range(a + 1, len(perm)) if perm[a] > perm[c])
                sgn *= (-1) ** inv
            new.append(tuple(sorted(shifted)))
        Td[idx[tuple(new)], i] = sgn
    w, v = np.linalg.eig(Td)
    phase = np.angle(w) / (2 * math.pi) * g            # eigenvalue exp(2 pi i m / g) -> m (up to a global offset)
    m = np.round(phase) % g
    keep = v[:, np.isclose(m, kx % g)]
    if keep.shape[1] == 0:
        return np.zeros(0)
    q, _ = np.linalg.qr(keep)
    return np.sort(np.linalg.eigvalsh(q.conj().T @ H @ q).real)


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


SU2_KEYS = {(0, 0): "UpUp", (1, 1): "DownDown", (0, 1): "UpDown"}
SU3_KEYS = {(0, 0): "11", (0, 1): "12", (0, 2): "13", (1, 1): "22", (1, 2): "23", (2, 2): "33"}
_SU4 = ["UpPlus", "UpMinus", "DownPlus", "DownMinus"]
SU4_KEYS = {(i, j): _SU4[i] + _SU4[j] for i in range(4) for j in range(i, 4)}
ONE_KEY = {(0, 0): ""}


def _kernel(spec):
    """(kernel, keep_q0) from a pair specification."""
    if spec[0] == "pp":
        return pseudopotential_kernel(spec[1]), True
    if spec[0] == "coulomb":
        return coulomb_kernel, False
    if spec[0] == "coulomb_ll":
        return landau_level_kernel(coulomb_kernel, spec[1]), False
    if spec[0] == "bilayer":
        return bilayer_coulomb_kernel(spec[1]), False
    raise ValueError(spec)


F8 = {(0, 0): ("pp", [0, 1, 0, 0.3]), (1, 1): ("pp", [0, 1, 0, 0.3]), (0, 1): ("pp", [1, 0.5, 0, 0.2])}
B6 = {(0, 0): ("pp", [1, 0, 0.2]), (1, 1): ("pp", [1, 0, 0.2]), (0, 1): ("pp", [0.7, 0, 0.3])}
S3 = {(0, 0): ("pp", [1, 0, 0.2]), (1, 1): ("pp", [1, 0, 0]), (2, 2): ("pp", [0.5, 0, 0.1]),
      (0, 1): ("pp", [0.7, 0, 0.3]), (0, 2): ("pp", [0.6, 0, 0]), (1, 2): ("pp", [0.4, 0, 0.1])}
_S4V = [[1, 0, 0], [1, 0, 0.1], [0.8, 0, 0], [0.6, 0, 0.2], [0.7, 0, 0.3], [0.5, 0, 0], [0.4, 0, 0.1], [0.3, 0, 0], [0.2, 0, 0.2], [0.9, 0, 0]]
S4 = {p: ("pp", v) for p, v in zip(sorted(SU4_KEYS), _S4V)}
V1V3 = {(0, 0): ("pp", [0, 1, 0, 0.3])}
BV0V2 = {(0, 0): ("pp", [1, 0, 0.3])}
C2 = {(0, 0): ("coulomb",), (1, 1): ("coulomb",), (0, 1): ("coulomb",)}


def _bilayer(d):
    return {(0, 0): ("coulomb",), (1, 1): ("coulomb",), (0, 1): ("bilayer", d)}


# name, numbers (per species), Nphi, fermion, ratio, ky (or "all"), kx (or None), pairs, keys, note
SPECIES_CASES = [
    ("fermions_n3_nphi9_coulomb_ll1_ky0", [3], 9, True, 1.0, 0, None, {(0, 0): ("coulomb_ll", 1)}, ONE_KEY, "Coulomb in the first excited Landau level"),
    ("fermions_n4_nphi8_v1v3_allky", [4], 8, True, 1.0, "all", None, V1V3, ONE_KEY, "V1 + V3, every Ky sector (the AllMomenta program writes them unresolved)"),
    ("fermions_n4_nphi8_v1v3_kx2_ky0", [4], 8, True, 1.0, 0, 2, V1V3, ONE_KEY, "V1 + V3, magnetic-translation sector (Kx, Ky) = (2, 0)"),
    ("fermions_n4_nphi8_v1v3_kx1_ky2", [4], 8, True, 1.0, 2, 1, V1V3, ONE_KEY, "V1 + V3, (Kx, Ky) = (1, 2)"),
    ("bosons_n4_nphi6_v0v2_kx0_ky0", [4], 6, False, 1.0, 0, 0, BV0V2, ONE_KEY, "bosons V0 + 0.3 V2, (Kx, Ky) = (0, 0)"),
    ("bosons_n4_nphi6_v0v2_kx1_ky2", [4], 6, False, 1.0, 2, 1, BV0V2, ONE_KEY, "bosons V0 + 0.3 V2, (Kx, Ky) = (1, 2)"),
    ("fermions_su2_n4_nphi8_generic_sz0_ky0", [2, 2], 8, True, 1.0, 0, None, F8, SU2_KEYS, "SU(2): spin-dependent pseudopotentials, 2 up 2 down, Ky = 0"),
    ("fermions_su2_n4_nphi8_generic_sz2_ky0", [3, 1], 8, True, 1.0, 0, None, F8, SU2_KEYS, "SU(2): 3 up 1 down (2Sz = 2), Ky = 0"),
    ("bosons_su2_n4_nphi6_generic_sz0_ky0", [2, 2], 6, False, 1.0, 0, None, B6, SU2_KEYS, "SU(2) bosons, 2 up 2 down, Ky = 0"),
    ("fermions_su2_n4_nphi8_coulomb_sz0_ky0", [2, 2], 8, True, 1.0, 0, None, C2, SU2_KEYS, "SU(2)-symmetric Coulomb (FQHETorusFermionsWithSpin, no layer separation), Ky = 0"),
    ("fermions_su2_n4_nphi8_bilayer_d1_sz0_ky0", [2, 2], 8, True, 1.0, 0, None, _bilayer(1.0), SU2_KEYS, "bilayer Coulomb, inter-layer 2 pi exp(-q d)/q at d = 1, Ky = 0"),
    ("fermions_su2_n4_nphi8_bilayer_d0.5_sz0_kx0_ky0", [2, 2], 8, True, 1.0, 0, 0, _bilayer(0.5), SU2_KEYS, "bilayer d = 0.5 with translations, (Kx, Ky) = (0, 0)"),
    ("fermions_su2_n4_nphi8_bilayer_d0.5_sz0_kx1_ky3", [2, 2], 8, True, 1.0, 3, 1, _bilayer(0.5), SU2_KEYS, "bilayer d = 0.5 with translations, (Kx, Ky) = (1, 3)"),
    ("bosons_su2_n4_nphi6_generic_sz0_kx1_ky1", [2, 2], 6, False, 1.0, 1, 1, B6, SU2_KEYS, "SU(2) bosons with translations, (Kx, Ky) = (1, 1)"),
    ("bosons_su3_n4_nphi6_generic_ky0", [2, 1, 1], 6, False, 1.0, 0, None, S3, SU3_KEYS, "SU(3) bosons (2,1,1), six pair channels, Ky = 0"),
    ("bosons_su3_n3_nphi6_generic_kx1_ky2", [1, 1, 1], 6, False, 1.0, 2, 1, S3, SU3_KEYS, "SU(3) bosons (1,1,1) with translations, (Kx, Ky) = (1, 2)"),
    ("bosons_su4_n4_nphi5_generic_ky0", [1, 1, 1, 1], 5, False, 1.0, 0, None, S4, SU4_KEYS, "SU(4) bosons, ten pair channels, Ky = 0"),
]


def render_species(files, only=None):
    for name, numbers, nphi, fermion, r, ky, kx, pairs, keys, note in SPECIES_CASES:
        if only and name != only:
            continue
        kernels, keep = {}, {}
        for p, spec in pairs.items():
            kernels[p], keep[p] = _kernel(spec)
        if any(spec[0] == "pp" for spec in pairs.values()):
            files[f"{name}_pp.dat"] = ("# generated by tests/oracles/torus_ed.py; do not edit\nName = ed\n"
                                       + "".join(f"Pseudopotentials{keys[p]} = " + " ".join(f"{x:.16g}" for x in pairs[p][1]) + "\n"
                                                 for p in sorted(pairs) if pairs[p][0] == "pp"))
        sectors = range(nphi) if ky == "all" else [ky]
        rows = []
        for k in sectors:
            e = spectrum_species(numbers, nphi, kernels, fermion, r, k, keep, kx=kx)
            rows += [(f"{kx} {k} {x:.15g}\n" if kx is not None else f"{k} {x:.15g}\n") for x in e]
        files[f"{name}_spectrum.dat"] = (f"# {'fermions' if fermion else 'bosons'} N={sum(numbers)} per species {numbers} Nphi={nphi} ratio={r} "
                                         f"{'Kx=%d ' % kx if kx is not None else ''}Ky={ky}: {note}; {len(rows)} eigenvalues from tests/oracles/torus_ed.py "
                                         f"({'Kx Ky E' if kx is not None else 'Ky E'})\n" + "".join(rows))


def render_all(only=None):
    files = {}
    render_species(files, only)
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
                av = [float(l.split()[-1]) for l in p.read_text().splitlines() if l and not l.startswith("#")]
                bv = [float(l.split()[-1]) for l in text.splitlines() if l and not l.startswith("#")]
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
