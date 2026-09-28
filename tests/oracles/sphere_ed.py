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

With SU(2) spin (spectrum_su2): N_up and N_down particles in the same
orbitals, three pseudopotential sets V^{uu}, V^{dd}, V^{ud}; the same-spin
pairs enter as above and an up-down pair of relative angular momentum m has
energy 2 V^{ud}_m too (established with one up and one down particle: the
L = 2S multiplet of V^{ud}_0 = 1 sits at 2). The interaction file then has
the keys PseudopotentialsUpUp, PseudopotentialsDownDown, PseudopotentialsUpDown.

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


def _remove(block, j, fermion):
    """a_j on a sorted occupation tuple: (amplitude, new tuple) or None."""
    if j not in block:
        return None
    l = list(block)
    i = l.index(j)
    amp = (-1) ** i if fermion else math.sqrt(block.count(j))
    del l[i]
    return amp, tuple(l)


def _add(block, j, fermion):
    """a+_j on a sorted occupation tuple: (amplitude, new tuple) or None."""
    if fermion and j in block:
        return None
    l = sorted(block + (j,))
    amp = (-1) ** l.index(j) if fermion else math.sqrt(l.count(j))
    return amp, tuple(l)


def spectrum_su2(n_up, n_down, two_s, pp_uu, pp_dd, pp_ud, fermion, two_lz=0, scale=2.0, sz_parity=0, lz_parity=0):
    """Every eigenvalue of the (Lz, Sz) sector, Sz = (n_up - n_down)/2, for
    H = 2 sum_m [ V^uu_m sum_{up pairs} P_m + V^dd_m sum_{down pairs} P_m + V^ud_m sum_{up-down pairs} P_m ]."""
    import numpy as np
    norb = two_s + 1
    gen = itertools.combinations if fermion else itertools.combinations_with_replacement
    ups = list(gen(range(norb), n_up))
    downs = list(gen(range(norb), n_down))
    lz = lambda c: sum(2 * j - two_s for j in c)
    states = [(u, d) for u in ups for d in downs if lz(u) + lz(d) == two_lz]
    idx = {s: i for i, s in enumerate(states)}
    cache = {}

    def V(pp, j1, j2, j3, j4):
        k = (id(pp), j1, j2, j3, j4)
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
    for s in states:
        u, d = s
        col = idx[s]
        # same-spin terms: 1/2 sum V(1234) a+_1 a+_2 a_4 a_3 on one block, the other block untouched
        for which, block, pp in ((0, u, pp_uu), (1, d, pp_dd)):
            for j3 in set(block):
                r3 = _remove(block, j3, fermion)
                for j4 in set(r3[1]):
                    r4 = _remove(r3[1], j4, fermion)
                    pre = r3[0] * r4[0]
                    M2 = j3 + j4          # j1 + j2 must equal j3 + j4 (Lz conservation)
                    for j1 in range(norb):
                        j2 = M2 - j1
                        if j2 < 0 or j2 >= norb:
                            continue
                        v = V(pp, j1, j2, j3, j4)
                        if v == 0:
                            continue
                        a2 = _add(r4[1], j2, fermion)
                        if a2 is None:
                            continue
                        a1 = _add(a2[1], j1, fermion)
                        if a1 is None:
                            continue
                        new = (a1[1], d) if which == 0 else (u, a1[1])
                        H[idx[new], col] += 0.5 * scale * pre * a2[0] * a1[0] * v
        # up-down term: sum V(1234) a+_{1 up} a+_{2 down} a_{4 down} a_{3 up}, each up-down pair once
        for j3 in set(u):
            r3 = _remove(u, j3, fermion)
            for j4 in set(d):
                r4 = _remove(d, j4, fermion)
                pre = r3[0] * r4[0]
                M2 = j3 + j4
                for j1 in range(norb):
                    j2 = M2 - j1
                    if j2 < 0 or j2 >= norb:
                        continue
                    v = V(pp_ud, j1, j2, j3, j4)
                    if v == 0:
                        continue
                    a2 = _add(r4[1], j2, fermion)
                    if a2 is None:
                        continue
                    a1 = _add(r3[1], j1, fermion)
                    if a1 is None:
                        continue
                    H[idx[(a1[1], a2[1])], col] += scale * pre * a2[0] * a1[0] * v
    if sz_parity == 0 and lz_parity == 0:
        return np.linalg.eigvalsh(H)
    # symmetry-resolved spectrum: restrict H to the eigenspace of the exchange operator(s)
    B = np.eye(len(states))
    if sz_parity:
        if n_up != n_down:
            raise ValueError("Sz parity needs n_up == n_down")
        X = np.zeros_like(H)          # up <-> down: (u, d) -> (d, u), blocks reordered: sign (-1)^(n_up n_down) for fermions
        sgn = (-1) ** (n_up * n_down) if fermion else 1
        for s, i in idx.items():
            X[idx[(s[1], s[0])], i] = sgn
        B = _parity_basis(B, X, sz_parity)
    if lz_parity:
        if two_lz != 0:
            raise ValueError("Lz parity needs Lz = 0")
        X = np.zeros_like(H)          # j -> 2S - j in both blocks; a sorted block of n fermions reverses: (-1)^(n(n-1)/2)
        for s, i in idx.items():
            sgn = 1
            new = []
            for b in s:
                new.append(tuple(sorted(two_s - j for j in b)))
                if fermion:
                    sgn *= (-1) ** (len(b) * (len(b) - 1) // 2)
            X[idx[tuple(new)], i] = sgn
        B = _parity_basis(B, X, lz_parity)
    if B.shape[1] == 0:
        return np.zeros(0)
    return np.linalg.eigvalsh(B.T @ H @ B)


def _parity_basis(B, X, sign):
    """Orthonormal basis of the +/- eigenspace of the involution X inside the span of the columns of B."""
    import numpy as np
    P = B.T @ (np.eye(X.shape[0]) + sign * X) @ B / 2      # projector in the B basis
    w, v = np.linalg.eigh(P)
    keep = v[:, w > 0.5]
    return B @ keep


def spectrum_sun(numbers, two_s, pp, fermion, two_lz=0, scale=2.0):
    """k species (SU(k) spin, k = len(numbers)) with numbers[i] particles of species i in the same
    orbitals; pp[(i, j)] for i <= j are the pseudopotentials of an (i, j) pair. Same-species pairs
    enter as in spectrum(), different-species pairs as in spectrum_su2(): a pair of relative
    angular momentum m has energy 2 pp[(i, j)][m]. Every eigenvalue of the Lz sector."""
    import numpy as np
    k = len(numbers)
    norb = two_s + 1
    gen = itertools.combinations if fermion else itertools.combinations_with_replacement
    blocks = [list(gen(range(norb), n)) for n in numbers]
    lz = lambda c: sum(2 * j - two_s for j in c)
    states = [s for s in itertools.product(*blocks) if sum(lz(b) for b in s) == two_lz]
    idx = {s: i for i, s in enumerate(states)}
    cache = {}

    def V(key, j1, j2, j3, j4):
        ck = (key, j1, j2, j3, j4)
        if ck in cache:
            return cache[ck]
        v = 0.0
        M = (2 * j1 - two_s) + (2 * j2 - two_s)
        for m, Vm in enumerate(pp[key]):
            J = 2 * two_s - 2 * m
            if abs(M) > J or Vm == 0:
                continue
            v += Vm * cg(two_s, 2 * j1 - two_s, two_s, 2 * j2 - two_s, J, M) * cg(two_s, 2 * j3 - two_s, two_s, 2 * j4 - two_s, J, M)
        cache[ck] = v
        return v

    H = np.zeros((len(states), len(states)))
    for s in states:
        col = idx[s]
        for i in range(k):
            for j in range(i, k):
                if (i, j) not in pp:
                    continue
                if i == j:      # 1/2 sum V(1234) a+_1 a+_2 a_4 a_3 within block i
                    block = s[i]
                    for j3 in set(block):
                        r3 = _remove(block, j3, fermion)
                        for j4 in set(r3[1]):
                            r4 = _remove(r3[1], j4, fermion)
                            for j1 in range(norb):
                                j2 = j3 + j4 - j1
                                if j2 < 0 or j2 >= norb:
                                    continue
                                v = V((i, i), j1, j2, j3, j4)
                                if v == 0:
                                    continue
                                a2 = _add(r4[1], j2, fermion)
                                if a2 is None:
                                    continue
                                a1 = _add(a2[1], j1, fermion)
                                if a1 is None:
                                    continue
                                new = s[:i] + (a1[1],) + s[i + 1:]
                                H[idx[new], col] += 0.5 * scale * r3[0] * r4[0] * a2[0] * a1[0] * v
                else:           # sum V(1234) a+_{1 i} a+_{2 j} a_{4 j} a_{3 i}, each (i, j) pair once
                    for j3 in set(s[i]):
                        r3 = _remove(s[i], j3, fermion)
                        for j4 in set(s[j]):
                            r4 = _remove(s[j], j4, fermion)
                            for j1 in range(norb):
                                j2 = j3 + j4 - j1
                                if j2 < 0 or j2 >= norb:
                                    continue
                                v = V((i, j), j1, j2, j3, j4)
                                if v == 0:
                                    continue
                                a2 = _add(r4[1], j2, fermion)
                                if a2 is None:
                                    continue
                                a1 = _add(r3[1], j1, fermion)
                                if a1 is None:
                                    continue
                                new = list(s)
                                new[i] = a1[1]
                                new[j] = a2[1]
                                H[idx[tuple(new)], col] += scale * r3[0] * r4[0] * a2[0] * a1[0] * v
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
    ("fermions_n4_2s8_ll1_coulomb", 4, 10, True, 0, lambda: fs.coulomb_pseudopotentials(8, 1),
     "Coulomb in the first excited Landau level at 2S=8 (11 orbitals, l=5): the program is run with -l 10, Lz=0"),
]


SU2_KEYS = {(0, 0): "UpUp", (1, 1): "DownDown", (0, 1): "UpDown"}
SU3_KEYS = {(0, 0): "11", (0, 1): "12", (0, 2): "13", (1, 1): "22", (1, 2): "23", (2, 2): "33"}
SU4_NAMES = ["UpPlus", "UpMinus", "DownPlus", "DownMinus"]
SU4_KEYS = {(i, j): SU4_NAMES[i] + SU4_NAMES[j] for i in range(4) for j in range(i, 4)}


def _same(keys, pp):
    return {k: pp for k in keys}


def _c(two_s):
    return coulomb(two_s)


# The spinful cases: numbers = particles per species, pp = pseudopotentials per species pair (i <= j),
# keys = the interaction-file key of each pair, parity = symmetry sector (Sz <-> -Sz, Lz <-> -Lz),
# two_lz = one sector or "all" (every Lz >= 0 sector, for a program that writes them all).
SPIN_CASES = [
    dict(name="fermions_su2_n4_2s6_coulomb_sz0", numbers=[2, 2], two_s=6, fermion=True, two_lz=0, keys=SU2_KEYS,
         pp=lambda: _same(SU2_KEYS, _c(6)), note="SU(2)-symmetric Coulomb, Sz=0, Lz=0, 47 states"),
    dict(name="fermions_su2_n4_2s6_aniso_sz0", numbers=[2, 2], two_s=6, fermion=True, two_lz=0, keys=SU2_KEYS,
         pp=lambda: {(0, 0): [0, 1, 0, 0.3, 0, 0, 0], (1, 1): [0, 1, 0, 0.3, 0, 0, 0], (0, 1): [1, 0.5, 0, 0.2, 0, 0, 0]},
         note="spin-dependent pseudopotentials (V^ud with even m), Sz=0"),
    dict(name="fermions_su2_n4_2s7_coulomb_sz2", numbers=[3, 1], two_s=7, fermion=True, two_lz=0, keys=SU2_KEYS,
         pp=lambda: _same(SU2_KEYS, _c(7)), note="Coulomb, 3 up 1 down (2Sz=2), Lz=0"),
    dict(name="bosons_su2_n4_2s4_v0_v2_sz0", numbers=[2, 2], two_s=4, fermion=False, two_lz=0, keys=SU2_KEYS,
         pp=lambda: {(0, 0): [1, 0, 0, 0, 0], (1, 1): [1, 0, 0, 0, 0], (0, 1): [0.7, 0, 0.3, 0, 0]},
         note="bosons, V^uu=V^dd=V0, V^ud=0.7 V0 + 0.3 V2, Sz=0"),
    dict(name="bosons_su2_n4_2s5_v0_v2_sz2", numbers=[3, 1], two_s=5, fermion=False, two_lz=0, keys=SU2_KEYS,
         pp=lambda: {(0, 0): [1, 0, 0.2, 0, 0, 0], (1, 1): [1, 0, 0.2, 0, 0, 0], (0, 1): [0.7, 0, 0.3, 0, 0, 0]},
         note="bosons, 3 up 1 down (2Sz=2)"),
    # the four symmetry sectors of the first case (--szsymmetrized-basis, --lzsymmetrized-basis, +/- parity)
    dict(name="fermions_su2_n4_2s6_coulomb_sz0_szplus", numbers=[2, 2], two_s=6, fermion=True, two_lz=0, keys=SU2_KEYS,
         pp=lambda: _same(SU2_KEYS, _c(6)), parity=dict(sz_parity=+1), note="Coulomb Sz=0 Lz=0, Sz <-> -Sz even sector (25 states)"),
    dict(name="fermions_su2_n4_2s6_coulomb_sz0_szminus", numbers=[2, 2], two_s=6, fermion=True, two_lz=0, keys=SU2_KEYS,
         pp=lambda: _same(SU2_KEYS, _c(6)), parity=dict(sz_parity=-1), note="Sz <-> -Sz odd sector (22 states)"),
    dict(name="fermions_su2_n4_2s6_coulomb_sz0_lzplus", numbers=[2, 2], two_s=6, fermion=True, two_lz=0, keys=SU2_KEYS,
         pp=lambda: _same(SU2_KEYS, _c(6)), parity=dict(lz_parity=+1), note="Lz <-> -Lz even sector (28 states)"),
    dict(name="fermions_su2_n4_2s6_coulomb_sz0_lzminus", numbers=[2, 2], two_s=6, fermion=True, two_lz=0, keys=SU2_KEYS,
         pp=lambda: _same(SU2_KEYS, _c(6)), parity=dict(lz_parity=-1), note="Lz <-> -Lz odd sector (19 states)"),
    dict(name="fermions_su2_n4_2s6_coulomb_sz0_szplus_lzplus", numbers=[2, 2], two_s=6, fermion=True, two_lz=0, keys=SU2_KEYS,
         pp=lambda: _same(SU2_KEYS, _c(6)), parity=dict(sz_parity=+1, lz_parity=+1),
         note="both symmetries, even-even sector (20 states); the program's basis with both flags is empty (U32)"),
    # the legacy program QHEFermionsSphereWithSpin -v V -w W: a pair in relative m has energy V_m (no factor 2),
    # V in the up-down m=0 channel, W in every m=1 channel; it writes every Lz >= 0 sector
    dict(name="fermions_su2_n4_2s6_legacy_v1_w0.5", numbers=[2, 2], two_s=6, fermion=True, two_lz=0, keys=SU2_KEYS,
         pp=lambda: {(0, 0): [0, 0.25, 0, 0, 0, 0, 0], (1, 1): [0, 0.25, 0, 0, 0, 0, 0], (0, 1): [0.5, 0.25, 0, 0, 0, 0, 0]},
         note="QHEFermionsSphereWithSpin -v 1 -w 0.5 in this normalisation (V^ud_0 = v/2, V_1 = w/2), Lz=0"),
    # SU(3): three species, keys Pseudopotentials11 .. 33 (FQHESphere*WithSU3Spin; fermions: -t 2Tz -y 3Y, bosons: --nbr-n1/2/3)
    dict(name="fermions_su3_n3_2s4_coulomb", numbers=[1, 1, 1], two_s=4, fermion=True, two_lz=0, keys=SU3_KEYS,
         pp=lambda: _same(SU3_KEYS, _c(4)), note="SU(3)-symmetric Coulomb, one particle per species, Lz=0, 19 states"),
    dict(name="fermions_su3_n4_2s5_generic_tz1_y1", numbers=[2, 1, 1], two_s=5, fermion=True, two_lz=0, keys=SU3_KEYS,
         pp=lambda: {(0, 0): [0, 1, 0, 0.3, 0, 0], (1, 1): [0, 1, 0, 0.3, 0, 0], (2, 2): [0, 0.8, 0, 0, 0, 0],
                     (0, 1): [1, 0.5, 0, 0, 0, 0], (0, 2): [0.6, 0, 0.2, 0, 0, 0], (1, 2): [0.4, 0.1, 0, 0, 0, 0]},
         note="species-dependent pseudopotentials, (2,1,1) particles: 2Tz=1, 3Y=1, Lz=0, 64 states"),
    dict(name="bosons_su3_n4_2s3_generic", numbers=[2, 1, 1], two_s=3, fermion=False, two_lz=0, keys=SU3_KEYS,
         pp=lambda: {(0, 0): [1, 0, 0, 0], (1, 1): [1, 0, 0, 0], (2, 2): [0.5, 0, 0, 0],
                     (0, 1): [0.7, 0, 0.3, 0], (0, 2): [0.6, 0, 0, 0], (1, 2): [0.4, 0, 0.1, 0]},
         note="bosons, (2,1,1) particles, Lz=0, 26 states"),
    dict(name="bosons_su3_n3_2s5_coulomb", numbers=[1, 1, 1], two_s=5, fermion=False, two_lz=1, keys=SU3_KEYS,
         pp=lambda: _same(SU3_KEYS, _c(5)), note="bosons, SU(3) Coulomb, 2Lz=1 (three half-integer Lz), 27 states"),
    # SU(4): four species up/down x plus/minus, ten keys (FQHESphereBosonsWithSU4Spin --nbr-n1..n4)
    dict(name="bosons_su4_n4_2s3_generic", numbers=[1, 1, 1, 1], two_s=3, fermion=False, two_lz=0, keys=SU4_KEYS,
         pp=lambda: {(0, 0): [1, 0, 0, 0], (1, 1): [1, 0, 0, 0], (2, 2): [0.8, 0, 0, 0], (3, 3): [0.6, 0, 0, 0],
                     (0, 1): [0.7, 0, 0.3, 0], (0, 2): [0.5, 0, 0, 0], (0, 3): [0.4, 0, 0.1, 0],
                     (1, 2): [0.3, 0, 0, 0], (1, 3): [0.2, 0, 0.2, 0], (2, 3): [0.9, 0, 0, 0]},
         note="bosons, one particle per species, all ten pair channels different, Lz=0, 44 states"),
]


def spin_case_spectrum(c, two_lz):
    pp = c["pp"]()
    if len(c["numbers"]) == 2:
        return spectrum_su2(c["numbers"][0], c["numbers"][1], c["two_s"], pp[(0, 0)], pp[(1, 1)], pp[(0, 1)],
                            c["fermion"], two_lz, **c.get("parity", {}))
    return spectrum_sun(c["numbers"], c["two_s"], pp, c["fermion"], two_lz)

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
    for c in SPIN_CASES:
        name, pp, keys = c["name"], c["pp"](), c["keys"]
        files[f"{name}_pp.dat"] = (f"# {c['note']}; generated by tests/oracles/sphere_ed.py, do not edit\n"
                                   + "".join(f"Pseudopotentials{keys[k]} = " + " ".join(f"{x:.16g}" for x in pp[k]) + "\n"
                                             for k in sorted(pp)))
        n = sum(c["numbers"])
        if c["two_lz"] == "all":
            max_two_lz = sum(c["two_s"] - 2 * i for i in range(n)) if c["fermion"] else n * c["two_s"]
            sectors = list(range(max_two_lz % 2, max_two_lz + 1, 2))
        else:
            sectors = [c["two_lz"]]
        rows = []
        for two_lz in sectors:
            for x in spin_case_spectrum(c, two_lz):
                rows.append(f"{two_lz // 2 if two_lz % 2 == 0 else two_lz / 2} {x:.15g}\n")
        files[f"{name}_spectrum.dat"] = (f"# {'fermions' if c['fermion'] else 'bosons'} N={n} per species {c['numbers']} 2S={c['two_s']} 2Lz={c['two_lz']}: {c['note']}; "
                                         f"{len(rows)} eigenvalues from the independent SU(k) ED in tests/oracles/sphere_ed.py\n" + "".join(rows))
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
