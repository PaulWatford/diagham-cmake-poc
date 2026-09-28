#!/usr/bin/env python3
"""Independent exact diagonalisation of k-body pseudopotential Hamiltonians on the sphere (numpy, no DiagHam).

    python3 tests/oracles/sphere_nbody.py --write tests/data/sphere_nbody
    python3 tests/oracles/sphere_nbody.py --check tests/data/sphere_nbody     (selftest.sphere_nbody_oracle)

A k-body pseudopotential Hamiltonian in the lowest Landau level is
    H = sum_m V^(k)_m sum_{i1<...<ik} P^(k)_m(i1..ik),
P^(k)_m the projector onto the k-particle states of relative angular momentum m,
i.e. of total angular momentum L = k l - m (l = S). For the m with a single
symmetric (bosons) or antisymmetric (fermions) relative state -- m = 0, 2, 3, 4, 5
for three bosons, m = 3, 5, 6, 7 for three fermions -- that projector is the
projector onto the L = k l - m multiplet of the k-particle space, which is
built here from L^2 (L_z, L_+, L_- summed over the particles) with no
Clebsch-Gordan algebra at all. The k-body "hard-core" interaction is the lowest
of them (m = 0 for bosons, m = k(k-1)/2 for fermions).

Second quantisation: with |b> the normalised k-particle Fock states and A_b
the operator that annihilates them (bosons: prod_j a_j^{b_j}/sqrt(b_j!);
fermions: a_{j_k} ... a_{j_1}, j_1 < ... < j_k), H = sum_ab <a|v|b> A+_a A_b.
DiagHam's normalisation (a k-particle state in the multiplet of V^(k)_m = 1
has energy SCALE[k]) is fixed with the k-particle case; nothing else is taken
from the program. The two-body term of the ThreeBodyGeneric programs goes
through the same machinery with the sphere_ed.py convention (energy 2 V_m).
"""
import argparse, itertools, math, sys
from collections import Counter
from pathlib import Path


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


def _annihilate(state, b, fermion):
    """A_b |state> for the sorted k-tuple b: (amplitude, remaining tuple) or None."""
    amp, block = 1.0, state
    for j in b:                       # a_{j_1} acts first (rightmost)
        r = _remove(block, j, fermion)
        if r is None:
            return None
        amp *= r[0]
        block = r[1]
    if not fermion:
        amp /= math.sqrt(math.prod(math.factorial(c) for c in Counter(b).values()))
    return amp, block


def _create(block, a, fermion):
    """A+_a |block> for the sorted k-tuple a: (amplitude, new tuple) or None."""
    amp = 1.0
    for j in reversed(a):             # A+_a = a+_{j_1} ... a+_{j_k}: a+_{j_k} acts first
        r = _add(block, j, fermion)
        if r is None:
            return None
        amp *= r[0]
        block = r[1]
    if not fermion:
        amp /= math.sqrt(math.prod(math.factorial(c) for c in Counter(a).values()))
    return amp, block


def kbody_projectors(k, two_s, fermion, ms):
    """{m: matrix of P^(k)_m} in the normalised k-particle Fock basis (all Lz), and that basis."""
    import numpy as np
    norb = two_s + 1
    gen = itertools.combinations if fermion else itertools.combinations_with_replacement
    basis = list(gen(range(norb), k))
    idx = {b: i for i, b in enumerate(basis)}
    n = len(basis)
    l = two_s / 2
    Lp, Lm = np.zeros((n, n)), np.zeros((n, n))
    Lz = np.zeros(n)
    for b in basis:
        Lz[idx[b]] = sum(j - l for j in b)
        for j in set(b):
            r = _remove(b, j, fermion)
            m = j - l
            for up, M in ((True, Lp), (False, Lm)):
                jj = j + 1 if up else j - 1
                if jj < 0 or jj > two_s:
                    continue
                c = math.sqrt(l * (l + 1) - m * (m + 1)) if up else math.sqrt(l * (l + 1) - m * (m - 1))
                a = _add(r[1], jj, fermion)
                if a is None:
                    continue
                M[idx[a[1]], idx[b]] += r[0] * a[0] * c
    L2 = np.diag(Lz * Lz) + (Lp @ Lm + Lm @ Lp) / 2
    w, v = np.linalg.eigh(L2)
    out = {}
    for m in ms:
        L = k * l - m
        if L < 0:
            continue
        sel = np.isclose(w, L * (L + 1), atol=1e-6)
        if sel.any():
            V = v[:, sel]
            out[m] = V @ V.T
    return out, basis


def multiplicities(k, two_s, fermion, mmax):
    """Number of k-particle multiplets with relative angular momentum m (dimension of the relative space):
    the projector P^(k)_m is the L-multiplet projector only where this is 1."""
    import numpy as np
    projs, basis = kbody_projectors(k, two_s, fermion, range(mmax + 1))
    return {m: int(round(np.trace(P) / (2 * (k * two_s / 2 - m) + 1))) for m, P in projs.items()}


# DiagHam's k-body pseudopotentials are not projector weights: for V^(k)_m = 1 a k-particle state in the
# L = k l - m multiplet has the energy below, which depends on m and on 2S (the program normalises its
# k-body basis states differently from the orthonormal multiplet used here). These numbers are read
# off the k-particle runs of the programs (the only thing taken from them); every N > k spectrum is
# then a prediction. The two-body convention (2 V_m, as in sphere_ed.py) is used for k = 2.
MULTIPLET_SCALE = {
    # (k, fermion, 2S): {m: energy}      programs FQHESphere{Bosons,Fermions}ThreeBodyGeneric
    (3, False, 6): {0: 1.5, 2: 2.9993111365882399, 3: 0.16391468575463, 4: 4.4401713055544567, 5: 0.31048185526808125},
    (3, False, 8): {0: 1.5, 2: 2.9996295838935980},
    (3, True, 10): {3: 0.22745249101270959, 5: 0.53751110975852878, 6: 0.0044512401606894292, 7: 0.94670491142186042},
}
HARDCORE_SCALE = {
    # (k, fermion, 2S): energy of the k-particle m_min multiplet   programs FQHESphereBosonsNBodyHardCore, QHEFermionsNBodyHardCore
    (3, False, 4): 0.091335201600665217, (3, False, 6): 0.076213083309916496, (3, False, 8): 0.069246696445559919, (3, False, 10): 0.065254110690944612,
    (4, False, 4): 0.055580539579273225, (4, False, 6): 0.043019963430062980,
    (3, True, 6): 0.25278087107535296, (3, True, 8): 0.23602486729577055, (3, True, 10): 0.22745249101270959,
}


def spectrum(n, two_s, kbody, fermion, two_lz=0, pp2=None, hardcore=False):
    """Every eigenvalue of the Lz sector for H = sum_k sum_m V^(k)_m P^(k)_m, kbody = {k: [V_0, V_1, ...]}
    (index m = relative angular momentum), plus an optional two-body pseudopotential list pp2;
    hardcore=True uses the NBodyHardCore programs' normalisation instead of the ThreeBodyGeneric one."""
    import numpy as np
    norb = two_s + 1
    gen = itertools.combinations if fermion else itertools.combinations_with_replacement
    states = [c for c in gen(range(norb), n) if sum(2 * j - two_s for j in c) == two_lz]
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    terms = dict(kbody)
    if pp2:
        terms[2] = list(pp2)
    for k, pps in terms.items():
        ms = [m for m, v in enumerate(pps) if v]
        if not ms:
            continue
        projs, kbasis = kbody_projectors(k, two_s, fermion, ms)
        kidx = {b: i for i, b in enumerate(kbasis)}
        if k == 2:
            scale = {m: 2.0 for m in ms}
        elif hardcore:
            scale = {m: HARDCORE_SCALE[(k, fermion, two_s)] for m in ms}
        else:
            scale = MULTIPLET_SCALE[(k, fermion, two_s)]
        if not projs:
            raise ValueError(f"no {k}-particle multiplet for m in {ms} at 2S={two_s} (fermions start at m={k * (k - 1) // 2})")
        V = sum(scale[m] * pps[m] * P for m, P in projs.items())
        nz = {j: [(i, V[i, j]) for i in range(len(kbasis)) if abs(V[i, j]) > 1e-14] for j in range(len(kbasis))}
        for s in states:
            col = idx[s]
            for b in set(itertools.combinations(s, k)):          # distinct k-sub-multisets of s
                r = _annihilate(s, b, fermion)
                if r is None:
                    continue
                for ia, v in nz[kidx[b]]:
                    c = _create(r[1], kbasis[ia], fermion)
                    if c is None:
                        continue
                    H[idx[c[1]], col] += r[0] * c[0] * v
    return np.linalg.eigvalsh(H)


CASES = [
    # name, N, 2S, fermion, 2Lz, {k: k-body pseudopotentials}, two-body pseudopotentials (or None), hardcore, note
    ("bosons_n5_2s6_hardcore3", 5, 6, False, 0, {3: [1]}, None, True, "three-body hard core, N=5 bosons at 2S=6 (two flux above the Moore-Read point), Lz=0: 32 states"),
    ("bosons_n6_2s4_hardcore4", 6, 4, False, 0, {4: [1]}, None, True, "four-body hard core, N=6 bosons at 2S=4, Lz=0: 18 states"),
    ("fermions_n5_2s10_hardcore3", 5, 10, True, 0, {3: [0, 0, 0, 1]}, None, True, "three-body hard core (m=3 projector), N=5 fermions at 2S=10, Lz=0: 32 states"),
    ("bosons_n5_2s6_v0_v2_v3", 5, 6, False, 0, {3: [1, 0, 0.5, 0.3]}, None, False, "three-body pseudopotentials V0=1, V2=0.5, V3=0.3, N=5 bosons 2S=6, Lz=0"),
    ("bosons_n5_2s6_v0_v2_v3_plus_twobody", 5, 6, False, 0, {3: [1, 0, 0.5, 0.3]}, [1, 0, 0.2], False, "the same plus a two-body V0=1, V2=0.2"),
    ("bosons_n5_2s8_v0_v2", 5, 8, False, 0, {3: [1, 0, 0.5]}, None, False, "three-body V0=1, V2=0.5 at 2S=8, N=5 bosons, Lz=0"),
    ("fermions_n5_2s10_v3_v5", 5, 10, True, 0, {3: [0, 0, 0, 1, 0, 0.5]}, None, False, "three-body V3=1, V5=0.5, N=5 fermions 2S=10, Lz=0"),
    ("fermions_n5_2s10_v3_v5_plus_twobody", 5, 10, True, 0, {3: [0, 0, 0, 1, 0, 0.5]}, [0, 1, 0, 0.3], False, "the same plus a two-body V1=1, V3=0.3"),
]


def render_all():
    files = {}
    for name, n, two_s, fermion, two_lz, kbody, pp2, hardcore, note in CASES:
        if not hardcore:
            lines = [f"# {note}; generated by tests/oracles/sphere_nbody.py, do not edit"]
            for k, pps in sorted(kbody.items()):
                lines.append("ThreebodyPseudopotentials = " + " ".join(f"{x:.16g}" for x in pps))
            if pp2:     # the program wants all 2S+1 two-body values
                padded = list(pp2) + [0] * (two_s + 1 - len(pp2))
                lines.append("Pseudopotentials = " + " ".join(f"{x:.16g}" for x in padded))
            files[f"{name}_pp.dat"] = "\n".join(lines) + "\n"
        e = spectrum(n, two_s, kbody, fermion, two_lz, pp2, hardcore)
        files[f"{name}_spectrum.dat"] = (f"# {'fermions' if fermion else 'bosons'} N={n} 2S={two_s} 2Lz={two_lz}: {note}; "
                                         f"{len(e)} eigenvalues from the independent k-body ED in tests/oracles/sphere_nbody.py\n"
                                         + "".join(f"{two_lz // 2 if two_lz % 2 == 0 else two_lz / 2} {x:.15g}\n" for x in e))
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
            if name.endswith("_pp.dat"):
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
