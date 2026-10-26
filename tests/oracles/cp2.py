#!/usr/bin/env python3
"""Independent two-body problem in the lowest Landau level of CP2 (numpy, nothing shared with DiagHam).

    python3 tests/oracles/cp2.py run PROGRAM CASE WORK_DIR TOLERANCE      (the tests; needs numpy)

A particle in the lowest Landau level of CP2 with N_phi = p flux quanta lives in the SU(3) representation (p, 0): the
degree-p monomials z1^a1 z2^a2 z3^a3. In the orthonormal basis |a> the generators act as E_ij |a> =
sqrt((a_i + 1) a_j) |a + e_i - e_j>. Two particles decompose as (p,0) x (p,0) = sum_k (2p - 2k, k), k = 0..p; the
quadratic Casimir of the pair, sum_ij E_ij E_ji - (sum_i E_ii)^2 / 3 with E = E(1) + E(2), takes the value
(l^2 + m^2 + l m + 3 l + 3 m) / 3 on (l, m), distinct for every k, so its eigenvectors split the ordered pair space
into the representations. A two-body interaction given by pseudopotentials is H = sum_k V_k sum_{i<j} P_k(ij).

FQHESphereCP2GenericInteraction reads the pair decomposition from Clebsch-Gordan files
CGCSU3_p1_p_q1_0_p2_p_q2_0_p_l_q_m.dat ("Q(M) Q(M1) Q(M2) CGC", 1-based): this script writes them from the Casimir
eigenvectors. Q(M1) is the program's index of a single-particle state, the Gelfand-Tsetlin pattern (m12, m11) of
(p, 0) in its order (m12 outer, m11 inner), which is the monomial a = (m11, m12 - m11, p - m12) with weights
tz = a1 - a3 and y = a1 + a3 - 2 a2; Q(M) is any orthonormal basis of the representation. Phases of the states do not
change the spectrum of any (tz, y) sector. The oracle diagonalises H on the explicit Fock basis of each sector. The
program's H is exactly sum_k V_k sum_{i<j} P_k(ij), the k-th pseudopotential on (2p - 2k, k) (checked on two
particles: factor 1). It writes the sectors with tz >= 0; each of them is compared.
"""
import itertools, math, sys
from pathlib import Path

def orbitals(p):
    """The program's single-particle order: Gelfand-Tsetlin (m12 outer, m11 inner) -> exponents (a1, a2, a3)."""
    return [(m11, m12 - m11, p - m12) for m12 in range(p + 1) for m11 in range(m12 + 1)]


def weight(a):
    return (a[0] - a[2], a[0] + a[2] - 2 * a[1])


def generators(p):
    import numpy as np
    orbs = orbitals(p)
    idx = {a: i for i, a in enumerate(orbs)}
    E = {}
    for i in range(3):
        for j in range(3):
            m = np.zeros((len(orbs), len(orbs)))
            for a in orbs:
                if i == j:
                    m[idx[a], idx[a]] = a[i]
                elif a[j] > 0:
                    b = list(a)
                    b[i] += 1
                    b[j] -= 1
                    m[idx[tuple(b)], idx[a]] = math.sqrt((a[i] + 1) * a[j])
            E[(i, j)] = m
    return E


def pair_decomposition(p):
    """{k: matrix whose columns are an orthonormal basis of (2p - 2k, k) in the ordered pair space (d^2 rows)}."""
    import numpy as np
    E = generators(p)
    d = len(orbitals(p))
    I = np.eye(d)
    T = {k: np.kron(v, I) + np.kron(I, v) for k, v in E.items()}
    C = sum(T[(i, j)] @ T[(j, i)] for i in range(3) for j in range(3))
    Hc = sum(T[(i, i)] for i in range(3))
    C = C - Hc @ Hc / 3
    # diagonalise the Casimir weight block by weight block, so that every basis vector has a definite weight
    orbs = orbitals(p)
    wts = [tuple(x + y for x, y in zip(weight(orbs[i]), weight(orbs[j]))) for i in range(d) for j in range(d)]
    out = {k: [] for k in range(p + 1)}
    for w in sorted(set(wts)):
        idx = [n for n, x in enumerate(wts) if x == w]
        e, v = np.linalg.eigh(C[np.ix_(idx, idx)])
        for k in range(p + 1):
            l, m = 2 * p - 2 * k, k
            c = 2 * (l * l + m * m + l * m + 3 * l + 3 * m) / 3     # this operator is twice the usual normalisation
            for col in v[:, np.abs(e - c) < 1e-8].T:
                full = np.zeros(d * d)
                full[idx] = col
                out[k].append(full)
    for k in range(p + 1):
        l, m = 2 * p - 2 * k, k
        out[k] = np.array(out[k]).T
        dim = (l + 1) * (m + 1) * (l + m + 2) // 2
        assert out[k].shape[1] == dim, (k, out[k].shape, dim)
    return out


def write_cg_files(p, work):
    """The program's CGCSU3 files, one per representation (2p - 2k, k)."""
    d = len(orbitals(p))
    for k, cols in pair_decomposition(p).items():
        lines = ["#Q(M) Q(M1) Q(M2) CGC"]
        for qm in range(cols.shape[1]):
            for q1 in range(d):
                for q2 in range(d):
                    c = cols[q1 * d + q2, qm]
                    if abs(c) > 1e-15:
                        lines.append(f"{qm + 1} {q1 + 1} {q2 + 1} {float(c)!r}")
        Path(work, f"CGCSU3_p1_{p}_q1_0_p2_{p}_q2_0_p_{2 * p - 2 * k}_q_{k}.dat").write_text("\n".join(lines) + "\n")


def spectra(p, n, pp, boson, factor):
    """{(tz, y): sorted eigenvalues} of H = factor * sum_k V_k sum_{i<j} P_k(ij), every sector, in second quantisation:
    H = (factor / 2) sum_{abcd} <ab|V|cd> c+_a c+_b c_d c_c on occupation-number states."""
    import numpy as np
    orbs = orbitals(p)
    d = len(orbs)
    V = sum(pp[k] * (cols @ cols.T) for k, cols in pair_decomposition(p).items() if k < len(pp) and pp[k])

    def annihilate(occ, j):
        if occ[j] == 0:
            return None, 0.0
        amp = math.sqrt(occ[j]) if boson else (-1.0) ** sum(occ[:j])
        o = list(occ); o[j] -= 1
        return tuple(o), amp

    def create(occ, j):
        if not boson and occ[j] == 1:
            return None, 0.0
        amp = math.sqrt(occ[j] + 1) if boson else (-1.0) ** sum(occ[:j])
        o = list(occ); o[j] += 1
        return tuple(o), amp

    pick = itertools.combinations_with_replacement if boson else itertools.combinations
    sectors = {}
    for s in pick(range(d), n):
        occ = tuple(s.count(i) for i in range(d))
        w = tuple(sum(weight(orbs[i])[c] for i in s) for c in (0, 1))
        sectors.setdefault(w, []).append(occ)
    nz = [(a, b, c, dd, V[a * d + b, c * d + dd]) for a in range(d) for b in range(d) for c in range(d) for dd in range(d)
          if abs(V[a * d + b, c * d + dd]) > 1e-15]
    out = {}
    for w, states in sectors.items():
        pos = {s: i for i, s in enumerate(states)}
        H = np.zeros((len(states), len(states)))
        for s in states:
            for a, b, c, dd, v in nz:
                o1, x1 = annihilate(s, c)
                if o1 is None:
                    continue
                o2, x2 = annihilate(o1, dd)
                if o2 is None:
                    continue
                o3, x3 = create(o2, b)
                if o3 is None:
                    continue
                o4, x4 = create(o3, a)
                if o4 is None:
                    continue
                H[pos[o4], pos[s]] += 0.5 * factor * v * x1 * x2 * x3 * x4
        out[w] = sorted(np.linalg.eigvalsh(H))
    return out


CASES = {  # name -> (p, N, bosons?, pseudopotentials V_0 .. V_p)
    "fermions_n3_nphi2": (2, 3, False, [0.0, 1.0, 0.3]),
    "fermions_n3_nphi3": (3, 3, False, [0.0, 1.0, 0.0, 0.6]),
    "bosons_n3_nphi3": (3, 3, True, [1.0, 0.0, 0.5, 0.0]),
    "bosons_n4_nphi2": (2, 4, True, [1.0, 0.3, 0.7]),
}


def run(argv):
    import os, shutil, subprocess
    prog, name, work, tol = os.path.abspath(argv[0]), argv[1], os.path.abspath(argv[2]), argv[3]
    rel = float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)
    p, n, boson, pp = CASES[name]
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    write_cg_files(p, work)
    Path(work, "pseudopotentials.dat").write_text("Pseudopotentials = " + " ".join(map(repr, pp)) + "\n")
    cmd = [prog, "-p", str(n), "-l", str(p), "--interaction-file", "pseudopotentials.dat", "--interaction-name", "oracle",
           "--full-diag", "5000"] + (["--boson"] if boson else [])
    r = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=600)
    outs = [f for f in Path(work).glob("*.dat") if f.name.startswith(("fermions", "bosons"))]
    if r.returncode != 0 or len(outs) != 1:
        print(r.stdout[-1500:], r.stderr[-1500:]); print(f"program exit {r.returncode}, outputs {outs}"); return 1
    got = {}
    for l in outs[0].read_text().splitlines():
        if l.strip() and not l.startswith("#"):
            w = l.split()
            got.setdefault((int(w[0]), int(w[1])), []).append(float(w[2]))
    want = spectra(p, n, pp, boson, 1.0)
    problems, worst = [], 0.0
    if not got:
        problems.append("no sector written")
    for k, g in got.items():
        if k not in want:
            problems.append(f"sector {k} is not a sector of the problem")
            continue
        if len(g) != len(want[k]):
            problems.append(f"sector {k}: {len(g)} levels, expected {len(want[k])}")
            continue
        for x, y in zip(sorted(g), want[k]):
            worst = max(worst, abs(x - y) / max(abs(x), abs(y), 1.0))
    missing = [k for k in want if k[0] >= 0 and k not in got]
    if missing:
        problems.append(f"sectors with tz >= 0 not written: {missing}")
    if worst > rel:
        problems.append(f"largest relative deviation {worst:.2e} above {rel:g}")
    print(f"{name}: {sum(map(len, got.values()))} levels in {len(got)} sectors, largest relative deviation {worst:.1e} (allowed {rel:g})")
    for p_ in problems[:10]:
        print("  " + p_)
    return 1 if problems else 0


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "run":
        sys.exit(run(sys.argv[2:]))
    print(__doc__)
    sys.exit(2)
