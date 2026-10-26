#!/usr/bin/env python3
"""Independent Laughlin 1/3 quasihole spaces, and the two-layer time-reversal-symmetric problem built on them
(numpy; nothing shared with DiagHam).

    python3 tests/oracles/quasiholes_trs.py matrixelements PROGRAM NPHI WORK_DIR TOLERANCE
    python3 tests/oracles/quasiholes_trs.py pairing QH_PROGRAM PROGRAM CASE WORK_DIR TOLERANCE

Both modes run the programs themselves and need numpy at test time (the answers are computed, not stored).

Quasihole space. For fermions in the 2S + 1 orbitals of the lowest Landau level, the Laughlin 1/3 quasihole states
of N particles are the zero modes of the V1 pseudopotential Hamiltonian; tests/oracles/sphere_ed.py builds that
Hamiltonian in each Lz sector, and its null space is computed here, with its own orthonormal basis.

FQHESphereQuasiholeMatrixElements writes, for every (N, Lz) and orbital m, the dimension of the quasihole space, the
matrix of the density c+_m c_m within it and the matrix of c_m from N to N - 1 particles. Its basis is its own, so
the comparison uses what does not depend on the basis: the dimensions, the eigenvalues of each density matrix, the
singular values of each c_m matrix, and sum_m c+_m c_m = N (which holds only in an orthonormal basis).

FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing works in the product of two quasihole spaces (up and down
layers) with
    H = sum_m [Vup_m n_{m,up} + Vdown_m n_{m,down}] + sum_m Delta_m (c+_{m,up} c+_{m',down} + h.c.) + Ec (N - N0)^2,
every operator projected onto the quasihole spaces, for a fixed 2 Sz = N_up - N_down. The spectrum does not depend on
the basis of each quasihole sector nor on the fermion sign between layers (both are changes of basis). The one
convention taken from the program is which down orbital m' pairs with m and how the down layer's Lz enters the total
Lz; it is fixed by the per-sector spectra (see PAIR_CONVENTION).
"""
import itertools, math, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


def qh_layer(two_s):
    """{(N, 2Lz): (Q, states)}: an orthonormal basis (columns) of the V1 zero modes on the Fock states of the sector."""
    import numpy as np
    import sphere_ed
    out = {(0, 0): (np.ones((1, 1)), [()])}
    for n in range(1, two_s + 2):
        for two_lz in range(-(n * (two_s + 1 - n)), n * (two_s + 1 - n) + 1, 2):
            if n == 1:
                states = [(j,) for j in range(two_s + 1) if 2 * j - two_s == two_lz]
                if states:
                    out[(n, two_lz)] = (np.ones((1, 1)), states)
                continue
            H, states = sphere_ed.hamiltonian(n, two_s, [0.0, 1.0], True, two_lz)
            if not states:
                continue
            e, v = np.linalg.eigh(H)
            Q = v[:, np.abs(e) < 1e-9]
            if Q.shape[1]:
                out[(n, two_lz)] = (Q, states)
    return out


def density(Q, states, j):
    import numpy as np
    return Q.T @ np.diag([1.0 if j in s else 0.0 for s in states]) @ Q


def annihilation(layer, n, two_lz, j, two_s):
    """Projected c_j from (N, 2Lz) to (N - 1, 2Lz - 2Lz_j), or None."""
    import numpy as np
    key = (n - 1, two_lz - (2 * j - two_s))
    if (n, two_lz) not in layer or key not in layer:
        return None
    Q1, st1 = layer[(n, two_lz)]
    Q0, st0 = layer[key]
    idx = {s: i for i, s in enumerate(st0)}
    m = np.zeros((len(st0), len(st1)))
    for i, s in enumerate(st1):
        if j in s:
            l = list(s)
            k = l.index(j)
            l.remove(j)
            m[idx[tuple(l)], i] = (-1) ** k
    return Q0.T @ m @ Q1


# The program's convention, fixed by the per-sector spectra (the three others are 0.37 off on the first case):
# the down layer's potential index m is orbital m, Delta_m pairs up-orbital m with down-orbital m, and the sectors are
# labelled by 2 Lz_up - 2 Lz_down (the opposite chirality of the down layer). Relabelling the down orbitals m -> 2S - m
# together with Lz_up + Lz_down describes the same Hamiltonian.
PAIR_CONVENTION = (lambda m, two_s: m, lambda m, two_s: m, -1)


def trs_spectra(two_s, two_sz, vup, vdown, delta, ec, n0, conv, fixed_n=None):
    """{2 Lz_total: sorted eigenvalues}."""
    import numpy as np
    layer = qh_layer(two_s)
    pot_dn, partner, lz_sign = conv
    basis = []
    for (nu, lu), (Qu, _) in sorted(layer.items()):
        nd = nu - two_sz
        if nd < 0 or (fixed_n is not None and nu + nd != fixed_n):
            continue
        for (nd2, ld), (Qd, _) in sorted(layer.items()):
            if nd2 != nd:
                continue
            for a in range(Qu.shape[1]):
                for b in range(Qd.shape[1]):
                    basis.append((nu, lu, a, nd, ld, b))
    index = {s: i for i, s in enumerate(basis)}
    dim = len(basis)
    H = np.zeros((dim, dim))
    dens_up = {(n, lz, j): density(*layer[(n, lz)], j) for (n, lz) in layer for j in range(two_s + 1)}
    for i, (nu, lu, a, nd, ld, b) in enumerate(basis):
        H[i, i] += ec * (nu + nd - n0) ** 2
        for a2 in range(layer[(nu, lu)][0].shape[1]):
            v = sum(vup[j] * dens_up[(nu, lu, j)][a2, a] for j in range(two_s + 1))
            H[index[(nu, lu, a2, nd, ld, b)], i] += v
        for b2 in range(layer[(nd, ld)][0].shape[1]):
            v = sum(vdown[m] * dens_up[(nd, ld, pot_dn(m, two_s))][b2, b] for m in range(two_s + 1))
            H[index[(nu, lu, a, nd, ld, b2)], i] += v
    if fixed_n is None:
        for m in range(two_s + 1):
            if delta[m] == 0:
                continue
            jd = partner(m, two_s)
            for i, (nu, lu, a, nd, ld, b) in enumerate(basis):
                # annihilate a pair: c_{m,up} c_{jd,down} (the Hermitian conjugate gives the creation part)
                cu = annihilation(layer, nu, lu, m, two_s)
                cd = annihilation(layer, nd, ld, jd, two_s)
                if cu is None or cd is None:
                    continue
                for a2 in range(cu.shape[0]):
                    for b2 in range(cd.shape[0]):
                        amp = delta[m] * cu[a2, a] * cd[b2, b] * (-1) ** (nu - 1)
                        key = (nu - 1, lu - (2 * m - two_s), a2, nd - 1, ld - (2 * jd - two_s), b2)
                        k = index.get(key)
                        if k is not None and amp != 0:
                            H[k, i] += amp
                            H[i, k] += amp
    sectors = {}
    for i, (nu, lu, a, nd, ld, b) in enumerate(basis):
        sectors.setdefault(lu + lz_sign * ld, []).append(i)
    out = {}
    for L, idx in sectors.items():
        blk = H[np.ix_(idx, idx)]
        out[L] = sorted(np.linalg.eigvalsh(blk))
        off = np.delete(H[idx, :], idx, axis=1)
        assert np.max(np.abs(off), initial=0) < 1e-12, "Hamiltonian couples Lz sectors under this convention"
    return out


def read_txt_matrix(path):
    import numpy as np
    rows = [l.split() for l in Path(path).read_text().splitlines() if l.strip()]
    m = np.zeros((1 + max(int(r[0]) for r in rows), 1 + max(int(r[1]) for r in rows)))
    for r in rows:
        m[int(r[0]), int(r[1])] = float(r[2])
    return m


def run_matrixelements(argv):
    """FQHESphereQuasiholeMatrixElements -l NPHI against the independent quasihole spaces (basis-free invariants)."""
    import os, re, shutil, subprocess
    import numpy as np
    prog, two_s, work, tol = os.path.abspath(argv[0]), int(argv[1]), os.path.abspath(argv[2]), argv[3]
    rel = float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    r = subprocess.run([prog, "-l", str(two_s)], cwd=work, capture_output=True, text=True, timeout=600)
    d = Path(work) / f"nphi_{two_s}_sphere"
    if r.returncode != 0 or not d.is_dir():
        print(r.stdout[-1500:], r.stderr[-1500:]); print(f"program exit {r.returncode}"); return 1
    layer = qh_layer(two_s)
    problems, worst, counts = [], 0.0, {"sectors": 0, "density": 0, "annihilation": 0}
    for l in (d / f"fermions_qh_states_k_1_r_2_nphi_{two_s}.dat").read_text().splitlines():
        if l.strip() and not l.startswith("#"):
            _, n, lz, dim = map(int, l.split())
            counts["sectors"] += 1
            want = layer[(n, lz)][0].shape[1] if (n, lz) in layer else 0
            if dim != want:
                problems.append(f"N = {n}, 2Lz = {lz}: dimension {dim}, expected {want}")
    if counts["sectors"] != len(layer):
        problems.append(f"{counts['sectors']} sectors listed, {len(layer)} expected")
    sums = {}
    for f in sorted(d.glob("*_cdc_*.mat.txt")):
        n, lz, x = map(int, re.search(r"_n_(\d+)_nphi_\d+_lz_(-?\d+)_cdc_(-?\d+)\.mat", f.name).groups())
        prog_m = read_txt_matrix(f)
        sums[(n, lz)] = sums.get((n, lz), 0) + prog_m
        if n == 0:
            continue
        ora = density(*layer[(n, lz)], (x + two_s) // 2)
        dev = float(np.max(np.abs(np.linalg.eigvalsh(prog_m) - np.linalg.eigvalsh(ora))))
        worst = max(worst, dev)
        counts["density"] += 1
    for (n, lz), sm in sums.items():
        worst = max(worst, float(np.max(np.abs(sm - n * np.eye(len(sm))))))
    for f in sorted(d.glob("*_c_*.mat.txt")):
        n, lz, x = map(int, re.search(r"_n_(\d+)_nphi_\d+_lz_(-?\d+)_c_(-?\d+)\.mat", f.name).groups())
        prog_m = read_txt_matrix(f)
        ora = annihilation(layer, n, lz, (x + two_s) // 2, two_s)
        sp = np.sort(np.linalg.svd(prog_m, compute_uv=False))[::-1]
        so = np.sort(np.linalg.svd(ora, compute_uv=False))[::-1]
        k = min(len(sp), len(so))
        dev = float(np.max(np.abs(sp[:k] - so[:k]))) if k else 0.0
        dev = max(dev, float(np.max(sp[k:], initial=0)), float(np.max(so[k:], initial=0)))
        worst = max(worst, dev)
        counts["annihilation"] += 1
    if worst > rel:
        problems.append(f"largest deviation {worst:.2e} above {rel:g}")
    if counts["density"] == 0 or counts["annihilation"] == 0:
        problems.append("no matrix element files found")
    print(f"2S = {two_s}: {counts}, largest deviation {worst:.1e} (sum_m c+c = N included; allowed {rel:g})")
    for p_ in problems[:10]:
        print("  " + p_)
    return 1 if problems else 0


PAIRING_CASES = {  # name -> (2S, 2Sz, Vup, Vdown, Delta, Ec, N0, fixed N, extra options)
    "pairing_2s6_sz0": (6, 0, [0.9, 0.1, 0.35, 0.0, 0.2, 0.05, 0.7], [0.15, 0.6, 0.0, 0.45, 0.3, 0.8, 0.25],
                        [0.3, 0.5, 0.1, 0.25, 0.0, 0.4, 0.15], 0.3, 2.5, None, ["--memory", "0"]),
    "pairing_2s5_sz1": (5, 1, [0.4, 0.8, 0.05, 0.6, 0.2, 0.9], [0.7, 0.0, 0.3, 0.55, 0.1, 0.35],
                        [0.2, 0.45, 0.6, 0.0, 0.3, 0.15], 0.5, 1.5, None, ["--memory", "0"]),
    "fixed_n2_2s6_sz0": (6, 0, [0.9, 0.1, 0.35, 0.0, 0.2, 0.05, 0.7], [0.15, 0.6, 0.0, 0.45, 0.3, 0.8, 0.25],
                         [0.0] * 7, 0.0, 0.0, 2, ["--memory", "0"]),
    # U49: with the default --memory the program crashes after the first sector that uses fast multiplication
    "pairing_2s6_sz0_default_memory": (6, 0, [0.9, 0.1, 0.35, 0.0, 0.2, 0.05, 0.7], [0.15, 0.6, 0.0, 0.45, 0.3, 0.8, 0.25],
                                       [0.3, 0.5, 0.1, 0.25, 0.0, 0.4, 0.15], 0.3, 2.5, None, []),
}


def run_pairing(argv):
    """Quasihole matrix elements, then the two-layer program, every Lz sector against the oracle."""
    import glob, os, shutil, subprocess
    qh_prog, prog, name, work, tol = os.path.abspath(argv[0]), os.path.abspath(argv[1]), argv[2], os.path.abspath(argv[3]), argv[4]
    rel = float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)
    two_s, two_sz, vup, vdn, dl, ec, n0, fixed_n, extra = PAIRING_CASES[name]
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    r = subprocess.run([qh_prog, "-l", str(two_s)], cwd=work, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        print(r.stdout[-1500:]); print(f"FQHESphereQuasiholeMatrixElements exit {r.returncode}"); return 1
    Path(work, "potentials.dat").write_text(
        "OneBodyPotentialUpUp = " + " ".join(map(repr, vup)) + "\nOneBodyPotentialDownDown = " + " ".join(map(repr, vdn))
        + "\nOneBodyPotentialPairing = " + " ".join(map(repr, dl)) + "\n")
    cmd = [prog, "-l", str(two_s), "-s", str(two_sz), "--interaction-file", "potentials.dat", "--interaction-name", "oracle",
           "--directory", f"nphi_{two_s}_sphere", "--force-negativelz", "--full-diag", "5000",
           "--charging-energy", repr(ec), "--average-nbrparticles", repr(n0)] + extra
    if fixed_n is not None:
        cmd += ["--nbr-particles", str(fixed_n)]
    r = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=600)
    outs = glob.glob(os.path.join(work, "fermions_su2_quasiholes_oracle_*_lz.dat"))
    if r.returncode != 0 or len(outs) != 1:
        print(r.stdout[-800:], r.stderr[-800:]); print(f"program exit {r.returncode}, outputs {outs}"); return 1
    got = {}
    for l in open(outs[0]):
        if l.strip() and not l.startswith("#"):
            w = l.split()
            got.setdefault(int(w[0]), []).append(float(w[1]))
    want = trs_spectra(two_s, two_sz, vup, vdn, dl, ec, n0, PAIR_CONVENTION, fixed_n)
    problems, worst = [], 0.0
    if sorted(got) != sorted(want):
        problems.append(f"sectors {sorted(got)}, expected {sorted(want)}")
    for L in want:
        g, w = sorted(got.get(L, [])), want[L]
        if len(g) != len(w):
            problems.append(f"2Lz = {L}: {len(g)} levels, expected {len(w)}")
            continue
        for x, y in zip(g, w):
            worst = max(worst, abs(x - y) / max(abs(x), abs(y), 1.0))
    if worst > rel:
        problems.append(f"largest relative deviation {worst:.2e} above {rel:g}")
    print(f"{name}: {sum(len(v) for v in got.values())} levels in {len(got)} sectors, largest relative deviation {worst:.1e} (allowed {rel:g})")
    for p_ in problems[:10]:
        print("  " + p_)
    return 1 if problems else 0


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "matrixelements":
        sys.exit(run_matrixelements(sys.argv[2:]))
    if len(sys.argv) > 2 and sys.argv[1] == "pairing":
        sys.exit(run_pairing(sys.argv[2:]))
    print(__doc__)
    sys.exit(2)
