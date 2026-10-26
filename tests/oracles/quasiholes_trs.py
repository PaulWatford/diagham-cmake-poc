#!/usr/bin/env python3
"""Independent Laughlin 1/3 quasihole spaces, and the two-layer time-reversal-symmetric problem built on them
(numpy; nothing shared with DiagHam).

    python3 tests/oracles/quasiholes_trs.py matrixelements PROGRAM NPHI WORK_DIR TOLERANCE
    python3 tests/oracles/quasiholes_trs.py pairing QH_PROGRAM PROGRAM CASE WORK_DIR TOLERANCE
    python3 tests/oracles/quasiholes_trs.py density QH_PROGRAM PROGRAM DENSITY_PROGRAM CASE WORK_DIR TOLERANCE [CONVERT_PROGRAM]

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


def trs_hamiltonian(two_s, two_sz, vup, vdown, delta, ec, n0, conv, fixed_n=None):
    """(H, basis, layer): the two-layer Hamiltonian on the basis states (N_up, 2Lz_up, a, N_down, 2Lz_down, b)."""
    import numpy as np
    layer = qh_layer(two_s)
    pot_dn, partner, _ = conv
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
    return H, basis, layer


def trs_sectors(basis, conv):
    """{2 Lz_total: basis indices}."""
    sectors = {}
    for i, (nu, lu, a, nd, ld, b) in enumerate(basis):
        sectors.setdefault(lu + conv[2] * ld, []).append(i)
    return sectors


def trs_spectra(two_s, two_sz, vup, vdown, delta, ec, n0, conv, fixed_n=None):
    """{2 Lz_total: sorted eigenvalues}."""
    import numpy as np
    H, basis, layer = trs_hamiltonian(two_s, two_sz, vup, vdown, delta, ec, n0, conv, fixed_n)
    sectors = trs_sectors(basis, conv)
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
    "fixed_n4_2s6_sz0": (6, 0, [0.9, 0.1, 0.35, 0.0, 0.2, 0.05, 0.7], [0.15, 0.6, 0.0, 0.45, 0.3, 0.8, 0.25],
                         [0.0] * 7, 0.0, 0.0, 4, ["--memory", "0"]),
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


DENSITY_CASES = {  # name -> (pairing case, 2Lz sector, number of states, check the charge imbalance as the wiki defines it, convert)
    "pairing_2s6_sz0_lz0_three_states": ("pairing_2s6_sz0", 0, 3, False, False),
    "pairing_2s5_sz1_lz1_two_states": ("pairing_2s5_sz1", 1, 2, True, False),
    "pairing_2s5_sz1_lz1_single_state": ("pairing_2s5_sz1", 1, 1, True, False),
    # U54: with an odd number of orbitals the program puts the central orbital into Q_R
    "pairing_2s6_sz0_lz0_imbalance_central_orbital": ("pairing_2s6_sz0", 0, 3, True, False),
    # FQHESphereQuasiholesWithSpinTimeReversalSymmetryConvertStates: fixed-N states to the space of all N, then the density
    "convert_fixed_n4_2s6_sz0_lz0_three_states": ("fixed_n4_2s6_sz0", 0, 3, False, True),
}


def run_density(argv):
    """The pairing program's lowest eigenstates of one sector, then FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity
    (after ConvertStates for the fixed-N case): <i|n_{m,sigma}|j>, the integrated charges and the charge imbalance."""
    import glob, os, re, shutil, subprocess
    import numpy as np
    qh_prog, prog, dens_prog, name, work, tol = (os.path.abspath(argv[0]), os.path.abspath(argv[1]), os.path.abspath(argv[2]), argv[3],
                                                 os.path.abspath(argv[4]), argv[5])
    conv_prog = os.path.abspath(argv[6]) if len(argv) > 6 else None
    rel = float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)
    case, sector, k, wiki_imbalance, convert = DENSITY_CASES[name]
    two_s, two_sz, vup, vdn, dl, ec, n0, fixed_n, extra = PAIRING_CASES[case]
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    run = lambda cmd: subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=600)
    r = run([qh_prog, "-l", str(two_s)])
    if r.returncode != 0:
        print(r.stdout[-1500:]); print(f"FQHESphereQuasiholeMatrixElements exit {r.returncode}"); return 1
    Path(work, "potentials.dat").write_text(
        "OneBodyPotentialUpUp = " + " ".join(map(repr, vup)) + "\nOneBodyPotentialDownDown = " + " ".join(map(repr, vdn))
        + "\nOneBodyPotentialPairing = " + " ".join(map(repr, dl)) + "\n")
    cmd = [prog, "-l", str(two_s), "-s", str(two_sz), "--interaction-file", "potentials.dat", "--interaction-name", "oracle",
           "--directory", f"nphi_{two_s}_sphere", "--full-diag", "5000", "--charging-energy", repr(ec),
           "--average-nbrparticles", repr(n0), "--initial-lz", str(sector), "--nbr-lz", "1", "--eigenstate", "--all-eigenstates"] + extra
    if fixed_n is not None:
        cmd += ["--nbr-particles", str(fixed_n)]
    r = run(cmd)
    vecs = [sorted(glob.glob(os.path.join(work, f"fermions_su2_quasiholes_oracle_*_lz_{sector}.{i}.vec"))) for i in range(k)]
    if r.returncode != 0 or any(len(v) != 1 for v in vecs):
        print(r.stdout[-800:], r.stderr[-800:]); print(f"pairing program exit {r.returncode}, eigenvectors {vecs}"); return 1
    vecs = [os.path.basename(v[0]) for v in vecs]
    if convert:
        Path(work, "fixed.dat").write_text("\n".join(vecs) + "\n")
        r = run([conv_prog, "--directory", f"nphi_{two_s}_sphere", "--degenerate-states", "fixed.dat"])
        conv = [v.replace(f"_n_{fixed_n}_", f"_fixedn_{fixed_n}_n_0_") for v in vecs]
        if r.returncode != 0 or not all(os.path.exists(os.path.join(work, c)) for c in conv):
            print(r.stdout[-800:], r.stderr[-800:]); print(f"ConvertStates exit {r.returncode}, expected {conv}"); return 1
        vecs = conv
    if k == 1:
        cmd = [dens_prog, "--input-state", vecs[0]]
    else:
        Path(work, "states.dat").write_text("\n".join(vecs) + "\n")
        cmd = [dens_prog, "--degenerate-states", "states.dat"]
    r = run(cmd + ["--directory", f"nphi_{two_s}_sphere"])
    dens = glob.glob(os.path.join(work, "*.density.dat"))
    if r.returncode != 0 or len(dens) != 1:
        print(r.stdout[-800:], r.stderr[-800:]); print(f"density program exit {r.returncode}, outputs {dens}"); return 1

    # the oracle: the same k lowest eigenstates of the sector, and every density matrix between them
    H, basis, layer = trs_hamiltonian(two_s, two_sz, vup, vdn, dl, ec, n0, PAIR_CONVENTION, fixed_n)
    idx = trs_sectors(basis, PAIR_CONVENTION)[sector]
    e, v = np.linalg.eigh(H[np.ix_(idx, idx)])
    problems = []
    if np.min(np.diff(e[:k + 1]), initial=1.0) < 1e-6:
        problems.append("the oracle's lowest levels are degenerate; the case cannot fix the states")
    V = v[:, :k]
    sub = [basis[i] for i in idx]
    pos = {s_: j for j, s_ in enumerate(sub)}
    A = np.zeros((2, two_s + 1, k, k))
    for m in range(two_s + 1):
        Nup, Ndn = np.zeros((len(sub), len(sub))), np.zeros((len(sub), len(sub)))
        for j, (nu, lu, a, nd, ld, b) in enumerate(sub):
            Du, Dd = density(*layer[(nu, lu)], m), density(*layer[(nd, ld)], m)
            for a2 in range(Du.shape[0]):
                Nup[pos[(nu, lu, a2, nd, ld, b)], j] += Du[a2, a]
            for b2 in range(Dd.shape[0]):
                Ndn[pos[(nu, lu, a, nd, ld, b2)], j] += Dd[b2, b]
        A[0, m], A[1, m] = V.T @ Nup @ V, V.T @ Ndn @ V
    worst = 0.0

    def dev(x, y):
        return abs(x - y) / max(abs(x), abs(y), 1.0)
    # <i|n|j>: the diagonal fixes each state's density; off the diagonal only |<i|n|j>| (each state's sign is free)
    got = {}
    for l in open(dens[0]):
        if l.strip() and not l.startswith("#"):
            t = l.split()
            got[tuple(int(x) for x in t[:4])] = float(t[4])
    for s_ in range(2):
        for m in range(two_s + 1):
            for i in range(k):
                for j in range(i, k):
                    g = got.get((s_, m, i, j))
                    if g is None:
                        problems.append(f"no <{i}|n_({m},{s_})|{j}> in the density file")
                        continue
                    w = A[s_, m, i, j] if i == j else abs(A[s_, m, i, j])
                    worst = max(worst, dev(g if i == j else abs(g), w))
    # integrated charges: the upper layer from orbital 0, the lower layer from orbital 2S (its time-reversed labels)
    lines = r.stdout.splitlines()
    start = next((n_ for n_, l in enumerate(lines) if l.startswith("Integrated charge")), None)
    if start is None:
        problems.append("no integrated charges printed")
    else:
        for m in range(two_s + 1):
            t = lines[start + 1 + m].split()
            vals = np.array([float(x) for x in t[1:]]).reshape(k, 3)
            up = np.linalg.eigvalsh(A[0, :m + 1].sum(axis=0))
            dn = np.linalg.eigvalsh(A[1, two_s - m:].sum(axis=0))
            tot = np.linalg.eigvalsh(A[0, :m + 1].sum(axis=0) + A[1, two_s - m:].sum(axis=0))
            if int(t[0]) != m:
                problems.append(f"integrated charge line {t[0]} where {m} was expected")
            for col, want in ((0, up), (1, dn), (2, tot)):
                for x, y in zip(sorted(vals[:, col]), want):
                    worst = max(worst, dev(x, y))
    if wiki_imbalance:
        # Q_L - Q_R = sum_{n < (N_phi + 1) / 2} sum_sigma (n_{n,sigma} - n_{N_phi - n,sigma}), the program's wiki page
        Q = sum(A[s_, n] - A[s_, two_s - n] for s_ in range(2) for n in range(two_s + 1) if n < (two_s + 1) / 2)
        want = np.linalg.eigvalsh(Q)
        imb = sorted(float(x) for x in re.findall(r"Charge imbalance between left and right \(orbital cut\) (\S+)", r.stdout))
        if len(imb) != k:
            problems.append(f"{len(imb)} charge imbalance values, expected {k}")
        else:
            worst = max([worst] + [dev(x, y) for x, y in zip(imb, want)])
    if worst > rel:
        problems.append(f"largest relative deviation {worst:.2e} above {rel:g}")
    print(f"{name}: {k} state(s) of 2Lz = {sector} ({'converted from fixed N' if convert else 'from the pairing program'}), "
          f"{len(got)} density matrix elements, integrated charges{', charge imbalance' if wiki_imbalance else ''}: "
          f"largest relative deviation {worst:.1e} (allowed {rel:g})")
    for p_ in problems[:10]:
        print("  " + p_)
    return 1 if problems else 0


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "matrixelements":
        sys.exit(run_matrixelements(sys.argv[2:]))
    if len(sys.argv) > 2 and sys.argv[1] == "pairing":
        sys.exit(run_pairing(sys.argv[2:]))
    if len(sys.argv) > 2 and sys.argv[1] == "density":
        sys.exit(run_density(sys.argv[2:]))
    print(__doc__)
    sys.exit(2)
