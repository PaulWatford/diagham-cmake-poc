#!/usr/bin/env python3
"""Independent Laughlin 1/3 quasihole spaces on the cylinder, and the two-layer time-reversal-symmetric problem with
Josephson pairing built on them (numpy; nothing shared with DiagHam).

    python3 tests/oracles/cylinder_quasiholes.py matrixelements PROGRAM NPHI PERIMETER WORK_DIR TOLERANCE
    python3 tests/oracles/cylinder_quasiholes.py josephson QH_PROGRAM PROGRAM CASE WORK_DIR TOLERANCE

Quasihole space. The V1 zero modes are the polynomials that vanish as (z_i - z_j)^3; the space of their coefficients
on Slater determinants does not depend on the geometry. tests/oracles/quasiholes_trs.py builds it on the sphere, in
normalised sphere orbitals; here each configuration's amplitude is carried to the cylinder's normalised orbitals,
multiplied by prod_i sqrt(C(2S, l_i)) exp(kappa^2 (l_i - S)^2 / 2) (kappa = 2 pi / L; the same normalisation that
tests/oracles/mps_oes.py checks against the MPS programs), and the result orthonormalised.

FQHECylinderQuasiholeMatrixElements (run with --offdiagonal-density 2 --enable-ascii) writes, per (N, 2Lz), the
dimension, the matrices of c+_m c_m, of c+_m c_n (|m - n| <= 2; the file's Lz is the sector the operator acts on) and
of c_m. The comparison uses what does not depend on the basis: the dimensions, the eigenvalues of each density, the
singular values of each other matrix.

FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing works in the product of two quasihole spaces with
    H = sum_{m,n,sigma} V_{m,n,sigma} c+_{m,sigma} c_{n,sigma} + sum_m (Delta_m e^{i pi phi_m} c+_{m,up} c+_{m,down} + h.c.)
        + Ec (N - N0)^2,
for a fixed 2 Sz. The confining file lists (m, n, V_up, V_down[, phase_up, phase_down]) and the superconducting file
(m, m, Delta, phase). A line with m < n stands for V e^{-i pi phase} c+_m c_n and its Hermitian conjugate (the program
adds the conjugate itself; the phase columns are not on the wiki page, and their sign relative to the pairing's
e^{+i pi phi} is fixed by the spectrum of the off-diagonal case, the other three sign choices being 1e-2 off). The labels are those of FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing (fixed
there by the spectra): the down layer's orbital m is orbital m, Delta_m pairs up-orbital m with down-orbital m, and
without off-diagonal terms the sectors are labelled 2Lz_up - 2Lz_down. With off-diagonal terms the program writes one
spectrum. U57: a line with m > n is not read as the conjugate partner of (n, m) but as a coupling of orbitals m and
2m - n; the open-finding case lists both triangles of a Hermitian V, as the wiki's notation suggests.
"""
import cmath, glob, math, os, re, shutil, subprocess, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import quasiholes_trs as qt  # noqa: E402


def cyl_layer(two_s, perimeter):
    """{(N, 2Lz): (Q, states)}: an orthonormal basis of the zero modes in normalised cylinder orbitals."""
    import numpy as np
    kappa = 2 * math.pi / perimeter
    out = {}
    for key, (Q, states) in qt.qh_layer(two_s).items():
        f = np.array([math.prod(math.sqrt(math.comb(two_s, k)) * math.exp(kappa ** 2 * (k - two_s / 2) ** 2 / 2) for k in s)
                      for s in states])
        U, _, _ = np.linalg.svd(f[:, None] * Q, full_matrices=False)
        out[key] = (U, states)
    return out


def hop(layer, n, two_lz, a, b):
    """Projected c+_a c_b from (N, 2Lz) to (N, 2Lz + 2(a - b)), or None."""
    import numpy as np
    dst = (n, two_lz + 2 * (a - b))
    if (n, two_lz) not in layer or dst not in layer:
        return None
    Q1, st1 = layer[(n, two_lz)]
    Q0, st0 = layer[dst]
    idx = {s: i for i, s in enumerate(st0)}
    m = np.zeros((len(st0), len(st1)))
    for i, s in enumerate(st1):
        if b in s and (a not in s or a == b):
            l = list(s)
            kb = l.index(b)
            l.pop(kb)
            l2 = sorted(l + [a])
            m[idx[tuple(l2)], i] += (-1) ** kb * (-1) ** l2.index(a)
    return Q0.T @ m @ Q1


def read_txt_matrix(path):
    return qt.read_txt_matrix(path)


def rel_of(tol):
    return float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)


def run_matrixelements(argv):
    import numpy as np
    prog, two_s, perimeter, work, tol = os.path.abspath(argv[0]), int(argv[1]), float(argv[2]), os.path.abspath(argv[3]), rel_of(argv[4])
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    r = subprocess.run([prog, "-l", str(two_s), "--cylinder-perimeter", repr(perimeter), "--enable-ascii", "--offdiagonal-density", "2"],
                       cwd=work, capture_output=True, text=True, timeout=600)
    d = Path(work) / f"nphi_{two_s}_cylinder_perimeter_{perimeter:.6f}"
    if r.returncode != 0 or not d.is_dir():
        print(r.stdout[-1500:], r.stderr[-800:]); print(f"program exit {r.returncode}"); return 1
    layer = cyl_layer(two_s, perimeter)
    problems, worst, counts = [], 0.0, {"sectors": 0, "density": 0, "hopping": 0, "annihilation": 0}

    def sv_dev(p, o):
        sp = np.sort(np.linalg.svd(p, compute_uv=False))[::-1] if p.size else np.zeros(0)
        so = np.sort(np.linalg.svd(o, compute_uv=False))[::-1] if o.size else np.zeros(0)
        k = min(len(sp), len(so))
        return max(float(np.max(np.abs(sp[:k] - so[:k]), initial=0)), float(np.max(sp[k:], initial=0)), float(np.max(so[k:], initial=0)))
    for l in (d / f"fermions_cylinder_perimeter_{perimeter:.6f}_qh_states_k_1_r_2_nphi_{two_s}.dat").read_text().splitlines():
        if l.strip() and not l.startswith("#"):
            _, n, lz, dim = map(int, l.split())
            counts["sectors"] += 1
            want = layer[(n, lz)][0].shape[1] if (n, lz) in layer else 0
            if dim != want:
                problems.append(f"N = {n}, 2Lz = {lz}: dimension {dim}, expected {want}")
    if counts["sectors"] != len(layer):
        problems.append(f"{counts['sectors']} sectors listed, {len(layer)} expected")
    sums = {}
    for f in sorted(d.glob("*.mat.txt")):
        m = re.search(r"_n_(\d+)_nphi_\d+_lz_(-?\d+)_(cdc|c)_(-?\d+)(?:_(-?\d+))?\.mat", f.name)
        if not m:
            continue
        n, lz, kind, x, y = int(m.group(1)), int(m.group(2)), m.group(3), int(m.group(4)), m.group(5)
        prog_m = read_txt_matrix(f)
        if kind == "c":
            worst = max(worst, sv_dev(prog_m, qt.annihilation(layer, n, lz, (x + two_s) // 2, two_s)))
            counts["annihilation"] += 1
        elif y is None:
            sums[(n, lz)] = sums.get((n, lz), 0) + prog_m
            if n == 0:
                continue
            ora = qt.density(*layer[(n, lz)], (x + two_s) // 2)
            worst = max(worst, float(np.max(np.abs(np.linalg.eigvalsh(prog_m) - np.linalg.eigvalsh(ora)))))
            counts["density"] += 1
        else:
            ora = hop(layer, n, lz, (x + two_s) // 2, (int(y) + two_s) // 2)
            if ora is None:
                problems.append(f"{f.name}: no such sectors")
                continue
            worst = max(worst, sv_dev(prog_m, ora))
            counts["hopping"] += 1
    for (n, lz), sm in sums.items():
        worst = max(worst, float(np.max(np.abs(sm - n * np.eye(len(sm))))))
    if worst > tol:
        problems.append(f"largest deviation {worst:.2e} above {tol:g}")
    if min(counts.values()) == 0:
        problems.append(f"missing matrix element files: {counts}")
    print(f"2S = {two_s}, perimeter {perimeter:g}: {counts}, largest deviation {worst:.1e} (sum_m c+c = N included; allowed {tol:g})")
    for p_ in problems[:10]:
        print("  " + p_)
    return 1 if problems else 0


_VU = [0.9, 0.1, 0.35, 0.0, 0.2, 0.05, 0.7]
_VD = [0.15, 0.6, 0.0, 0.45, 0.3, 0.8, 0.25]
_DL = [0.3, 0.5, 0.1, 0.25, 0.0, 0.4, 0.15]
_PH = [0.0, 0.5, 0.25, 1.0, 0.0, 0.3, 0.7]
_DIAG = [(m, m, _VU[m], _VD[m], 0.0, 0.0) for m in range(7)]
JOSEPHSON_CASES = {  # name -> (2S, perimeter, 2Sz, confining lines, superconducting lines, Ec, N0, fixed N)
    "diagonal_2s6_perimeter6_sz0": (6, 6.0, 0, _DIAG, [(m, m, _DL[m], 0.0) for m in range(7)], 0.3, 2.5, None),
    "diagonal_2s5_perimeter5.5_sz1": (5, 5.5, 1, [(m, m, _VU[m], _VD[m], 0.0, 0.0) for m in range(6)],
                                      [(m, m, _DL[m], 0.0) for m in range(6)], 0.5, 1.5, None),
    "offdiagonal_phases_2s6_perimeter6_sz0": (6, 6.0, 0, _DIAG + [(0, 1, 0.2, 0.1, 0.3, -0.4), (2, 4, 0.15, 0.05, 0.0, 0.6), (3, 5, 0.0, 0.12, 0.0, 0.25)],
                                              [(m, m, _DL[m], _PH[m]) for m in range(7)], 0.3, 2.5, None),
    "offdiagonal_fixed_n4_2s6_perimeter6_sz0": (6, 6.0, 0, _DIAG + [(0, 1, 0.2, 0.1, 0.3, -0.4), (1, 3, 0.25, 0.0, 0.5, 0.0)],
                                                [(m, m, 0.0, 0.0) for m in range(7)], 0.0, 0.0, 4),
    # U57: both triangles of a Hermitian V, as the wiki's notation suggests
    "both_triangles_2s6_perimeter6_sz0": (6, 6.0, 0, _DIAG + [(0, 1, 0.2, 0.1, 0.0, 0.0), (1, 0, 0.2, 0.1, 0.0, 0.0)],
                                          [(m, m, _DL[m], 0.0) for m in range(7)], 0.3, 2.5, None),
}


def josephson_hamiltonian(two_s, perimeter, two_sz, conf, sc, ec, n0, fixed_n=None):
    """(H, basis): the Hamiltonian on the basis states (N_up, 2Lz_up, a, N_down, 2Lz_down, b)."""
    import numpy as np
    layer = cyl_layer(two_s, perimeter)
    vu, vd = {}, {}
    for m, n, a, b, pa, pb in conf:
        lo, hi = min(m, n), max(m, n)
        for table, v, ph in ((vu, a, pa), (vd, b, pb)):
            if v == 0:
                continue
            amp = v * cmath.exp(-1j * math.pi * ph)     # the program's sign, relative to the pairing's e^{+i pi phi}
            if m == n:
                table[(m, m)] = table.get((m, m), 0) + v
            else:
                # V_{m,n} c+_m c_n, and its conjugate; a pair listed twice (both triangles) is one Hermitian term
                if (lo, hi) not in table:
                    table[(m, n)] = amp
                    table[(n, m)] = amp.conjugate()
    basis = []
    for (nu, lu), (Qu, _) in sorted(layer.items()):
        nd = nu - two_sz
        if nd < 0 or (fixed_n is not None and nu + nd != fixed_n):
            continue
        for (nd2, ld), (Qd, _) in sorted(layer.items()):
            if nd2 == nd:
                basis += [(nu, lu, a, nd, ld, b) for a in range(Qu.shape[1]) for b in range(Qd.shape[1])]
    index = {s: i for i, s in enumerate(basis)}
    H = np.zeros((len(basis), len(basis)), dtype=complex)
    for i, (nu, lu, a, nd, ld, b) in enumerate(basis):
        H[i, i] += ec * (nu + nd - n0) ** 2
        for (m1, m2), v in vu.items():
            M = hop(layer, nu, lu, m1, m2)
            if M is not None:
                for a2 in range(M.shape[0]):
                    k = index.get((nu, lu + 2 * (m1 - m2), a2, nd, ld, b))
                    if k is not None:
                        H[k, i] += v * M[a2, a]
        for (m1, m2), v in vd.items():
            M = hop(layer, nd, ld, m1, m2)
            if M is not None:
                for b2 in range(M.shape[0]):
                    k = index.get((nu, lu, a, nd, ld + 2 * (m1 - m2), b2))
                    if k is not None:
                        H[k, i] += v * M[b2, b]
    if fixed_n is None:
        for m, _, dlt, ph in sc:
            if dlt == 0:
                continue
            amp_d = dlt * cmath.exp(1j * math.pi * ph)
            for i, (nu, lu, a, nd, ld, b) in enumerate(basis):
                cu = qt.annihilation(layer, nu, lu, m, two_s)
                cd = qt.annihilation(layer, nd, ld, m, two_s)
                if cu is None or cd is None:
                    continue
                for a2 in range(cu.shape[0]):
                    for b2 in range(cd.shape[0]):
                        amp = cu[a2, a] * cd[b2, b] * (-1) ** (nu - 1)
                        k = index.get((nu - 1, lu - (2 * m - two_s), a2, nd - 1, ld - (2 * m - two_s), b2))
                        if k is not None and amp != 0:
                            H[k, i] += amp_d.conjugate() * amp
                            H[i, k] += amp_d * amp
    return H, basis


def run_josephson(argv):
    import numpy as np
    qh_prog, prog, name, work, tol = os.path.abspath(argv[0]), os.path.abspath(argv[1]), argv[2], os.path.abspath(argv[3]), rel_of(argv[4])
    two_s, perimeter, two_sz, conf, sc, ec, n0, fixed_n = JOSEPHSON_CASES[name]
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    maxt = max([abs(m - n) for m, n, *_ in conf] + [1])
    r = subprocess.run([qh_prog, "-l", str(two_s), "--cylinder-perimeter", repr(perimeter), "--offdiagonal-density", str(maxt)],
                       cwd=work, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        print(r.stdout[-1500:]); print(f"FQHECylinderQuasiholeMatrixElements exit {r.returncode}"); return 1
    Path(work, "confining.dat").write_text("".join(f"{m} {n} {a!r} {b!r} {pa!r} {pb!r}\n" for m, n, a, b, pa, pb in conf))
    Path(work, "superconducting.dat").write_text("".join(f"{m} {n} {d!r} {p!r}\n" for m, n, d, p in sc))
    cmd = [prog, "-l", str(two_s), "-s", str(two_sz), "--cylinder-perimeter", repr(perimeter),
           "--directory", f"nphi_{two_s}_cylinder_perimeter_{perimeter:.6f}", "--interaction-name", "oracle",
           "--confining-file", "confining.dat", "--superconducting-file", "superconducting.dat",
           "--charging-energy", repr(ec), "--average-nbrparticles", repr(n0), "--full-diag", "5000", "--force-negativeky"]
    if fixed_n is not None:
        cmd += ["--nbr-particles", str(fixed_n)]
    r = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=600)
    outs = glob.glob(os.path.join(work, "fermions_cylinder_*su2_quasiholes_oracle_*.dat"))
    if r.returncode != 0 or len(outs) != 1:
        print(r.stdout[-800:], r.stderr[-800:]); print(f"program exit {r.returncode}, outputs {outs}"); return 1
    got = {}
    for l in open(outs[0]):
        if l.strip() and not l.startswith("#"):
            w = l.split()
            got.setdefault(int(w[0]) if len(w) > 1 else None, []).append(float(w[-1]))
    H, basis = josephson_hamiltonian(two_s, perimeter, two_sz, conf, sc, ec, n0, fixed_n)
    problems = []
    if not np.allclose(H, H.conj().T, atol=1e-13):
        problems.append("the oracle's Hamiltonian is not Hermitian")
    if None in got:
        want = {None: sorted(np.linalg.eigvalsh(H))}
    else:
        sectors = {}
        for i, (nu, lu, a, nd, ld, b) in enumerate(basis):
            sectors.setdefault(lu - ld, []).append(i)
        want = {}
        for L, idx in sectors.items():
            want[L] = sorted(np.linalg.eigvalsh(H[np.ix_(idx, idx)]))
            if np.max(np.abs(np.delete(H[idx, :], idx, axis=1)), initial=0) > 1e-12:
                problems.append("the Hamiltonian couples the Lz sectors, but the program wrote them separately")
                break
    worst = 0.0
    if sorted(got, key=str) != sorted(want, key=str):
        problems.append(f"sectors {sorted(got, key=str)}, expected {sorted(want, key=str)}")
    for L in want:
        g, w = sorted(got.get(L, [])), want[L]
        if len(g) != len(w):
            problems.append(f"sector {L}: {len(g)} levels, expected {len(w)}")
            continue
        worst = max([worst] + [abs(x - y) / max(abs(x), abs(y), 1.0) for x, y in zip(g, w)])
    if worst > tol:
        problems.append(f"largest relative deviation {worst:.2e} above {tol:g}")
    print(f"{name}: {sum(len(v) for v in got.values())} levels in {len(got)} sector(s), largest relative deviation {worst:.1e} (allowed {tol:g})")
    for p_ in problems[:10]:
        print("  " + p_)
    return 1 if problems else 0


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "matrixelements":
        sys.exit(run_matrixelements(sys.argv[2:]))
    if len(sys.argv) > 2 and sys.argv[1] == "josephson":
        sys.exit(run_josephson(sys.argv[2:]))
    print(__doc__)
    sys.exit(2)
