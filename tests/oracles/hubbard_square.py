#!/usr/bin/env python3
"""Independent answer for HubbardSquareLatticeModel (numpy; nothing shared with DiagHam).

    python3 tests/oracles/hubbard_square.py run PROGRAM WORK_DIR TOLERANCE LX LY N 2SZ U T TP [gutzwiller|none] [PROGRAM OPTIONS...]
    python3 tests/oracles/hubbard_square.py singlets PROGRAM PROJECTION WORK_DIR TOLERANCE LX LY N 2SZ U T TP KX KY NBR_STATES MIN [OPTIONS...]

The singlets mode takes the lowest eigenstates of one momentum sector from the program and checks
FTIMinNbrSingletProjection's weight of each on the configurations with at least MIN doubly occupied sites (on-site
singlets), for every eigenstate that is not degenerate within the sector (the weight does not depend on the basis, so
the program's 2D-translation basis is never needed).

Spin-1/2 fermions on the Lx x Ly square lattice with periodic boundaries, every momentum sector:
    H = sum_{k,sigma} eps(k) n_{k,sigma} + U sum_i n_{i,up} n_{i,down},
    eps(k) = 2 t (cos kx + cos ky) + 2 t' (cos(kx + ky) + cos(kx - ky)),
with k = 2 pi (kx / Lx, ky / Ly). The one-body part is built in real space from eps(k) (T_ij = (1/N_s) sum_k eps(k)
e^{i k (r_i - r_j)}, which also covers Lx or Ly = 2, where a neighbour is reached both ways). The Fock basis of the
2 Sz sector is built here, H is diagonalised in each joint eigenspace of the two translations (projector
(1 / N_s) sum_{a,b} e^{-i (Kx a + Ky b)} Tx^a Ty^b), and every (kx, ky) spectrum is compared with the program's.
With "gutzwiller" doubly occupied sites are excluded (and U plays no role). The signs in eps are the program's, fixed
on the first cases (the other three sign choices are 4 to 5 off); the sign of the momentum labels cannot be told apart
(the lattice is inversion symmetric). Options after the eleventh argument go to the program: --use-lapack in the
LAPACK builds, without which the built-in full diagonalisation is off by up to 4e-8 inside degenerate levels (U46).
"""
import cmath, itertools, math, os, re, shutil, subprocess, sys

# one BLAS thread: ctest runs many oracles at once, and a thread pool per process oversubscribes the machine
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")


def fock(ns, n, two_sz, gutz):
    nup = (n + two_sz) // 2
    ndn = n - nup
    out = []
    for up in itertools.combinations(range(ns), nup):
        for dn in itertools.combinations(range(ns), ndn):
            if gutz and set(up) & set(dn):
                continue
            out.append(frozenset([2 * i + 1 for i in up] + [2 * i for i in dn]))
    return out


def apply_cdc(state, a, b):
    """c+_a c_b |state>: (sign, new state) or None"""
    if b not in state or (a in state and a != b):
        return None
    s1 = sum(1 for x in state if x < b)
    st = state - {b}
    if a in st:
        return None
    s2 = sum(1 for x in st if x < a)
    return (-1) ** (s1 + s2), st | {a}


def spectra(lx, ly, n, two_sz, u, t, tp, gutz, ksign=1, min_doublons=None):
    """{(kx, ky): sorted eigenvalues}, or with min_doublons {(kx, ky): [(E, weight of the configurations with at least
    min_doublons doubly occupied sites)]}."""
    import numpy as np
    ns = lx * ly
    site = lambda x, y: (x % lx) + lx * (y % ly)
    ks = [(2 * math.pi * a / lx, 2 * math.pi * b / ly) for a in range(lx) for b in range(ly)]
    eps = lambda kx, ky: -2 * t * (math.cos(kx) + math.cos(ky)) - 2 * tp * (math.cos(kx + ky) + math.cos(kx - ky))
    T = np.zeros((ns, ns))
    for x1, y1, x2, y2 in itertools.product(range(lx), range(ly), range(lx), range(ly)):
        T[site(x1, y1), site(x2, y2)] = sum(eps(kx, ky) * math.cos(kx * (x1 - x2) + ky * (y1 - y2)) for kx, ky in ks) / ns
    basis = fock(ns, n, two_sz, gutz)
    index = {s: k for k, s in enumerate(basis)}
    dim = len(basis)
    H = np.zeros((dim, dim))
    for k, s in enumerate(basis):
        H[k, k] += u * sum(1 for i in range(ns) if 2 * i in s and 2 * i + 1 in s)
        for i, j in itertools.product(range(ns), range(ns)):
            if abs(T[i, j]) < 1e-14:
                continue
            for sg in (0, 1):
                r = apply_cdc(s, 2 * i + sg, 2 * j + sg)
                if r is not None and r[1] in index:
                    H[index[r[1]], k] += T[i, j] * r[0]
    def translation(dx, dy):
        M = np.zeros((dim, dim))
        perm = {}
        for x, y in itertools.product(range(lx), range(ly)):
            for sg in (0, 1):
                perm[2 * site(x, y) + sg] = 2 * site(x + dx, y + dy) + sg
        for k, s in enumerate(basis):
            img = [perm[m] for m in sorted(s)]
            # sign of the permutation that sorts the image (fermion reordering)
            inv = sum(1 for a in range(len(img)) for b in range(a + 1, len(img)) if img[a] > img[b])
            M[index[frozenset(img)], k] = (-1) ** inv
        return M
    Tx, Ty = translation(1, 0), translation(0, 1)
    out = {}
    Txp = [np.linalg.matrix_power(Tx, a) for a in range(lx)]
    Typ = [np.linalg.matrix_power(Ty, b) for b in range(ly)]
    for a, b in itertools.product(range(lx), range(ly)):
        Kx, Ky = 2 * math.pi * a / lx, 2 * math.pi * b / ly
        P = sum(cmath.exp(-1j * ksign * (Kx * p + Ky * q)) * (Txp[p] @ Typ[q]) for p in range(lx) for q in range(ly)) / ns
        e, v = np.linalg.eigh((P + P.conj().T) / 2)
        Q = v[:, e > 0.5]
        if Q.shape[1] and min_doublons is None:
            out[(a, b)] = sorted(np.linalg.eigvalsh(Q.conj().T @ H @ Q))
        elif Q.shape[1]:
            e, w = np.linalg.eigh(Q.conj().T @ H @ Q)
            V = Q @ w
            mask = np.array([sum(1 for i in range(ns) if 2 * i in s and 2 * i + 1 in s) >= min_doublons for s in basis])
            out[(a, b)] = [(float(e[i]), float(np.sum(np.abs(V[mask, i]) ** 2))) for i in range(len(e))]
    return out


def run_singlets(argv):
    """HubbardSquareLatticeModel eigenstates of one sector, then FTIMinNbrSingletProjection: the weight of each
    non-degenerate eigenstate on the configurations with at least MIN doubly occupied sites (on-site singlets)."""
    import glob
    prog, proj, work, tol = os.path.abspath(argv[0]), os.path.abspath(argv[1]), os.path.abspath(argv[2]), rel_of(argv[3])
    lx, ly, n, two_sz = int(argv[4]), int(argv[5]), int(argv[6]), int(argv[7])
    u, t, tp = float(argv[8]), float(argv[9]), float(argv[10])
    kx, ky, nst, mink = int(argv[11]), int(argv[12]), int(argv[13]), int(argv[14])
    extra = argv[15:]
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    r = subprocess.run([prog, "-p", str(n), "-x", str(lx), "-y", str(ly), "--sz-value", str(two_sz), "--u-potential", repr(u),
                        "--nn-t", repr(t), "--nnn-t", repr(tp), "--full-diag", "5000", "--only-kx", str(kx), "--only-ky", str(ky),
                        "--eigenstate", "-n", str(nst)] + extra, cwd=work, capture_output=True, text=True, timeout=600)
    vecs = sorted(glob.glob(os.path.join(work, f"*_kx_{kx}_ky_{ky}_sz_{two_sz}.*.vec")), key=lambda p: int(p.rsplit(".", 2)[1]))
    if r.returncode != 0 or len(vecs) != nst:
        print(r.stdout[-800:], r.stderr[-400:]); print(f"HubbardSquareLatticeModel exit {r.returncode}, {len(vecs)} vectors"); return 1
    with open(os.path.join(work, "states.dat"), "w") as f:
        f.write("\n".join(os.path.basename(v) for v in vecs) + "\n")
    r = subprocess.run([proj, "--degenerate-states", "states.dat", "--min-nbrsinglets", str(mink), "--only-weights"], cwd=work,
                       capture_output=True, text=True, timeout=600)
    got = {}
    for line in r.stdout.splitlines():
        m = re.match(r"\s*(\S+\.(\d+)\.vec) : (\S+)\s*$", line)
        if m:
            got[int(m.group(2))] = float(m.group(3))
    if r.returncode != 0 or len(got) != nst:
        print(r.stdout[-800:], r.stderr[-400:]); print(f"FTIMinNbrSingletProjection exit {r.returncode}, {len(got)} weights"); return 1
    levels = spectra(lx, ly, n, two_sz, u, -t, -tp, False, min_doublons=mink)
    sector = sorted(levels.get((kx, ky), levels.get(((-kx) % lx, (-ky) % ly), [])))
    problems, worst, checked = [], 0.0, 0
    for i in range(nst):
        e = sector[i][0]
        gap = min([abs(e - x[0]) for j, x in enumerate(sector) if j != i] + [1.0])
        if gap < 1e-8:
            continue                       # inside a degenerate level the weight depends on the chosen vector
        checked += 1
        worst = max(worst, abs(got[i] - sector[i][1]) / max(abs(got[i]), abs(sector[i][1]), 1.0))
    if checked == 0:
        problems.append("no non-degenerate state among the ones computed")
    if worst > tol:
        problems.append(f"largest relative deviation {worst:.2e} above {tol:g}")
    print(f"{lx} x {ly}, N = {n}, 2Sz = {two_sz}, sector ({kx}, {ky}), at least {mink} on-site singlet(s): {checked} of {nst} "
          f"weights checked (non-degenerate states), largest relative deviation {worst:.1e} (allowed {tol:g})")
    for p_ in problems[:10]:
        print("  " + p_)
    return 1 if problems else 0


def rel_of(tol):
    return float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)


def run(argv):
    prog, work, tol = os.path.abspath(argv[0]), os.path.abspath(argv[1]), rel_of(argv[2])
    lx, ly, n, two_sz = int(argv[3]), int(argv[4]), int(argv[5]), int(argv[6])
    u, t, tp = float(argv[7]), float(argv[8]), float(argv[9])
    gutz = len(argv) > 10 and argv[10] == "gutzwiller"
    argv = argv if len(argv) > 10 else argv + ["none"]
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    cmd = [prog, "-p", str(n), "-x", str(lx), "-y", str(ly), "--sz-value", str(two_sz), "--u-potential", repr(u),
           "--nn-t", repr(t), "--nnn-t", repr(tp), "--full-diag", "5000"] + (["--gutzwiller"] if gutz else []) + argv[11:]
    r = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=600)
    outs = [p for p in os.listdir(work) if p.endswith(".dat") and "tightbinding" not in p]
    if r.returncode != 0 or len(outs) != 1:
        print(r.stdout[-800:], r.stderr[-400:])
        print(f"program exit {r.returncode}, outputs {outs}")
        return 1
    got = {}
    for line in open(os.path.join(work, outs[0])):
        if line.strip() and not line.startswith("#"):
            w = line.split()
            got.setdefault((int(w[0]), int(w[1])), []).append(float(w[-1]))
    want = spectra(lx, ly, n, two_sz, u, -t, -tp, gutz)   # spectra() takes the hopping as -t, -t'
    problems, worst = [], 0.0
    if sorted(got) != sorted(want):
        problems.append(f"sectors {sorted(got)}, expected {sorted(want)}")
    for key, w in want.items():
        g = sorted(got.get(key, []))
        if len(g) != len(w):
            problems.append(f"sector {key}: {len(g)} levels, expected {len(w)}")
            continue
        worst = max([worst] + [abs(a - b) / max(abs(a), abs(b), 1.0) for a, b in zip(g, w)])
    if worst > tol:
        problems.append(f"largest relative deviation {worst:.2e} above {tol:g}")
    print(f"{lx} x {ly}, N = {n}, 2Sz = {two_sz}, U = {u:g}, t = {t:g}, t' = {tp:g}{' (Gutzwiller)' if gutz else ''}: "
          f"{sum(len(v) for v in got.values())} levels in {len(got)} sectors, largest relative deviation {worst:.1e} (allowed {tol:g})")
    for p_ in problems[:10]:
        print("  " + p_)
    return 1 if problems else 0


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "run":
        sys.exit(run(sys.argv[2:]))
    if len(sys.argv) > 2 and sys.argv[1] == "singlets":
        sys.exit(run_singlets(sys.argv[2:]))
    print(__doc__)
    sys.exit(2)
