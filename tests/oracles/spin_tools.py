#!/usr/bin/env python3
"""Independent answers for spin-chain post-processing tools (numpy; nothing shared with DiagHam).

    python3 tests/oracles/spin_tools.py s2 SHOWBASIS PROGRAM WORK_DIR TOLERANCE N 2S 2SZ [NBR_DEGENERATE]
    python3 tests/oracles/spin_tools.py multient SHOWBASIS PROGRAM WORK_DIR TOLERANCE N 2SZ NBR_STATES

The input states are written here, with generic amplitudes, in DiagHam's basis order of an open chain in one Sz
sector, read from GenericSpinChainShowBasis at test time (spin 1 prints m; spin 1/2 prints "-" for up and "+" for
down: the 2Sz = +2 sector of two sites prints "- -").
- SpinChainComputeS2: <S^2> of one state, built from the single-site spin operators; with --degenerate-states the
  eigenvalues of S^2 in the span of several orthonormal states (here mixtures of different total spins), which are
  exact S(S + 1) values.
- SpinChainMultipleEntanglementSpectra: the entanglement spectrum of every column of a binary matrix of states (DiagHam's
  matrix layout), cut into the first l_a and the last N - l_a sites in GenericSpinChainShowBasis's printed order
  (fixed on the first cases with 2Sz != 0), block by block in Sz_A. The
  binary output is read as its wiki page describes. Without LAPACK the program crashes (U51).
"""
import math, re, struct, subprocess, sys
from pathlib import Path


def basis(showbasis, n, two_s, two_sz):
    """The states of the sector as tuples of 2m per printed site."""
    out = subprocess.run([showbasis, "-p", str(n), "-s", str(two_s), "-z", str(two_sz)], capture_output=True, text=True, timeout=120).stdout
    rows = []
    for l in out.splitlines():
        t = l.split()
        if len(t) == n and all(re.fullmatch(r"[+-]|-?\d+", x) for x in t):
            rows.append(tuple((-1 if x == "+" else 1) if two_s == 1 else 2 * int(x) for x in t))   # spin 1/2: "+" is down
    return rows


def s2_matrix(states, n, two_s):
    import numpy as np
    idx = {s: i for i, s in enumerate(states)}
    S = two_s / 2
    H = np.zeros((len(states), len(states)))
    for k, s in enumerate(states):
        m = [x / 2 for x in s]
        H[k, k] += sum(m) ** 2 + sum(S * (S + 1) - x * x for x in m)
        for i in range(n):
            for j in range(n):
                if i != j and m[i] < S and m[j] > -S:       # S+_i S-_j for every ordered pair i != j
                    t = list(s)
                    t[i] += 2
                    t[j] -= 2
                    H[idx[tuple(t)], k] += math.sqrt(S * (S + 1) - m[i] * (m[i] + 1)) * math.sqrt(S * (S + 1) - m[j] * (m[j] - 1))
    return H


def write_vec(path, v):
    with open(path, "wb") as f:
        f.write(struct.pack("<i", len(v)) + b"".join(struct.pack("<d", float(x)) for x in v))


def rel_of(tol):
    return float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)


def run_s2(argv):
    import os, shutil
    import numpy as np
    showbasis, prog, work, tol = os.path.abspath(argv[0]), os.path.abspath(argv[1]), os.path.abspath(argv[2]), rel_of(argv[3])
    n, two_s, two_sz = int(argv[4]), int(argv[5]), int(argv[6])
    ndeg = int(argv[7]) if len(argv) > 7 else 0
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    states = basis(showbasis, n, two_s, two_sz)
    H = s2_matrix(states, n, two_s)
    spin = "1_2" if two_s == 1 else str(two_s // 2)
    if ndeg == 0:
        v = np.array([math.cos(1.1 * i + 0.2) + 0.3 * math.sin(0.4 * i) for i in range(len(states))])
        v /= np.linalg.norm(v)
        name = f"spin_{spin}_oracle_n_{n}_sz_{two_sz}.0.vec"
        write_vec(Path(work) / name, v)
        cmd = [prog, "-i", name]
        want = [float(v @ H @ v)]
    else:
        # orthonormal mixtures of S^2 eigenstates with different S: the program must recover the S(S+1) values
        e, U = np.linalg.eigh(H)
        values = sorted(set(np.round(e, 6)))
        picks = [U[:, int(np.argmin(np.abs(e - values[k])))] for k in range(ndeg)]
        mix = np.array([[math.cos(0.7 * a + 0.3 * b + 0.1) for b in range(ndeg)] for a in range(ndeg)])
        Q, _ = np.linalg.qr(mix)
        names = []
        for k in range(ndeg):
            v = sum(Q[a, k] * picks[a] for a in range(ndeg))
            name = f"spin_{spin}_oracle_n_{n}_sz_{two_sz}.{k}.vec"
            write_vec(Path(work) / name, v)
            names.append(name)
        Path(work, "states.dat").write_text("\n".join(names) + "\n")
        cmd = [prog, "--degenerate-states", "states.dat"]
        want = sorted(float(values[k]) for k in range(ndeg))
    r = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=300)
    got = sorted(float(x) for x in re.findall(r"S\^2 = ([-0-9.e+]+)", r.stdout))
    if r.returncode != 0 or len(got) != len(want):
        print(r.stdout[-1500:], r.stderr[-800:]); print(f"program exit {r.returncode}, S^2 values {got}, expected {len(want)}"); return 1
    worst = max(abs(x - y) / max(abs(x), abs(y), 1.0) for x, y in zip(got, want))
    print(f"N = {n}, 2s = {two_s}, 2Sz = {two_sz}, {len(want)} value(s): S^2 = {got}, expected {want}; largest relative deviation {worst:.1e} (allowed {tol:g})")
    return 0 if worst <= tol else 1


def run_multient(argv):
    import os, shutil
    import numpy as np
    showbasis, prog, work, tol = os.path.abspath(argv[0]), os.path.abspath(argv[1]), os.path.abspath(argv[2]), rel_of(argv[3])
    n, two_sz, nstates = int(argv[4]), int(argv[5]), int(argv[6])
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    states = basis(showbasis, n, 1, two_sz)
    dim = len(states)
    M = np.array([[math.cos(0.7 * i + 1.3 * j + 0.1) + 0.2 * math.sin(0.3 * i * (j + 1)) for j in range(nstates)] for i in range(dim)])
    M, _ = np.linalg.qr(M)
    name = f"spin_1_2_oracle_n_{n}_sz_{two_sz}.eigenvec.mat"
    with open(Path(work) / name, "wb") as f:
        f.write(struct.pack("<iii", 1, dim, nstates))
        for i in range(dim):
            f.write(b"".join(struct.pack("<d", float(M[i, j])) for j in range(nstates)))
    r = subprocess.run([prog, "--multiple-states", name], cwd=work, capture_output=True, text=True, timeout=300)
    ents = list(Path(work).glob("*.full.ent"))
    if r.returncode != 0 or len(ents) != 1:
        print(r.stdout[-800:], r.stderr[-800:]); print(f"program exit {r.returncode}, outputs {ents}"); return 1
    data = ents[0].read_bytes()
    pos = 0

    def ints(k):
        nonlocal pos
        v = struct.unpack_from(f"<{k}i", data, pos)
        pos += 4 * k
        return v
    nspec, nsec = ints(2)
    secs = ints(nsec)
    got = []
    for _ in range(nspec):
        spec = {}
        for sec in secs:
            (nl,) = ints(1)
            spec[sec] = list(struct.unpack_from(f"<{nl}d", data, pos))
            pos += 8 * nl
        got.append(spec)
    la = int(re.search(r"_la_(\d+)", ents[0].name).group(1))
    problems, worst = [], 0.0
    if nspec != nstates:
        problems.append(f"{nspec} spectra, expected {nstates}")
    for s in range(min(nspec, nstates)):
        blocks = {}
        for amp, st in zip(M[:, s], states):
            blocks.setdefault(sum(st[:la]), {}).setdefault(st[la:], {})[st[:la]] = amp
        for sza in secs:
            byB = blocks.get(sza, {})
            rows = sorted({a for d in byB.values() for a in d})
            idx = {a: i for i, a in enumerate(rows)}
            m = np.zeros((len(rows), len(byB)))
            for j, d in enumerate(byB.values()):
                for a, v in d.items():
                    m[idx[a], j] = v
            lam = sorted((float(x) ** 2 for x in np.linalg.svd(m, compute_uv=False)), reverse=True) if m.size else []
            g = got[s][sza]
            k = min(len(g), len(lam))
            dev = max([abs(x - y) for x, y in zip(g[:k], lam[:k])] + [abs(x) for x in g[k:]] + [abs(x) for x in lam[k:]] + [0.0])
            worst = max(worst, dev)
    if worst > tol:
        problems.append(f"largest deviation {worst:.2e} above {tol:g}")
    print(f"N = {n}, 2Sz = {two_sz}, l_a = {la}, {nspec} spectra, Sz_A sectors {list(secs)}: largest deviation {worst:.1e} (allowed {tol:g})")
    for p_ in problems:
        print("  " + p_)
    return 1 if problems else 0


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "s2":
        sys.exit(run_s2(sys.argv[2:]))
    if len(sys.argv) > 2 and sys.argv[1] == "multient":
        sys.exit(run_multient(sys.argv[2:]))
    print(__doc__)
    sys.exit(2)
