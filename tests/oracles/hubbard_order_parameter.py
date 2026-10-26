#!/usr/bin/env python3
"""Independent answer for HubbardSuperconductorOrderParameter (numpy; nothing shared with DiagHam).

    python3 tests/oracles/hubbard_order_parameter.py run SHOWBASIS PROGRAM WORK_DIR TOLERANCE NS N_RIGHT N_LEFTSTATES N_RIGHTSTATES [gutzwiller] [output-option]

Generic complex states (the program reads complex vectors) of N_RIGHT and N_RIGHT + 2 spin-1/2 fermions on NS sites, every
Sz, are written here in DiagHam's basis order, read from HubbardShowBasis --su2-spin at test time (each site printed as
0, u, d or X; site 0 first; an extra trailing token is ignored), with --gutzwiller in the Gutzwiller-projected basis.
Several states on one side are orthonormalised and passed with --degenerate-leftstates / --degenerate-rightstates.

The program writes, for every pair of sites i <= j, the four pair operators
    O = c+_{i,up} c+_{j,up} + c+_{i,down} c+_{j,down},   c+_{i,up} c+_{j,up} - c+_{i,down} c+_{j,down},
        c+_{i,up} c+_{j,down} + c+_{i,down} c+_{j,up},   c+_{i,up} c+_{j,down} - c+_{i,down} c+_{j,up}.
For one state on each side it gives |<L|O|R>|^2, |<L|O|R>| and the argument of <L|O|R>; with several states, the
min(n_L, n_R) largest eigenvalues of (1 / (n_L n_R)) M M+, M_ab = <L_a|O|R_b> (the order parameter matrix of the wiki
page, P_L O P_R O+ with P = (1/n) sum |Psi><Psi|). Fermion modes are ordered site by site, down before up (mode
2 i + sigma, sigma = 0 down); that and the overall sign of <L|O|R> (the arguments, which are fixed up to a common
sign by the choice of operator order) are taken from the program, fixed on the first single-state case: any other
mode order is at least 0.1 off. The site labels are those of HubbardShowBasis.
The output is read from the file the program names after the first left state (.vec -> .orderparam.dat); with
"output-option" the run passes -o order.dat and reads that file instead (U59: the program allocates the name but never
copies it, so the file is written under an uninitialised name).
"""
import cmath, math, os, re, shutil, struct, subprocess, sys


def basis(showbasis, n, ns, gutzwiller):
    cmd = [showbasis, "-p", str(n), "-x", str(ns), "--su2-spin"] + (["--gutzwiller"] if gutzwiller else [])
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=120).stdout
    rows = []
    for line in out.splitlines():
        t = line.split()
        if len(t) in (ns, ns + 1) and all(x in ("0", "u", "d", "X") for x in t[:ns]) and (len(t) == ns or t[ns] == "0"):
            modes = set()
            for s, x in enumerate(t[:ns]):
                if x in ("d", "X"):
                    modes.add(2 * s)
                if x in ("u", "X"):
                    modes.add(2 * s + 1)
            assert len(modes) == n
            rows.append(frozenset(modes))
    return rows


def states(dim, k, seed):
    import numpy as np
    M = np.array([[complex(math.cos(0.7 * i + 0.3 * a + seed), math.sin(0.4 * i * (a + 1) + 1.1 + seed)) for a in range(k)]
                  for i in range(dim)])
    Q, _ = np.linalg.qr(M)
    return Q


def pair_amplitudes(bl, br, vl, vr, ns):
    """{(a, b): <L| c+_a c+_b |R>} for every pair of modes, with the program's overall sign."""
    import numpy as np
    index = {m: k for k, m in enumerate(bl)}
    amp = {}
    for k, m in enumerate(br):
        if vr[k] == 0:
            continue
        for b in range(2 * ns):
            if b in m:
                continue
            s1 = -1 if sum(1 for x in m if x < b) % 2 else 1
            m2 = m | {b}
            for a in range(2 * ns):
                if a in m2:
                    continue
                s2 = -1 if sum(1 for x in m2 if x < a) % 2 else 1
                kk = index.get(m2 | {a})
                if kk is not None:
                    amp[(a, b)] = amp.get((a, b), 0) - np.conj(vl[kk]) * s1 * s2 * vr[k]
    return amp


def operators(amp, i, j):
    A = lambda si, sj: amp.get((2 * i + si, 2 * j + sj), 0)
    return [A(1, 1) + A(0, 0), A(1, 1) - A(0, 0), A(1, 0) + A(0, 1), A(1, 0) - A(0, 1)]


def rel_of(tol):
    return float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)


def write_cvec(path, v):
    with open(path, "wb") as f:
        f.write(struct.pack("<i", len(v)) + b"".join(struct.pack("<dd", x.real, x.imag) for x in v))


def run(argv):
    import numpy as np
    showbasis, prog, work, tol = os.path.abspath(argv[0]), os.path.abspath(argv[1]), os.path.abspath(argv[2]), rel_of(argv[3])
    ns, nr, kl, kr = int(argv[4]), int(argv[5]), int(argv[6]), int(argv[7])
    gutz = "gutzwiller" in argv[8:]
    use_o = "output-option" in argv[8:]
    nl = nr + 2
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    bl, br = basis(showbasis, nl, ns, gutz), basis(showbasis, nr, ns, gutz)
    VL, VR = states(len(bl), kl, 0.2), states(len(br), kr, 0.9)
    tag = "_gutzwiller" if gutz else ""
    names = {}
    for side, n, V, k in (("left", nl, VL, kl), ("right", nr, VR, kr)):
        names[side] = [f"fermions_hubbard{tag}_n_{n}_ns_{ns}_u_0.{a}.vec" for a in range(k)]
        for a, name in enumerate(names[side]):
            write_cvec(os.path.join(work, name), VL[:, a] if side == "left" else VR[:, a])
    cmd = [prog]
    for side, k in (("left", kl), ("right", kr)):
        if k == 1:
            cmd += [f"--{side}-state", names[side][0]]
        else:
            with open(os.path.join(work, f"{side}.dat"), "w") as f:
                f.write("\n".join(names[side]) + "\n")
            cmd += [f"--degenerate-{side}states", f"{side}.dat"]
    r = subprocess.run(cmd + (["-o", "order.dat"] if use_o else []), cwd=work, capture_output=True, text=True, timeout=600)
    path = os.path.join(work, "order.dat" if use_o else names["left"][0].replace(".vec", ".orderparam.dat"))
    if r.returncode != 0 or not os.path.exists(path):
        print(r.stdout[-800:], r.stderr[-400:])
        print(f"program exit {r.returncode}, no output")
        return 1
    single = kl == 1 and kr == 1
    kmin = min(kl, kr)
    amps = [[pair_amplitudes(bl, br, VL[:, a], VR[:, b], ns) for b in range(kr)] for a in range(kl)]
    problems, worst, nlines = [], 0.0, 0
    seen = set()
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        t = line.split()
        i, j = int(t[0]), int(t[1])
        seen.add((i, j))
        nlines += 1
        ops = [[operators(amps[a][b], i, j) for b in range(kr)] for a in range(kl)]
        rest = t[2:]
        for c in range(4):
            if single:
                w = ops[0][0][c]
                chunk, rest = rest[:3], rest[3:]
                n2 = 0.0 if chunk[0] == "(0,0)" else float(chunk[0])
                nrm, arg = float(chunk[1]), float(chunk[2])
                worst = max(worst, abs(n2 - abs(w) ** 2), abs(nrm - abs(w)))
                if abs(w) > 1e-6:
                    worst = max(worst, abs(cmath.phase(cmath.exp(1j * (arg - cmath.phase(w))))) * abs(w))
            else:
                M = np.array([[ops[a][b][c] for b in range(kr)] for a in range(kl)])
                want = sorted(np.linalg.eigvalsh(M @ M.conj().T / (kl * kr)), reverse=True)[:kmin]
                if rest and rest[0] == "(0,0)":      # i = j, the same-spin pairs vanish: "(0,0) 0 0"
                    chunk, rest = [], rest[3:]
                    worst = max([worst] + [abs(x) for x in want])
                    continue
                chunk, rest = rest[:kmin], rest[kmin:]
                worst = max([worst] + [abs(float(x) - y) / max(abs(y), 1.0) for x, y in zip(chunk, want)])
    expected = {(i, j) for i in range(ns) for j in range(i + (1 if gutz else 0), ns)}   # no on-site pairs with Gutzwiller
    if seen != expected:
        problems.append(f"{len(seen)} site pairs written, expected {len(expected)}")
    if worst > tol:
        problems.append(f"largest deviation {worst:.2e} above {tol:g}")
    print(f"{ns} sites, N = {nr} -> {nl}{' (Gutzwiller)' if gutz else ''}, {kl} left and {kr} right state(s): {nlines} site pairs, "
          f"largest deviation {worst:.1e} (allowed {tol:g})")
    for p_ in problems[:10]:
        print("  " + p_)
    return 1 if problems else 0


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "run":
        sys.exit(run(sys.argv[2:]))
    print(__doc__)
    sys.exit(2)
