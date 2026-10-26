#!/usr/bin/env python3
"""Independent answer for FQHECylinderDensity (pure Python; nothing shared with DiagHam).

    python3 tests/oracles/cylinder_density.py run SHOWBASIS PROGRAM WORK_DIR TOLERANCE N 2S 2LZ PERIMETER NBR_STATES [ROOT]

States with generic amplitudes are written here in DiagHam's basis order, read from FQHESphereShowBasis at test time
(occupations printed from orbital 0; every printed configuration is checked to have the requested Lz), in the
n-body basis or, with ROOT (a 0/1 string), in the Haldane basis of that root configuration (--haldane --reference-file).
The vector file name carries the perimeter, as the program expects. With several states they are orthonormalised here
and passed with --degenerate-states.
- Orbital occupations: <i| n_k |j> for every pair of states, printed in the header as "# k i j (value,0)" (the
  operator is diagonal in the basis, so <i| n_k |j> = sum_c v_i(c) v_j(c) n_k(c)).
- With one state, the density integrated over the perimeter, rho(x) = sum_k n_k exp(-(x - X_k)^2) / sqrt(pi), and the
  integrated charge sum_k n_k erfc(X_k - x) / 2, with X_k = 2 pi (k - N_phi / 2) / L (orbital 0 at the left end, the
  program's orientation), at every x the program prints; and the header's cylinder length 2 pi (N_phi + 1) / L.
"""
import math, os, re, shutil, struct, subprocess, sys


def basis(showbasis, n, two_s, two_lz, root, work):
    cmd = [showbasis, "--fermion", "-p", str(n), "-l", str(two_s), "-z", str(two_lz)]
    if root:
        with open(os.path.join(work, "root.dat"), "w") as f:
            f.write(f"NbrParticles = {n}\nLzMax = {two_s}\nReferenceState = {' '.join(root)}\n")
        cmd += ["--haldane", "--reference-file", "root.dat"]
    out = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=120).stdout
    rows = []
    for line in out.splitlines():
        t = line.split()
        if len(t) == two_s + 1 and all(x in ("0", "1") for x in t):
            rows.append([int(x) for x in t])
    for r in rows:
        assert sum(r) == n and sum(2 * k - two_s for k, o in enumerate(r) if o) == two_lz, r
    return rows


def rel_of(tol):
    return float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)


def run(argv):
    showbasis, prog, work, tol = os.path.abspath(argv[0]), os.path.abspath(argv[1]), os.path.abspath(argv[2]), rel_of(argv[3])
    n, two_s, two_lz, perimeter, nstates = int(argv[4]), int(argv[5]), int(argv[6]), float(argv[7]), int(argv[8])
    root = argv[9] if len(argv) > 9 else ""
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    rows = basis(showbasis, n, two_s, two_lz, root, work)
    dim = len(rows)
    # generic states, orthonormalised by Gram-Schmidt
    states = []
    for s in range(nstates):
        v = [math.cos(0.9 * i + 0.3 + 1.1 * s) + 0.4 * math.sin(1.7 * i * (s + 1)) for i in range(dim)]
        for u in states:
            p = sum(a * b for a, b in zip(u, v))
            v = [a - p * b for a, b in zip(v, u)]
        norm = math.sqrt(sum(a * a for a in v))
        states.append([a / norm for a in v])
    prefix = f"fermions_{'haldane_' if root else ''}cylinder_perimeter_{perimeter:.6f}_n_{n}_2s_{two_s}_lz_{two_lz}"
    names = [f"{prefix}.{s}.vec" for s in range(nstates)]
    for name, v in zip(names, states):
        with open(os.path.join(work, name), "wb") as f:
            f.write(struct.pack("<i", dim) + b"".join(struct.pack("<d", a) for a in v))
    cmd = [prog]
    if nstates == 1:
        cmd += ["--input-state", names[0]]
    else:
        with open(os.path.join(work, "states.dat"), "w") as f:
            f.write("\n".join(names) + "\n")
        cmd += ["--degenerate-states", "states.dat"]
    if root:
        cmd += ["--haldane", "--reference-file", "root.dat"]
    r = subprocess.run(cmd + ["--nbr-points", "60", "-o", "density.dat"], cwd=work, capture_output=True, text=True, timeout=300)
    path = os.path.join(work, "density.dat")
    if r.returncode != 0 or not os.path.exists(path):
        print(r.stdout[-800:], r.stderr[-400:])
        print(f"program exit {r.returncode}, no density file")
        return 1
    occ, pts, length = {}, [], None
    for line in open(path):
        m = re.match(r"#\s*(\d+)\s+(\d+)\s+(\d+)\s+\(([^,]+),([^)]+)\)", line)
        if m:
            occ[(int(m.group(1)), int(m.group(2)), int(m.group(3)))] = (float(m.group(4)), float(m.group(5)))
            continue
        m = re.match(r"#\s*Length = (\S+),", line)
        if m:
            length = float(m.group(1))
            continue
        if line.strip() and not line.startswith("#"):
            pts.append([float(x) for x in line.split()])
    problems, worst = [], 0.0

    def dev(a, b):
        return abs(a - b) / max(abs(a), abs(b), 1.0)
    want_occ = {}
    for k in range(two_s + 1):
        for i in range(nstates):
            for j in range(nstates):
                want_occ[(k, i, j)] = sum(states[i][c] * states[j][c] * rows[c][k] for c in range(dim))
    for key, w in want_occ.items():
        if key[1] > key[2] and key not in occ:
            continue                       # the program may print only i <= j
        if key not in occ:
            problems.append(f"no occupation for orbital {key[0]}, states {key[1]}, {key[2]}")
            continue
        worst = max(worst, dev(occ[key][0], w), abs(occ[key][1]))
    if nstates == 1:
        lx = 2 * math.pi * (two_s + 1) / perimeter
        if length is None or dev(length, lx) > tol:
            problems.append(f"cylinder length {length}, expected {lx}")
        nk = [want_occ[(k, 0, 0)] for k in range(two_s + 1)]
        xs = [2 * math.pi * (k - two_s / 2) / perimeter for k in range(two_s + 1)]
        if len(pts) < 10 or any(len(p) != 3 for p in pts):
            problems.append(f"{len(pts)} density lines")
        for x, rho, q in pts:
            worst = max(worst, dev(rho, sum(a * math.exp(-(x - b) ** 2) / math.sqrt(math.pi) for a, b in zip(nk, xs))),
                        dev(q, sum(a * math.erfc(b - x) / 2 for a, b in zip(nk, xs))))
    if worst > tol:
        problems.append(f"largest relative deviation {worst:.2e} above {tol:g}")
    print(f"N = {n}, 2S = {two_s}, 2Lz = {two_lz}, perimeter {perimeter:g}, {'Haldane basis of ' + root if root else 'n-body basis'}, "
          f"{dim} states in the sector, {nstates} input state(s): {len(occ)} occupations, {len(pts)} density points; "
          f"largest relative deviation {worst:.1e} (allowed {tol:g})")
    for p_ in problems[:10]:
        print("  " + p_)
    return 1 if problems else 0


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "run":
        sys.exit(run(sys.argv[2:]))
    print(__doc__)
    sys.exit(2)
