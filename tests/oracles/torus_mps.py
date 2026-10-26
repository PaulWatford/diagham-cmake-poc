#!/usr/bin/env python3
"""Independent answer for FQHETorusMPSCreateState with the Laughlin MPS (numpy; nothing shared with DiagHam).

    python3 tests/oracles/torus_mps.py run PROGRAM WORK_DIR TOLERANCE NPHI RATIO SECTOR P_TRUNCATION [M]
    python3 tests/oracles/torus_mps.py sectors PROGRAM WORK_DIR TOLERANCE NPHI RATIO P1 P2 [M]

On the torus the Laughlin 1/m states are the zero modes of the pseudopotentials V_1, V_3, ... V_{m-2} (fermions, m odd)
or V_0, ..., V_{m-2} (bosons, m even); there is one in each of m momentum sectors Ky. The Hamiltonian is built by
tests/oracles/torus_ed.py (Landau gauge, aspect ratio Lx / Ly, the q = 0 term kept so that the zero modes are at 0).
run: the program's text output (-t: a coefficient and the occupation of every orbital per line) must be a normalised
vector equal, up to its sign, to the zero mode of the Ky sector in its file name, every coefficient within the
tolerance, and no zero mode may be missing a configuration.
sectors: --topological-sector 0, 1, ... m - 1 at two P truncations must give the same Ky each (U61: the program selects
the internal charge index modulo m, and that index is shifted by the truncation, so the sector moves with P mod m).
"""
import glob, math, os, re, shutil, subprocess, sys
from pathlib import Path

# one BLAS thread: ctest runs many oracles at once, and a thread pool per process oversubscribes the machine
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import torus_ed  # noqa: E402


def rel_of(tol):
    return float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)


def create(prog, work, nphi, ratio, sector, ptrunc, m):
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    cmd = [prog, "--p-truncation", str(ptrunc), "--nbr-fluxquanta", str(nphi), "--trim-qsector", "--aspect-ratio", repr(ratio),
           "--topological-sector", str(sector), "--laughlin-index", str(m), "-t"] + (["--boson", "--boson-truncation", str(nphi)] if m % 2 == 0 else [])
    r = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=600)
    txt = glob.glob(os.path.join(work, "*.vec.txt"))
    if r.returncode != 0 or len(txt) != 1:
        return None, r
    ky = int(re.search(r"_ky_(-?\d+)\.", txt[0]).group(1))
    coef = {}
    for line in open(txt[0]):
        t = line.split()
        if t:
            occ = [int(x) for x in t[1:]]
            coef[tuple(j for j, o in enumerate(occ) for _ in range(o))] = float(t[0])
    return (ky, coef), r


def zero_modes(n, nphi, ratio, ky, m):
    import numpy as np
    fermion = m % 2 == 1
    pp = [1.0 if (k < m - 1 and k % 2 == (1 if fermion else 0)) else 0.0 for k in range(m - 1)]
    H, states = torus_ed.spectrum(n, nphi, torus_ed.pseudopotential_kernel(pp), fermion, ratio, ky % nphi, True, matrix=True)
    e, v = np.linalg.eigh(H)
    return v[:, np.abs(e) < 1e-9], states


def run(argv):
    import numpy as np
    prog, work, tol = os.path.abspath(argv[0]), os.path.abspath(argv[1]), rel_of(argv[2])
    nphi, ratio, sector, ptrunc = int(argv[3]), float(argv[4]), int(argv[5]), int(argv[6])
    m = int(argv[7]) if len(argv) > 7 else 3
    res, r = create(prog, work, nphi, ratio, sector, ptrunc, m)
    if res is None:
        print(r.stdout[-800:], r.stderr[-400:]); print(f"program exit {r.returncode}, no text output"); return 1
    ky, coef = res
    n = len(next(iter(coef)))
    Z, states = zero_modes(n, nphi, ratio, ky, m)
    problems = []
    if Z.shape[1] != 1:
        problems.append(f"{Z.shape[1]} zero modes in Ky = {ky} (one expected)")
        print("\n".join(problems)); return 1
    z = Z[:, 0]
    index = {s: i for i, s in enumerate(states)}
    extra = [s for s in coef if s not in index]
    g = np.array([coef.get(s, 0.0) for s in states])
    norm = float(np.linalg.norm(g))
    sign = 1.0 if float(g @ z) >= 0 else -1.0
    worst = max(abs(norm - 1.0), float(np.max(np.abs(g - sign * z))))
    if extra:
        problems.append(f"{len(extra)} configurations outside the Ky = {ky} sector")
    if worst > tol:
        problems.append(f"largest deviation {worst:.2e} above {tol:g}")
    print(f"N = {n}, N_phi = {nphi}, aspect ratio {ratio:g}, sector {sector}, P {ptrunc}: Ky = {ky}, {len(coef)} coefficients "
          f"against the zero mode of {len(states)} states; norm {norm:.15f}, largest deviation {worst:.1e} (allowed {tol:g})")
    for p_ in problems[:10]:
        print("  " + p_)
    return 1 if problems else 0


def run_sectors(argv):
    prog, work, tol = os.path.abspath(argv[0]), os.path.abspath(argv[1]), rel_of(argv[2])
    nphi, ratio, p1, p2 = int(argv[3]), float(argv[4]), int(argv[5]), int(argv[6])
    m = int(argv[7]) if len(argv) > 7 else 3
    kys = {}
    for p in (p1, p2):
        for s in range(m):
            res, r = create(prog, os.path.join(work, f"p{p}_s{s}"), nphi, ratio, s, p, m)
            if res is None:
                print(r.stdout[-400:]); print(f"program exit {r.returncode} at P = {p}, sector {s}"); return 1
            kys[(p, s)] = res[0]
    same = all(kys[(p1, s)] == kys[(p2, s)] for s in range(m))
    print(f"N_phi = {nphi}: Ky of --topological-sector 0..{m - 1} at P = {p1}: {[kys[(p1, s)] for s in range(m)]}, "
          f"at P = {p2}: {[kys[(p2, s)] for s in range(m)]}")
    if not same:
        print("  the sector labels depend on the truncation")
    return 0 if same else 1


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "run":
        sys.exit(run(sys.argv[2:]))
    if len(sys.argv) > 2 and sys.argv[1] == "sectors":
        sys.exit(run_sectors(sys.argv[2:]))
    print(__doc__)
    sys.exit(2)
