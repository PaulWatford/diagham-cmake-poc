#!/usr/bin/env python3
"""Independent answer for FQHESphereMPSEntanglementSpectrumParticlePartition --orbital-es with the bosonic Laughlin MPS
on a finite cylinder (numpy for the singular values; the state itself from tests/oracles/mps_laughlin.py).

    python3 tests/oracles/mps_oes.py run PROGRAM WORK_DIR TOLERANCE N M P_TRUNCATION PERIMETER

The Laughlin state prod_{i<j} (z_i - z_j)^m with z = exp(2 pi (x + i y) / L) is expanded exactly; on the cylinder the
orbital l carries exp(kappa^2 (l - N_phi / 2)^2 / 2) relative to z^l (kappa = 2 pi / L), and a monomial coefficient
becomes a normalised boson Fock amplitude after division by sqrt(prod_l n_l!). The state is normalised, cut after the
first (N_phi + 1) / 2 orbitals (the program's default cut), and its Schmidt values are found block by block in
(N_A, Lz_A). The program's .full.ent file lists every eigenvalue lambda with its auxiliary labels (x, Q, P); the
sorted eigenvalues, and the entropy in the .ent file, are compared with the oracle's.
"""
import math, os, re, shutil, subprocess, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mps_laughlin import laughlin


def rel_of(tol):
    return float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)


def spectrum(n, m, perimeter):
    import numpy as np
    nphi = m * (n - 1)
    kappa = 2 * math.pi / perimeter
    amp = {}
    for lam, c in laughlin(n, m).items():
        o = [0] * (nphi + 1)
        for l in lam:
            o[l] += 1
        f = 1.0
        for l, k in enumerate(o):
            f *= math.exp(k * kappa ** 2 * (l - nphi / 2) ** 2 / 2) / math.sqrt(math.factorial(k))
        amp[tuple(o)] = c * f
    norm = math.sqrt(sum(a * a for a in amp.values()))
    la = (nphi + 1) // 2
    blocks = {}
    for o, a in amp.items():
        A, B = o[:la], o[la:]
        blocks.setdefault((sum(A), sum(l * k for l, k in enumerate(A))), {}).setdefault(B, {})[A] = a / norm
    ev = []
    for d in blocks.values():
        rows = sorted({a for x in d.values() for a in x})
        idx = {a: i for i, a in enumerate(rows)}
        M = np.zeros((len(rows), len(d)))
        for j, x in enumerate(d.values()):
            for a, v in x.items():
                M[idx[a], j] = v
        ev += [float(s) ** 2 for s in np.linalg.svd(M, compute_uv=False)]
    return la, sorted(ev, reverse=True)


def run(argv):
    prog, work, tol = os.path.abspath(argv[0]), os.path.abspath(argv[1]), rel_of(argv[2])
    n, m, ptrunc, perimeter = int(argv[3]), int(argv[4]), int(argv[5]), float(argv[6])
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    lzmax = m * (n - 1)
    root = " ".join("1" if o % m == 0 else "0" for o in range(lzmax + 1))
    Path(work, "root.dat").write_text(f"NbrParticles = {n}\nLzMax = {lzmax}\nReferenceState = {root}\n")
    cmd = [prog, "--boson", "--boson-truncation", str(n), "--laughlin-index", str(m), "--reference-file", "root.dat",
           "--p-truncation", str(ptrunc), "--orbital-es", "--normalize-cylinder", "--cylinder-perimeter", argv[6]]
    r = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=1200)
    full = list(Path(work).glob("*.full.ent"))
    ent = [p for p in Path(work).glob("*.ent") if not p.name.endswith(".full.ent")]
    got = []
    if full:
        for line in full[0].read_text().splitlines():
            t = line.split()
            if t and not line.startswith("#"):
                got.append(float(t[3]))
    if r.returncode != 0 or not got or len(ent) != 1:
        print(r.stdout[-800:], r.stderr[-600:])
        print(f"program exit {r.returncode}, {len(got)} eigenvalues written")
        return 1
    la, want = spectrum(n, m, perimeter)
    got.sort(reverse=True)
    k = min(len(got), len(want))
    # the program lists the non-zero part of the spectrum; anything beyond it in either list must be below tolerance
    worst = max([abs(a - b) for a, b in zip(got[:k], want[:k])] + [abs(x) for x in got[k:]] + [abs(x) for x in want[k:]] + [0.0])
    t = ent[0].read_text().split()
    s_got, s_want = float(t[1]), -sum(x * math.log(x) for x in want if x > 0)
    s_dev = abs(s_got - s_want) / max(abs(s_got), abs(s_want), 1.0)
    print(f"N = {n}, m = {m}, perimeter {perimeter:g}, P truncation {ptrunc}, cut after {la} orbitals: {len(got)} eigenvalues "
          f"(oracle {sum(1 for x in want if x > tol)} above tolerance); largest deviation {worst:.1e}; S_A = {s_got!r}, "
          f"expected {s_want!r} (relative deviation {s_dev:.1e}, allowed {tol:g})")
    return 0 if worst <= tol and s_dev <= tol and int(t[0]) == la else 1


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "run":
        sys.exit(run(sys.argv[2:]))
    print(__doc__)
    sys.exit(2)
