#!/usr/bin/env python3
"""Independent answer for FQHESphereMPSCreateState with the Laughlin MPS (pure Python; nothing shared with DiagHam).

    python3 tests/oracles/mps_laughlin.py run PROGRAM WORK_DIR TOLERANCE N M P_TRUNCATION [sphere|fermion]
    python3 tests/oracles/mps_laughlin.py --check

The Laughlin state prod_{i<j} (z_i - z_j)^m is expanded here exactly, in integers, and its coefficients on the
monomials (bosons, m even) or Slater determinants (fermions, m odd) are compared with the program's text output
(-t: one "coefficient [partition]" line per basis state, the partition printed without its zero parts). The root
configuration (one particle every m orbitals, 2S = m (N - 1)) is written here. The program's unnormalised basis is
the monomial one, with the root coefficient 1, so the comparison needs no rescaling; a basis state the expansion does
not reach must carry 0. Bosons need --boson-truncation N (its default 1 keeps only partitions with distinct parts) and
a P truncation large enough for the whole squeezed space. Fermions (no --boson): the program exits 0 without writing
anything (U53), so the fermion case is an open-finding test; there the overall sign is left free.
With "sphere" the program's -n output (the normalised sphere basis) is compared with the same expansion carried to
normalised orbitals here: a monomial's coefficient times prod_l C(2S, l)^(-n_l / 2) / sqrt(n_l!) for occupations n_l,
then the whole state normalised, its sign fixed by the root.
--check (a self-test) evaluates the expansion, symmetrised (m even) or antisymmetrised (m odd), at random complex
points and compares it with the product itself.
"""
import math, os, re, shutil, subprocess, sys
from math import comb


def laughlin(n, m):
    """{decreasing partition with its zero parts: integer coefficient of z_1^l_1 ... z_N^l_N}."""
    poly = {tuple([0] * n): 1}
    for i in range(n):
        for j in range(i + 1, n):
            new = {}
            for e, c in poly.items():
                for k in range(m + 1):          # (z_i - z_j)^m = sum_k C(m, k) z_i^(m - k) (-z_j)^k
                    t = list(e)
                    t[i] += m - k
                    t[j] += k
                    t = tuple(t)
                    new[t] = new.get(t, 0) + c * comb(m, k) * (-1) ** k
            poly = {k: v for k, v in new.items() if v}
    return {e: c for e, c in poly.items() if list(e) == sorted(e, reverse=True) and (m % 2 == 0 or len(set(e)) == n)}


def rel_of(tol):
    return float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)


def run(argv):
    prog, work, tol = os.path.abspath(argv[0]), os.path.abspath(argv[1]), rel_of(argv[2])
    n, m, ptrunc = int(argv[3]), int(argv[4]), int(argv[5])
    mode = argv[6] if len(argv) > 6 else "monomial"
    fermion = mode == "fermion"
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    lzmax = m * (n - 1)
    root = " ".join("1" if o % m == 0 else "0" for o in range(lzmax + 1))
    with open(os.path.join(work, "root.dat"), "w") as f:
        f.write(f"NbrParticles = {n}\nLzMax = {lzmax}\nReferenceState = {root}\n")
    cmd = [prog, "--laughlin-index", str(m), "--reference-file", "root.dat", "--p-truncation", str(ptrunc), "-t", "state.txt"]
    if not fermion:
        cmd += ["--boson", "--boson-truncation", str(n)]
    if mode == "sphere":
        cmd += ["-n"]
    r = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=1200)
    got = {}
    path = os.path.join(work, "state.txt")
    if os.path.exists(path):
        for line in open(path):
            mm = re.match(r"\s*(\S+)\s+\[([0-9,]*)\]", line)
            if mm:
                got[tuple(int(x) for x in mm.group(2).split(",") if x and int(x))] = float(mm.group(1))
    exact = laughlin(n, m)
    if mode == "sphere":
        def factor(lam):
            f = 1.0
            for l in set(lam):
                k = lam.count(l)
                f *= comb(lzmax, l) ** (-k / 2) / math.sqrt(math.factorial(k))
            return f
        exact = {k: v * factor(k) for k, v in exact.items()}
        norm = math.sqrt(sum(v * v for v in exact.values()))
        exact = {k: v / norm for k, v in exact.items()}
    want = {tuple(x for x in k if x): v for k, v in exact.items()}
    if r.returncode != 0 or not got:
        print(r.stdout[-800:], r.stderr[-400:])
        print(f"program exit {r.returncode}, {len(got)} coefficients written (expected {len(want)})")
        return 1
    sign = 1.0
    if mode != "monomial":
        root_part = max(want)
        sign = 1.0 if got.get(root_part, 0) * want[root_part] >= 0 else -1.0
    worst = 0.0
    for k in set(want) | set(got):
        a, b = sign * got.get(k, 0.0), float(want.get(k, 0))
        worst = max(worst, abs(a - b) / max(abs(a), abs(b), 1.0))
    missing = [k for k in want if k not in got]
    print(f"N = {n}, m = {m}, P truncation {ptrunc}: {len(got)} basis states written, {len(want)} non-zero exact coefficients, "
          f"{len(missing)} missing; largest relative deviation {worst:.1e} (allowed {tol:g})")
    return 0 if worst <= tol and not missing else 1


def check():
    import cmath, itertools, random
    rng = random.Random(4493)
    worst = 0.0
    for n, m in ((2, 2), (3, 2), (3, 3), (4, 2), (4, 3), (3, 6), (5, 2)):
        coeffs = laughlin(n, m)
        for _ in range(3):
            z = [cmath.rect(rng.uniform(0.5, 1.5), rng.uniform(0, 2 * math.pi)) for _ in range(n)]
            direct = 1
            for i in range(n):
                for j in range(i + 1, n):
                    direct *= (z[i] - z[j]) ** m
            total = 0
            for lam, c in coeffs.items():
                if m % 2:   # Slater determinant: signed sum over every ordering of the distinct parts
                    for perm in itertools.permutations(range(n)):
                        inv = sum(1 for x in range(n) for y in range(x + 1, n) if perm[x] > perm[y])
                        term = -c if inv % 2 else c
                        for i in range(n):
                            term *= z[i] ** lam[perm[i]]
                        total += term
                else:       # monomial symmetric function: each distinct arrangement once
                    for arr in set(itertools.permutations(lam)):
                        term = c
                        for i in range(n):
                            term *= z[i] ** arr[i]
                        total += term
            worst = max(worst, abs(total - direct) / max(abs(direct), 1.0))
    print(f"exact Laughlin expansions (N = 2 to 5, m = 2, 3, 6) against the product at random points: largest relative deviation {worst:.1e}")
    return 0 if worst < 1e-10 else 1


if __name__ == "__main__":
    if sys.argv[1:] == ["--check"]:
        sys.exit(check())
    if len(sys.argv) > 2 and sys.argv[1] == "run":
        sys.exit(run(sys.argv[2:]))
    print(__doc__)
    sys.exit(2)
