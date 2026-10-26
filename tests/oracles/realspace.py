#!/usr/bin/env python3
"""Closed-form orbital weights of a real-space region in the lowest Landau level (nothing shared with DiagHam).

    python3 tests/oracles/realspace.py run PROGRAM CASE WORK_DIR TOLERANCE      (the tests)

The *RealSpacePartitionCoefficients programs write, for every orbital, the weight of its density inside region A
("OrbitalSquareWeights = ..."). With the magnetic length set to 1:
- sphere, cap theta < theta_c (theta_c in units of pi): orbital j, counted from the north pole, has density
  proportional to cos(theta/2)^(2(2S - j)) sin(theta/2)^(2j); with t = sin^2(theta/2) this is a Beta(j + 1, 2S - j + 1)
  law, so the weight is the regularised incomplete beta function I_x(j + 1, 2S - j + 1), x = sin^2(theta_c / 2);
- cylinder of perimeter L, cut at x = 0: orbital j is a Gaussian exp(-(x - x_j)^2) centred at
  x_j = (2 pi / L) (j - N_phi / 2 + f), f the flux insertion; the weight on the x < 0 side is erfc(x_j) / 2;
- disk of radius R around the origin: orbital m has radial density r^(2m) exp(-r^2 / 2), so the weight is the
  regularised lower incomplete gamma function P(m + 1, R^2 / 2).
Each is summed with positive terms only. The comparison uses the suite's rel: rule (|a - b| <= tol max(|a|, |b|, 1)):
the cylinder program's smallest weights (1e-11 and below) carry about 1e-17 absolute, not relative, precision.
"""
import math, sys
from pathlib import Path


def sphere_cap(two_s, theta):
    x = math.sin(math.pi * theta / 2) ** 2
    n = two_s
    # I_x(j + 1, n - j + 1) = sum_{k = j + 1}^{n + 1} C(n + 1, k) x^k (1 - x)^(n + 1 - k)
    return [math.fsum(math.comb(n + 1, k) * x ** k * (1 - x) ** (n + 1 - k) for k in range(j + 1, n + 2)) for j in range(n + 1)]


def cylinder_cut(two_s, perimeter, flux=0.0):
    return [0.5 * math.erfc(2 * math.pi / perimeter * (j - two_s / 2 + flux)) for j in range(two_s + 1)]


def disk(two_s, radius):
    x = radius * radius / 2
    out = []
    for m in range(two_s + 1):
        # P(m + 1, x) = exp(-x) x^(m + 1) sum_k x^k / (m + 1 + k)!, positive terms
        terms, k, t = [], 0, x ** (m + 1) / math.factorial(m + 1)
        while True:
            terms.append(t)
            k += 1
            t *= x / (m + 1 + k)
            if t < 1e-18 * math.fsum(terms):
                break
        out.append(math.exp(-x) * math.fsum(terms))
    return out


CASES = {  # name -> (geometry, arguments for the program, oracle)
    "sphere_2s15_hemisphere": ("sphere", ["-s", "15", "--theta", "0.5"], lambda: sphere_cap(15, 0.5)),
    "sphere_2s10_theta0.3": ("sphere", ["-s", "10", "--theta", "0.3"], lambda: sphere_cap(10, 0.3)),
    "sphere_2s7_theta0.8": ("sphere", ["-s", "7", "--theta", "0.8"], lambda: sphere_cap(7, 0.8)),
    "cylinder_2s17_perimeter10": ("cylinder", ["-s", "17", "--cylinder-perimeter", "10"], lambda: cylinder_cut(17, 10.0)),
    "cylinder_2s12_perimeter8": ("cylinder", ["-s", "12", "--cylinder-perimeter", "8"], lambda: cylinder_cut(12, 8.0)),
    "disk_2s27_radius5": ("disk", ["-s", "27", "--radius", "5"], lambda: disk(27, 5.0)),
    "disk_2s10_radius2.5": ("disk", ["-s", "10", "--radius", "2.5"], lambda: disk(10, 2.5)),
}


def run(argv):
    import os, shutil, subprocess
    prog, name, work, tol = os.path.abspath(argv[0]), argv[1], os.path.abspath(argv[2]), argv[3]
    rel = float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)
    geometry, args, oracle = CASES[name]
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    r = subprocess.run([prog] + args + ["-o", "weights.dat"], cwd=work, capture_output=True, text=True, timeout=120)
    f = Path(work, "weights.dat")
    if r.returncode != 0 or not f.exists():
        print(r.stdout[-1500:], r.stderr[-1500:]); print(f"program exit {r.returncode}"); return 1
    got = None
    for l in f.read_text().splitlines():
        if l.strip().startswith("OrbitalSquareWeights"):
            got = [float(x) for x in l.split("=", 1)[1].split()]
    want = oracle()
    if got is None or len(got) != len(want):
        print(f"{name}: {None if got is None else len(got)} weights, expected {len(want)}"); return 1
    worst = max(abs(x - y) / max(abs(x), abs(y), 1.0) for x, y in zip(got, want))   # the suite's rel: definition
    print(f"{name}: {len(got)} weights, largest relative deviation {worst:.1e} (allowed {rel:g}); smallest weight {min(want):.1e}")
    return 0 if worst <= rel else 1


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "run":
        sys.exit(run(sys.argv[2:]))
    print(__doc__)
    sys.exit(2)
