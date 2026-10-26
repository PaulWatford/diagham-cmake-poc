#!/usr/bin/env python3
"""Independent answer for FQHECylinderConfiningPotentialCoefficients (pure Python; nothing shared with DiagHam).

    python3 tests/oracles/confining.py run PROGRAM WORK_DIR TOLERANCE [PROGRAM OPTIONS...]
    python3 tests/oracles/confining.py --check

The lowest-Landau-level orbital k of a cylinder of perimeter L with N_phi flux quanta is a Gaussian strip centred at
X_k = 2 pi (k - N_phi / 2) / L, |phi_k(x)|^2 = exp(-(x - X_k)^2) / sqrt(pi) along the axis (magnetic length 1). For a
potential that depends on x only, the one-body matrix is diagonal and its k-th element is the integral of V against
that Gaussian. For V_R(x) = V0 (x - x0)^alpha (x > x0) it is V0 K_alpha(x0 - X_k) / sqrt(pi), and for the left
potential V_L(x) = V0 (x0 - x)^alpha (x < x0) it is V0 K_alpha(X_k - x0) / sqrt(pi), with
    K_alpha(c) = int_0^inf t^alpha exp(-(t + c)^2) dt,
    K_0 = sqrt(pi) erfc(c) / 2,  K_1 = exp(-c^2) / 2 - c K_0,  2 K_alpha = (alpha - 1) K_(alpha-2) - 2 c K_(alpha-1)
(integration by parts). The options the program is given are read here with the program's defaults (right: power 1,
strength 1, offset 0; left: power 1, strength 0, offset 0); the perimeter is --cylinder-perimeter, or
sqrt(2 pi (N_phi + 1) r) from the aspect ratio. The output's "OneBodyPotentials" line (or, with --spinful, both
"OneBodyPotentialUpUp" and "OneBodyPotentialDownDown") is compared element by element.
--check (a self-test) compares K_alpha with a direct Simpson-rule integration, for alpha = 0 to 4 and c of either sign.
"""
import math, os, re, shutil, subprocess, sys


def k_alpha(alpha, c):
    k0 = math.sqrt(math.pi) * math.erfc(c) / 2
    if alpha == 0:
        return k0
    k1 = math.exp(-c * c) / 2 - c * k0
    ks = [k0, k1]
    for a in range(2, alpha + 1):
        ks.append(((a - 1) * ks[a - 2] - 2 * c * ks[a - 1]) / 2)
    return ks[alpha]


def options(argv):
    o = {"s": 0, "r": 1.0, "perimeter": 0.0, "ra": 1, "rv": 1.0, "rx": 0.0, "la": 1, "lv": 0.0, "lx": 0.0, "spinful": False}
    names = {"-s": "s", "--nbr-flux": "s", "-r": "r", "--aspect-ratio": "r", "--cylinder-perimeter": "perimeter",
             "--confining-rightpower": "ra", "--confining-rightstrength": "rv", "--confining-rightoffset": "rx",
             "--confining-leftpower": "la", "--confining-leftstrength": "lv", "--confining-leftoffset": "lx"}
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--spinful":
            o["spinful"] = True
        elif a in names:
            key = names[a]
            o[key] = int(argv[i + 1]) if key in ("s", "ra", "la") else float(argv[i + 1])
            i += 1
        else:
            raise SystemExit(f"option {a} is not modelled by this oracle")
        i += 1
    if o["perimeter"] == 0.0:
        o["perimeter"] = math.sqrt(2 * math.pi * (o["s"] + 1) * o["r"])
    return o


def coefficients(o):
    out = []
    for k in range(o["s"] + 1):
        x = 2 * math.pi * (k - o["s"] / 2) / o["perimeter"]
        v = 0.0
        if o["rv"] != 0.0:
            v += o["rv"] * k_alpha(o["ra"], o["rx"] - x) / math.sqrt(math.pi)
        if o["lv"] != 0.0:
            v += o["lv"] * k_alpha(o["la"], x - o["lx"]) / math.sqrt(math.pi)
        out.append(v)
    return out


def rel_of(tol):
    return float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)


def run(argv):
    prog, work, tol, opts = os.path.abspath(argv[0]), os.path.abspath(argv[1]), rel_of(argv[2]), argv[3:]
    o = options(opts)
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    r = subprocess.run([prog] + opts + ["-o", "potential.dat"], cwd=work, capture_output=True, text=True, timeout=300)
    path = os.path.join(work, "potential.dat")
    lines = {}
    if os.path.exists(path):
        for line in open(path):
            m = re.match(r"\s*(OneBodyPotentials|OneBodyPotentialUpUp|OneBodyPotentialDownDown)\s*=(.*)", line)
            if m:
                lines[m.group(1)] = [float(x) for x in m.group(2).split()]
    keys = ["OneBodyPotentialUpUp", "OneBodyPotentialDownDown"] if o["spinful"] else ["OneBodyPotentials"]
    if r.returncode != 0 or any(k not in lines for k in keys):
        print(r.stdout[-800:], r.stderr[-400:])
        print(f"program exit {r.returncode}, lines {sorted(lines)} (expected {keys})")
        return 1
    want = coefficients(o)
    worst, where = 0.0, None
    for key in keys:
        got = lines[key]
        if len(got) != len(want):
            print(f"{key}: {len(got)} coefficients, expected {len(want)}")
            return 1
        for k, (a, b) in enumerate(zip(got, want)):
            d = abs(a - b) / max(abs(a), abs(b), 1.0)
            if d > worst:
                worst, where = d, (key, k, a, b)
    print(f"N_phi = {o['s']}, perimeter {o['perimeter']:.6g}, right ({o['ra']}, {o['rx']:g}, {o['rv']:g}), "
          f"left ({o['la']}, {o['lx']:g}, {o['lv']:g}): {len(want)} coefficients, largest relative deviation {worst:.1e} "
          f"(allowed {tol:g})" + (f"; worst at {where[0]}[{where[1]}]: {where[2]!r}, expected {where[3]!r}" if worst > tol else ""))
    return 0 if worst <= tol else 1


def check():
    worst = 0.0
    for alpha in range(5):
        for c in (-3.0, -1.2, -0.3, 0.0, 0.4, 1.5, 3.5):
            n, top = 40000, 12.0 + abs(c)       # composite Simpson rule on [0, top]
            h = top / n
            f = [(i * h) ** alpha * math.exp(-(i * h + c) ** 2) for i in range(n + 1)]
            direct = h / 3 * (f[0] + f[n] + 4 * sum(f[1:n:2]) + 2 * sum(f[2:n:2]))
            worst = max(worst, abs(direct - k_alpha(alpha, c)) / max(abs(direct), 1e-300))
    print(f"K_alpha(c), alpha = 0 to 4, seven values of c: largest relative deviation from direct integration {worst:.1e}")
    return 0 if worst < 1e-10 else 1


if __name__ == "__main__":
    if sys.argv[1:] == ["--check"]:
        sys.exit(check())
    if len(sys.argv) > 4 and sys.argv[1] == "run":
        sys.exit(run(sys.argv[2:]))
    print(__doc__)
    sys.exit(2)
