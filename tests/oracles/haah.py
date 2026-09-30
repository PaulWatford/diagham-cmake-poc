#!/usr/bin/env python3
"""Independent oracle for HaahCodeEntropy (Spin/src/Programs/HaahCodeEntropy.cc).

    python3 tests/oracles/haah.py --write tests/data/haah
    python3 tests/oracles/haah.py --check tests/data/haah      (selftest.haah_oracle)

The program builds the state sum_{v in V} |v> / sqrt|V|, V = GF(2) span of the translates of the
eight-spin "Z term" mask (spin s at site (x,y,z) = qubit s + 2 (z + L (y + L x))). A uniform
superposition over a linear subspace is a stabilizer state, so for a region A of qubits
    S(A) / ln 2 = rank P_A V + rank P_B V - rank V        (GF(2) linear algebra, nothing shared with DiagHam).
The check mode also
  * verifies the closed form S = 6 a^2 - 6 a + 2 for an a x a x a cube on the L^3 torus (a < L),
  * (if mpmath and numpy are installed) builds the state explicitly for L = 2 (and L = 3 with the
    environment variable HAAH_FULL set, about 2 GB) with a one-site region, forms the reduced density matrix and diagonalizes it with mpmath at 40 digits,
    which shares neither the rank formula nor DiagHam's matrix code.
"""
import os
import sys
from pathlib import Path

S0 = [(1, 0, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1)]   # offsets of the spin-0 terms of one mask
S1 = [(0, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1)]   # offsets of the spin-1 terms


def qubit(x, y, z, s, n):
    nx, ny, nz = n
    return s + 2 * ((z % nz) + nz * ((y % ny) + ny * (x % nx)))


def masks(n):
    out = []
    for x in range(n[0]):
        for y in range(n[1]):
            for z in range(n[2]):
                m = 0
                for a, b, c in S0:
                    m ^= 1 << qubit(x + a, y + b, z + c, 0, n)
                for a, b, c in S1:
                    m ^= 1 << qubit(x + a, y + b, z + c, 1, n)
                out.append(m)
    return out


def rank(vals):
    basis, r = {}, 0
    for v in vals:
        while v:
            p = v.bit_length() - 1
            if p in basis:
                v ^= basis[p]
            else:
                basis[p] = v
                r += 1
                break
    return r


def box(n, ax, ay, az):
    m = 0
    for x in range(ax):
        for y in range(ay):
            for z in range(az):
                for s in (0, 1):
                    m |= 1 << qubit(x, y, z, s, n)
    return m


def entropy(n, ax, ay, az):
    vs = masks(n)
    a = box(n, ax, ay, az)
    b = ((1 << (2 * n[0] * n[1] * n[2])) - 1) & ~a
    return rank([v & a for v in vs]) + rank([v & b for v in vs]) - rank(vs)


# (nx, ny, nz, ax, ay, az): system and region sizes
CASES = [(2, 2, 2, 1, 1, 1), (2, 2, 2, 1, 1, 2), (2, 2, 2, 1, 2, 2), (2, 2, 2, 2, 1, 1), (4, 4, 4, 2, 2, 2),
         (2, 2, 3, 1, 1, 2), (2, 3, 2, 1, 2, 1), (3, 2, 2, 2, 1, 1), (4, 2, 2, 2, 2, 1)]


def name(c):
    return "haah_%dx%dx%d_A%dx%dx%d" % c


def render_all():
    out = {}
    for c in CASES:
        out[name(c) + ".txt"] = "# Entangement entropy = %d\n" % entropy(c[:3], *c[3:])
    return out


def mpmath_check():
    try:
        import numpy as np
        import mpmath as mp
    except ImportError:
        print("mpmath/numpy not installed: explicit-state check skipped")
        return 0
    mp.mp.dps = 40
    bad = 0
    for L in ((2, 3) if os.environ.get("HAAH_FULL") else (2,)):
        vs = masks((L, L, L))
        basis = []
        for v in vs:                       # reduced basis of V
            for p in basis:
                v = min(v, v ^ p)
            if v:
                basis.append(v)
        r = len(basis)
        elems = np.zeros(1 << r, dtype=np.uint64)
        for i, v in enumerate(basis):
            elems[1 << i:2 << i] = elems[:1 << i] ^ np.uint64(v)
        elems.sort()                       # region A = qubits 0, 1 = site (0,0,0), the two lowest bits
        b = elems >> np.uint64(2)
        a = (elems & np.uint64(3)).astype(np.uint8)
        del elems
        start = np.concatenate([[0], np.nonzero(b[1:] != b[:-1])[0] + 1])
        del b
        sets = np.bitwise_or.reduceat(np.left_shift(np.uint8(1), a), start)      # bit a set: (a, b) in V
        cnt = np.bincount(sets, minlength=16)
        rho = mp.zeros(4, 4)
        for m in range(1, 16):
            if cnt[m]:
                s = [1 if (m >> i) & 1 else 0 for i in range(4)]
                for i in range(4):
                    for j in range(4):
                        rho[i, j] += int(cnt[m]) * s[i] * s[j]
        rho = rho / mp.mpf(2) ** r
        tr = sum(rho[i, i] for i in range(4))
        ev = mp.eigsy(rho, eigvals_only=True)
        ent = -sum(x * mp.log(x) for x in ev if x > mp.mpf(10) ** -30)
        want = 2 * mp.log(2)
        ok = abs(tr - 1) < mp.mpf(10) ** -35 and abs(ent - want) < mp.mpf(10) ** -35
        print("mpmath L=%d rank=%d: trace-1 = %s, S - 2 ln 2 = %s, levels = %s" % (
            L, r, mp.nstr(tr - 1, 3), mp.nstr(ent - want, 3), [mp.nstr(x, 6) for x in ev]))
        if not ok:
            bad += 1
    return bad


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in ("--write", "--check"):
        print(__doc__)
        return 2
    d = Path(sys.argv[2])
    files = render_all()
    if sys.argv[1] == "--write":
        d.mkdir(parents=True, exist_ok=True)
        for n, t in files.items():
            (d / n).write_text(t)
        return 0
    bad = 0
    for n, t in files.items():
        if not (d / n).is_file():
            print("missing:", n)
            bad += 1
        elif (d / n).read_text().replace(chr(13) + chr(10), chr(10)) != t:
            print("differs:", n)
            bad += 1
    for L in range(2, 7):
        for a in range(1, L):
            if entropy((L, L, L), a, a, a) != 6 * a * a - 6 * a + 2:
                print("cube L=%d a=%d: entropy %d, closed form %d" % (L, a, entropy((L, L, L), a, a, a), 6 * a * a - 6 * a + 2))
                bad += 1
    bad += mpmath_check()
    print("ok" if not bad else "%d problems" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
