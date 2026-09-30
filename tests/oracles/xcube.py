#!/usr/bin/env python3
"""Independent oracle for XCubeEntanglementEntropy (Spin/src/Programs/XCubeEntanglementEntropy.cc).

    python3 tests/oracles/xcube.py --write tests/data/xcube
    python3 tests/oracles/xcube.py --check tests/data/xcube      (selftest.xcube_oracle)

The program builds the uniform superposition over the GF(2) span V of the Z terms of the X-cube
model: three qubits per site (qubit s at site (x,y,z) = s + 3 (z + nz (y + ny x))), one four-spin
term of each of two types per site,
    type 1: (0,0,0)s0 (0,0,0)s1 (-1,0,0)s0 (0,-1,0)s1
    type 2: (0,0,0)s1 (0,0,0)s2 (0,-1,0)s1 (0,0,-1)s2
(the relations among the terms make the program use nz (nx ny - 1) + nx (ny nz - 1) of them; the
span is the same). A uniform superposition over a linear subspace is a stabilizer state, so for a
region A of qubits
    S(A) / ln 2 = rank P_A V + rank P_B V - rank V        (GF(2) linear algebra, nothing shared with DiagHam).
The oracle was validated against the program on 8 cases in the default and LAPACK builds.
The check mode also verifies, with mpmath at 40 digits, the entropy of a one-site region of the
2x2x2 system from an explicitly built state (needs mpmath and numpy; skipped if they are missing).
"""
import sys
from pathlib import Path

T1 = [(0, 0, 0, 0), (0, 0, 0, 1), (-1, 0, 0, 0), (0, -1, 0, 1)]
T2 = [(0, 0, 0, 1), (0, 0, 0, 2), (0, -1, 0, 1), (0, 0, -1, 2)]


def qubit(x, y, z, s, n):
    nx, ny, nz = n
    return s + 3 * ((z % nz) + nz * ((y % ny) + ny * (x % nx)))


def masks(n):
    t1, t2 = [], []
    for x in range(n[0]):
        for y in range(n[1]):
            for z in range(n[2]):
                for terms, out in ((T1, t1), (T2, t2)):
                    m = 0
                    for a, b, c, s in terms:
                        m ^= 1 << qubit(x + a, y + b, z + c, s, n)
                    out.append(m)
    return t1 + t2


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
                for s in (0, 1, 2):
                    m |= 1 << qubit(x, y, z, s, n)
    return m


def entropy(n, ax, ay, az):
    vs = masks(n)
    a = box(n, ax, ay, az)
    b = ((1 << (3 * n[0] * n[1] * n[2])) - 1) & ~a
    return rank([v & a for v in vs]) + rank([v & b for v in vs]) - rank(vs)


# (nx, ny, nz, ax, ay, az): system and region sizes. 3x3x3 is the U43 case (the program would size
# its array 2^48; the rank is 46).
CASES = [(2, 2, 2, 1, 1, 1), (2, 2, 2, 1, 1, 2), (2, 2, 2, 1, 2, 2), (2, 2, 2, 2, 1, 1),
         (2, 2, 3, 1, 1, 2), (2, 3, 2, 1, 2, 1), (3, 2, 2, 2, 1, 1), (3, 2, 2, 1, 1, 1),
         (3, 3, 3, 2, 2, 2)]


def name(c):
    return "xcube_%dx%dx%d_A%dx%dx%d" % c


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
    n = (2, 2, 2)
    basis = []
    for v in masks(n):
        for p in basis:
            v = min(v, v ^ p)
        if v:
            basis.append(v)
    r = len(basis)
    elems = np.zeros(1 << r, dtype=np.uint64)
    for i, v in enumerate(basis):
        elems[1 << i:2 << i] = elems[:1 << i] ^ np.uint64(v)
    elems.sort()                           # region A = qubits 0..2 = site (0,0,0), the three lowest bits
    b = elems >> np.uint64(3)
    a = (elems & np.uint64(7)).astype(np.uint8)
    start = np.concatenate([[0], np.nonzero(b[1:] != b[:-1])[0] + 1])
    sets = np.bitwise_or.reduceat(np.left_shift(np.uint16(1), a.astype(np.uint16)), start)   # bit a set: (a, b) in V
    cnt = {}
    for m in sets.tolist():
        cnt[m] = cnt.get(m, 0) + 1
    rho = mp.zeros(8, 8)
    for m, c in cnt.items():
        s = [1 if (m >> i) & 1 else 0 for i in range(8)]
        for i in range(8):
            for j in range(8):
                rho[i, j] += c * s[i] * s[j]
    rho = rho / mp.mpf(2) ** r
    tr = sum(rho[i, i] for i in range(8))
    ev = mp.eigsy(rho, eigvals_only=True)
    ent = -sum(x * mp.log(x) for x in ev if x > mp.mpf(10) ** -30)
    want = entropy(n, 1, 1, 1) * mp.log(2)
    ok = abs(tr - 1) < mp.mpf(10) ** -35 and abs(ent - want) < mp.mpf(10) ** -35
    print("mpmath 2x2x2 rank=%d: trace-1 = %s, S - %d ln 2 = %s" % (r, mp.nstr(tr - 1, 3), entropy(n, 1, 1, 1), mp.nstr(ent - want, 3)))
    return 0 if ok else 1


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
    if rank(masks((3, 3, 3))) != 46:
        print("rank of the 3x3x3 span is not 46")
        bad += 1
    bad += mpmath_check()
    print("ok" if not bad else "%d problems" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
