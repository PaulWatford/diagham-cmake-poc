#!/usr/bin/env python3
"""Independent entanglement-entropy oracles.

    python3 tests/oracles/entanglement.py --write tests/data/entanglement
    python3 tests/oracles/entanglement.py --check tests/data/entanglement      (selftest.entanglement_oracle)

Toric code (2DToricCodeEntanglementEntropy). The program's ground state is the uniform superposition
over the coset of the GF(2) span V of the vertex masks
    A_{x,y} = {(x,y,0), (x,y,1), (x-1,y,0), (x,y-1,1)}     (edge (x,y,s) -> qubit s + 2 (y + Ny x))
on the torus. A uniform superposition over an affine subspace is a stabilizer state, so for any
bipartition A|B its entropy is (rank P_A V + rank P_B V - rank V) ln 2, P_A the restriction of the
vectors to the qubits of A (linear algebra over GF(2), nothing shared with DiagHam; the offset of the
coset, i.e. --gs-parity, does not enter). A region of n_A x n_B sites with both its edges is
a rectangle; a kept-sites list picks single qubits.

Particle-partition counting on the torus (FQHETorusEntanglementEntropyParticlePartition): for the
Laughlin state at filling 1/m (N = 4, Nphi = 12, m = 3, the three degenerate ground states) the
Ky-resolved number of non-zero levels of the reduced density matrix of N_A particles equals the
number of cyclically admissible N_A-particle configurations (no two particles within m - 1
orbitals) with total momentum Ky mod Nphi: 1 per sector for N_A = 1 and 3,4,3,4,... (42 in all)
for N_A = 2.
"""
import itertools, math, struct, sys
from pathlib import Path


def qubit(x, y, s, nx, ny):
    return s + 2 * ((y % ny) + ny * (x % nx))


def toric_entropy(nx, ny, region):
    """entropy / ln 2 of the qubit set `region`: rank P_A V + rank P_B V - rank V over GF(2)"""
    vs = []
    for x in range(nx):
        for y in range(ny):
            v = 0
            for q in (qubit(x, y, 0, nx, ny), qubit(x, y, 1, nx, ny), qubit(x - 1, y, 0, nx, ny), qubit(x, y - 1, 1, nx, ny)):
                v ^= 1 << q
            vs.append(v)
    mask_a = sum(1 << q for q in region)
    mask_b = ((1 << (2 * nx * ny)) - 1) & ~mask_a
    return gf2_rank_int([v & mask_a for v in vs]) + gf2_rank_int([v & mask_b for v in vs]) - gf2_rank_int(vs)


def gf2_rank_int(vals):
    basis = {}
    rank = 0
    for v in vals:
        while v:
            p = v.bit_length() - 1
            if p in basis:
                v ^= basis[p]
            else:
                basis[p] = v
                rank += 1
                break
    return rank


def rectangle(nx, ny, xa, ya):
    return {qubit(x, y, s, nx, ny) for x in range(xa) for y in range(ya) for s in (0, 1)}


TORIC = [
    # name, nx, ny, kind, spec
    ("3x3_A2x2", 3, 3, "rect", (2, 2)),
    ("3x3_A1x1", 3, 3, "rect", (1, 1)),
    ("4x3_A2x2", 4, 3, "rect", (2, 2)),
    ("3x4_A2x3", 3, 4, "rect", (2, 3)),
    ("3x3_A1x3", 3, 3, "rect", (1, 3)),
]


def toric_cases():
    for name, nx, ny, kind, spec in TORIC:
        if kind == "rect":
            reg = rectangle(nx, ny, *spec)
        else:
            reg = {qubit(x, y, s, nx, ny) for x, y, s in spec}
        yield name, nx, ny, kind, spec, toric_entropy(nx, ny, reg)


def admissible_counts(nphi, n_a, m):
    """number of n_a-subsets of the ring Z_nphi with no m consecutive orbitals holding two particles
    (the Laughlin 1/m generalized Pauli rule, cyclically), per total momentum mod nphi"""
    c = [0] * nphi
    for occ in itertools.combinations(range(nphi), n_a):
        ok = all(((b - a) % nphi) >= m and ((a - b) % nphi) >= m for a, b in itertools.combinations(occ, 2))
        if ok:
            c[sum(occ) % nphi] += 1
    return c


def render_all():
    out = {}
    for name, nx, ny, kind, spec, s in toric_cases():
        out[f"toric_{name}.txt"] = f"# Entangement entropy = {s}\n"
        if kind == "kept":
            out[f"toric_{name}_kept.txt"] = "".join(f"{x} {y} {i}\n" for x, y, i in spec)
    for n_a in (1, 2):
        pc = admissible_counts(12, n_a, 3)
        out[f"torus_pes_laughlin_n4_2s12_na{n_a}.txt"] = "".join(f"{k} {v}\n" for k, v in enumerate(pc))
    # spin-1/2 chains (SpinChainEntanglementEntropy, LAPACK build): binary vector = int32 dimension + doubles
    r = 1 / math.sqrt(2)
    out["spin_1_2_n_2_sz_0.0.vec"] = struct.pack("<i2d", 2, r, -r)
    out["spin_1_2_n_4_sz_4.0.vec"] = struct.pack("<id", 1, 1.0)
    return out


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in ("--write", "--check"):
        print(__doc__)
        return 2
    d = Path(sys.argv[2])
    files = render_all()
    if sys.argv[1] == "--write":
        d.mkdir(parents=True, exist_ok=True)
        for n, t in files.items():
            if isinstance(t, bytes):
                (d / n).write_bytes(t)
            else:
                (d / n).write_text(t)
        return 0
    bad = 0
    for n, t in files.items():
        if not (d / n).is_file():
            print("missing:", n)
            bad += 1
        elif isinstance(t, bytes):
            if (d / n).read_bytes() != t:
                print("differs:", n)
                bad += 1
        elif (d / n).read_text().replace(chr(13) + chr(10), chr(10)) != t:
            print("differs:", n)
            bad += 1
    # closed-form cross-checks independent of the linear algebra above
    if toric_entropy(3, 3, rectangle(3, 3, 3, 3)) != 0:
        print("whole system must have zero entropy"); bad += 1
    if toric_entropy(3, 3, set()) != 0:
        print("empty region must have zero entropy"); bad += 1
    if toric_entropy(3, 3, rectangle(3, 3, 2, 2)) != 6:
        print("2x2 on 3x3 must be 6"); bad += 1
    print("ok" if not bad else f"{bad} problems")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
