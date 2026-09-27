#!/usr/bin/env python3
"""Zero-mode counts for the torus, the cylinder and the disk, without DiagHam.

    python3 tests/oracles/fqhe_geometries.py --write tests/data/fqhe_geometries
    python3 tests/oracles/fqhe_geometries.py --check tests/data/fqhe_geometries  (selftest.fqhe_geometries_oracle)

Same principle as tests/oracles/fqhe_sphere.py: the zero-energy states of
the (k+1)-body hard-core interaction are counted by the (k,r)-admissible
occupation configurations (at most k particles in r consecutive orbitals;
r = 2 for bosons, r = k+2 for fermions), the geometry only changing which
configurations exist:
- torus: N_phi orbitals on a ring, the admissibility window wraps around;
  the count is the total over all momentum sectors (a program run with
  --redundant-kymomenta, or with magnetic translations, lists every sector);
- cylinder: N_phi+1 orbitals on a line (same counting as the sphere); the
  programs list the sectors 2Ky >= 0 only, so the count is over those;
- disk: orbitals m = 0, 1, 2, ... with a fixed total Lz; the zero modes at
  Lz = Lz_0 + dL are the edge excitations, p(dL) of them for dL <= N
  (the partitions of dL), which the admissible count reproduces.
"""
import argparse, sys
from pathlib import Path


def _rec_count(n, norb, k, r, fermion, accept, weight=None):
    """Enumerate occupations of n particles on norb orbitals with the window rule
    applied along the line; `accept(occ)` filters at the end; returns {weight(occ): count}."""
    max_occ = 1 if fermion else k
    out = {}

    def rec(j, left, occ):
        if j == norb:
            if left == 0 and accept(occ):
                w = weight(occ) if weight else 0
                out[w] = out.get(w, 0) + 1
            return
        for m in range(0, min(max_occ, left) + 1):
            if sum(occ[max(0, j - r + 1):j]) + m <= k:
                rec(j + 1, left - m, occ + [m])
    rec(0, n, [])
    return out


def torus_total(n, nphi, k, r, fermion):
    """All (k,r)-admissible configurations on the ring of nphi orbitals."""
    def wraps_ok(occ):
        return all(sum(occ[(i + t) % nphi] for t in range(r)) <= k for i in range(nphi))
    return sum(_rec_count(n, nphi, k, r, fermion, wraps_ok).values())


def cylinder_nonnegative_ky(n, nphi, k, r, fermion):
    """Configurations on nphi+1 orbitals (2m = -nphi..nphi) with 2Ky >= 0."""
    norb = nphi + 1
    counts = _rec_count(n, norb, k, r, fermion, lambda occ: True,
                        weight=lambda occ: sum(o * (2 * i - (norb - 1)) for i, o in enumerate(occ)))
    return sum(c for w, c in counts.items() if w >= 0)


def disk_count(n, lz, k, r, fermion):
    """Configurations on orbitals m = 0..lz with total angular momentum lz."""
    counts = _rec_count(n, lz + 1, k, r, fermion, lambda occ: True,
                        weight=lambda occ: sum(o * i for i, o in enumerate(occ)))
    return counts.get(lz, 0)


CASES = [
    # name, geometry, N, size, k, r, fermion, note
    ("torus_bosons_laughlin_1_2_n4_nphi8", "torus", 4, 8, 1, 2, False, "bosonic Laughlin nu=1/2 at Nphi=2N: 2-fold topological degeneracy"),
    ("torus_bosons_laughlin_1_2_n5_nphi10", "torus", 5, 10, 1, 2, False, "same at N=5"),
    ("torus_bosons_laughlin_1_2_n4_nphi9", "torus", 4, 9, 1, 2, False, "one quasihole: Nphi zero modes"),
    ("torus_fermions_laughlin_1_3_n4_nphi12", "torus", 4, 12, 1, 3, True, "Laughlin nu=1/3 at Nphi=3N: 3-fold topological degeneracy"),
    ("torus_fermions_laughlin_1_3_n4_nphi13", "torus", 4, 13, 1, 3, True, "one quasihole: Nphi zero modes"),
    ("torus_fermions_laughlin_1_3_n4_nphi14", "torus", 4, 14, 1, 3, True, "two quasiholes"),
    ("torus_bosons_mooreread_n4_nphi4", "torus", 4, 4, 2, 2, False, "bosonic Moore-Read nu=1 at Nphi=N: 3-fold degeneracy"),
    ("torus_bosons_mooreread_n6_nphi6", "torus", 6, 6, 2, 2, False, "same at N=6"),
    ("torus_fermions_mooreread_n4_nphi8", "torus", 4, 8, 2, 4, True, "fermionic Moore-Read nu=1/2 at Nphi=2N: 6-fold degeneracy"),
    ("cylinder_bosons_laughlin_1_2_n4_l6", "cylinder", 4, 6, 1, 2, False, "unique zero mode at Nphi=2(N-1) (sectors 2Ky>=0 listed)"),
    ("cylinder_bosons_laughlin_1_2_n5_l8", "cylinder", 5, 8, 1, 2, False, "same at N=5"),
    ("cylinder_bosons_laughlin_1_2_n4_l7", "cylinder", 4, 7, 1, 2, False, "one quasihole"),
    ("cylinder_bosons_laughlin_1_2_n4_l8", "cylinder", 4, 8, 1, 2, False, "two quasiholes"),
    ("cylinder_bosons_mooreread_n4_l2", "cylinder", 4, 2, 2, 2, False, "bosonic Moore-Read, three-body delta, unique zero mode"),
    ("cylinder_bosons_mooreread_n6_l4", "cylinder", 6, 4, 2, 2, False, "same at N=6"),
    ("cylinder_bosons_mooreread_n6_l5", "cylinder", 6, 5, 2, 2, False, "one extra flux"),
    ("disk_bosons_laughlin_1_2_n4_lz12", "disk", 4, 12, 1, 2, False, "Lz_0 = N(N-1): the Laughlin state"),
    ("disk_bosons_laughlin_1_2_n4_lz13", "disk", 4, 13, 1, 2, False, "dL=1: p(1)=1 edge mode"),
    ("disk_bosons_laughlin_1_2_n4_lz14", "disk", 4, 14, 1, 2, False, "dL=2: p(2)=2"),
    ("disk_bosons_laughlin_1_2_n4_lz15", "disk", 4, 15, 1, 2, False, "dL=3: p(3)=3"),
    ("disk_bosons_laughlin_1_2_n4_lz16", "disk", 4, 16, 1, 2, False, "dL=4: p(4)=5"),
    ("disk_bosons_laughlin_1_2_n5_lz20", "disk", 5, 20, 1, 2, False, "N=5 Laughlin state"),
    ("disk_bosons_laughlin_1_2_n5_lz22", "disk", 5, 22, 1, 2, False, "N=5, dL=2: p(2)=2"),
]


def count(case):
    _, geom, n, size, k, r, fermion, _ = case
    if geom == "torus":
        return torus_total(n, size, k, r, fermion)
    if geom == "cylinder":
        return cylinder_nonnegative_ky(n, size, k, r, fermion)
    return disk_count(n, size, k, r, fermion)


def render():
    lines = ["# zero-mode counts: name  count  -- (k,r)-admissible configurations; torus = all sectors, cylinder = sectors 2Ky>=0, disk = one Lz",
             "# generated by tests/oracles/fqhe_geometries.py; do not edit"]
    for case in CASES:
        name, geom, n, size, k, r, fermion, note = case
        lines.append(f"{name} {count(case)}  # {geom} N={n} {'Nphi' if geom != 'disk' else 'Lz'}={size} ({k},{r}) {'fermions' if fermion else 'bosons'}: {note}")
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", metavar="DIR")
    ap.add_argument("--check", metavar="DIR")
    a = ap.parse_args()
    text = render()
    if a.write:
        Path(a.write).mkdir(parents=True, exist_ok=True)
        (Path(a.write) / "zero_modes.txt").write_text(text)
        print(f"wrote {len(CASES)} cases to {a.write}/zero_modes.txt")
    if a.check:
        p = Path(a.check) / "zero_modes.txt"
        ok = p.exists() and p.read_text() == text
        print("zero_modes.txt", "matches" if ok else "MISMATCH")
        sys.exit(0 if ok else 1)
    if not a.write and not a.check:
        print(text)


if __name__ == "__main__":
    main()
