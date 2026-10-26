#!/usr/bin/env python3
"""Independent answer for FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis (pure Python).

    python3 tests/oracles/trs_generate_basis.py run PROGRAM WORK_DIR TOLERANCE CASE [degeneracies]

The program combines single-layer spectra (one "lz E" file per particle number N_layer, as the pairing programs write
them) into the spectrum of two decoupled layers with N_up = N_down (2 Sz = 0) and 2Lz_up - 2Lz_down = 0 (2 Ky = 0):
    E = e_up + e_down + Ec (N_up + N_down - N0)^2,
and keeps an effective basis of the lowest levels. The single-layer spectra are written here (made-up levels with
exact degeneracies), so everything is known:
- --write-fullspectrum: every two-layer level (N, 2Sz, 2Lz, N_up, N_down, 2Lz_up, 2Lz_down, E), as a multiset;
- the .basis file: at most --nbr-states states, the lowest levels, each with E - Ec and the indices of its two
  single-layer states in the .vectors file, whose names give N_layer, 2Lz and the level index; the energies those
  indices point to must add up to E - Ec.
- with "degeneracies", also the wiki's rule that the cut is moved down so that no degenerate level is split ("the
  actual number of kept states can be lower to preserve the spectrum degeneracies"). U60: the program tests the
  degeneracy on the unsorted energy array, so the cut can split a level.
"""
import os, re, shutil, subprocess, sys

PREFIX = "fermions_cylinder_perimeter_6.000000_su2_quasiholes_oracle_cenergy_0.000000_n0_0.000000_pairing"

CASES = {  # name -> (2S, largest N_layer, Ec, N0, nbr-states, level rule)
    "2s6_n3_ec0.3_n0_3_keep20": (6, 3, 0.3, 3.0, 20, 0),
    "2s6_n3_ec0.05_n0_4_keep35": (6, 3, 0.05, 4.0, 35, 1),
    "2s8_n3_ec0.2_n0_4_keep50": (8, 3, 0.2, 4.0, 50, 0),
}


def single_layer(two_s, n, rule):
    """{2Lz: [energies, ascending]}: made-up levels, with the same energy in mirror sectors and some equal gaps."""
    if n == 0:
        return {0: [0.0]}
    lzmax = n * (two_s + 1 - n)
    out = {}
    for lz in range(-lzmax, lzmax + 1, 2):
        k = 1 + (abs(lz) // 2 + rule) % 3
        out[lz] = sorted(round(0.37 * n * n + 0.11 * abs(lz) + 0.5 * a + 0.013 * lz * lz + (0.2 if a == 1 else 0.0) + 0.07 * rule * n, 10)
                         for a in range(k))
    return out


def rel_of(tol):
    return float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)


def run(argv):
    prog, work, tol, name = os.path.abspath(argv[0]), os.path.abspath(argv[1]), rel_of(argv[2]), argv[3]
    strict = len(argv) > 4 and argv[4] == "degeneracies"
    two_s, nmax, ec, n0, keep, rule = CASES[name]
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    layers = {n: single_layer(two_s, n, rule) for n in range(nmax + 1)}
    files = []
    for n, spec in layers.items():
        f = f"{PREFIX}_n_{n}_2s_{two_s}_sz_{n}_lz.dat"
        with open(os.path.join(work, f), "w") as h:
            h.write("# lz E\n" + "".join(f"{lz} {e!r}\n" for lz, es in spec.items() for e in es))
        files.append(f)
    with open(os.path.join(work, "spectra.dat"), "w") as h:
        h.write("\n".join(files) + "\n")
    r = subprocess.run([prog, "--singlelayer-spectra", "spectra.dat", "--charging-energy", repr(ec), "--average-nbrparticles", repr(n0),
                        "--nbr-states", str(keep), "--write-fullspectrum"], cwd=work, capture_output=True, text=True, timeout=300)
    names = os.listdir(work)
    full = [x for x in names if x.endswith("_fullspectrum.dat")]
    basis = [x for x in names if x.endswith(".basis")]
    vectors = [x for x in names if x.endswith(".vectors")]
    if r.returncode != 0 or len(full) != 1 or len(basis) != 1 or len(vectors) != 1:
        print(r.stdout[-800:], r.stderr[-400:])
        print(f"program exit {r.returncode}, outputs {sorted(names)}")
        return 1
    # the two-layer levels
    want = []
    for n, spec in layers.items():
        for lz, es in spec.items():
            for a, ea in enumerate(es):
                for b, eb in enumerate(es):
                    want.append((ea + eb + ec * (2 * n - n0) ** 2, ea + eb, n, lz, a, b))
    want.sort()
    problems, worst = [], 0.0

    def dev(x, y):
        return abs(x - y) / max(abs(x), abs(y), 1.0)
    got_full = []
    for line in open(os.path.join(work, full[0])):
        if line.strip() and not line.startswith("#"):
            t = line.split()
            got_full.append((float(t[7]), int(t[0]), int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5]), int(t[6])))
    if len(got_full) != len(want):
        problems.append(f"full spectrum: {len(got_full)} levels, expected {len(want)}")
    for (e, N, sz, lz, nu, nd, lu, ld), w in zip(sorted(got_full), want):
        worst = max(worst, dev(e, w[0]))
        if (N, sz, lz, nu, nd, lu, ld) != (2 * w[2], 0, 0, w[2], w[2], w[3], w[3]) and abs(e - w[0]) < 1e-12:
            pass   # levels of equal energy may come in any order; the labels are checked as multisets below
    lab_got = sorted((round(e, 9), N, lu) for e, N, sz, lz, nu, nd, lu, ld in got_full)
    lab_want = sorted((round(w[0], 9), 2 * w[2], w[3]) for w in want)
    if lab_got != lab_want:
        problems.append("full spectrum: the (E, N, 2Lz_up) labels differ")
    # the effective basis: the lowest levels, the cut not splitting a degenerate level
    k = min(keep, len(want))
    while 0 < k < len(want) and abs(want[k][0] - want[k - 1][0]) < 1e-12:
        k -= 1
    vec = {}
    for line in open(os.path.join(work, vectors[0])):
        t = line.split()
        if len(t) == 2:
            m = re.search(r"_n_(\d+)_2s_\d+_sz_\d+_lz_(-?\d+)\.(\d+)\.vec$", t[1])
            vec[int(t[0])] = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
    got_basis = []
    for line in open(os.path.join(work, basis[0])):
        if line.strip() and not line.startswith("#"):
            t = line.split()
            got_basis.append((float(t[7]), float(t[8]), int(t[0]), int(t[3]), int(t[5]), int(t[6]), int(t[9]), int(t[10])))
    if strict and len(got_basis) != k:
        problems.append(f"effective basis: {len(got_basis)} states, expected {k} (--nbr-states {keep}, degeneracies kept whole)")
    if len(got_basis) > keep or len(got_basis) == 0:
        problems.append(f"effective basis: {len(got_basis)} states for --nbr-states {keep}")
    eb = sorted(round(x[0], 9) for x in got_basis)
    ew = sorted(round(w[0], 9) for w in want[:len(got_basis)])
    if eb != ew:
        problems.append("effective basis: not the lowest levels")
    for e, e0, N, nu, lu, ld, iu, idn in got_basis:
        worst = max(worst, dev(e - e0, ec * (N - n0) ** 2))
        for idx, lz in ((iu, lu), (idn, ld)):
            if idx not in vec or vec[idx][0] != nu or vec[idx][1] != lz:
                problems.append(f"effective basis: index {idx} does not name a state with N = {nu}, 2Lz = {lz}")
        if iu in vec and idn in vec:
            ea = layers[vec[iu][0]][vec[iu][1]][vec[iu][2]]
            eb_ = layers[vec[idn][0]][vec[idn][1]][vec[idn][2]]
            worst = max(worst, dev(ea + eb_, e0))
    if worst > tol:
        problems.append(f"largest relative deviation {worst:.2e} above {tol:g}")
    print(f"{name}: full spectrum {len(got_full)} levels (expected {len(want)}), effective basis {len(got_basis)} states "
          f"(expected {k}); largest relative deviation {worst:.1e} (allowed {tol:g})")
    for p_ in problems[:10]:
        print("  " + p_)
    return 1 if problems else 0


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "run":
        sys.exit(run(sys.argv[2:]))
    print(__doc__)
    sys.exit(2)
