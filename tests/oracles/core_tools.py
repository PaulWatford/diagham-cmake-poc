#!/usr/bin/env python3
"""Independent answers for core tools that read or list without a physical model of their own.

    python3 tests/oracles/core_tools.py --write tests/data/core_tools
    python3 tests/oracles/core_tools.py --check tests/data/core_tools    (selftest.core_tools_oracle, needs numpy)

- GenericHamiltonianDiagonalization reads a matrix as a text list of elements, one per line, every element
  listed (both triangles): "value column row" for a real matrix, "Re Im column row" for a complex one
  (--data-column 0, the default; 0-based indices). The matrices here are given by closed formulas and written
  with 17 significant digits, so the program and numpy diagonalise exactly the same doubles; the reference is
  numpy.linalg.eigvalsh.
- FTIShowBasis lists the n-body basis of a sector of total momentum: one line per state, the occupied
  single-particle orbitals (kx,ky) or, with several bands, (kx,ky,band), where the program writes the band as
  + and - (fermions) or up and down (bosons) for two bands. The independent answer enumerates every set (fermions) or multiset (bosons) of
  orbitals whose momenta add up to (Kx, Ky) modulo the lattice, and the two sets of states must be equal;
  no tolerance. Run as: core_tools.py showbasis PROGRAM ARGS... (the test runs the program itself).
- FQHESphereBosonEntanglementEntropy cuts a state on the sphere into the first l_a orbitals and the rest and
  writes, for each l_a, the von Neumann entropy (natural logarithm) and, with --density-matrix, every
  eigenvalue of the reduced density matrix. The independent answer: the Hamiltonian of the case is built by
  tests/oracles/sphere_ed.py (nothing shared with DiagHam), its ground state in the Lz = 0 sector must be
  unique (checked), and the Schmidt decomposition of that vector over (occupations of the first l_a orbitals,
  occupations of the rest) gives the eigenvalues and the entropy. An L = 0 state is unchanged by Lz -> -Lz,
  so it does not matter which pole the program starts the cut from. The test computes the ground state with
  FQHESphereBosonsTwoBodyGeneric and the entanglement with the program under test, and compares both outputs
  with the stored answer: core_tools.py entanglement GROUND_PROGRAM ENT_PROGRAM DATA_DIR CASE WORK_DIR TOLERANCE.
- FQHESphereUnnormalizeState converts a sphere state between the normalised Fock basis and the unnormalised
  (monomial) basis of the polynomial. The bosonic Laughlin state is prod_{i<j} (z_i - z_j)^2, whose monomial
  coefficients are exact integers, expanded here. Two checks, neither assuming the program's normalisation
  conventions: (a) --normalize applied to those integers must give a normalised state that the V0 Hamiltonian
  of tests/oracles/sphere_ed.py annihilates (an exact zero mode), and (b) unnormalising that zero mode, as
  computed by numpy, must give a vector proportional to the integers.
- FQHESphereWithSU2SpinRotation rotates a spinful boson state about the y axis. A fully polarised input state
  (any amplitudes; here a fixed formula) has an exact image: each boson's spin turns to cos(t/2) up - sin(t/2)
  down, so an orbital with n up-bosons goes to amplitude sqrt(C(n, k)) cos(t/2)^(n-k) (-sin(t/2))^k on (n-k up,
  k down). The sign is the program's sense of rotation, exp(+i t S_y) (the active rotation by -t), the only
  convention taken from the program, fixed by these tests; the help text does not state it. The rotation
  acts orbital by orbital, so which pole is orbital 0 does not matter.
Both run the programs on vectors this script writes in DiagHam's basis order, read from FQHESphereShowBasis at
test time, so no DiagHam diagonalisation enters: core_tools.py unnormalize|su2rotation SHOWBASIS PROGRAM ...
- HubbardKitaevHeisenbergModel on the single hexagon of its wiki page (geometry file: site, site, bond type
  0 = x, 1 = y, 2 = z), two cases that separate its terms: (a) isotropic hopping t = 1 and U = 2, no Gutzwiller
  projection, 4 electrons: the Hubbard model on a 6-site ring, every (N_up, N_down) sector, built by
  tests/oracles/lattice_fermions.py; (b) the Gutzwiller projection at half filling, where no hopping is possible:
  the spin model J1 sum sigma_i.sigma_j + J2 sum (sigma^a sigma^a - sigma^b sigma^b - sigma^c sigma^c) on each
  bond of type a, built by tests/oracles/spin_models.py. The couplings multiply Pauli matrices (4 S.S), not the
  S.S of the wiki's formula (U47); that factor is taken from the program, everything else is independent. The
  anisotropic hopping t' is not covered.
- FTIEntanglementEntropyParticlePartition: the particle-partition entanglement of momentum-space fermion states
  (one band). The input states are written here, with generic complex amplitudes, in DiagHam's basis order read
  from FTIShowBasis at test time (the program always reads complex vectors). For N_A particles the reduced density
  matrix is rho_A = Tr_B |psi><psi| / C(N, N_A), block by block in (N_A, total momentum of A), with the fermion sign
  of splitting the occupied set; with --degenerated-groundstate it is the equal-weight average over the states.
  Entropies (natural logarithm) and every eigenvalue are compared. Run as: core_tools.py ftientanglement
  SHOWBASIS PROGRAM WORK_DIR TOLERANCE N NX NY KX,KY [KX,KY ...]  (several sectors = a degenerate ground state).
- MatrixExtractColumns: a binary matrix (DiagHam's layout: a 4-byte type, 1 real or 2 complex, the numbers of rows and
  columns, then the elements row by row, little-endian doubles) written here with values that need every bit; the
  extracted columns must come back bit-exact as binary vectors (a 4-byte dimension, then the doubles, real and
  imaginary parts interleaved for a complex vector). Run as: core_tools.py matrixcolumns PROGRAM WORK_DIR real|complex
  NROW NCOL FIRST NBR.
"""
import argparse, math, sys
from pathlib import Path


def real_dense(n):
    """Real symmetric, every element non-zero, no degeneracy."""
    return [[(math.cos(0.37 * i + 0.11) * 3.0 if i == j else math.sin(0.5 * (i + 1) * (j + 1)) / (1 + abs(i - j)))
             for j in range(n)] for i in range(n)]


def real_sparse(n):
    """Real symmetric, banded (Lanczos-sized): a chain in a linear potential with a longer-range coupling, so that
    the lowest eigenvalues are well separated (a Wannier-Stark ladder) and Lanczos converges on each of them."""
    h = {}
    for i in range(n):
        h[(i, i)] = 0.05 * i + 0.2 * math.sin(1.3 * i)
        for d, a in ((1, -1.0), (7, 0.3)):
            if i + d < n:
                v = a * (1 + 0.1 * math.sin(0.9 * i + d))
                h[(i, i + d)] = h[(i + d, i)] = v
    return h


def complex_dense(n):
    """Hermitian, complex off-diagonal elements."""
    m = [[0j] * n for _ in range(n)]
    for i in range(n):
        m[i][i] = complex(2.0 * math.cos(0.45 * i + 0.3), 0.0)
        for j in range(i + 1, n):
            v = complex(math.cos(0.3 * i * j + 0.2), math.sin(0.17 * (i + 2 * j))) / (1 + 0.5 * (j - i))
            m[i][j] = v
            m[j][i] = v.conjugate()
    return m


CASES = {  # name -> (builder, kind)
    "real_dense_40": (lambda: real_dense(40), "real"),
    "real_sparse_600": (lambda: real_sparse(600), "sparse"),
    "complex_dense_30": (lambda: complex_dense(30), "complex"),
}


def matrix_lines(name):
    build, kind = CASES[name]
    m = build()
    if kind == "real":
        return [f"{m[i][j]!r} {j} {i}" for i in range(len(m)) for j in range(len(m))]
    if kind == "sparse":
        return [f"{v!r} {j} {i}" for (i, j), v in sorted(m.items())]
    return [f"{m[i][j].real!r} {m[i][j].imag!r} {j} {i}" for i in range(len(m)) for j in range(len(m))]


def spectrum(name):
    import numpy as np
    build, kind = CASES[name]
    m = build()
    if kind == "sparse":
        n = 1 + max(i for i, _ in m)
        a = np.zeros((n, n))
        for (i, j), v in m.items():
            a[i, j] = v
    else:
        a = np.array(m)
    return np.linalg.eigvalsh(a)


def render(name):
    return {f"{name}.txt": "\n".join(matrix_lines(name)) + "\n",
            f"{name}_spectrum.dat": f"# eigenvalues of {name}.txt (numpy.linalg.eigvalsh); generated by tests/oracles/core_tools.py, do not edit\n"
                                    + "".join(f"{e:.17g}\n" for e in spectrum(name))}


KH_BONDS = [(0, 1, 2), (0, 2, 1), (1, 3, 0), (2, 4, 0), (3, 5, 1), (4, 5, 2)]     # the wiki's single hexagon


def kh_spectra():
    """Every eigenvalue of the two HubbardKitaevHeisenbergModel cases."""
    import numpy as np
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import lattice_fermions, spin_models
    hop = [[0.0] * 6 for _ in range(6)]
    for i, j, _ in KH_BONDS:
        hop[i][j] = hop[j][i] = -1.0
    ring = np.sort(np.concatenate([lattice_fermions.hubbard_spectrum(hop, 2.0, nu, 4 - nu) for nu in range(5)]))
    terms = []
    for i, j, t in KH_BONDS:
        for b in "xyz":
            terms.append((4 * (0.3 + 0.5 * (1 if b == "xyz"[t] else -1)), [(i, b), (j, b)]))
    spin = np.sort(np.linalg.eigvalsh(spin_models.build(6, 1, terms)))
    return {"kh_hubbard_ring_n4_u2": ring, "kh_spin_j1_0.3_j2_0.5": spin}


def render_kh():
    out = {"kh_hexagon.dat": "".join(f"{i} {j} {t}\n" for i, j, t in KH_BONDS)}
    for name, ev in kh_spectra().items():
        out[f"{name}_spectrum.dat"] = (f"# every eigenvalue of {name} (see tests/oracles/core_tools.py); generated, do not edit\n"
                                       + "".join(f"{e:.17g}\n" for e in ev))
    return out


def ftishowbasis_states(n, nx, ny, kx, ky, bands=1, boson=False):
    """Every n-particle state of total momentum (kx, ky) mod (nx, ny), as a sorted tuple of orbitals."""
    import itertools
    orbitals = [(x, y, b) for x in range(nx) for y in range(ny) for b in range(bands)]
    pick = itertools.combinations_with_replacement if boson else itertools.combinations
    return {s for s in pick(orbitals, n) if sum(o[0] for o in s) % nx == kx % nx and sum(o[1] for o in s) % ny == ky % ny}


def parse_ftishowbasis(text):
    """The states FTIShowBasis prints: lines such as [(0,1)(1,1)] or [(1,1,+)(1,1,-)]."""
    import re
    sign = {"+": 0, "-": 1, "up": 0, "down": 1}
    states = []
    for line in text.splitlines():
        line = line.strip()
        if not (line.startswith("[") and line.endswith("]")):
            continue
        orbs = []
        for t in re.findall(r"\(([^)]*)\)", line):
            f = t.split(",")
            band = (sign[f[2]] if f[2] in sign else int(f[2])) if len(f) > 2 else 0
            orbs.append((int(f[0]), int(f[1]), band))
        states.append(tuple(sorted(orbs)))
    return states


def run_ftishowbasis(argv):
    """Run FTIShowBasis with the given options and compare its basis with the enumeration."""
    import subprocess
    prog, args = argv[0], argv[1:]
    opt = {"-p": 4, "-x": 3, "-y": 3, "--kx": 0, "--ky": 0, "-s": 1}
    i = 0
    while i < len(args):
        if args[i] in opt:
            opt[args[i]] = int(args[i + 1])
            i += 2
        else:
            i += 1
    boson = "--boson" in args
    out = subprocess.run([prog] + args, capture_output=True, text=True, timeout=120)
    if out.returncode != 0:
        print(out.stdout[-2000:], out.stderr[-2000:])
        print(f"FTIShowBasis exited with {out.returncode}")
        return 1
    got = parse_ftishowbasis(out.stdout)
    want = ftishowbasis_states(opt["-p"], opt["-x"], opt["-y"], opt["--kx"], opt["--ky"], opt["-s"], boson)
    dims = [int(w) for w in out.stdout.split() if w.isdigit()]
    problems = []
    if len(got) != len(set(got)):
        problems.append("the program lists a state twice")
    if set(got) != want:
        problems.append(f"{len(set(got) - want)} listed states are not in the sector, {len(want - set(got))} states of the sector are missing")
    if len(want) not in dims:
        problems.append(f"the printed dimension is not {len(want)}")
    print(f"FTIShowBasis {' '.join(args)}: {len(got)} states listed, {len(want)} expected")
    for p in problems:
        print("  " + p)
    return 1 if problems else 0


def _coulomb(two_s):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import sphere_ed
    return [float(x) for x in sphere_ed.coulomb(two_s)]


ENT_CASES = {  # name -> (N, 2S, pseudopotentials V_0, V_1, ...)
    "bosons_laughlin_n4_2s6": (4, 6, [1.0]),
    "bosons_v0_v2_n5_2s6": (5, 6, [1.0, 0.0, 0.5]),
    "bosons_coulomb_n4_2s6": (4, 6, None),          # the Coulomb pseudopotentials of tests/oracles/sphere_ed.py
}


def ent_case(name):
    n, two_s, pp = ENT_CASES[name]
    return n, two_s, (_coulomb(two_s) if pp is None else pp)


def entanglement_reference(name):
    """Ground energy, gap, and per l_a the entropy and the sorted non-zero reduced density matrix eigenvalues."""
    import numpy as np
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import sphere_ed
    n, two_s, pp = ent_case(name)
    H, states = sphere_ed.hamiltonian(n, two_s, pp, False, 0)
    e, v = np.linalg.eigh(H)
    psi = v[:, 0]
    out = {"e0": float(e[0]), "gap": float(e[1] - e[0]), "cuts": {}}
    for la in range(1, (two_s + 1) // 2 + 1):
        rows, cols = {}, {}
        entries = []
        for amp, st in zip(psi, states):
            a = tuple(j for j in st if j < la)
            b = tuple(j for j in st if j >= la)
            entries.append((rows.setdefault(a, len(rows)), cols.setdefault(b, len(cols)), amp))
        m = np.zeros((len(rows), len(cols)))
        for i, j, amp in entries:
            m[i, j] += amp
        lam = np.linalg.svd(m, compute_uv=False) ** 2
        lam = sorted(float(x) for x in lam if x > 1e-14)
        out["cuts"][la] = (float(-sum(x * math.log(x) for x in lam)), lam)
    return out


def render_entanglement(name):
    n, two_s, pp = ent_case(name)
    r = entanglement_reference(name)
    head = (f"# {name}: N = {n}, 2S = {two_s}, Lz = 0, V = {pp}; ground energy {r['e0']:.6g}, gap {r['gap']:.6g}\n"
            "# generated by tests/oracles/core_tools.py, do not edit\n")
    return {f"{name}_pp.dat": head + "Pseudopotentials = " + " ".join(repr(x) for x in pp) + "\n",
            f"{name}_entropy.dat": head + "# l_a entropy\n" + "".join(f"{la} {S:.17g}\n" for la, (S, _) in r["cuts"].items()),
            f"{name}_dm.dat": head + "# l_a lambda (every eigenvalue above 1e-14, sorted)\n"
                              + "".join(f"{la} {x:.17g}\n" for la, (_, lam) in r["cuts"].items() for x in lam)}


def run_entanglement(argv):
    """Ground state with DiagHam, entanglement with the program under test, both compared with the stored answer."""
    import glob, os, subprocess
    ground, ent, data, name, work, tol = argv
    ground, ent, data, work = (os.path.abspath(x) for x in (ground, ent, data, work))
    rel = float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)
    n, two_s, pp = ent_case(name)
    os.makedirs(work, exist_ok=True)
    for f in glob.glob(os.path.join(work, "*")):
        os.remove(f)
    cmd = [ground, "-p", str(n), "-l", str(two_s), "--interaction-file", os.path.join(data, f"{name}_pp.dat"),
           "--interaction-name", name, "--nbr-lz", "1", "--eigenstate", "-n", "1", "--full-diag", "5000"]
    r = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=300)
    vec = glob.glob(os.path.join(work, "*_lz_0.0.vec"))
    if r.returncode != 0 or len(vec) != 1:
        print(r.stdout[-2000:], r.stderr[-2000:])
        print(f"ground state: exit {r.returncode}, vectors {vec}")
        return 1
    r = subprocess.run([ent, os.path.basename(vec[0]), "--density-matrix", "dm.full.ent"], cwd=work, capture_output=True, text=True, timeout=300)
    entf = vec[0][:-4] + ".ent"
    if r.returncode != 0 or not os.path.exists(entf):
        print(r.stdout[-2000:], r.stderr[-2000:])
        print(f"entanglement program: exit {r.returncode}")
        return 1

    def rows(path):
        return [l.split() for l in open(path) if l.strip() and not l.startswith("#")]
    want_s = {int(a): float(b) for a, b in rows(os.path.join(data, f"{name}_entropy.dat"))}
    got_s = {int(f[0]): float(f[1]) for f in rows(entf)}
    want_l, got_l = {}, {}
    for a, x in rows(os.path.join(data, f"{name}_dm.dat")):
        want_l.setdefault(int(a), []).append(float(x))
    for f in rows(os.path.join(work, "dm.full.ent")):
        if float(f[-1]) > 1e-12:
            got_l.setdefault(int(f[0]), []).append(float(f[-1]))
    problems = []
    worst = 0.0
    if sorted(got_s) != sorted(want_s):
        problems.append(f"subsystem sizes {sorted(got_s)}, expected {sorted(want_s)}")
    for la in want_s:
        if la in got_s:
            d = abs(got_s[la] - want_s[la]) / max(abs(want_s[la]), 1.0)
            worst = max(worst, d)
            if d > rel:
                problems.append(f"l_a = {la}: entropy {got_s[la]!r}, expected {want_s[la]!r}")
        g = sorted(got_l.get(la, []))
        w = [x for x in want_l.get(la, []) if x > 1e-12]
        if len(g) != len(w):
            problems.append(f"l_a = {la}: {len(g)} non-zero density-matrix eigenvalues, expected {len(w)}")
        else:
            for x, y in zip(g, w):
                worst = max(worst, abs(x - y))
                if abs(x - y) > rel:
                    problems.append(f"l_a = {la}: eigenvalue {x!r}, expected {y!r}")
    print(f"{name}: entropies for l_a = {sorted(got_s)}, largest deviation {worst:.2e} (allowed {rel:g})")
    for p in problems:
        print("  " + p)
    return 1 if problems else 0


def showbasis_states(showbasis, args):
    """DiagHam's basis order: one occupation list per state (spinless: n_m; spinful: (n_up, n_down) per orbital)."""
    import re, subprocess
    out = subprocess.run([showbasis] + args, capture_output=True, text=True, timeout=120).stdout
    states = []
    for line in out.splitlines():
        if "|" in line and "(" in line:
            states.append(tuple(tuple(int(x) for x in t.split(",")) for t in re.findall(r"\((\d+,\d+)\)", line)))
        elif line.strip() and all(w.isdigit() for w in line.split()) and len(line.split()) > 1:
            states.append(tuple(int(w) for w in line.split()))
    return states


def laughlin_monomials(n):
    """prod_{i<j} (z_i - z_j)^2 expanded: {exponent tuple: integer coefficient}."""
    poly = {tuple([0] * n): 1}
    for i in range(n):
        for j in range(i + 1, n):
            new = {}
            for e, c in poly.items():
                for di, dj, f in ((2, 0, 1), (1, 1, -2), (0, 2, 1)):
                    k = list(e)
                    k[i] += di
                    k[j] += dj
                    k = tuple(k)
                    new[k] = new.get(k, 0) + c * f
            poly = {k: v for k, v in new.items() if v}
    return poly


def run_unnormalize(argv):
    """Laughlin 1/2 bosons: the integer monomial coefficients, normalised by the program, must be an exact zero
    mode of V0; the numpy zero mode, unnormalised by the program, must be proportional to the integers."""
    import os, subprocess
    import numpy as np
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import sphere_ed, vector_io
    showbasis, prog, n, work, tol = os.path.abspath(argv[0]), os.path.abspath(argv[1]), int(argv[2]), os.path.abspath(argv[3]), argv[4]
    rel = float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)
    two_s = 2 * (n - 1)
    os.makedirs(work, exist_ok=True)
    occ = showbasis_states(showbasis, ["-p", str(n), "-l", str(two_s), "-z", "0", "--boson"])
    mono = laughlin_monomials(n)
    ints = np.array([float(mono.get(tuple(sorted((m for m in range(two_s + 1) for _ in range(o[m])), reverse=True)), 0)) for o in occ])
    H, states = sphere_ed.hamiltonian(n, two_s, [1.0], False, 0)
    idx = {st: i for i, st in enumerate(states)}
    perm = [idx[tuple(m for m in range(two_s + 1) for _ in range(o[m]))] for o in occ]
    Hd = H[np.ix_(perm, perm)]                       # the V0 Hamiltonian in DiagHam's basis order
    problems = []
    # (a) normalise the integers
    name = f"bosons_unnormalized_laughlin_n_{n}_2s_{two_s}_lz_0.0.vec"
    vector_io.write_vec(Path(work) / name, list(ints))
    r = subprocess.run([prog, "-i", name, "-o", "normalized.vec", "--normalize"], cwd=work, capture_output=True, text=True, timeout=120)
    if r.returncode != 0 or not os.path.exists(os.path.join(work, "normalized.vec")):
        print(r.stdout[-1500:], r.stderr[-1500:]); print(f"--normalize: exit {r.returncode}"); return 1
    a = np.array(vector_io.read_vec(Path(work) / "normalized.vec"))
    norm, res = float(np.linalg.norm(a)), float(np.linalg.norm(Hd @ a))
    if abs(norm - 1) > rel or res > rel * np.abs(H).max():
        problems.append(f"--normalize: norm {norm!r}, |H v| = {res:.2e} (a zero mode needs 0)")
    # (b) unnormalise the numpy zero mode
    e, v = np.linalg.eigh(Hd)
    if e[1] - e[0] < 1e-3 or abs(e[0]) > 1e-10:
        problems.append("the V0 zero mode is not unique")
    zero = v[:, 0]
    name = f"bosons_laughlin_n_{n}_2s_{two_s}_lz_0.0.vec"
    vector_io.write_vec(Path(work) / name, list(zero))
    ref = int(np.argmax(np.abs(ints)))
    r = subprocess.run([prog, "-i", name, "-o", "unnormalized.vec", "--normalization", str(ref)], cwd=work, capture_output=True, text=True, timeout=120)
    if r.returncode != 0 or not os.path.exists(os.path.join(work, "unnormalized.vec")):
        print(r.stdout[-1500:], r.stderr[-1500:]); print(f"unnormalize: exit {r.returncode}"); return 1
    u = np.array(vector_io.read_vec(Path(work) / "unnormalized.vec"))
    want = ints / ints[ref]
    dev = float(np.max(np.abs(u - want)) / np.max(np.abs(want)))
    if dev > rel:
        problems.append(f"unnormalised zero mode differs from the integer coefficients by {dev:.2e} (relative)")
    print(f"Laughlin N = {n}, 2S = {two_s}: {len(occ)} states; --normalize gives norm {norm:.15f}, |H v| = {res:.1e}; "
          f"unnormalised vs integers: {dev:.1e} (allowed {rel:g})")
    for p_ in problems:
        print("  " + p_)
    return 1 if problems else 0


def run_su2rotation(argv):
    """Rotate a fully polarised spinful boson state about y and compare with the exact image."""
    import os, subprocess, glob
    import numpy as np
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import vector_io
    showbasis, prog, n, two_s, angle, work, tol = os.path.abspath(argv[0]), os.path.abspath(argv[1]), int(argv[2]), int(argv[3]), float(argv[4]), os.path.abspath(argv[5]), argv[6]
    rel = float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)
    os.makedirs(work, exist_ok=True)
    for f in glob.glob(os.path.join(work, "*.vec")):
        os.remove(f)
    base = ["-p", str(n), "-l", str(two_s), "-z", "0", "--boson", "--su2-spin", "--use-alt"]
    src = showbasis_states(showbasis, base + ["-s", str(n)])
    dst = showbasis_states(showbasis, base + ["--all-sz"])
    amp = np.array([math.cos(1.3 * i + 0.4) + 0.1 for i in range(len(src))])
    amp /= np.linalg.norm(amp)
    name = f"bosons_sphere_su2_test_n_{n}_2s_{two_s}_sz_{n}_lz_0.0.vec"
    vector_io.write_vec(Path(work) / name, list(amp))
    r = subprocess.run([prog, "-i", name, "-a", repr(angle)], cwd=work, capture_output=True, text=True, timeout=120)
    out = [f for f in glob.glob(os.path.join(work, "*syrotation*.vec"))]
    if r.returncode != 0 or len(out) != 1:
        print(r.stdout[-1500:], r.stderr[-1500:]); print(f"rotation: exit {r.returncode}, outputs {out}"); return 1
    got = np.array(vector_io.read_vec(Path(out[0])))
    index = {st: i for i, st in enumerate(dst)}
    c, s_ = math.cos(math.pi * angle / 2), -math.sin(math.pi * angle / 2)    # the program's sense: exp(+i t S_y)
    want = np.zeros(len(dst))
    import itertools
    for a_, st in zip(amp, src):
        ups = [o[0] for o in st]
        for ks in itertools.product(*[range(u + 1) for u in ups]):
            w = a_
            for u, k in zip(ups, ks):
                w *= math.sqrt(math.comb(u, k)) * c ** (u - k) * s_ ** k
            want[index[tuple((u - k, k) for u, k in zip(ups, ks))]] += w
    dev = float(np.max(np.abs(got - want))) if len(got) == len(want) else float("inf")
    print(f"rotation by {angle} pi, N = {n}, 2S = {two_s}: {len(src)} -> {len(dst)} states; norm {np.linalg.norm(got):.15f}; "
          f"largest deviation from the exact image {dev:.1e} (allowed {rel:g})")
    if dev > rel:
        flip = np.array([w * (-1) ** sum(t[1] for t in st) for w, st in zip(want, dst)])
        print(f"  differs; with the opposite rotation sense the deviation would be {float(np.max(np.abs(got - flip))):.1e}")
        return 1
    return 0


def write_cvec(path, values):
    """DiagHam's binary complex vector: a 4-byte dimension, then (real, imaginary) little-endian doubles."""
    import struct
    with open(path, "wb") as f:
        f.write(struct.pack("<i", len(values)))
        for z in values:
            f.write(struct.pack("<dd", z.real, z.imag))


def particle_partition(states_amps, n, nx, ny):
    """{N_A: {K_A: sorted eigenvalues}} of the averaged particle-partition reduced density matrix."""
    import itertools
    import numpy as np
    order = {(x, y): x * ny + y for x in range(nx) for y in range(ny)}
    res = {}
    for na in range(1, n):
        rho = {}
        for states, amp in states_amps:
            blocks = {}
            for s, a in zip(states, amp):
                for A in itertools.combinations(s, na):
                    Bp = tuple(m for m in s if m not in A)
                    sign = (-1) ** sum(1 for x in A for y in Bp if order[y] < order[x])
                    kA = (sum(m[0] for m in A) % nx, sum(m[1] for m in A) % ny)
                    blocks.setdefault(kA, {}).setdefault(Bp, {})[A] = sign * a
            for kA, byB in blocks.items():
                r = rho.setdefault(kA, {})
                for d in byB.values():
                    for A1, v1 in d.items():
                        for A2, v2 in d.items():
                            r[(A1, A2)] = r.get((A1, A2), 0) + v1 * np.conj(v2) / (math.comb(n, na) * len(states_amps))
        lam = {}
        for kA, r in rho.items():
            rows = sorted({A for A, _ in r})
            idx = {A: i for i, A in enumerate(rows)}
            m = np.zeros((len(rows), len(rows)), dtype=complex)
            for (A1, A2), v in r.items():
                m[idx[A1], idx[A2]] += v
            lam[kA] = sorted(float(x) for x in np.linalg.eigvalsh(m) if x > 1e-14)
        res[na] = lam
    return res


def run_ftientanglement(argv):
    import os, glob, subprocess
    import numpy as np
    showbasis, prog, work, tol = os.path.abspath(argv[0]), os.path.abspath(argv[1]), os.path.abspath(argv[2]), argv[3]
    n, nx, ny = int(argv[4]), int(argv[5]), int(argv[6])
    sectors = [tuple(int(v) for v in s.split(",")) for s in argv[7:]]
    rel = float(tol.split(":", 1)[1]) if tol.startswith("rel:") else float(tol)
    os.makedirs(work, exist_ok=True)
    for f in glob.glob(os.path.join(work, "*")):
        os.remove(f)
    states_amps, names = [], []
    for j, (kx, ky) in enumerate(sectors):
        out = subprocess.run([showbasis, "-p", str(n), "-x", str(nx), "-y", str(ny), "--kx", str(kx), "--ky", str(ky)],
                             capture_output=True, text=True, timeout=120).stdout
        states = [tuple(sorted((o[0], o[1]) for o in s)) for s in parse_ftishowbasis(out)]
        amp = np.array([complex(math.cos(0.9 * i + 0.3 + j), 0.6 * math.sin(1.7 * i + 0.2 - j)) for i in range(len(states))])
        amp /= np.linalg.norm(amp)
        name = f"fermions_singleband_oracle_n_{n}_x_{nx}_y_{ny}_kx_{kx}_ky_{ky}.0.vec"
        write_cvec(Path(work) / name, list(amp))
        states_amps.append((states, amp))
        names.append(name)
    if len(names) == 1:
        cmd = [prog, names[0], "--density-matrix", "dm.full.parent"]
    else:
        Path(work, "ground.dat").write_text("\n".join(names) + "\n")
        cmd = [prog, "--degenerated-groundstate", "ground.dat", "--density-matrix", "dm.full.parent"]
    r = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=600)
    ents = glob.glob(os.path.join(work, "*.partent"))
    if r.returncode != 0 or len(ents) != 1 or not os.path.exists(os.path.join(work, "dm.full.parent")):
        print(r.stdout[-1500:], r.stderr[-1500:]); print(f"program exit {r.returncode}, {ents}"); return 1
    ent = {int(l.split()[0]): float(l.split()[1]) for l in open(ents[0]) if l.strip() and not l.startswith("#")}
    got = {}
    for l in open(os.path.join(work, "dm.full.parent")):
        w = l.split()
        if w and not l.startswith("#") and float(w[3]) > 1e-12:
            got.setdefault(int(w[0]), {}).setdefault((int(w[1]), int(w[2])), []).append(float(w[3]))
    ref = particle_partition(states_amps, n, nx, ny)
    problems, worst = [], 0.0
    if not ent:
        problems.append("no entropy written")
    for na in ent:
        S = -sum(x * math.log(x) for v in ref[na].values() for x in v)
        d = abs(S - ent[na]) / max(abs(S), 1.0)
        worst = max(worst, d)
        if d > rel:
            problems.append(f"N_A = {na}: entropy {ent[na]!r}, expected {S!r}")
        for kA in set(ref[na]) | set(got.get(na, {})):
            g, w = sorted(got.get(na, {}).get(kA, [])), [x for x in ref[na].get(kA, []) if x > 1e-12]
            if len(g) != len(w):
                problems.append(f"N_A = {na}, K_A = {kA}: {len(g)} non-zero eigenvalues, expected {len(w)}")
                continue
            for x, y in zip(g, w):
                worst = max(worst, abs(x - y))
                if abs(x - y) > rel:
                    problems.append(f"N_A = {na}, K_A = {kA}: eigenvalue {x!r}, expected {y!r}")
    print(f"N = {n}, {nx} x {ny}, sectors {sectors}: N_A = {sorted(ent)}, largest deviation {worst:.1e} (allowed {rel:g})")
    for p_ in problems[:10]:
        print("  " + p_)
    return 1 if problems else 0


def run_matrixcolumns(argv):
    import os, shutil, struct, subprocess
    prog, work, kind = os.path.abspath(argv[0]), os.path.abspath(argv[1]), argv[2]
    nrow, ncol, first, nbr = (int(x) for x in argv[3:7])
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    special = [1 / 3, -1e-300, 2.0 ** 53 + 2, -0.1, math.pi, 5e-324, 1 - 2 ** -53, -7.25]
    val = lambda i, j, part: special[(3 * i + 5 * j + part) % len(special)] * (1 + i) - (j + 0.5) * 2 ** -40 * part
    with open(os.path.join(work, "m.mat"), "wb") as f:
        f.write(struct.pack("<iii", 1 if kind == "real" else 2, nrow, ncol))
        for i in range(nrow):
            for j in range(ncol):
                f.write(struct.pack("<d", val(i, j, 0)))
                if kind == "complex":
                    f.write(struct.pack("<d", val(i, j, 1)))
    cmd = [prog, "-i", "m.mat", "-o", "col", "--first-column", str(first), "--nbr-columns", str(nbr)] + (["-c"] if kind == "complex" else [])
    r = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        print(r.stdout[-1500:], r.stderr[-1500:]); print(f"program exit {r.returncode}"); return 1
    problems = []
    for j in range(first, first + nbr):
        path = os.path.join(work, f"col.{j}.vec")
        if not os.path.exists(path):
            problems.append(f"col.{j}.vec not written ({sorted(os.listdir(work))})")
            continue
        data = open(path, "rb").read()
        want = struct.pack("<i", nrow) + b"".join(
            struct.pack("<d", val(i, j, 0)) + (struct.pack("<d", val(i, j, 1)) if kind == "complex" else b"") for i in range(nrow))
        if data != want:
            problems.append(f"col.{j}.vec differs from column {j} ({len(data)} bytes, expected {len(want)})")
    print(f"{kind} {nrow} x {ncol}, columns {first}..{first + nbr - 1}: {'bit-exact' if not problems else 'DIFFERENT'}")
    for p_ in problems:
        print("  " + p_)
    return 1 if problems else 0


def main():
    if len(sys.argv) > 2 and sys.argv[1] == "matrixcolumns":
        sys.exit(run_matrixcolumns(sys.argv[2:]))
    if len(sys.argv) > 2 and sys.argv[1] == "ftientanglement":
        sys.exit(run_ftientanglement(sys.argv[2:]))
    if len(sys.argv) > 2 and sys.argv[1] == "unnormalize":
        sys.exit(run_unnormalize(sys.argv[2:]))
    if len(sys.argv) > 2 and sys.argv[1] == "su2rotation":
        sys.exit(run_su2rotation(sys.argv[2:]))
    if len(sys.argv) > 2 and sys.argv[1] == "entanglement":
        sys.exit(run_entanglement(sys.argv[2:]))
    if len(sys.argv) > 2 and sys.argv[1] == "showbasis":
        sys.exit(run_ftishowbasis(sys.argv[2:]))
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", metavar="DIR")
    ap.add_argument("--check", metavar="DIR")
    a = ap.parse_args()
    if a.write:
        d = Path(a.write)
        d.mkdir(parents=True, exist_ok=True)
        for name in CASES:
            for f, text in render(name).items():
                (d / f).write_text(text)
        for name in ENT_CASES:
            for f, text in render_entanglement(name).items():
                (d / f).write_text(text)
        for f, text in render_kh().items():
            (d / f).write_text(text)
        print(f"wrote {2 * len(CASES) + 3 * len(ENT_CASES) + 3} files to {d}")
    if a.check:
        import numpy as np
        d = Path(a.check)
        bad = 0
        for name in CASES:
            files = render(name)
            if (d / f"{name}.txt").read_text() != files[f"{name}.txt"]:
                print(f"{name}.txt differs from the oracle's matrix"); bad += 1
            ref = [float(l) for l in (d / f"{name}_spectrum.dat").read_text().split("\n") if l and not l.startswith("#")]
            new = spectrum(name)
            dev = max(abs(x - y) / max(abs(x), abs(y), 1.0) for x, y in zip(ref, new))
            if len(ref) != len(new) or dev > 1e-13:
                print(f"{name}: stored spectrum differs (relative {dev:.2e})"); bad += 1
        for name in ENT_CASES:
            r = entanglement_reference(name)
            if r["gap"] < 1e-3:
                print(f"{name}: ground state not unique (gap {r['gap']:.2e})"); bad += 1
            for f, text in render_entanglement(name).items():
                stored = (d / f).read_text()
                if f.endswith("_pp.dat"):
                    if stored != text:
                        print(f"{f} differs from the oracle's"); bad += 1
                    continue
                vals = lambda t: [float(l.split()[1]) for l in t.splitlines() if l and not l.startswith("#")]
                a_, b_ = vals(stored), vals(text)
                if len(a_) != len(b_) or max(abs(x - y) for x, y in zip(a_, b_)) > 1e-13:
                    print(f"{f}: stored values differ from the oracle's"); bad += 1
        for f, text in render_kh().items():
            stored = (d / f).read_text()
            if f == "kh_hexagon.dat":
                if stored != text:
                    print(f"{f} differs from the oracle's"); bad += 1
                continue
            vals = lambda t: [float(l) for l in t.splitlines() if l and not l.startswith("#")]
            a_, b_ = vals(stored), vals(text)
            if len(a_) != len(b_) or max(abs(x - y) for x, y in zip(a_, b_)) > 1e-12:
                print(f"{f}: stored values differ from the oracle's"); bad += 1
        print("core_tools oracle:", "OK" if not bad else f"{bad} problems")
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
