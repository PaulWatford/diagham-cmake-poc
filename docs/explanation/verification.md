# What is verified, what is only watched, what is untested

Purpose: one honest page on how far the test suite goes — what a green run proves, what it merely detects, and what it says nothing about.
Source: the r4493 build of 2026-09-28, `docs/reference/test-coverage.md` (generated), `docs/reference/tests.md`, `docs/reference/known-defects.md`. Numbers are those of the `default` preset; the `full` and `lapack` presets add one cross-check and one regression each.

## The four levels

| Level | Tests | What a pass means | What it does not mean |
|---|---:|---|---|
| **Physics golden** (`physics`) | 191 (+5 with GSL) | the program reproduced a number known *without* DiagHam: exact counting, a closed form, or an independent solver in this repository | nothing about options or regimes the golden does not exercise |
| **Cross-check** (`crosscheck`) | 3 (+1 with LAPACK) | two DiagHam implementations, or two algorithms (Lanczos vs full diagonalisation, LAPACK vs internal), agree | both could be wrong the same way |
| **Regression** (`regression`) | 11 (+1 with GSL) | the program reproduces what this build produced on 2026-09-28, to 10⁻⁹ | *nothing* about correctness — a wrong number reproduced is a pass |
| **Smoke** (`smoke`) | 599 | the program links, parses `--help` and exits 0 | nothing about physics |

Plus 20 self-tests of the harness (the manifest, twelve oracle re-derivations, the checker's and runner's rejection of wrong values, missing files, NaN, wrong-length spectra and missing output), 6 `known-bug` reproducers (`WILL_FAIL`), 1 install test and the Python cross-checks.

## Programs, by module

| Module directory | Programs | Physics golden | Cross-check | Regression | Any of the three |
|---|---:|---:|---:|---:|---:|
| Analysis | 19 | 0 | 0 | 0 | 0 |
| FCI | 56 | 1 | 0 | 2 | 2 |
| FQHEOnCylinder | 35 | 3 | 0 | 0 | 3 |
| FQHEOnDisk | 27 | 3 | 0 | 0 | 3 |
| FQHEOnLattice | 22 | 0 | 0 | 1 | 1 |
| FQHEOnSphere | 193 | 18 | 1 | 1 | 18 |
| FQHEOnTorus | 64 | 19 | 3 | 1 | 19 |
| FTI | 36 | 1 | 0 | 1 | 2 |
| HubbardModels | 25 | 5 | 0 | 2 | 6 |
| QuantumDots | 19 | 0 | 0 | 1 | 1 |
| Spin | 71 | 15 | 0 | 2 | 15 |
| core (`src/Programs`) | 36 | 1 | 0 | 0 | 1 |
| **all** | **603** | **66** | **4** | **11** | **71** |

So: **66 programs have their physics checked** (66 against independent answers), **5 more are watched for change**, and **528 are smoke-tested only**. A green suite says the 71 are right or unchanged; it says nothing about the other 528.

## What the physics goldens establish

- **Hilbert spaces** (27 tests): the sphere, torus (magnetic-translation sectors), disk, single-band lattice, Hubbard and spin-chain classes enumerate the right number of states in the right sectors — exact combinatorics, including the fermionic sign structure of torus translations.
- **Model states** (43 tests): the Laughlin, Moore–Read and Read–Rezayi states are the unique zero modes of their parent interactions on the sphere, torus, cylinder and disk, with the exact quasihole and edge-mode countings ((k,r)-admissible configurations); the Jack-polynomial construction and exact diagonalisation give the same state (overlap 1); the particle entanglement spectrum of the Laughlin state has exactly the quasihole counting.
- **Interactions** (14 tests): the sphere Coulomb pseudopotentials equal the closed form from Wigner algebra to 10⁻¹⁵, in the lowest Landau level and (monopole-harmonic form, l = S + n) in the first and second; for any pseudopotentials, fermions or bosons, the two-body sphere programs reproduce the full spectrum of an independent exact diagonalisation (7 cases).
- **Chains and lattices** (49 tests): XX, AKLT and Haldane–Shastry rings against closed forms; Heisenberg (spin ½ and 1), XXZ, J₁–J₂ and transverse-Ising chains against independent dense diagonalisation; Hubbard at U = 0 against the tight-binding sum and at U = 4 on 2×2 against −4√2; the plain-Lanczos path (`-n 1`) gives every momentum sector of the 2×4 lattice its tight-binding energy (the U30 case, fixed here); the XYZ chain (open, and periodic with its parity-twisted boundary bond), the fully generic open chain with three fields, the 2D Heisenberg and transverse-field Ising models on periodic clusters, the J₁–J₂ chain in three momentum sectors, the double-triangle chain, the spin-2 generalised AKLT model (P₃ + P₄) in two momentum sectors, the O'Brien–Fendley program with its default factors, the three-state Potts chain, and (with GSL) a disordered chain from a fixed field file, all whole sectors against `tests/oracles/spin_models.py` (models given as explicit operator lists, dense ED, Sz and momentum projectors); the Hubbard model on a 2×3 cluster at U = 2, the Haldane-Hubbard model on a 3×3 honeycomb cluster (U = 0, U = 2, with a staggered potential), the SSH chain and the Kitaev–Heisenberg honeycomb spin model against `tests/oracles/lattice_fermions.py` and `spin_models.py`.
- **Spinful sphere interactions** (15 tests): the same sphere oracle with two, three and four species (SU(2), SU(3), SU(4) spin, species-dependent pseudopotentials, unpolarised and polarised sectors), including the Sz ↔ −Sz and Lz ↔ −Lz parity sectors of the symmetrised bases and the legacy `QHEFermionsSphereWithSpin` program in its own normalisation; whole spectra to 10⁻⁹.
- **Cylinder and disk interactions** (4 tests, +4 with GSL): the delta and Laplacian-delta cylinder programs, and with GSL the cylinder pseudopotential program at three aspect ratios and two momentum sectors, reproduce an independent Landau-gauge exact diagonalisation (`tests/oracles/geometry_ed.py`) to 10⁻¹²; the disk Coulomb pseudopotentials equal Γ(m+½)/(2 m!) to 10⁻¹⁵. The disk pseudopotential Hamiltonian itself is **wrong** for N ≥ 3 (U34 below): its three known-bug tests carry the exact spectra.
- **k-body interactions** (8 tests): three- and four-body pseudopotential Hamiltonians on the sphere (hard core, generic three-body sets, with and without a two-body term, bosons and fermions, N = 5 and 6) against an independent ED built from the L-multiplet projectors of the k-particle space (`tests/oracles/sphere_nbody.py`); DiagHam's per-multiplet normalisation of its k-body pseudopotentials is read off the k-particle runs, every N > k spectrum is then predicted to 10⁻⁹.
- **Torus interactions** (25 tests): the Coulomb and pseudopotential spectra of fermions and bosons on the torus, whole K_y sectors at aspect ratio 1 and 2, against an independent Landau-gauge exact diagonalisation (`tests/oracles/torus_ed.py`), to 10⁻⁹; with SU(2), SU(3) and SU(4) spin, bilayer Coulomb at a layer separation, the first excited Landau level, every K_y sector at once, and the magnetic-translation (K_x, K_y) sectors of the translation programs, spinless and with spin.
- **FCI bands** (6 tests): both bands of the checkerboard lattice model at every lattice momentum, at and away from the flat-band hoppings, against the published Bloch Hamiltonian (`tests/oracles/fci_bands.py`), to 10⁻¹⁵.

What they do **not** establish: Coulomb energies on the cylinder (the matrix elements need a regularised q = 0 term, not yet in the oracle); any interaction energy on the disk beyond zero modes (the program is wrong, U34, and the fermionic one hangs, U29); spin on the cylinder and disk beyond dimension counting (the sphere and torus with spin are covered); the two-Landau-level programs on the sphere and torus, `FQHESphereFermionsWithSpinFull` (four-index pseudopotential keys, convention not matched), the twisted torus and the product geometries (T²×T², T²×S², T²×cylinder); the FCI/FTI many-body programs (regression only; the checkerboard single-particle bands are verified); Landau-level mixing, real-space entanglement, Monte Carlo, MPS/DMRG, time evolution; every MPI and ScaLAPACK code path (compiled and linked in the `hpc` preset, never run by a test); QuantumDots (regression only); the Analysis tools (need input spectra).

## Defects the suite found

| | Program | Status |
|---|---|---|
| U28 | cylinder pseudopotential Hamiltonian crashed without GSL (undefined behaviour) | fixed here; report in `docs/upstream-reports/` |
| U29 | `FQHEDiskFermionsTwoBodyGeneric` hangs (N ≥ 3; segfaults at N = 2 for Lz ≥ 5) | open; reproducer report |
| U33 | `FQHECylinderFermionsCoulomb` segfaulted in builds without GSL | fixed here (fails clearly, as U28); report 17 |
| U34 | the disk pseudopotential Hamiltonian gives wrong energies for N ≥ 3 (zero modes right) | open; three `known-bug` tests with the exact spectra; report `open-U34.md` |
| U35 | `TwoDimensionalTransverseFieldIsingModel` adds the wrap-around Ising bonds whether or not `--use-periodic` is given (its "open" spectrum is the periodic one) | open; `known-bug` test; reproducer `open-U35.md` |
| U30 | plain Lanczos (`-n 1`) returned an energy below the ground state on a degenerate spectrum: the Krylov space had closed and the iteration went on with round-off | **fixed here** (closure detected, iteration stops); report 16 in `docs/upstream-reports/`; test `physics.hubbard.2x4.U0.plain_lanczos_krylov_closure` |
| U31 | `FCICheckerboardLatticeModel` two-band model segfaults | open; `known-bug` test; use `--single-band` |
| U32 | `FQHESphereFermionsWithSpin` with both `--lzsymmetrized-basis` and `--szsymmetrized-basis` builds an empty basis | open; `known-bug` test; use one symmetry at a time |

The full register (31 entries, most of them build-level) is `docs/reference/known-defects.md`.

## How to read a failure

A red `physics` test after a change means the change altered a number that is known independently — read the failure before touching anything; a golden is never edited to make a test pass. A red `regression` test means the output changed — which may be a fix (regenerate the reference and say why in the commit) or a break. A red `selftest` means the harness itself is broken. A `known-bug` test turning red is good news: the defect is gone; promote it.

## How this will move

The coverage page is generated from the build and CI refuses a stale one, so these numbers cannot silently drift. Every prompt of the verification cycle raised them; the next candidates are listed in `docs/reference/test-coverage.md` (the smoke-only programs by module) — the sphere programs with spin and the FCI programs are the largest uncovered families.
