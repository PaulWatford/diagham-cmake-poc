# Test suite reference

Source: `TESTING.md` (2026-09-24), moved 2026-09-26; labels `crosscheck`, `regression`, `selftest` and the coverage page added 2026-09-27; the "Adding a golden" section is now `docs/how-to/develop/add-a-golden-test.md`.

`ctest` is the test runner. The suite lives in `tests/` and is registered
whenever programs are built (`DIAGHAM_BUILD_TESTS`, default ON).

```bash
ctest --preset default                 # everything
ctest --preset physics                 # physics goldens only (seconds)
ctest --preset default -L smoke        # the per-program --help checks
ctest --preset default -R hubbard      # by name
```

## What is tested

**What the suite is, and is not.** Only a `physics` test says a program
computes the right numbers. A `crosscheck` says two DiagHam implementations
agree (both could be wrong the same way). A `regression` test says the
output has not changed since it was saved (right or wrong). A `smoke` test
says the program links and parses `--help`. The per-program picture is in
[test-coverage.md](test-coverage.md), generated from the build; the
headline today is **23 of 603 programs with a physics or cross-check test** (20 of them against an independently known answer).

| Label | Tests | What passes means |
|---|---|---|
| `physics` | goldens with independently known answers (table below) | a DiagHam program, run from a clean directory, writes a spectrum that matches an answer known without DiagHam (analytic, exact counting, an independent solver) |
| `crosscheck` | `crosscheck.fqhe.torus.su2_coulomb.*` | two DiagHam programs (or two algorithms) agree on the same problem; no independent answer is involved |
| `regression` | none yet | a spectrum saved from the r4493 build is reproduced; detects change, not correctness |
| `selftest` | `selftest.manifest` | the registered tests equal `tests/manifest.txt` (a green run that silently lost tests fails); regenerate the manifest with `python3 tests/coverage.py <build-dir> --manifest tests/manifest.txt` when a test is deliberately added or removed |
| `python` | `python.hubbard_ed_cross_check` | DiagHam's Hubbard ground state agrees with the independent Python ED in `benchmarks/hubbard_ed.py` to 1e-10, on four lattices/couplings (needs numpy; skipped at configure time without it) |
| `smoke` | `smoke.<target>`, one per program | the program starts, parses options and exits 0 on `--help` (catches link and static-initialisation breakage in the programs no golden reaches). Four programs are left out, with the reason in `tests/CMakeLists.txt`: `QHEBosons`, `MultipleSpinChain` and `TestDiagHamVectors` have no option parser; `EvaluateBroadening` checks its required `--input` before `--help` |
| `install` | `install.find_package_consumer` | the `Development` install component, used through `find_package(DiagHam)` alone from a separate CMake project, compiles, links, and diagonalises a tight-binding ring correctly |
| `known-bug` | reproducers for documented upstream defects | the defect is still present (`WILL_FAIL`); see below |

### Hilbert-space dimension goldens (`-L dimension`)

27 tests, `physics.dimension.*`, one per (program, case). The number of
basis states in a sector is exact combinatorics: `tests/oracles/dimensions.py`
enumerates the Fock states in pure Python (sphere Lz and L, torus (k_x,k_y)
sectors with magnetic translations by Burnside's lemma including the
fermion sign, single-band lattice (kx,ky), Hubbard with/without Sz, spin
chains, disk Lz) and writes `tests/data/dimensions/*.txt`; each file states
its formula. `tests/RunAndMatch.cmake` runs the program and requires either
every expected line verbatim in its output (`GetDimension` programs) or
exactly N printed basis states (`ShowBasis` programs).
`selftest.dimension_oracle` regenerates all 27 files and fails if a
committed one differs. Programs covered: `FQHESphereGetDimension`,
`FQHESphereShowBasis`, `FQHETorusGetDimension`, `FQHETorusShowBasis`,
`FQHEDiskShowBasis`, `FTIGetDimension`, `HubbardGetDimension`,
`GenericSpinChainShowBasis`. What it proves: the Hilbert-space classes
enumerate the right states in the right sectors — nothing about matrix
elements or eigenvalues.

Two things learnt writing them: `FTIGetDimension` assumes inversion
symmetry by default and prints only the kx ≤ Nx/2 sectors (`--no-inversion`
prints all); `FQHETorusShowBasis` always uses magnetic translations, so
`-x`/`-y` select one (k_x, k_y) sector, not a plain Ky sector.

### FQHE sphere goldens (`-L sphere`)

20 tests, `physics.fqhe.sphere.*`, answers from `tests/oracles/fqhe_sphere.py`
(pure Python; `selftest.fqhe_sphere_oracle` re-derives the committed files in
`tests/data/fqhe_sphere/`):

- **Zero-mode counts** (11): the zero-energy states of the (k+1)-body
  hard-core interaction are counted by (k,r)-admissible configurations
  (Read–Rezayi; Bernevig–Haldane): bosonic Laughlin ν=1/2 (V₀), bosonic
  and fermionic Moore–Read (three-body), bosonic Read–Rezayi k=3
  (four-body), fermionic Laughlin with three quasiholes — ground states,
  quasihole multiplets, a half-integer-Lz sector, N up to 9. Programs:
  `FQHESphereBosonsTwoBodyGeneric`, `FQHESphereBosonsNBodyHardCore`,
  `QHEFermionsNBodyHardCore`, `FQHESphereFermionsTwoBodyGeneric`.
- **Coulomb pseudopotentials** (3): `CoulombPseudopotentials` at 2S = 6,
  10, 15 against the closed form from Wigner 3j/6j algebra (Fano, Ortolani
  & Colombo 1986), every V_m to 1e-13 (the oracle agrees with the program
  to 8e-16).
- **Jack polynomial = exact diagonalisation** (3): `FQHESphereJackGenerator`
  (α = −2 Laughlin fermions and bosons, α = −3 bosonic Moore–Read) →
  `FQHESphereConvertHaldaneBasis` → `GenericOverlap` against the ED zero
  mode: squared overlap 1 within 1e-12. Two independent constructions of
  the same state.
- **Particle entanglement spectrum** (3):
  `FQHESphereFermionEntanglementEntropyParticlePartition` on the N=6,
  2S=15 Laughlin state: the number of non-zero levels for N_A = 2 (91),
  N_A = 3 (220) and in the (N_A = 3, 2Lz_A = −1) sector (15) equals the
  quasihole counting of N_A particles (Sterdyniak, Regnault & Bernevig
  2011).

Chains (eigenstate → Jack → overlap; eigenstate → PES) use the `PRE_STEPS`
option of `diagham_physics_test`; the checker gained the modes `line`
(pseudopotential lines) and `nonzero` (level counts per sector).

A convention worth knowing: options that `--help` lists with no short
letter and that name an input file (`--ground-file`, `--input-file`) are
registered with option code `'\0'`, which the option manager treats as
*the first positional argument* — `program state.vec …`, not
`program --ground-file state.vec`. The help text is misleading; the
manuals adapted from the wiki show the positional form.

### Physics goldens and cross-checks

| Test | Program | Check | Why the answer is known |
|---|---|---|---|
| `physics.hubbard.2x2.U4.ground_state_is_minus_4sqrt2` | `HubbardSquareLatticeModel` | lowest eigenvalue within 4 ulp of -5.6568542494923806 | analytic -4√2 (the locked golden of `AGENTS.md`); 0 ulp with GCC 13 |
| `physics.hubbard.2x2.U0.ground_state_is_minus_8` | `HubbardSquareLatticeModel` | lowest eigenvalue = -8 ± 1e-12 | free fermions, filled band |
| `physics.hubbard.2x2.U4.full_spectrum` | `HubbardSquareLatticeModel` | all 36 eigenvalues vs `benchmarks/` reference, 1e-11 | saved, independently cross-checked output |
| `physics.hubbard.2x4.U4.sector_ground_states` | `HubbardSquareLatticeModel` | lowest eigenvalue of each of 8 momentum sectors vs `benchmarks/` reference, 1e-10 | saved, independently cross-checked output |
| `physics.fqhe.sphere.laughlin_1_3.unique_zero_mode` | `FQHESphereFermionsTwoBodyGeneric` | exactly 1 zero eigenvalue at N=6, 2S=15, Lz=0 | Laughlin state is the unique V1 zero mode at 2S=3(N-1) |
| `physics.fqhe.sphere.laughlin_1_3.two_quasihole_count` | `FQHESphereFermionsTwoBodyGeneric` | exactly 4 zero eigenvalues at N=6, 2S=17, Lz=0 | two bosonic quasiholes of angular momentum N/2=3 give L=6,4,2,0 |
| `physics.fqhe.torus.laughlin_1_3.threefold_degeneracy` | `FQHETorusFermionsTwoBodyGeneric` | exactly 3 zero eigenvalues over all Ky at N=4, Nphi=12 | topological degeneracy of Laughlin 1/3 on the torus |
| `crosscheck.fqhe.torus.su2_coulomb.polarized_N3_matches_spinless` | `FQHETorusFermionsWithSpin` vs `FQHETorusFermionsCoulomb` | full 84-state spectra agree to 1e-10 | full polarisation must reproduce the spinless problem; regression for patch 14 (fails by 0.25 without it) |
| `crosscheck.fqhe.torus.su2_coulomb.unpolarized_N3_matches_translations` | `FQHETorusFermionsWithSpin` vs `FQHETorusFermionsWithSpinAndTranslations` | full 324-state Sz=1 spectra agree to 1e-8 | two independent implementations; regression for patch 14 (fails by 0.43 without it) |
| `physics.spin.heisenberg_ring_L4` | `GenericPeriodicSpinChain` | lowest eigenvalue = -2 ± 1e-12 | S=1/2 Heisenberg ring, exact |
| `physics.spin.heisenberg_ring_L6` | `GenericPeriodicSpinChain` | lowest eigenvalue = -(2+√13)/2 ± 1e-12 | exact (closed form), cross-checked with numpy |
| `physics.spin.heisenberg_ring_L10` | `GenericPeriodicSpinChain` | lowest eigenvalue = -4.5154463544920365 ± 1e-10 | independent numpy ED (Lanczos path; dimension 252) |

Tolerances are chosen so that compiler and FMA differences pass and a real
change in the physics (which moves the 12th digit or earlier) fails.
`tests/check_spectrum.cc` does the comparison; it needs no Python and no
DiagHam library, so the goldens run anywhere the programs do.

### Known upstream defects

There are no `known-bug` tests at the moment: the one defect that was
going to need one (the spinful torus Coulomb bug) is fixed by patch 14 and
covered by ordinary physics tests. The mechanism stays for the next one.
A `known-bug` test reproduces a defect documented in this repository and is
registered with `WILL_FAIL`: it passes while the defect is present, and the
day someone fixes the defect it fails, so it can't be forgotten -- turn it
into a `physics` test (drop `WILL_FAIL`) in the same change as the fix.

## Adding a golden

See [../how-to/develop/add-a-golden-test.md](../how-to/develop/add-a-golden-test.md).

## Relationship to `cmake/verify_build.sh`

`verify_build.sh` compares a CMake build against an autotools build of the
same tree (library coverage and `nm` symbol counts). It needs an autotools
build to compare against and is not part of `ctest`; the ctest suite does
not replace it, it checks different things (behaviour, not build parity).
It reports 79 passed / 0 failed in the default configuration (FQHE and FTI
on); with the default preset the ctest suite has 663 tests at r4493.
