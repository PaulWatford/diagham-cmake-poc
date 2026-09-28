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

**What the suite is, and is not** (the one-page version is [../explanation/verification.md](../explanation/verification.md)). Only a `physics` test says a program
computes the right numbers. A `crosscheck` says two DiagHam implementations
agree (both could be wrong the same way). A `regression` test says the
output has not changed since it was saved (right or wrong). A `smoke` test
says the program links and parses `--help`. The per-program picture is in
[test-coverage.md](test-coverage.md), generated from the build; the
headline today is **55 of 603 programs with a physics or cross-check test** (55 of them against an independently known answer).

| Label | Tests | What passes means |
|---|---|---|
| `physics` | goldens with independently known answers (table below) | a DiagHam program, run from a clean directory, writes a spectrum that matches an answer known without DiagHam (analytic, exact counting, an independent solver) |
| `crosscheck` | `crosscheck.fqhe.torus.su2_coulomb.*`, `crosscheck.fqhe.sphere.*` | two DiagHam programs, or two algorithms (Lanczos vs full diagonalisation, LAPACK vs internal), agree on the same problem; no independent answer is involved |
| `regression` | `regression.*` (11 programs) | a spectrum saved from the r4493 build is reproduced; detects change, not correctness |
| `selftest` | `selftest.manifest`, `selftest.*_oracle`, `selftest.checker.*`, `selftest.runner.*` | the harness works: the registered tests equal `tests/manifest.txt` (a green run that silently lost tests fails); regenerate the manifest with `python3 tests/coverage.py <build-dir> --manifest tests/manifest.txt` when a test is deliberately added or removed; tests that exist only with an optional library are listed in `tests/manifest-optional.txt` and may be present or absent |
| `python` | `python.hubbard_ed_cross_check` | DiagHam's Hubbard ground state agrees with the independent Python ED in `benchmarks/hubbard_ed.py` to 1e-10, on four lattices/couplings (needs numpy; skipped at configure time without it) |
| `smoke` | `smoke.<target>`, one per program | the program starts, parses options and exits 0 on `--help` (catches link and static-initialisation breakage in the programs no golden reaches). Four programs are left out, with the reason in `tests/CMakeLists.txt`: `QHEBosons`, `MultipleSpinChain` and `TestDiagHamVectors` have no option parser; `EvaluateBroadening` checks its required `--input` before `--help` |
| `install` | `install.find_package_consumer` | the `Development` install component, used through `find_package(DiagHam)` alone from a separate CMake project, compiles, links, and diagonalises a tight-binding ring correctly |
| `known-bug` | `knownbug.fci.checkerboard.two_band_model_segfault` (U31), `knownbug.fqhe.sphere.fermions_with_spin.lz_sz_symmetrized_basis_empty` (U32), `knownbug.fqhe.disk.bosons_two_body_generic.*` (U34, three cases) | the defect is still present (`WILL_FAIL`); see below |

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

### Torus, cylinder and disk goldens (`-L torus`, `-L cylinder`, `-L disk`)

22 tests, answers from `tests/oracles/fqhe_geometries.py` (re-derived by
`selftest.fqhe_geometries_oracle`): the same (k,r)-admissible counting as
on the sphere, with the geometry changing which configurations exist —
periodic windows on the torus (total over every momentum sector), a line
of Nφ+1 orbitals on the cylinder (sectors 2Ky ≥ 0), orbitals m ≥ 0 at
fixed total Lz on the disk.

- **Torus** (8): bosonic Laughlin ν=½ two-fold degeneracy (N=4, 5) and one
  quasihole; fermionic Laughlin ⅓ with one and two quasiholes (13 and 35
  zero modes: (Nφ/(N+n))·C(N+n,N)); bosonic Moore–Read three-fold (N=4, 6)
  and fermionic Moore–Read six-fold degeneracy with the three-body
  hard-core / hollow-core programs (`--all-points`: by default they compute
  only a reduced set of momentum sectors).
- **Cylinder** (7): `FQHECylinderBosonsDeltaInteraction` and
  `…ThreeBodyDeltaInteraction`: Laughlin and Moore–Read zero modes,
  quasiholes.
- **Disk** (7): `FQHEDiskBosonsTwoBodyGeneric` with V₀: the Laughlin state
  at Lz₀ = N(N−1) and its edge excitations at Lz₀+ΔL, p(ΔL) = 1, 1, 2, 3, 5
  of them for ΔL = 0…4 (N=4) and N=5.

Two defects met on the way, both in the register: `FQHECylinderFermionsTwoBodyGeneric`
crashed with "stack smashing detected" in every run of a build without
GSL (U28, fixed — it now says it needs GSL); `FQHEDiskFermionsTwoBodyGeneric`
hangs after printing `start` for every case tried, so there is no fermionic
disk golden yet (U29, open, reproducer in the register).

### Spin chains and Hubbard goldens (`-L spin`, `-L hubbard`)

16 physics tests, values and tolerances from `tests/oracles/spin_hubbard.py`
(`tests/data/spin_hubbard/values.txt`, re-derived by
`selftest.spin_hubbard_oracle`, which needs numpy):

- **Closed forms**: XX rings L = 6, 8, 10 (free fermions, −Σ of the lowest
  cos k with the right periodic/antiperiodic momenta); AKLT rings L = 6, 8
  (E₀ = −2L/3 exactly); Haldane–Shastry rings L = 6, 8, 10
  (E₀ = −(π²/24)(L + 5/L)); the open AKLT chain as a sum of projectors has
  exactly four zero-energy edge states (three listed, the 2Sz ≥ 0 sectors);
  Hubbard at U = 0 on 2×4, 4×2 and 3×2 by full diagonalisation equals the
  tight-binding sum (−16, −16, −12).
- **Independent dense diagonalisation in numpy** (no symmetry, nothing
  shared with DiagHam): spin-½ Heisenberg ring L = 8, spin-1 Heisenberg
  ring L = 6, transverse-field Ising ring L = 8 (`SpinChainXYZ`, whose
  convention — Pauli matrices for the couplings, S = σ/2 for the field —
  was established by matching and is stated in the oracle).

Programs: `GenericPeriodicSpinChain`, `PeriodicSpinChainAKLT`,
`SpinChainAKLT`, `HaldaneShastrySpinChain`, `SpinChainXYZ`,
`HubbardSquareLatticeModel`.

A defect found on the way (**U30**) and since fixed here: plain Lanczos
(`-n 1`, no reorthogonalisation) on the U = 0 2×4 Hubbard case returned
−18.41 and −17.87 in two momentum sectors, *below* the exact ground energy
−16. On this degenerate spectrum the Krylov space of the start vector
closes after 7 to 9 steps; the recurrence went on with the normalised
round-off and produced Ritz values outside the spectrum. `BasicLanczosAlgorithm`
and `ComplexBasicLanczosAlgorithm` now detect the closure and stop
(`docs/upstream-reports/`, report 16). The former `known-bug` test is
`physics.hubbard.2x4.U0.plain_lanczos_krylov_closure`: every one of the
eight momentum sectors must give the tight-binding minimum of that sector,
from `tight_binding_square_sectors` in `tests/oracles/spin_hubbard.py`
(`tests/data/spin_hubbard/hubbard_2x4_N8_U0_sectors.dat`).

**QuantumDots** has no golden yet: its programs are continuum
finite-difference solvers parameterised in Ångström and effective masses,
with no exactly solvable case exposed at the command line; they get
`regression` tests (prompt 7) and the coverage page says so.

### Independent solvers (`-L ed`), algorithm consensus, and self-tests

- **Sphere pseudopotential ED** (8 tests, `physics.fqhe.sphere.ed.*`, one of them with the first-Landau-level Coulomb pseudopotentials):
  `tests/oracles/sphere_ed.py` builds the lowest-Landau-level Hamiltonian
  for any Haldane pseudopotentials from scratch — Clebsch–Gordan pair
  projectors from exact Wigner 3j symbols, an explicit Fock basis of one
  Lz sector, dense numpy diagonalisation, nothing shared with DiagHam —
  and the program's **whole spectrum** of that sector must match it to
  10⁻⁹: Coulomb (the closed-form pseudopotentials of this repository, not
  the program's) and generic pseudopotential sets, fermions and bosons,
  N = 4–6, Lz = 0 and 2. The one convention taken from the program is its
  normalisation H = 2 Σ_{i<j} Σ_m V_m P_m, fixed with the two-particle case.
  While writing it the oracle disagreed with the program twice; both times
  the oracle was wrong (operator order; summing over particles instead of
  orbitals for bosons) — the disagreement→check loop the plan asks for.
- **Cylinder and disk ED** (4 tests in every build, 4 more with GSL,
  `physics.fqhe.cylinder.ed.*`, `physics.fqhe.disk.coulomb_pseudopotentials.*`):
  `tests/oracles/geometry_ed.py` carries the torus matrix element over to
  an infinite cylinder (no periodic images, the q_x sum an integral; the
  programs' circumference is L = √(2π r N_orb) for their `-r`, and they
  label momentum sectors by twice the momentum, `-y 2k --nbr-ky 1`) and
  builds the disk pair projectors from (z₁−z₂)^m (z₁+z₂)^M expanded with
  binomials. `FQHECylinderBosonsDeltaInteraction` (delta = V₀, two
  ratios and sectors), `FQHECylinderFermionsLaplacianDelta` (= V₁) and,
  with GSL, `FQHECylinderFermionsTwoBodyGeneric` (V₁ + V₃ at ratio 1, 2
  and ½ and in the Ky = 2 sector) reproduce whole sectors to 10⁻¹²;
  `FQHEDiskCoulombPseudopotentials` equals Γ(m+½)/(2 m!) to 10⁻¹⁵. The
  disk Hamiltonian is another story (**U34**): `FQHEDiskBosonsTwoBodyGeneric`
  gives the right zero modes and wrong nonzero energies from N = 3 on
  (three bosons at Lz = 6 with V₀: 1.5951, 1.9825, 2.0000 … where the
  exact values are 3/2, 15/8, 33/16 …). Two independent computations in
  the oracle agree to 10⁻¹⁵: the second-quantised projector assembly, and
  an exact rational first-quantised projection in the basis of monomial
  symmetric polynomials (the relative-angular-momentum decomposition is
  orthogonal, so no metric is needed); `selftest.geometry_ed_oracle`
  checks them against each other. The three `knownbug.fqhe.disk.*` tests
  carry the exact spectra and turn red when the program is fixed; the
  fermionic program hangs (U29) and could not be tested at all.
- **k-body sphere ED** (8 tests, `physics.fqhe.sphere.ed.nbody.*`):
  `tests/oracles/sphere_nbody.py` builds three- and four-body
  pseudopotential Hamiltonians from the projectors onto the L = kl − m
  multiplets of the k-particle space (L² from the summed ladder operators,
  no Clebsch–Gordan algebra; for the m with a single relative state that
  multiplet projector is the pseudopotential projector) and assembles the
  N-body matrix as Σ ⟨a|v|b⟩ A†_a A_b over normalised k-particle Fock
  states. DiagHam's k-body pseudopotentials are not projector weights: for
  V^(k)_m = 1 the k-particle multiplet has an energy that depends on m and
  2S (1.5 for three bosons at m = 0 and any 2S, 2.9993 at m = 2 and 2S = 6,
  0.1639 at m = 3 …; the hard-core programs use yet another constant).
  Those constants are read off the k-particle runs and are the only thing
  taken from the programs; the N = 5 and 6 spectra of
  `FQHESphereBosonsNBodyHardCore` (3- and 4-body), `QHEFermionsNBodyHardCore`,
  `FQHESphereBosonsThreeBodyGeneric` and `FQHESphereFermionsThreeBodyGeneric`
  (three-body sets with and without a two-body term) are then predicted and
  match to 10⁻¹³. The two-body term through the same machinery reproduces
  `sphere_ed.py`. Not covered: the spinful three-body programs, the torus,
  cylinder and disk n-body programs (their contact matrix elements need
  the geometry; sessions 4 and later).
- **Spinful sphere ED** (15 tests, `physics.fqhe.sphere.ed.su2.*`,
  `su3.*`, `su4.*`): the same oracle with k species (`spectrum_su2`,
  `spectrum_sun`), one Fock block per species and pseudopotentials per
  species pair; a different-species pair of relative angular momentum m
  costs 2 V_m like a same-species one (fixed with one up and one down
  particle). `FQHESphereFermionsWithSpin`, `FQHESphereBosonsWithSpin`
  (SU(2)-symmetric Coulomb and spin-dependent sets, Sz = 0 and 1),
  `FQHESphereFermionsWithSU3Spin`, `FQHESphereBosonsWithSU3Spin` (keys
  `Pseudopotentials11` … `33`, sectors by 2Tz/3Y or particle numbers),
  `FQHESphereBosonsWithSU4Spin` (ten channels) reproduce whole spectra to
  10⁻⁹. The `--szsymmetrized-basis` and `--lzsymmetrized-basis` sectors
  (± parity) are checked one by one against the oracle's projection onto
  the eigenspaces of the exchange operators (25 + 22 and 28 + 19 of the 47
  Sz = 0, Lz = 0 states); both flags together give an empty fermionic basis
  (**U32**, `known-bug`). The legacy `QHEFermionsSphereWithSpin -v -w` has
  its own normalisation (a pair costs V_m, V in the up-down m = 0 channel,
  W in every m = 1 channel) and is checked in it. Not matched:
  `FQHESphereFermionsWithSpinFull` (four-index keys, its convention is not
  established) and the two-Landau-level programs (next session).
- **Spin-chain ED** (6 tests, `physics.spin.ed.*`): `tests/oracles/spin_ed.py`,
  dense numpy diagonalisation of XXZ / J₁–J₂ chains with a field, open
  and periodic, spin ½ and 1, L up to 10, against `GenericOpenSpinChain`
  and `GenericPeriodicSpinChain`.
- **Algorithm consensus** (`crosscheck.fqhe.sphere.*`): reorthogonalised
  Lanczos (`-n 4`) against full diagonalisation on the same Coulomb
  Hamiltonian (338 states, 4 lowest eigenvalues to 10⁻¹⁰); and, in builds
  with LAPACK, `--use-lapack` against DiagHam's own diagonaliser (full
  spectrum to 10⁻¹⁰).
- **Self-tests of the harness** (`selftest.checker.*`, `selftest.runner.*`,
  all but one `WILL_FAIL`): the checker must reject a wrong value, a wrong
  count, a missing file, NaN and a spectrum of the wrong length, and accept
  the right value; the runner must fail when a program writes no output
  file (stale output can never satisfy a check). If any of these ever
  "passes" the harness itself is broken.

### Torus exact diagonalisation and FCI band goldens (`-L ed`, `-L bands`)

`tests/oracles/torus_ed.py` is an independent lowest-Landau-level exact
diagonalisation on the torus in Landau gauge (Yoshioka's two-body matrix
element with the momentum sums converged to 10⁻¹⁴; Coulomb with the q = 0
term dropped, or Haldane pseudopotentials 4π Σ_m V_m L_m(q²) e^{−q²/2}
with q = 0 kept, normalised so that a pair of relative angular momentum m
costs 2 V_m as in `sphere_ed.py`; fermions or bosons; one K_y sector).
Eight tests `physics.fqhe.torus.ed.*`: Coulomb fermions N = 3, N_φ = 9 at
aspect ratio 1 and 2, N = 4 at N_φ = 8 and 12, Coulomb bosons N = 3,
N_φ = 6, V₁ and V₁ + V₃ fermions, V₀ + V₂ bosons. `FQHETorusFermionsCoulomb`,
`FQHETorusBosonsCoulomb` and the two `TwoBodyGeneric` programs reproduce
every eigenvalue of the sector to 10⁻⁹ (the tolerance is 10⁻⁸: DiagHam
truncates its own momentum sums at about 10⁻¹⁰, visible in its degenerate
pairs; the oracle is converged further). The reference spectra and the
pseudopotential files are in `tests/data/torus_ed/`,
`selftest.torus_ed_oracle` re-derives them.

`spectrum_species` in the same oracle takes several species (one Fock block
each, kernels per species pair, a different-species pair of relative
angular momentum m costing 2 V_m like a same-species one), the Landau
level index (form factor [L_n(q²/2)]² e^{−q²/2}), the inter-layer Coulomb
kernel 2π e^{−qd}/q, and the magnetic-translation sectors: inside a K_y
sector the translation of every particle by N_φ/gcd(N, N_φ) orbitals is a
signed permutation whose eigenvalues e^{2πi K_x/g} label the K_x sectors
(the program's labelling of K_x, fixed once on N = 4, N_φ = 8, is the only
convention taken). Seventeen further tests: `FQHETorusFermionsCoulomb` in
Landau level 1; `FQHETorusFermionsTwoBodyGenericAllMomenta` (every K_y at
once); `FQHETorusFermionsWithTranslations` and `BosonsWithTranslations`
((K_x, K_y) sectors); `FQHETorusFermionsWithSpinTwoBodyGeneric`,
`BosonsWithSpinTwoBodyGeneric` (spin-dependent pseudopotentials, keys
`PseudopotentialsUpUp/DownDown/UpDown`, Sz = 0 and 1);
`FQHETorusFermionsWithSpin` (SU(2) Coulomb and bilayer Coulomb at d = 1);
`FQHETorusFermionsWithSpinAndTranslations` (bilayer d = 0.5, two (K_x, K_y)
sectors); `FQHETorusBosonsWithSpinAndTranslations`;
`FQHETorusBosonsWithSU3SpinTwoBodyGeneric`, `WithSU3SpinAndTranslations`
(keys `Pseudopotentials11` … `33`) and `FQHETorusBosonsWithSU4SpinTwoBodyGeneric`
(ten channels). All whole sectors to 10⁻⁹ (tolerance 10⁻⁸). While writing
it the oracle disagreed with the program once: the fermionic operator
order of the two-species term was wrong in the oracle (the momentum
transfer must stay within a species); the bosonic cases had hidden it.

`tests/oracles/fci_bands.py` writes the two bands of the checkerboard
lattice model from the published Bloch Hamiltonian (Sun, Gu, Katsura and
Das Sarma, PRL 106, 236803) at every lattice momentum. Six tests
`physics.fci.checkerboard.bands.*` (3×3 and 4×4 at the flat-band hoppings,
4×3 away from them, one test per band) compare the E₀ and E₁ columns of
`FCICheckerboardLatticeModel --singleparticle-spectrum --export-onebodytext`
with `tests/data/fci_bands/` to 10⁻¹²; `selftest.fci_bands_oracle`
re-derives the files. It is the first independent number for an FCI
program; the many-body FCI spectra remain regression-only.

### Regression spectra (`-L regression`)

12 tests, `regression.*`, one per module directory that has no physics
golden yet (FCI, FTI, FQHEOnLattice, QuantumDots, more of HubbardModels
and Spin) plus a few otherwise-untested families (spinful sphere, torus
Coulomb, the GSL-only cylinder fermions). Each reproduces, to 10⁻⁹, the
spectrum saved from this repository's r4493 build in
`tests/data/regression/<case>.dat`, whose header records program,
arguments, commit, preset, compiler and date. **A regression reference is
not a golden**: it is what the build produced, right or wrong; the test
detects change. The cases are in `tests/data/regression/cases.txt` and
the references are regenerated with `python3 tests/regression_reference.py
build/default --write` — only for an understood, intended change, said so
in the commit. `check_spectrum numbers` compares a program's standard
output when it writes no file (`PeriodicQuantumDot2D`).

A regression case must give the same spectrum in builds with and without
LAPACK (it is checked before being added). `FQSHCheckerboardModelTwoBands`
does not: its bands are degenerate, so the projected model depends on the
eigensolver's choice of basis (differences of 0.14–0.17 between the
internal diagonaliser and LAPACK). That is a property of the model set-up,
not a defect, but it makes the program unsuitable for a portable reference;
`FTI3DSimpleTI` (identical to 10⁻¹⁵ in both builds) stands in for the FTI
module.

A second `known-bug` test came out of this batch (**U31**):
`FCICheckerboardLatticeModel` without `--single-band` (the two-band model)
segfaults after `start` for every size tried; the regression case uses
`--single-band`, and `knownbug.fci.checkerboard.two_band_model_segfault`
(`WILL_FAIL`) holds the reproducer.

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

Five `known-bug` tests exist (U31, the two-band checkerboard segfault; U32,
the empty fermionic basis with both symmetrisation flags; U34, the disk
pseudopotential Hamiltonian, three cases; see above. U30,
plain Lanczos below the ground state on a degenerate spectrum, was fixed
here and its test promoted to `physics`). A `known-bug` test reproduces a defect documented in this repository and is
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
on); with the default preset the ctest suite has 738 tests at r4493.
