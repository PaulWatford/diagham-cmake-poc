# Test suite reference

Source: `docs/reference/tests.md` (2026-09-24), moved 2026-09-26; the "Adding a golden" section is now `docs/how-to/develop/add-a-golden-test.md`; counts updated to the r4493 build.

`ctest` is the test runner. The suite lives in `tests/` and is registered
whenever programs are built (`DIAGHAM_BUILD_TESTS`, default ON).

```bash
ctest --preset default                 # everything
ctest --preset physics                 # physics goldens only (seconds)
ctest --preset default -L smoke        # the per-program --help checks
ctest --preset default -R hubbard      # by name
```

## What is tested

| Label | Tests | What passes means |
|---|---|---|
| `physics` | analytic / reference goldens and cross-program comparisons (table below) | a DiagHam program, run from a clean directory, writes a spectrum that matches a known answer or an independent implementation |
| `python` | `python.hubbard_ed_cross_check` | DiagHam's Hubbard ground state agrees with the independent Python ED in `benchmarks/hubbard_ed.py` to 1e-10, on four lattices/couplings (needs numpy; skipped at configure time without it) |
| `smoke` | `smoke.<target>`, one per program | the program starts, parses options and exits 0 on `--help` (catches link and static-initialisation breakage in the programs no golden reaches). Four programs are left out, with the reason in `tests/CMakeLists.txt`: `QHEBosons`, `MultipleSpinChain` and `TestDiagHamVectors` have no option parser; `EvaluateBroadening` checks its required `--input` before `--help` |
| `install` | `install.find_package_consumer` | the `Development` install component, used through `find_package(DiagHam)` alone from a separate CMake project, compiles, links, and diagonalises a tight-binding ring correctly |
| `known-bug` | reproducers for documented upstream defects | the defect is still present (`WILL_FAIL`); see below |

### Physics goldens

| Test | Program | Check | Why the answer is known |
|---|---|---|---|
| `physics.hubbard.2x2.U4.ground_state_is_minus_4sqrt2` | `HubbardSquareLatticeModel` | lowest eigenvalue within 4 ulp of -5.6568542494923806 | analytic -4√2 (the locked golden of `AGENTS.md`); 0 ulp with GCC 13 |
| `physics.hubbard.2x2.U0.ground_state_is_minus_8` | `HubbardSquareLatticeModel` | lowest eigenvalue = -8 ± 1e-12 | free fermions, filled band |
| `physics.hubbard.2x2.U4.full_spectrum` | `HubbardSquareLatticeModel` | all 36 eigenvalues vs `benchmarks/` reference, 1e-11 | saved, independently cross-checked output |
| `physics.hubbard.2x4.U4.sector_ground_states` | `HubbardSquareLatticeModel` | lowest eigenvalue of each of 8 momentum sectors vs `benchmarks/` reference, 1e-10 | saved, independently cross-checked output |
| `physics.fqhe.sphere.laughlin_1_3.unique_zero_mode` | `FQHESphereFermionsTwoBodyGeneric` | exactly 1 zero eigenvalue at N=6, 2S=15, Lz=0 | Laughlin state is the unique V1 zero mode at 2S=3(N-1) |
| `physics.fqhe.sphere.laughlin_1_3.two_quasihole_count` | `FQHESphereFermionsTwoBodyGeneric` | exactly 4 zero eigenvalues at N=6, 2S=17, Lz=0 | two bosonic quasiholes of angular momentum N/2=3 give L=6,4,2,0 |
| `physics.fqhe.torus.laughlin_1_3.threefold_degeneracy` | `FQHETorusFermionsTwoBodyGeneric` | exactly 3 zero eigenvalues over all Ky at N=4, Nphi=12 | topological degeneracy of Laughlin 1/3 on the torus |
| `physics.fqhe.torus.su2_coulomb.polarized_N3_matches_spinless` | `FQHETorusFermionsWithSpin` vs `FQHETorusFermionsCoulomb` | full 84-state spectra agree to 1e-10 | full polarisation must reproduce the spinless problem; regression for patch 14 (fails by 0.25 without it) |
| `physics.fqhe.torus.su2_coulomb.unpolarized_N3_matches_translations` | `FQHETorusFermionsWithSpin` vs `FQHETorusFermionsWithSpinAndTranslations` | full 324-state Sz=1 spectra agree to 1e-8 | two independent implementations; regression for patch 14 (fails by 0.43 without it) |
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
on); with the default preset the ctest suite has 613 tests at r4493.
