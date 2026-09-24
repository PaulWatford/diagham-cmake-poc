# DiagHam: CMake Migration Proof-of-Concept

[DiagHam](http://www.nick-ux.org/diagham/wiki) is an exact-diagonalization
toolkit for strongly correlated quantum systems, fractional quantum Hall,
fractional Chern insulators, Hubbard models, spin chains. It is ~573,000
lines of C++ across ~2,350 classes and 565 executables, built on a 24-year-old
GNU autotools system that the upstream community has identified as a
modernisation target.

This repository holds initial work on replacing that autotools build with
a clean CMake build, following discussion with Gunnar Möller (University
of Kent). It also exercises one of DiagHam's Hubbard-model executables
end-to-end to confirm the build chain produces correct physics output.

A note on scope: most of this is assembly rather than invention. The
build fixes restore patterns already present in DiagHam's own working
files, and the physics, the algorithms, and the bulk of the engineering
remain the work of the DiagHam authors. This proof of concept collects
those pieces into a CMake build and verifies that the result is correct.

---

## What this repository contains

```
.
|---- CMakeLists.txt              top-level build: options, dependencies, modules
|---- CMakePresets.json           default / core / lapack / full / hpc / mkl / debug
|---- cmake/
|   |---- DiagHamHelpers.cmake          helper macros (diagham_add_library, ...)
|   |---- DiagHamInstall.cmake          install rules + find_package(DiagHam) export
|   |---- DiagHamConfig.cmake.in        package config template
|   |---- ApplyUpstreamPatches.cmake    applies upstream-bug patches at configure time
|   |---- config_ac.h.in                autoconf-equivalent header template
|   |---- Find*.cmake                   GMP, FFTW3, GSL, MKL, ... find modules
|   +---- verify_build.sh               CMake-vs-autotools build-parity script
|---- scripts_cmake/
|   |---- overlay.py               lays this build over an upstream DiagHam tree
|   |---- extract_autotools.py     reads upstream Makefile.am, emits CMakeLists.txt
|   +---- ci.sh                    the CI pipeline (GitHub Actions and GitLab CI call it)
|---- tests/                       ctest suite: physics goldens, smoke, install (TESTING.md)
|---- patches/
|   |---- PATCHES.md                15 upstream-bug patches, with audit trail
|   +---- 01-...15-*.patch
|---- benchmarks/
|   |---- BENCHMARK.md             physics verification log + independent Python ED
|   |---- hubbard_ed.py            100-line independent reference implementation
|   +---- fermions_hubbard_square_x_2_y_{2,4}*.dat.txt   (saved DiagHam output)
|---- .github/workflows/ci.yml     GitHub Actions
|---- .gitlab-ci.yml               same pipeline for the Kent GitLab move
|---- CONFIGURE_FLAGS.md           every configure.ac flag -> its CMake option
|---- HPC.md                       cluster build recipe
|---- TESTING.md                   what ctest checks and how to add a golden
|---- DEFERRED.md                  1 excluded file: what it does, what's missing
|---- BUG_torus_su2_coulomb.md     spinful torus Coulomb defect report
|---- HUBBARD_BENCHMARK.md         Hubbard ED demonstration with physics output
+---- README.md                    this file
```

`scripts_cmake/overlay.py` copies this repository's build files into an
upstream DiagHam checkout and runs `extract_autotools.py`, which generates
the per-subdirectory `CMakeLists.txt` files under `Base/src/`, `src/`,
`FQHE/src/`, `FTI/src/`, `Spin/src/` and `QuantumDots/src/` from
upstream's own `Makefile.am` files (88 files against the pinned upstream
revision; the number tracks upstream). Those generated files live in the
DiagHam tree, not in this repository.

## What works

```
git clone https://github.com/guysoft/DiagHam.git
python3 scripts_cmake/overlay.py DiagHam
cd DiagHam
cmake --preset default            # patches apply at configure time
cmake --build --preset default -j4
ctest --preset default            # 572 tests: physics goldens, smoke, install
```

`scripts_cmake/ci.sh` does exactly this against a pinned upstream
revision (`ed78a30`, the mirror's head) and is what CI runs.

Current outcome with the `default` preset (every module on, SMP, no
optional libraries), GCC 13, Ubuntu 24.04:

- **61 of 61 static libraries** build (Base/src, src, FQHE, FTI, Spin,
  QuantumDots).
- **561 programs** build: 473 core/FQHE/FTI + 69 Spin + 19 QuantumDots.
  1 legacy program is explicitly skipped (`DEFERRED.md`). Bringing Spin
  and QuantumDots in surfaced two more upstream compile defects and one
  defect in this PoC's own patch 08, now fixed (patches 12, 13; see
  `patches/PATCHES.md`, "Class H").
- **572 of 572 ctest tests pass**, including 13 physics tests with
  analytic, independently computed, or cross-implementation answers
  across Hubbard, FQHE (sphere and torus Laughlin zero modes, spinful
  torus Coulomb) and Spin (Heisenberg rings), an independent Python ED
  cross-check, `--help` on every program, and a `find_package(DiagHam)`
  downstream-consumer build. See `TESTING.md`.
- **The spinful torus Coulomb defect is fixed** (patch 14, pending
  maintainer review): two independent upstream bugs, one in the basis
  ordering behind `FindStateIndex` and one that dropped same-orbital
  up-down interaction terms, so every `FQHETorusFermionsWithSpin` result
  was wrong for N >= 3 and every unpolarised one even for N = 2. After
  the fix the program agrees with the spinless program and with the
  independent `FQHETorusFermionsWithSpinAndTranslations` to 1e-10-1e-12
  on full spectra. The ctest suite also caught an uninitialised-member
  bug in `SpinChainHamiltonianWithTranslations` (patch 15).
- The same suite passes with the `full` (LAPACK, GSL, GMP, FFTW3, bzip2)
  and `hpc` (MPI, LAPACK, ScaLAPACK, GSL) presets, and the `core`
  preset builds and passes with Clang.
- Hubbard 2x2 U=4 benchmark: ground state `-5.6568542494923806`, a
  bit-identical (0 ULP) match to the analytical answer
  `-4√2 = -5.656854249492381...`, now enforced by
  `physics.hubbard.2x2.U4.ground_state_is_minus_4sqrt2`.
- `cmake --install` gives `bin/`, `lib/diagham/`, `include/diagham/`
  and `lib/cmake/DiagHam/`; downstream projects use
  `find_package(DiagHam)` and link `DiagHam::DiagHam`.

`cmake --preset core` reduces the build to the original core-only scope.
`cmake/verify_build.sh` (CMake vs autotools build parity: library
existence and `nm` symbol counts) is separate from ctest and still needs
an autotools build to compare against; the "79 checks" figure quoted in
`AGENTS.md` is still not independently confirmed (see the note there).

## Why CMake (and what the autotools build does today)

The upstream `configure.ac` has **36 `AC_ARG_ENABLE` / `AC_ARG_WITH` flags**:
`--enable-lapack`, `--enable-gsl`, `--enable-mpi`, `--enable-bz2`,
`--enable-gmp`, `--enable-mpack`, `--enable-fftw`, `--enable-scalapack`,
`--with-blas-libs=...`, and so on. Each flag is some cluster admin's
hard-won dependency chain. Any CMake migration must preserve every one of
them; `CONFIGURE_FLAGS.md` maps each of the 34 distinct flags to its CMake
option or standard CMake variable, and `HPC.md` gives the cluster recipe.

The autotools build also relies on two custom Perl scripts
(`scripts/genmake.pl` and `scripts/genam.pl`) that auto-generate
`Makefile.am` entries from directory contents. Under CMake, that whole
layer goes away: adding a new program is just dropping a `.cc` file in the
relevant Programs directory. (Spin and QuantumDots follow their `bin_PROGRAMS`
list instead, because those directories contain sources upstream
deliberately doesn't build.)

## Methodology

The per-subdirectory CMakeLists.txt files are auto-generated by
`scripts_cmake/extract_autotools.py`, which parses every upstream `Makefile.am`
to extract its `libFOO_a_SOURCES` lists and emits the equivalent CMake.
The Python script is itself part of the PoC, a real migration
benefits from a *reproducible* derivation of CMake from autotools, not
just a one-shot hand-port. Re-running the script after any upstream
`Makefile.am` change gives an updated CMake snapshot.

The generated files are deliberately thin (each one a few `diagham_add_library`
calls); all complexity lives in `cmake/DiagHamHelpers.cmake`.

## Upstream issues surfaced during the migration

The CMake port surfaced **thirteen** long-standing issues in the upstream
codebase that the autotools build either silently hid or fenced off
behind broken include paths. Three are pre-existing notes from the
core-only iteration of this PoC; ten emerged after FQHE and FTI were
brought into scope. More (classes H and I in the table below,
24/09) came from building Spin and QuantumDots and from the new ctest
suite, two of them wrong-physics bugs rather than compile errors.

### Three from the core iteration

- **`Fermions.cc` orphan from May 2001**, still listed in
  `libQHEHilbertSpace_a_SOURCES` but its include path is missing from
  the upstream `Makefile.am`. A 24-year-old orphan in the build list.
- **`DelocalizedRealVector.cc` dormant dead code**, entirely wrapped
  in `#ifdef USE_CLUSTER_ARCHITECTURE`. The file doesn't
  `#include "config.h"`, so the macro never reaches it and the build
  produces a 1456-byte object file with no symbols.
- **`make -C FTI/src` fails out of the box**, FTI's `Makefile.am`
  only sets `-I src -I Base/src`, but FTI sources include from
  `FQHE/src` too. A clean `./configure && make -C FTI/src` fails.

### Thirteen more surfaced once FQHE + FTI came into scope, plus deeper inspection

| Class | Count | Issue | Resolution |
|---|---:|---|---|
| A | 5 | Missing `#ifdef __LAPACK__` gating around `LapackDiagonalize` | Patched (canonical pattern lift) |
| B | 3 | Stale constructor / API mismatches | 1 patched; 2 since resolved (now Classes F, G via patches 10, 11) |
| C | 1 | `#include` of a header that never existed (`*New.h`) | Patched (author left the working version commented out) |
| D | 1 | Undeclared variable `SubsystemSize` (copy-paste from sibling) | Patched (matched to working twin's convention) |
| E | 2 | `precision(14)` hardcoded at result-output sites, silently truncating output below IEEE-754 double precision | Patched in 2 parts (patches 08, 09) across 448 result-writing files; uses `std::numeric_limits<double>::max_digits10` |
| F | 1 | `FCIWannierConstruction.cc` (Wannier construction for fractional Chern insulators) references public getter methods that were never defined on parent Hilbert-space classes | Patched (patch 10); adds 4 inline getters to existing classes |
| G | 1 | `FCIDiceLatticeModel.cc` (Dice-lattice \|C\|=2 FCI) had unwired tight-binding parameters, depended on an abandoned stub Hilbert-space class, and on 3 classes orphaned from every Makefile.am | Patched (patch 11); wires t1/t2/l1/l2 from the Kagome sibling, routes to the completed SU2 boson class, recovers the 3 orphaned classes via allowlisted build-system recovery |
| - | 1 | `HAVE_FTI` macro not defined in our CMake config | Fixed in `cmake/config_ac.h.in` |
| H | 2 (+1) | Found once Spin and QuantumDots were built (24/09): a C++11-invalid `ostream << ostream&` in `ThreeDTwoParticles.cc` that broke every QuantumDots program; 22 more ungated `LapackDiagonalize` calls in 9 Spin programs; and a misplaced `#include <limits>` in this PoC's own patch 08 | Patches 12, 13; patch 08 corrected |
| I | 2 | Wrong physics / undefined behaviour caught by the ctest suite (24/09): the spinful torus Coulomb defect (basis order + up-down interaction) and an uninitialised `SpinChainHamiltonianWithTranslations` | Patches 14, 15 (need maintainer physics review) |

15 are patched in `patches/` with full audit trails in
[`patches/PATCHES.md`](patches/PATCHES.md). 1 is deferred to the
maintainer's review (a legacy duplicate with a canonical replacement), documented in [`DEFERRED.md`](DEFERRED.md) with
grep-verified evidence: file sizes, inline TODO/FIXME markers,
git-history snapshots, and the specific missing-input physics
parameters where applicable.

## Demonstration: Hubbard ED with machine-precision physics output

With the build working, the next check is that the produced binaries
give correct physics on problems with known answers. Two documents cover
this:

- **`HUBBARD_BENCHMARK.md`** (repo root), original U-sweep across
  `{0,1,2,3,4,5,6,7,8}` at 2x2 half-filling.
- **`benchmarks/BENCHMARK.md`**, physics verification log for this
  iteration of the PoC, including an independent **100-line Python
  exact diagonalisation** (`benchmarks/hubbard_ed.py`) that builds
  the Hubbard Hamiltonian from scratch in the 2nd-quantised Fock basis
  and reproduces DiagHam's output to machine precision on both 2x2
  (basis dim 36) and 2x4 (basis dim 4900) test cases.

Highlights:

- **2x2 Hubbard, half-filling, U=4**: DiagHam gives ground-state energy
  `-5.6568542494923806`. Analytical value: `-4*sqrt(2) = -5.656854249492381...`.
  This is a **bit-identical (0 ULP) match** to the double-precision rounding
  of the exact analytic value (verified with `mpmath`, 50-digit precision).
  *(Corrected 22/09, audited by a second Claude session: this file
  previously quoted `-5.6568542494924`, described as agreeing to
  "1.95×10⁻¹⁴, the limit of double-precision arithmetic." That figure was
  real output but from a run predating patch 08's precision fix — it is
  actually 22 ULP (≈1.95×10⁻¹⁴) away from the exact double-rounding, not
  at the precision limit. A freshly rebuilt, freshly run binary this
  session reproduced the corrected value above, matching what
  `benchmarks/BENCHMARK.md`'s own Test 1 section already independently
  stated.)*
- **2x4 Hubbard, half-filling, U=4**: DiagHam gives `E_0 = -10.252952955264`.
  Independent Python ED gives `-10.2529529552636`. Match to 4×10⁻¹³.
- **Strong-coupling U-scaling test (U ∈ {50, 100, 200, 500}):** E₀ · U
  converges to a constant (≈ −48.14) confirming the expected 1/U scaling
  in the Heisenberg limit.

## Finding your way around

The directory tree above is what's on disk; [`MODULE_MAP.md`](MODULE_MAP.md)
is organised by task instead ("I want to add a dependency", "I want to
verify the physics", ...) and points at the specific file for each.

## AI-assisted contribution rules

This repository accepts AI-assisted contributions under a strict
review policy: see [`AGENTS.md`](AGENTS.md) for the full rules. In
short — an AI agent may describe a diff and judge its usefulness, but
it never sets pass/fail, never edits a golden file, and never merges;
every suggested change is held for a named human maintainer to
approve, reject, or ask for tests to be restored. The automated
mailer/hold pipeline described there is not yet built; `AGENTS.md` is
the specification for it.

## Production migration context

This repository is a personal proof-of-concept on GitHub, developed
against a read-only git mirror for ground truth. It is not the
production repository. See
[`MIGRATION_ROADMAP.md`](MIGRATION_ROADMAP.md) for how this work
reconciles with what was actually agreed for the production migration
(Kent GitLab hosting, single-repo, full SVN history preservation, and
static linking) and what, if anything, changes as a result.

## Beyond this iteration

Still open, roughly in order of importance:

- **Maintainer review of patches 14 and 15**, which change upstream
  physics code (the spinful torus Coulomb fix and the spin-chain
  initialisation fix). Both are verified against independent
  implementations, but the sign-off is the maintainers'.
- **Physics sign-off** on the Dice-lattice defaults and the leftover
  `// IS IT CORRECT?` interaction mapping (patch 11).
- **Broader goldens**: overlaps, entanglement spectra, FTI/FCI
  (Chern-insulator) spectra; today's goldens cover Hubbard, Laughlin
  zero modes and Heisenberg rings.
- **MPI runtime tests**: the `hpc` preset compiles and links the MPI and
  ScaLAPACK paths, but no test launches a program under `mpirun`.
- **Production hosting and history** (Kent GitLab, history-preserving
  SVN->Git conversion; see `MIGRATION_ROADMAP.md`) and maintainer
  adoption of the patch series.
- **Packaging**: Spack / EasyBuild recipes, a container image, a binary
  cache.
- **Library API / bindings**: `find_package(DiagHam)` exposes the C++
  libraries, but there is no curated stable API or Python/Julia binding.
- **Doxygen documentation extraction.**

## Build dependencies

- CMake >= 3.16 (>= 3.21 for `CMakePresets.json`)
- C++11 compiler (GCC, Clang)
- pthread, `patch`, Python 3 (overlay and extraction scripts)

Optional: LAPACK/BLAS or MKL, MPI, ScaLAPACK, GSL, GMP, MPACK, FFTW3,
bzip2 (see `CONFIGURE_FLAGS.md`); numpy for the Python cross-check test.

See [`CHANGELOG.md`](CHANGELOG.md) for what's changed and when.

## License

DiagHam is licensed under the GNU General Public License, version 2 or
later (see `COPYING` in the upstream repository). The contributions in
this proof of concept are released under the same terms. See the
`LICENSE` file for the contribution copyright notice, a statement of what
was changed, and the full license text.

## Acknowledgements

The upstream codebase is the work of Nicolas Regnault, Gunnar Möller,
Duc Phuong Nguyen, and contributors over 24 years.
