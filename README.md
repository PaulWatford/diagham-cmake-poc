# DiagHam: CMake migration

[DiagHam](https://www.nick-ux.org/diagham/) is an exact-diagonalization
toolkit for strongly correlated quantum systems: fractional quantum Hall,
fractional Chern insulators, Hubbard models, spin chains, quantum dots. It is
some 1.5 million lines of C++ (headers included) across ~1,450 classes and
603 programs at SVN r4493, developed since
2003 (CVS, then Subversion) by Nicolas Regnault, Gunnar Möller, Zlatko Papić
and collaborators, and built with a GNU autotools system that the upstream
community has identified as a modernisation target.

This repository is the migration of that code base to Git and CMake, done in
discussion with Gunnar Möller (University of Kent). It contains the complete
DiagHam source with its full history, a CMake build for every module, a test
suite with physics goldens, and a small set of upstream fixes. The physics,
the algorithms and the bulk of the engineering remain the work of the DiagHam
authors; this repository collects them into a modern build and verifies that
the result is correct.

---

## Repository layout

Two branches matter:

- **`upstream`** — a faithful `git svn` conversion of the canonical DiagHam
  Subversion repository (trunk, branches and tag; 4,477 commits, 2003-05-03
  to 2026-09-18, currently SVN r4493). It is never edited by hand; it is
  refreshed from SVN with `git svn fetch`.
- **`main`** — this CMake work, merged with `upstream`, plus the accepted
  fixes as ordinary commits. `main` is where the project is developed and
  what CI builds.

```
.
|---- CMakeLists.txt              top-level build: options, modules, dependencies
|---- CMakePresets.json           default / core / lapack / full / hpc / mkl / debug
|---- cmake/                      helpers, install/export, config header template, Find modules
|---- scripts_cmake/
|   |---- extract_autotools.py    reads upstream Makefile.am, emits the per-directory CMakeLists.txt
|   +---- ci.sh                   the CI pipeline (GitHub Actions and GitLab CI call it)
|---- tests/                      ctest suite: physics goldens, smoke, install (TESTING.md)
|---- patches/                    audit trail of the upstream fixes (PATCHES.md); the fixes themselves are commits
|---- benchmarks/                 Hubbard verification log + independent Python exact diagonalisation
|---- Base/ src/ FQHE/ FTI/ Spin/ QuantumDots/   DiagHam itself, as in upstream, with generated CMakeLists.txt in each directory
|---- docs/                     knowledge base: how-to, reference, explanation, tutorials, history (docs/README.md)
+---- TESTING.md, DEFERRED.md, MIGRATION_ROADMAP.md, AGENTS.md, CONTRIBUTING.md
```

The per-directory `CMakeLists.txt` files are generated from upstream's own
`Makefile.am` files by `scripts_cmake/extract_autotools.py` (88 files at
r4493) and are committed, so the repository builds without running the
generator. Re-run it after an upstream `Makefile.am` change.

## Build and test

```
git clone <this repository> DiagHam && cd DiagHam
cmake --preset default            # every module, SMP, no optional libraries
cmake --build --preset default -j
ctest --preset default            # physics goldens, --help smoke test, install test
```

Current outcome with the `default` preset on DiagHam r4493 (Ubuntu, GCC 15):

- **61 of 61 static libraries** build (Base/src, src, FQHE, FTI, Spin, QuantumDots).
- **603 of 603 programs** build (including the 19 QuantumDots analysis
  tools). Five sources are deliberately not built: one
  legacy duplicate (`QHEFermionsTorusWithSpin`, see `DEFERRED.md`) and four
  that upstream itself never lists in `bin_PROGRAMS`.
- **613 of 613 ctest tests pass**: 13 physics tests with analytic,
  independently computed or cross-implementation answers (Hubbard, Laughlin
  zero modes on sphere and torus, Heisenberg rings, spinful torus Coulomb),
  an independent Python exact-diagonalisation cross-check, `--help` on every
  program, and a `find_package(DiagHam)` consumer build. See `TESTING.md`.
- The Hubbard 2x2 U=4 ground state is `-5.6568542494923806`, a bit-identical
  (0 ULP) match to the analytic `-4√2` (see `benchmarks/BENCHMARK.md`).
- The `lapack` preset (system LAPACK/BLAS) builds with no errors and passes
  the same suite; `full` and `hpc` (MPI + ScaLAPACK) were verified in CI.

`cmake/verify_build.sh` compares a CMake build against an autotools build of
the same tree (library coverage and `nm` symbol counts). It reports
79 passed / 0 failed on the default configuration (FQHE and FTI on); it needs
an autotools build to compare against and is not part of `ctest`.

## Why CMake (and what the autotools build does)

The upstream `configure.ac` has **34 `AC_ARG_ENABLE` / `AC_ARG_WITH` flags**:
`--enable-lapack`, `--enable-gsl`, `--enable-mpi`, `--enable-bz2`,
`--enable-gmp`, `--enable-mpack`, `--enable-fftw`, `--enable-scalapack`,
`--with-blas-libs=...`, and so on. Each flag is some cluster admin's hard-won
dependency chain. Every one of them has a CMake equivalent; the mapping is in
`docs/reference/configure-flag-map.md`.

The autotools build also relies on two custom Perl scripts
(`scripts/genmake.pl` and `scripts/genam.pl`) that generate `Makefile.am`
entries from directory contents. Under CMake that layer is replaced by the
generator above, which reads the `Makefile.am` files rather than the
directories: a program is built if upstream lists it in `bin_PROGRAMS`, plus
a short, documented allow-list of sources upstream never listed.

## Methodology

The per-subdirectory CMakeLists.txt files are generated, not hand-written,
because a real migration benefits from a *reproducible* derivation of CMake
from autotools rather than a one-shot port: after any upstream change, one
script run gives the updated build. The generated files are deliberately
thin (a few `diagham_add_library` / `diagham_add_programs` calls); all
complexity lives in `cmake/DiagHamHelpers.cmake`.

Correctness is checked at three levels: the build (every library and program
upstream builds, we build), build parity with autotools (`verify_build.sh`),
and physics (the ctest goldens and the independent Python diagonaliser).

## Upstream issues found during the migration

Moving to a stricter build surfaced a number of long-standing issues in the
upstream code base. Each fix is an ordinary commit on `main` with an
`Upstream-Patch:` trailer; `patches/PATCHES.md` is the audit trail (root
cause, what was compared, how it was verified) and `DEFERRED.md` records the
one file excluded instead of fixed.

| Class | Count | Issue | Status |
|---|---:|---|---|
| A | 5 + 9 | Missing `#ifdef __LAPACK__` gating around `LapackDiagonalize` (FQHE/FTI, then Spin) | fixed; one Spin file was fixed upstream after 2020 |
| C | 1 | `#include` of a header that never existed | fixed |
| D | 1 | Undeclared variable (copy-paste from a sibling) | fixed |
| E | 1 | `precision(14)` at result-output sites, truncating output below double precision | fixed at the 14-digit sites only (978 sites, 544 files at r4493); deliberate display widths left alone |
| F, G | 2 | `FCIWannierConstruction` and `FCIDiceLatticeModel` never compiled or linked upstream | enabled; the Dice-lattice defaults need a maintainer's sign-off |
| H | 1 | Invalid C++ (`ostream << ostream&`) that broke every QuantumDots program | fixed |
| I | 1 | **Physics defect**: spinful torus Coulomb wrong for N ≥ 3 and for every unpolarised case (basis order vs `FindStateIndex`; dropped same-orbital up-down terms) | fixed, verified against two independent programs to 1e-12; **pending maintainer review** |

Several of the "hidden" failures turned out to be programs upstream had
dropped from its build lists rather than failures autotools concealed; the
distinction is recorded per file in `PATCHES.md`. One earlier fix (an
uninitialised member in `SpinChainHamiltonianWithTranslations`) was dropped
because upstream fixed it after 2020.

## Demonstration: Hubbard ED at machine precision

`benchmarks/BENCHMARK.md` is the physics verification log: the 2x2 Hubbard
model at U=4 reproduces the analytic `-4√2` bit for bit; the 2x4 model
(basis dimension 4,900) and the 3x3 model agree with a from-scratch
100-line Python exact diagonalisation (`benchmarks/hubbard_ed.py`) on full
spectra; and the strong-coupling limit E₀·U → −48 is recovered. The same
checks run as ctest goldens on every build.

## Finding your way around

`MODULE_MAP.md` is organised by task ("I want to add a dependency", "I want
to verify the physics", ...) and points at the file for each.
`docs/reference/configure-flag-map.md` maps every configure flag; `docs/how-to/build/` has the per-platform build guides and the cluster
recipe; `TESTING.md` explains the test suite. The DiagHam user manual —
program-by-program pages — is on the upstream wiki at nick-ux.org and is
being brought into `docs/` (see `MIGRATION_ROADMAP.md`).

## AI-assisted contribution rules

AI tools were used heavily in this migration and the rules for their use are
explicit and versioned with the code: see `AGENTS.md`. In short, the build
and the tests decide pass or fail, an AI agent may describe a change and
judge its usefulness but never sets pass/fail, never edits a golden and
never merges, and DiagHam's physics code is out of bounds for automated
changes; every proposed change is held for a named human maintainer.
`CONTRIBUTING.md` has the human workflow.

## Production migration context

The agreed plan (see `MIGRATION_ROADMAP.md`): Git first, then CMake on top;
a single repository; the full SVN history preserved (done — the `upstream`
branch); static linking retained; the production home to be the University
of Kent's GitLab, with this GitHub repository as the development home in the
meantime. Until the maintainers move to Git themselves, the canonical source
remains the Subversion repository at nick-ux.org and `upstream` mirrors it.

## Still open

- Maintainer review of the physics fix (class I) and the Dice-lattice defaults.
- The Kent cluster recipe (`cmake/KentDefaults.cmake` is a scaffold; the
  site paths are not yet known) and a tested Intel MKL configuration.
- Broader goldens: overlaps, entanglement spectra, FCI spectra; an `mpirun` test.
- The "best copy of each" review of duplicated and orphaned upstream classes
  (`DEFERRED.md` lists the candidates).
- Spack/EasyBuild packaging; a curated library API.

## Build dependencies

- CMake >= 3.21 (presets), a C++11 compiler (GCC or Clang), pthreads,
  Python 3 (generator and cross-check; numpy for the Python test).
- Optional: LAPACK/BLAS or Intel MKL, MPI, ScaLAPACK, GSL, GMP, MPACK,
  FFTW3, bzip2 — see `docs/reference/configure-flag-map.md`.

See `CHANGELOG.md` for what changed and when.

## License

DiagHam is licensed under the GNU General Public License, version 2 or
later; the full text is in `COPYING`. The migration work is released under
the same terms; `LICENSE` carries the contribution notice.

## Citing

If DiagHam contributes to a publication, please cite it as software (see
`CITATION.cff`) and the relevant papers in the DiagHam publication list on
the upstream wiki.

## Acknowledgements

The upstream code base is the work of Nicolas Regnault, Gunnar Möller,
Zlatko Papić, Cécile Repellin, Antoine Sterdyniak, Duc Phuong Nguyen,
Niall Moran, Yang-Le Wu and other contributors since 2003 (`AUTHORS`).
