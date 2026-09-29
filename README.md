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
|---- tests/                      ctest suite: physics goldens, cross-checks, smoke, install, self-test (docs/reference/tests.md)
|---- patches/                    audit trail of the upstream fixes (PATCHES.md); the fixes themselves are commits
|---- benchmarks/                 Hubbard verification log + independent Python exact diagonalisation
|---- Base/ src/ FQHE/ FTI/ Spin/ QuantumDots/   DiagHam itself, as in upstream, with generated CMakeLists.txt in each directory
|---- docs/                     knowledge base: how-to, reference, explanation, tutorials, history (docs/README.md)
+---- docs/reference/tests.md, AGENTS.md, CONTRIBUTING.md, GOVERNANCE.md, SUPPORT.md, CODE_OF_CONDUCT.md, SECURITY.md
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
ctest --preset default            # goldens, cross-checks, --help smoke test, install test, manifest
```

Current outcome with the `default` preset on DiagHam r4493 (Ubuntu, GCC 15):

- **61 of 61 static libraries** build (Base/src, src, FQHE, FTI, Spin, QuantumDots).
- **603 of 603 programs** build (including the 19 QuantumDots analysis
  tools). Five sources are deliberately not built: one
  legacy duplicate (`QHEFermionsTorusWithSpin`, see `docs/explanation/deferred-code.md`) and four
  that upstream itself never lists in `bin_PROGRAMS`.
- **884 of 884 ctest tests pass** — but read that number honestly: 234 are
  physics goldens with independently known answers (Hilbert-space dimensions
  against exact counting; Laughlin, Moore–Read and Read–Rezayi zero modes
  and quasihole counts; Coulomb pseudopotentials against the closed form;
  Jack polynomials against exact diagonalisation; entanglement-spectrum
  counting; Hubbard tight binding; Heisenberg, XX, AKLT, Haldane–Shastry and transverse-Ising chains) plus independent Python exact
  diagonalisations (Hubbard, sphere and torus Coulomb and pseudopotential Hamiltonians with and without SU(2), SU(3) and SU(4) spin, spin chains) and the band structures of seven lattice models with the Haldane Chern number, checkerboard many-body spectra in the flat-band limit against a band-projected exact diagonalisation and the exactly solvable atomic limit, 4 are cross-checks between DiagHam
  programs or algorithms (spinful torus Coulomb; Lanczos vs full diagonalisation; LAPACK vs internal), 1 is a `find_package(DiagHam)` consumer
  build, 12 reproduce a spectrum saved from this build (regression: change detection, not correctness), 22 are self-tests (manifest, fourteen oracles, the checker and runner), 14 are known-bug reproducers, and **599 are `--help` smoke tests**,
  which prove that a program links and parses options and nothing about its
  physics. `docs/reference/test-coverage.md` gives the per-program truth:
  74 of 603 programs have a physics or cross-check test today (74 with an independently known answer) and 3 more a regression test (77 with any test beyond --help). Expanding that is the
  current work. See `docs/reference/tests.md` and, for what a green run does
  and does not prove, `docs/explanation/verification.md`.
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
`Upstream-Patch:` trailer; `docs/explanation/upstream-fixes.md` is the index, `patches/PATCHES.md`
the audit trail (root cause, what was compared, how it was verified) and
`docs/explanation/deferred-code.md` records the one file excluded instead of fixed.

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
distinction is recorded per file in `patches/PATCHES.md`. One earlier fix (an
uninitialised member in `SpinChainHamiltonianWithTranslations`) was dropped
because upstream fixed it after 2020.

## Demonstration: Hubbard ED at machine precision

`docs/tutorials/hubbard-walkthrough.md` walks through it; `benchmarks/BENCHMARK.md` is the verification log: the 2x2 Hubbard
model at U=4 reproduces the analytic `-4√2` bit for bit; the 2x4 model
(basis dimension 4,900) and the 3x3 model agree with a from-scratch
100-line Python exact diagonalisation (`benchmarks/hubbard_ed.py`) on full
spectra; and the strong-coupling limit E₀·U → −48 is recovered. The same
checks run as ctest goldens on every build.

## Finding your way around

`docs/README.md` is the index of the knowledge base, organised by what you
want to do; the same pages are published as a searchable site at
<https://paulwatford.github.io/diagham-cmake-poc/> (Material for MkDocs,
`.github/workflows/pages.yml`). `docs/reference/build-system.md` describes every option, file
and mechanism of the build; `docs/reference/configure-flag-map.md` maps every
configure flag; `docs/how-to/build/` has the per-platform build guides and
the cluster recipe; `docs/reference/tests.md` describes the test suite. The
program-by-program reference (`docs/reference/programs/`, generated from
`--help` for all 603 programs) carries the 132 manuals adapted from the
upstream wiki at nick-ux.org; what is still unfinished is listed in
`docs/drafts/INDEX.md`.

## Contributing, support, governance

`CONTRIBUTING.md` says how to propose a change, `SUPPORT.md` where to ask,
`GOVERNANCE.md` who decides what (the DiagHam authors own the physics; this
repository owns build, tests, docs and releases), `CODE_OF_CONDUCT.md` the
Contributor Covenant, `SECURITY.md` how to report privately. Issue and pull
request templates are in `.github/`.

## AI-assisted contribution rules

AI tools were used heavily in this migration and the rules for their use are
explicit and versioned with the code: see `AGENTS.md`. In short, the build
and the tests decide pass or fail, an AI agent may describe a change and
judge its usefulness but never sets pass/fail, never edits a golden and
never merges, and DiagHam's physics code is out of bounds for automated
changes; every proposed change is held for a named human maintainer.
`CONTRIBUTING.md` has the human workflow.

## Production migration context

The agreed plan (see `docs/explanation/migration-history.md`): Git first, then CMake on top;
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
  (`docs/explanation/deferred-code.md` lists the candidates).
- Spack/EasyBuild packaging; a curated library API.

## Build dependencies

- CMake >= 3.21 (presets), a C++11 compiler (GCC or Clang), pthreads,
  Python 3 (generator and cross-check; numpy for the Python test).
- Optional: LAPACK/BLAS or Intel MKL, MPI, ScaLAPACK, GSL, GMP, MPACK,
  FFTW3, bzip2 — see `docs/reference/configure-flag-map.md`.

See `CHANGELOG.md` for what changed and when.

## License

DiagHam is licensed under the GNU General Public License, version 2 or
later; the full text is in `LICENSE`. The migration work is released under
the same terms; `NOTICE` carries the contribution notice.

## Citing

If DiagHam contributes to a publication, please cite it as software (see
`CITATION.cff`; GitHub's "Cite this repository" button reads it) and the
relevant papers in `docs/reference/publications.md`. This is release
**2026.9.0** (DiagHam r4493); releases are listed on GitHub and described
in `CHANGELOG.md`. No DOI is minted for this repository: DiagHam is the
DiagHam authors' work, and this repository is only the build, tests and
documentation around it — cite them, not this.

## Acknowledgements

The upstream code base is the work of Nicolas Regnault, Gunnar Möller,
Zlatko Papić, Cécile Repellin, Antoine Sterdyniak, Duc Phuong Nguyen,
Niall Moran, Yang-Le Wu and other contributors since 2003 (`AUTHORS`).
