# DiagHam documentation

Purpose: the index of everything under `docs/`, organised by what you want to do; every page is at most one click from here. Published as a site at <https://paulwatford.github.io/diagham-cmake-poc/>; the repository's root files ([README](https://github.com/PaulWatford/diagham-cmake-poc/blob/main/README.md), [CONTRIBUTING](https://github.com/PaulWatford/diagham-cmake-poc/blob/main/CONTRIBUTING.md), [GOVERNANCE](https://github.com/PaulWatford/diagham-cmake-poc/blob/main/GOVERNANCE.md), [CHANGELOG](https://github.com/PaulWatford/diagham-cmake-poc/blob/main/CHANGELOG.md)) are on GitHub.
Source: new (2026-09-27), replacing the upstream wiki's hand-kept Main_Page. Layout follows the Diátaxis convention: tutorials (learning), how-to guides (doing), reference (looking up), explanation (understanding), plus `history/` (upstream files kept as found) and `drafts/` (pages that are not finished, with what is missing).

Every page starts with a `Purpose:` and a `Source:` line saying where its content came from (upstream wiki, an earlier repository document, the build itself) and what was changed. Program pages are generated from `--help` and are not edited by hand.

## I want to…

### …build DiagHam

| Situation | Page |
|---|---|
| I have never built it; take me from clone to a verified build | [tutorials/first-build-and-verify.md](tutorials/first-build-and-verify.md) |
| Which build guide is for my platform? | [how-to/build/README.md](how-to/build/README.md) |
| Linux with GCC (the default) | [how-to/build/linux-gcc.md](how-to/build/linux-gcc.md) |
| Linux with Clang | [how-to/build/linux-clang.md](how-to/build/linux-clang.md) |
| Intel compilers and MKL | [how-to/build/intel-mkl.md](how-to/build/intel-mkl.md) |
| Apple silicon Mac | [how-to/build/macos.md](how-to/build/macos.md) |
| A cluster with a module system, MPI and ScaLAPACK | [how-to/build/hpc-cluster.md](how-to/build/hpc-cluster.md) |
| What does each preset turn on? | [how-to/build/presets.md](how-to/build/presets.md) |
| I used to type `./configure --enable-…`; what is the CMake equivalent? | [reference/configure-flag-map.md](reference/configure-flag-map.md) |
| Every option, file and mechanism of the build | [reference/build-system.md](reference/build-system.md) |
| Install it and link my own code against it (`find_package(DiagHam)`) | [how-to/install-and-use-as-library.md](how-to/install-and-use-as-library.md) |
| Build with `long double` precision | [how-to/develop/long-double.md](how-to/develop/long-double.md) |
| Run DiagHam from a container, or record exactly what a build was made with | [how-to/build/container.md](how-to/build/container.md) |
| Why CMake, and why generated? | [explanation/why-cmake.md](explanation/why-cmake.md) |

### …run a calculation

| Situation | Page |
|---|---|
| My first physics run: the Laughlin state on the sphere | [tutorials/first-fqhe-run.md](tutorials/first-fqhe-run.md) |
| The Hubbard model, with results checked against known answers | [tutorials/hubbard-walkthrough.md](tutorials/hubbard-walkthrough.md) |
| Which program does what? (all 603, by module, from `--help`) | [reference/programs/README.md](reference/programs/README.md) |
| The longer manuals adapted from the wiki (132 programs) | [reference/programs/manuals/](reference/programs/README.md) — linked as "manual" in the program index |
| What do the output and input files look like; convert binary vectors and matrices | [reference/data-formats.md](reference/data-formats.md) |
| Run under MPI, or mixed MPI/SMP, and profile it | [how-to/run-mpi.md](how-to/run-mpi.md) |
| Use ScaLAPACK for full diagonalisation on a parallel machine | [how-to/run-scalapack.md](how-to/run-scalapack.md) |
| Full diagonalisation vs Lanczos; excited states with projectors | [drafts/lanczos-diagonalisation-schemes.md](drafts/lanczos-diagonalisation-schemes.md) (draft: two sections empty on the wiki) |
| Is the result I got affected by a known defect? | [reference/known-defects.md](reference/known-defects.md) |
| Papers that used DiagHam; how to cite it | [reference/publications.md](reference/publications.md), `CITATION.cff` at the root |

### …change or extend the code

| Situation | Page |
|---|---|
| Bring new DiagHam Subversion revisions in (update `upstream` and `main`) | [how-to/develop/sync-upstream.md](how-to/develop/sync-upstream.md) |
| Add my own program, in the tree or as a separate project | [how-to/develop/create-a-program.md](how-to/develop/create-a-program.md) |
| Add a spin Hilbert space with a discrete symmetry | [how-to/develop/new-spin-hilbert-space-with-symmetry.md](how-to/develop/new-spin-hilbert-space-with-symmetry.md) |
| Preview or publish this documentation as a site | [how-to/develop/docs-site.md](how-to/develop/docs-site.md) |
| Add a test that checks a spectrum against a known answer | [how-to/develop/add-a-golden-test.md](how-to/develop/add-a-golden-test.md) |
| What the test suite checks, and what it does not | [reference/tests.md](reference/tests.md) |
| Which programs actually have their physics tested (honest count) | [reference/test-coverage.md](reference/test-coverage.md) |
| Fix something in DiagHam's own sources: the rules and the index of fixes so far | `CONTRIBUTING.md` at the root, then [explanation/upstream-fixes.md](explanation/upstream-fixes.md) |
| Offer a fix or a finding to the DiagHam authors (reports + `svn patch` diffs) | [how-to/develop/send-a-fix-upstream.md](how-to/develop/send-a-fix-upstream.md), [upstream-reports/](upstream-reports/README.md) |
| Code that is excluded from the build, and why | [explanation/deferred-code.md](explanation/deferred-code.md) |
| The torus SU(2) Coulomb defect in detail | [explanation/torus-su2-coulomb-defect.md](explanation/torus-su2-coulomb-defect.md) |
| How changes are reviewed (the guard pipeline) and the rules for AI-assisted work | [explanation/guard-design.md](explanation/guard-design.md), `AGENTS.md` at the root |

### …understand where this repository comes from

| Situation | Page |
|---|---|
| How this repository relates to the canonical DiagHam (SVN), what was agreed, what is done | [explanation/migration-history.md](explanation/migration-history.md) |
| What changed, release by release | `CHANGELOG.md` at the root; release notes in [history/release-notes/](history/release-notes/2026.9.0.md) |
| Make a release (version, tag, GitHub Release, Zenodo) | [how-to/develop/release.md](how-to/develop/release.md) |
| The upstream INSTALL, ChangeLog, TODO and NEWS, kept as found | [history/README.md](history/README.md) |
| Who wrote DiagHam | `AUTHORS` at the root |
| Who decides what, and how a fix reaches the DiagHam authors | `GOVERNANCE.md` at the root; `SUPPORT.md`, `CODE_OF_CONDUCT.md`, `SECURITY.md` |

### …finish something

| Situation | Page |
|---|---|
| What is incomplete, what is missing, who can complete it; where every retired page went | [drafts/INDEX.md](drafts/INDEX.md) |
| The Kent cluster recipe (site values unknown) | [drafts/kent-cluster-build.md](drafts/kent-cluster-build.md) |
| Wiki program pages too thin to attach as manuals | [drafts/manual-stubs/INDEX.md](drafts/manual-stubs/INDEX.md) |
| Programs that have a `--help` page but no manual | [drafts/programs-without-manual.md](drafts/programs-without-manual.md) |

## Layout

```
docs/
  tutorials/     learning by doing: first build, first FQHE run, Hubbard walkthrough
  how-to/        build/ (per platform, presets), develop/ (programs, tests, symmetries), MPI, ScaLAPACK, install
  reference/     build system, configure-flag map, data formats, tests, defects, publications, programs/ (generated) + manuals/
  explanation/   migration history, upstream fixes, deferred code, the torus defect, guard design, why CMake
  history/       upstream INSTALL, ChangeLog, TODO, NEWS — unchanged
  drafts/        incomplete pages and the index of what went where
```

## Keeping it truthful

- Counts quoted in the docs (61 libraries, 603 programs, 738 tests, 107 physics goldens) are those of the r4493 build; [reference/build-system.md](reference/build-system.md) is where they are defined and the other pages follow it.
- `docs/reference/programs/` is regenerated with `python3 scripts_cmake/gen_program_reference.py build/lapack docs/reference/programs` after any change to a program's options; the manuals are re-attached with `scripts_cmake/attach_wiki_manuals.py`.
- A page that cannot be verified against the build goes to `drafts/` with a note on what is missing, not into the sections above.
