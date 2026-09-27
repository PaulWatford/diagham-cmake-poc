# Changelog

All notable changes to the DiagHam CMake migration, newest first, in the
[Keep a Changelog](https://keepachangelog.com/) style. Releases carry
calendar versions `YEAR.MONTH.PATCH` (tag `vYEAR.MONTH.PATCH`; see
`docs/how-to/develop/release.md`); each synchronisation with upstream is
also tagged `migration-r<SVN revision>`. Earlier entries are dated only.

For the history of DiagHam itself: upstream's own `ChangeLog` covers 2003 to
April 2005 only and is preserved as `docs/history/ChangeLog-upstream`; after
that the Subversion commit log is the record (4,477 trunk commits to r4493,
all present on the `upstream` branch: `git log upstream`).

## [Unreleased]

## [2026.9.0] — 2026-09-28

The first numbered release: DiagHam r4493 with the CMake build (61
libraries, 603 programs, presets default/lapack/full/hpc/core/mkl), the
ctest suite (738 tests; 107 physics goldens with independently known
answers, 4 cross-checks, 12 regression spectra, 14 self-tests, 2 known-bug
reproducers, 599 smokes; 34 of 603 programs physics-checked, 11 more with
a regression test), the documentation set and site, the upstream
synchronisation procedure, the environment record and container recipes,
and the community/governance files. Release notes:
`docs/history/release-notes/2026.9.0.md`. Everything below this heading
down to the 2026-09-25 entry is part of it.

### Added
- `docs/explanation/verification.md` (what is proven, what is only watched,
  what is untested, by module) and `docs/explanation/checklists.md` (scored
  line by line against JOSS, FAIR4RS, The Turing Way and CERTIFY-ED, with
  evidence and the gaps).
- The upstream collaboration pack: `scripts_cmake/upstream_pack.py` generates,
  from the commits with an `Upstream-Patch:` trailer, one report
  (`docs/upstream-reports/`) and one `svn patch`-applicable diff
  (`patches/upstream/`, upstream file names, line endings kept) per fix,
  15 in all, verified to apply in series to the pristine `upstream` tree;
  hand-written reproducer reports for the open defects U29–U31;
  `docs/how-to/develop/send-a-fix-upstream.md`.
- Reproducibility: `scripts_cmake/environment.sh` records commit, revision,
  CMake options, toolchain and library versions and `TestDiagHamConf`; CI
  keeps it with the ctest log as an artefact of every run;
  `container/Dockerfile` (full preset, whole test suite, installed runtime
  image) built and run by CI; `container/diagham.def` for Apptainer;
  `docs/how-to/build/container.md`.
- `scripts_cmake/sync_upstream.sh`: fetch new Subversion revisions into
  the git-svn clone, rebase them onto `upstream` with tree verification,
  merge into `main`, regenerate the CMake files, record the revision,
  build and test; `scripts_cmake/svn-authors.txt` (names only, `yangle`
  corrected) in the repository; the weekly drift job opens an
  `upstream-drift` issue; `docs/how-to/develop/sync-upstream.md`.
- The documentation site: `mkdocs.yml` (Material for MkDocs, MathJax,
  search), navigation generated from the tree by
  `scripts_cmake/gen_mkdocs_nav.py`, published to GitHub Pages by
  `.github/workflows/pages.yml` on every push to `main` and built with
  `--strict` on pull requests; `docs/reference/programs/manuals/README.md`
  (index of the attached manuals); `docs/how-to/develop/docs-site.md`.
- Community and governance files: `GOVERNANCE.md` (two owners: the DiagHam
  authors for the physics, the maintainer for build/tests/docs/releases; how
  a fix travels upstream; the GitHub settings kept), `SUPPORT.md`,
  `SECURITY.md`, `CODE_OF_CONDUCT.md` (Contributor Covenant 2.1), issue
  templates (bug / feature or test / documentation), a pull-request
  checklist, `CODEOWNERS`.
- 12 regression tests (`regression.*`), one spectrum saved from the r4493
  build per module directory without a physics golden and for a few more
  families; `tests/regression_reference.py` and `cases.txt`; the references
  say what they are (change detection, not correctness). `check_spectrum
  numbers` for programs that print rather than write. Second `known-bug`
  test: the two-band checkerboard segfault (U31). CI fails when the coverage
  page is stale.
- Independent solvers: `tests/oracles/sphere_ed.py` (lowest-Landau-level
  pseudopotential Hamiltonians from Clebsch–Gordan projectors, 7 full-spectrum
  tests, fermions and bosons, Coulomb and generic pseudopotentials) and
  `tests/oracles/spin_ed.py` (XXZ / J₁–J₂ chains with a field, 6 tests).
  Algorithm-consensus cross-checks (Lanczos vs full diagonalisation; LAPACK
  vs internal). Harness self-tests: the checker rejects wrong values, counts,
  missing files, NaN and wrong-length spectra; the runner rejects missing
  output. `check_spectrum lowest` mode.
- 15 spin-chain and Hubbard goldens (`physics.spin.*`, `physics.hubbard.*.U0.*`):
  XX, AKLT and Haldane–Shastry rings against closed forms, an open AKLT
  chain's four edge states, Hubbard at U = 0 against the tight-binding sum,
  Heisenberg (spin ½ and 1) and transverse-field Ising rings against an
  independent numpy diagonalisation (`tests/oracles/spin_hubbard.py`).
- First `known-bug` test: plain Lanczos below the ground state on a
  degenerate spectrum (U30).
- 22 torus / cylinder / disk goldens (`physics.fqhe.{torus,cylinder,disk}.*`):
  Laughlin and Moore–Read topological degeneracies and quasihole counts on
  the torus, cylinder zero modes, disk edge-mode counting p(ΔL) — answers
  from `tests/oracles/fqhe_geometries.py`, re-derived by
  `selftest.fqhe_geometries_oracle`.
- 20 FQHE-sphere goldens (`physics.fqhe.sphere.*`): zero-mode counts of
  the Laughlin, Moore–Read and Read–Rezayi states from (k,r)-admissible
  counting, Coulomb pseudopotentials against the Wigner-3j/6j closed form,
  Jack polynomials against exact diagonalisation (overlap 1), and the
  particle entanglement spectrum of the Laughlin state against quasihole
  counting — all answers from `tests/oracles/fqhe_sphere.py`, re-derived
  by `selftest.fqhe_sphere_oracle`. `check_spectrum` gained `line` and
  `nonzero`; `diagham_physics_test` gained `PRE_STEPS` for chains.
  23 of 603 programs now have a physics or cross-check test.
- 27 Hilbert-space dimension goldens (`physics.dimension.*`) across sphere,
  torus (magnetic-translation sectors), disk, single-band lattice, Hubbard
  and spin chains, checked against exact counting in
  `tests/oracles/dimensions.py`; `selftest.dimension_oracle` re-derives the
  expected files.

### Removed
- The Zenodo deposit (`.zenodo.json`, the DOI step of the release checklist):
  DiagHam is the DiagHam authors' work and a citable archive of it is theirs
  to create; `CITATION.cff` now says to cite them, not this repository.

### Fixed
- Mid-cycle audit (fresh clone, presets core/full/default): the `core`
  preset failed to configure (`set_tests_properties` on tests not
  registered when a module is off); the manifest self-test failed in the
  `full` preset (option-gated tests, now `tests/manifest-optional.txt`);
  the FQSH two-band checkerboard regression case depended on the
  eigensolver (degenerate bands) and was replaced by `FTI3DSimpleTI`.
- `ParticleOnCylinderPseudopotentialHamiltonian` without GSL: undefined
  behaviour (every `FQHECylinderFermionsTwoBodyGeneric` run crashed with
  "stack smashing detected") replaced by a clear message and exit 1 (U28).
- **CI built the wrong thing.** `scripts_cmake/ci.sh` still fetched the
  guysoft mirror (r4114) and overlaid the CMake files on it — the
  proof-of-concept flow — so both CI runs after the r4493 push failed. It
  now configures, builds and tests the repository in place; `overlay.py` is
  retired. The weekly job now checks whether the DiagHam Subversion trunk
  has moved past `scripts_cmake/upstream-revision.txt`.
- `upstream`: 30 commits carried the raw SVN username `yangle`; now Yang-Le
  Wu, as in the authors file. Trees unchanged; every SHA on `upstream` and
  the merge on `main` changed (re-clone if you had the earlier push).

### Changed
- Test labels tell the truth: `physics` (independently known answer),
  `crosscheck` (two DiagHam implementations agree), `regression` (saved
  spectrum reproduced), `smoke`, `selftest`. The two spinful-torus tests are
  `crosscheck.*`. `tests/coverage.py` generates
  `docs/reference/test-coverage.md` (7 of 603 programs have a physics or
  cross-check test, 4 against an independently known answer); README, `tests.md` and `AGENTS.md` quote that instead
  of "613 tests pass".
- `tests/manifest.txt` and `selftest.manifest`: the registered tests must
  equal the manifest (the guard design's missing artefact).
- `COPYING` → `LICENSE` (GPL-2 text) and the contribution notice →
  `NOTICE`, so the licence is detected by hosting sites.

### Documentation
- Root files brought up to date for the migration (README, LICENSE, AUTHORS,
  CITATION.cff, CONTRIBUTING.md); this changelog restructured; upstream's
  INSTALL, ChangeLog, TODO and NEWS moved to `docs/history/`.

## [2026-09-25] — main assembled on DiagHam r4493

### Added
- The complete DiagHam history: `upstream` branch = `git svn` conversion of
  the canonical Subversion repository (trunk, branches `testing_nr`,
  `testing_zp`, tag `first-release`; 4,477 commits, 2003-05-03 to
  2026-09-18, r4493), merged into `main`. Authors by name; the history
  stores no e-mail addresses. The converted r4114 tree is byte-identical to
  the GitHub mirror `guysoft/DiagHam` that earlier work was based on.
- Spin and QuantumDots modules, every configure flag
  (`CONFIGURE_FLAGS.md`), `diagham_external_deps`, install/export with
  `find_package(DiagHam)`, `CMakePresets.json`, the ctest suite
  (`tests/`, `TESTING.md`), `scripts_cmake/ci.sh`, GitHub and GitLab CI
  files, `HPC.md` — imported from an audited external contribution
  (PR #1, commit b414895; audit in the project notes) and corrected below.
- Generated per-directory `CMakeLists.txt` files committed in-tree (88 at
  r4493), including the 19 `QuantumDots/src/Tools/Analysis` programs the
  previous generator had silently dropped.
- Upstream fixes as individual commits with `Upstream-Patch:` /
  `Upstream-Base:` trailers: patches 01–07, 09–12, 14; 13 regenerated for
  r4493 (one file was fixed upstream after 2020); 08 regenerated and
  re-scoped (978 `precision(14)` sites in 544 files; deliberate display
  widths left alone). Patch 15 dropped: upstream fixed it after 2020.
- `.gitattributes` keeping patches and goldens LF on every platform.

### Changed
- `extract_autotools.py`: one discovery policy for every module — a
  directory emits programs iff its `Makefile.am` has `bin_PROGRAMS`, plus an
  explicit `EXTRA_PROGRAMS` allow-list; works on Windows (`as_posix`);
  duplicate `bin_PROGRAMS` entries collapsed.
- Top-level CMake: in-source builds refused; `CMAKE_CXX_EXTENSIONS ON`
  (gnu++11, as upstream's g++); Release `-O2` and RelWithDebInfo /
  RelWithChecks `-O2 -g` with no `NDEBUG` (autotools' defaults); the
  `__64_BITS__`, `__SMP__` and `__DEBUG__` injections removed (upstream's
  `config.h` defines them); static libraries no longer forced PIC.
- Patches are no longer applied at configure time
  (`cmake/ApplyUpstreamPatches.cmake` removed); the one program exclusion
  (`QHEFermionsTorusWithSpin`) is a cache variable.
- Two FQHEOnDisk programs renamed because their names collided with
  different FQHEOnSphere programs and the install clobbered one of each
  pair: `QHEBosonsDelta` → `FQHEDiskBosonsDelta`, `QHEFermionsOverlap` →
  `FQHEDiskLaughlinMonteCarloOverlap`.

### Fixed
- `FQHESphereApplyCreationOperatorToVector.cc`: `<limits>` placed after the
  first include (the file has includes mid-file).
- `Spin/src/Programs/PairHoppingModelFSA.cc`: one `LapackDiagonalize` call
  upstream left unguarded.
- `EvaluateBroadening` excluded from the `--help` smoke test: it checks for
  its required `--input` before `--help` (upstream quirk).

### Result
- Default preset on r4493: 61/61 libraries, 603/603 programs, ctest
  613/613 (13 physics goldens, Python cross-check, install consumer,
  `--help` on every program). `lapack` preset: 0 errors, same suite green.

## [2026-09-24] — external contribution (PR #1), audited, closed

- A draft PR adding Spin/QuantumDots, all flags, ctest, install, presets and
  CI was opened against this repository by an automated session without
  instruction. It was audited in ten rounds; the tests and the spinful torus
  Coulomb fix were verified correct, several defects were found (19
  programs dropped, generator broken on Windows, install collisions,
  unfixed top-level regressions, no review gates in CI). The PR was closed
  and its branch deleted; the verified parts were re-landed above, on the
  current upstream, as separate commits.

## [2026-09-22] — Audit-and-correction pass (findings of a second independent audit, re-verified)

A second independent audit examined this repo's claims read-only
and reported findings; per instruction, none of those findings were taken
on trust — every claim checked here was independently re-derived with
real tools (mpmath at 50-digit precision, fresh non-shallow git clones,
real `cmake`/`patch`/compile/run cycles), not just re-asserted from the
audit text. This surfaced one place where my own prior fix (the 22/09
`d49f816` commit, see below) was itself broken, and one place where I
found a more precise root cause than either my original diagnosis or the
auditor's.

- **Fixed `cmake/ApplyUpstreamPatches.cmake` for real this time.** The
  `d49f816` fix (patch, then fall back to `git apply`) was verified
  against a Windows checkout's CRLF copies of the patches, not what git
  actually has stored — against the real LF-committed blobs, `git apply`
  fails on patch 09 for the same reason plain `patch` does, so the
  previous fix did not actually work. Root cause, confirmed directly:
  all 11 patches are LF as committed; patch 09's one target file
  (`FTI/src/Programs/FCI/FCIHofstadterCorrelation.cc`) is genuinely CRLF
  in canonical upstream DiagHam, so an LF patch can never match it, under
  any tool. New strategy: `patch --binary -p1` first (gets 10 of 11
  through), and on failure, byte-convert a scratch copy of that one patch
  to CRLF and retry. *(Superseded 2026-09-25: patches are commits now.)*
- **Corrected a false "one commit" claim in `MIGRATION_ROADMAP.md` and
  `DEFERRED.md`.** Both files claimed the `guysoft/DiagHam` mirror "has
  exactly one commit"; that was an artifact of always using
  `git clone --depth 1`. A full clone shows 4,098 real commits spanning
  2003-05-03 to 2020-10-23 (= SVN r4114; the mirror has not moved since).
- **Corrected stale Hubbard 2x2 U=4 precision claims** in `README.md` and
  `HUBBARD_BENCHMARK.md` (`-5.6568542494924`, "1.95×10⁻¹⁴"): that figure
  predates patch 08. A rebuilt binary gives `-5.6568542494923806`, a
  bit-identical (0 ULP) match.
- **Added an argument-validation guard to `cmake/verify_build.sh`** (it
  printed "0 passed, 0 failed" and exited 0 on nonexistent directories).

## [2026-09-22] — Patch-application fragility fix, first real configure+build (superseded above)

- **`d49f816`**: found and fixed patch 09 failing to apply against a fresh
  upstream clone (previous entries had only reasoned about patch
  application, not run `cmake` for real). Added a `git apply` fallback;
  verified a real `cmake -S . -B build` and a real `cmake --build build
  --target Vector`. **This fix turned out to be incomplete** — see the
  audit-and-correction entry above.

## [2026-09-18] — Documentation and rules pass (KT plan steps 10-12)

Pure documentation/reconciliation; no CMake or patch content changed.

- **Added `MIGRATION_ROADMAP.md`.** Reconciles this PoC's scope and
  reference material against what was agreed on the call with Gunnar
  Möller: Git-before/with-CMake sequencing, Kent GitLab as the eventual
  host, single-repo, full SVN history preservation, static linking.
- **Added `AGENTS.md`.** First cut of the CMake AI Guard specification,
  versioned with the code it governs.
- **Added `MODULE_MAP.md`.** Task-oriented navigation.
- README updated with short pointers to all three new documents.

## [2026-09-18] — LAPACK/MPI fix, MKL, build-type validation, Kent defaults, GSL

- **Fixed a root-cause CMake translation bug**: the LAPACK and MPI option
  blocks were injecting the downstream macro (`__LAPACK__`, `__MPI__`)
  directly via compile definitions, instead of setting the upstream-named
  CMake variable (`HAVE_LAPACK`, `HAVE_MPI`) that `configure_file()`
  substitutes into the generated config header. This bypassed
  `src/config.h`'s own gating (the `doublecomplex` typedef), which was the
  real cause of the "5000 errors with LAPACK on" that earlier notes blamed
  on upstream.
- **Added MKL as an optional LAPACK provider** (`DIAGHAM_USE_MKL`), using
  `cmake/FindMKL.cmake` adapted from Gunnar Möller's `BDMC_UFL` project.
- **Added `cmake/CompilerDefaults.cmake`**: build-type validation and
  per-compiler warning flags (adapted from BDMC_UFL).
- **Added `cmake/KentDefaults.cmake`**: a scaffold for Kent-specific HPC
  defaults, modelled on the `TCMDefaults.cmake` pattern; all paths are
  placeholders pending real values.
- **Wired GSL as an optional dependency** (`DIAGHAM_USE_GSL`); staged the
  remaining `Find*.cmake` modules from `BDMC_UFL`.

## Baseline (initial proof of concept, 2026-05 to 2026-09)

- Core CMake build for `Base/src`, `src`; later extended to `FQHE/src`
  and `FTI/src`.
- `scripts_cmake/extract_autotools.py`: generates per-subdirectory
  `CMakeLists.txt` from upstream `Makefile.am`, rather than a
  hand-maintained source list.
- 11 patches for upstream bugs surfaced during the migration, with full
  audit trails in `patches/PATCHES.md`.
- `cmake/verify_build.sh`: build-parity verification against the autotools
  build (79 checks in the default FQHE+FTI configuration).
- Hubbard ED end-to-end demonstration matching the analytic ground state
  to machine precision, cross-checked against an independent 100-line
  Python exact-diagonalisation implementation.
- One file (`QHEFermionsTorusWithSpin.cc`) documented and excluded from
  the build, with a canonical replacement identified; see `DEFERRED.md`.
