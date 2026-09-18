# Changelog

All notable changes to this PoC are recorded here, newest first. This
repo doesn't cut version tags (it's a proof-of-concept, not a release
train), so entries are grouped by work session rather than by version
number.

## 2026-09-18 — Documentation and rules pass (KT plan steps 10-12)

Pure documentation/reconciliation; no CMake or patch content changed.

- **Added `MIGRATION_ROADMAP.md`.** Reconciles this PoC's scope and
  reference material against what was actually agreed on the call with
  Gunnar Möller: Git-before/with-CMake sequencing, Kent GitLab as the
  eventual host (this repo stays a GitHub scratch space), single-repo,
  full SVN history preservation in the production conversion (the
  `guysoft/DiagHam` mirror used as ground truth here has only one
  commit — flagged and explained), and static linking (already the
  only mode this build produces).
- **Added `AGENTS.md`.** First cut of the CMake AI Guard specification,
  versioned with the code it governs: locked facts (79/79 checks, the
  Hubbard machine-precision golden), scope, authority (ctest owns
  pass/fail, never a model), the intake -> gates -> damage-scan ->
  description -> hold pipeline, what an AI agent must never do, the
  review-mail field layout, and the approval rule. The automated
  mailer/status-check pipeline itself is not implemented yet; this file
  is its specification.
- **Added `MODULE_MAP.md`.** Task-oriented navigation ("I want to add a
  dependency," "I want to verify the physics," ...) pointing at
  specific files, on top of the existing directory-tree listing in
  `README.md`.
- README updated with short pointers to all three new documents.

## Earlier — LAPACK/MPI fix, MKL, build-type validation, Kent defaults, GSL

- **Fixed a root-cause CMake translation bug**: the LAPACK and MPI
  option blocks were injecting the downstream macro (`__LAPACK__`,
  `__MPI__`) directly via compile definitions, instead of setting the
  upstream-named CMake variable (`HAVE_LAPACK`, `HAVE_MPI`) that
  `configure_file()` substitutes into the generated config header. This
  bypassed `src/config.h`'s own gating logic (the
  `#ifdef HAVE_LAPACK -> #define __LAPACK__ -> typedef doublecomplex`
  chain). Two instances found and fixed: LAPACK (active) and MPI
  (dormant — option defaults `OFF`, so this was a latent bug waiting
  for anyone who turned MPI on).
- **Added MKL as an optional LAPACK provider** (`DIAGHAM_USE_MKL`),
  using `cmake/FindMKL.cmake` adapted from Gunnar Möller's `BDMC_UFL`
  project, feeding `DIAGHAM_LAPACK_LIBRARIES` so the rest of the build
  doesn't need to know which provider is active.
- **Added `cmake/CompilerDefaults.cmake`**: build-type validation
  (rejects unknown `CMAKE_BUILD_TYPE` values) and per-compiler
  (Intel/Clang/GNU) warning flags. Deliberately does not set a default
  build type or C++ standard, and does not duplicate the existing
  config-summary footer.
- **Added `cmake/KentDefaults.cmake`**: an explicit scaffold for
  Kent-specific HPC defaults (`DIAGHAM_KENT_DEFAULTS` option), modelled
  on the `TCMDefaults.cmake` pattern from Gunnar's prior project. All
  paths are placeholders pending real values from Gunnar or Kent's
  e-Research/RCS team; a `FATAL_ERROR` guard prevents silently building
  with empty paths.
- **Wired GSL as an optional dependency** (`DIAGHAM_USE_GSL`), and
  staged (but did not wire) the remaining `Find*.cmake` modules from
  `BDMC_UFL` (`FindNAG`, `FindMPC`, `FindMPFR`, `FindGMP`, `FindFFTW3`,
  `FindDLR`) and the `gsl_external_cblas.patch`, since DiagHam's own GSL
  usage was verified not to need the cblas coexistence patch BDMC_UFL
  needed for its different use of GSL.

## Baseline (initial PoC)

- Core CMake build for `Base/src`, `src`; later extended to `FQHE/src`
  and `FTI/src`.
- `scripts_cmake/extract_autotools.py`: generates per-subdirectory
  `CMakeLists.txt` from upstream `Makefile.am`, rather than a
  hand-maintained source list.
- 11 patches for upstream-bug surfaced during the migration, with full
  audit trails in `patches/PATCHES.md`.
- `cmake/verify_build.sh`: 79-check build-equivalence verification
  against the autotools build (core-only scope).
- Hubbard ED end-to-end demonstration matching the analytic ground
  state to machine precision, cross-checked against an independent
  100-line Python exact-diagonalisation implementation.
- One file (`QHEFermionsTorusWithSpin.cc`) documented and excluded from
  the build, with a canonical replacement identified; see `DEFERRED.md`.
