# Changelog

All notable changes to this PoC are recorded here, newest first. This
repo doesn't cut version tags (it's a proof-of-concept, not a release
train), so entries are grouped by work session rather than by version
number.

## 2026-09-22 — Audit-and-correction pass (second Claude session's findings, independently re-verified)

A second, independent Claude session audited this repo's claims read-only
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
  through byte-for-byte), and on failure, byte-convert a scratch copy of
  that one patch to CRLF and retry with `patch --binary -p1` again.
  Verified three ways: standalone against all 11 real patches, the CRLF
  retry alone, and end-to-end as an actual `cmake -S . -B build` run
  against a fresh `guysoft/DiagHam` clone overlaid with this repo's real
  (git-archived) `cmake/`/`CMakeLists.txt`/`patches/` — reaches
  `Configuring done`/`Generating done` and a real target (`Vector`)
  compiles and links. `patches/PATCHES.md` corrected to match (the
  previous entry's "GNU patch 2.7.6 CRLF-detection bug" diagnosis was
  wrong; also flags a genuine, separate, still-open issue this testing
  surfaced: re-configuring against an already-patched source tree fails
  with "Reversed (or previously applied) patch detected").
- **Corrected a false "one commit" claim in `MIGRATION_ROADMAP.md` and
  `DEFERRED.md`.** Both files claimed the `guysoft/DiagHam` mirror "has
  exactly one commit," used throughout this engagement to argue the
  mirror lacked real history. That was an artifact of always using
  `git clone --depth 1` (shallow clone); a full clone shows 4,098 real
  commits spanning 2003-05-03 to 2020-10-23. `DEFERRED.md`'s per-file
  "last touched" claims, built the same way, are consequently understated
  — spot-checked for `QHEFermionsTorusWithSpin.cc`, whose real last touch
  is 2009-11-27 by `moller`, not the 2020-10-23 import-sweep date
  originally recorded. Not every entry in `DEFERRED.md` was re-derived
  this pass; the remainder are flagged there as unverified pending a
  `git log --follow` re-check against full history.
- **Corrected stale Hubbard 2x2 U=4 precision claims** in `README.md`
  and `HUBBARD_BENCHMARK.md` (`-5.6568542494924`, described as agreeing
  to "1.95×10⁻¹⁴, the limit of double-precision arithmetic"). That figure
  predates patch 08's precision fix and is actually 22 ULP off from the
  exact double-rounding of `-4*sqrt(2)` (confirmed with mpmath). A
  LAPACK-enabled binary rebuilt and rerun fresh this session gives
  `-5.6568542494923806`, a bit-identical (0 ULP) match — matching what
  `benchmarks/BENCHMARK.md`'s own Test 1 section already independently
  stated, so that file was right the whole time.
- **Added an argument-validation guard to `cmake/verify_build.sh`.**
  Confirmed (by running the real git-tracked script against nonexistent
  directories) that it silently printed "RESULTS: 0 passed, 0 failed"
  and exited 0 — a false PASS with zero checks run. Now hard-fails with
  a clear message if either build directory argument doesn't exist.

## 2026-09-22 — Patch-application fragility fix, first real configure+build (superseded above)

- **`d49f816`**: found and fixed patch 09
  (`09-FCIHofstadterCorrelationPrecision.patch`) failing to apply against
  a fresh upstream clone (previous entries below had only reasoned about
  patch application, not run `cmake` for real — no `cmake` binary was
  available in the working environment until this session; installed a
  user-local one via `pip install --user cmake`). Added a `git apply`
  fallback in `cmake/ApplyUpstreamPatches.cmake`; verified a real
  `cmake -S . -B build` reaching `Configuring done`/`Generating done` and
  a real `cmake --build build --target Vector` compiling and linking —
  the first real compiler invocation this PoC's CMake build had in this
  environment, everything before this was configure-only or reasoned
  about from source. **This fix turned out to be incomplete** — see the
  2026-09-22 audit-and-correction entry above, which replaces it. No
  CHANGELOG entry was added for this commit at the time it was made;
  this entry backfills that gap.

## 2026-09-18 — Documentation and rules pass (KT plan steps 10-12)

Pure documentation/reconciliation; no CMake or patch content changed.

- **Added `MIGRATION_ROADMAP.md`.** Reconciles this PoC's scope and
  reference material against what was actually agreed on the call with
  Gunnar Möller: Git-before/with-CMake sequencing, Kent GitLab as the
  eventual host (this repo stays a GitHub scratch space), single-repo,
  full SVN history preservation in the production conversion, and static
  linking (already the only mode this build produces). *(This entry
  originally said the `guysoft/DiagHam` mirror "has only one commit" —
  that was wrong, see the 2026-09-22 audit-and-correction entry above.)*
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
