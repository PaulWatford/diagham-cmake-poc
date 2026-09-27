# The review pipeline ("CMake AI Guard") — design

Purpose: why the repository's review process is shaped the way it is, and how each part maps onto CI.
Source: the "CMake AI Guard" specification (Paul Watford, September 2026) and the pipeline section of `AGENTS.md`. Changed: the spec's "79-name manifest" is expressed in terms of the ctest suite that now exists; implementation status stated honestly.

## The idea in one sentence

Deterministic checks gate the merge, a language model may only *describe* a
change, and a named human approves — the same split serious CI already
uses; the rules exist because AI tools were used heavily in this migration
and one automated session pushed an unapproved change (2026-09-24).

## Authority

Pass/fail of the suite is owned by CMake and ctest. A model may write the
change description and an opinion on usefulness (at most eight lines,
labelled ADVICE). It cannot set a pass bit, edit a golden file, or merge.
Silence is not approval: a change stays held until a named maintainer
replies `approve`, `reject` or `restore-tests`. Physics code is out of
bounds for automated changes (`AGENTS.md`).

## Pipeline

1. **Intake.** A proposed change (pull/merge request, patch or dropped
   archive) is built in an isolated working tree; `main` is never built in
   place.
2. **Machine gates (blocking).** Out-of-source configure; build; `ctest`
   (physics goldens, `--help` smoke test, install test); the Hubbard demo
   against its locked analytic value (a ctest golden); the **manifest**
   check — the sorted list of test names (`ctest -N`) must equal the
   committed `tests/manifest.txt`, so a green run with silently dropped
   tests cannot pass.
3. **Damage scan (blocking for the mail class; never an auto-fix).**
   Flag even when ctest is green: `add_test`/`enable_testing` removed; a
   test name gone; the demo target unregistered; golden bits or `EXPECT`
   values changed; an option's default flipped to OFF; an in-source build
   reintroduced; anything under `patches/` changed.
4. **Description.** From the diff plus the gate output: files touched,
   options and cache variables added/removed/renamed, new `find_package`
   or `FetchContent` pins, compiler or generator assumptions, install and
   preset changes, tests added, skipped or renamed. It must not invent
   files that are not in the diff.
5. **Hold.** Status `HELD_FOR_APPROVAL` until a maintainer's reply.

## Locked artefacts

The test manifest (`tests/manifest.txt`, enforced by `selftest.manifest`), the golden values in
`tests/CMakeLists.txt` and the input data in `tests/data/`, the saved reference spectra in
`benchmarks/`, and the default values of the module options. Changing any
of them is Damage even if the build stays green: it is flagged and held,
never merged on a green tick.

## Status (2026-09-26)

| Part | State |
|---|---|
| Gates 1–2 without the manifest | implemented: `scripts_cmake/ci.sh` (build + ctest) via `.github/workflows/ci.yml` and `.gitlab-ci.yml` |
| Manifest check | done 2026-09-27: `tests/manifest.txt` + `selftest.manifest` (CheckManifest.cmake) |
| Damage scan | designed (a git-and-grep script with the seven classes above); to be added as a CI job |
| Description (model) | not built; the system prompt is in `AGENTS.md` |
| Hold / approval | to be a manual CI job plus branch protection on `main` (a named maintainer plays it); until then, human review of every merge request is the rule |

The earlier wording "79 checks" referred to `cmake/verify_build.sh`, the
autotools-parity script, which reports 79 passes in the default
configuration; it is a nightly parity job, not the per-change gate. The
per-change gate is the ctest suite (613 tests at r4493).
