# The review pipeline ("CMake AI Guard") — design

Purpose: why the repository's review process is shaped the way it is, and how each part maps onto CI.
Source: the "CMake AI Guard" specification (Paul Watford, September 2026) and the pipeline section of `AGENTS.md`. Changed: the spec's "79-name manifest" is expressed in terms of the ctest suite that now exists; implementation status stated plainly. Status table updated 2 November 2026 for the Guard workflow (`how-to/develop/guard.md`).

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

## Who may submit, and who must look

The specification speaks of a maintainer who approves. The Guard workflow makes the two sides of that
explicit as files on `main`: `.github/approved-developers.txt` (whose pull requests are not held for being
external) and `.github/guard-reviewers.txt` (who is asked to review a held pull request and who alone can
release it). Every pull request from anyone else, every pull request by a bot, and every pull request that
declares or shows AI assistance is held and the reviewers are asked. AI-assisted changes to DiagHam's physics
sources are reported "for maintainer review", the `GOVERNANCE.md` wording for physics-affecting changes.
Disclosure is required: the pull request template's "AI assistance" section must have a box ticked, or the
gate fails for everyone. No AI agent is on either list, now or ever.

## Status (2026-11-02)

| Part | State |
|---|---|
| Gates 1–2 without the manifest | implemented: `scripts_cmake/ci.sh` (build + ctest) via `.github/workflows/ci.yml` and `.gitlab-ci.yml` |
| Manifest check | done 2026-09-27: `tests/manifest.txt` + `selftest.manifest` (CheckManifest.cmake) |
| Damage scan | built 2026-11-02: `scripts_cmake/damage_scan.py` (the seven classes above, read with git only, exit 2 on damage; `damage_scan_selftest.sh` checks each class) run by `.github/workflows/guard.yml` as the check `guard / damage scan`, reported as one comment per pull request, held with the labels `damage` and `held-for-approval` |
| Author gate | built 2026-11-02: `.github/scripts/author_gate.py` as the check `guard / author gate`: approved-developer list, AI disclosure and detection, physics-source check, review requests |
| Description (model) | not built; the system prompt is in `AGENTS.md`. The damage scan's report already carries the deterministic parts of the mail (Damage, test-name add/remove list, diffstat) |
| Hold / approval | built 2026-11-02 as labels and required checks: `held-for-approval` until a login in `.github/guard-reviewers.txt` adds `approved-by-maintainer` or approves the head commit; `.github/rulesets/main.json` requires the two Guard checks, the five CI presets, the container job and one code-owner review on `main`. A release through a label is the `approve` reply of the specification; `reject` is closing the pull request; `restore-tests` is a review comment asking for the names and goldens back, after which the scan re-runs on the next push |
| Mail | not built; the comment on the pull request and the notification GitHub sends for the review request stand in for it |
| Gates for AI work (`CLAUDE.md`) | written 2026-11-02: gates G1 to G7, the stop-and-ask list and the vocabulary (`docs/how-to/develop/ai-gates.md`). G1 to G3, G6 and G7 exist; the numbers page of G5 (`tests/numbers_of_record.py`) is built. Built the same day: G4 (`selftest.claims`, `tests/check_claims.py`, with `tests/claims_selftest.sh`) and the claim index of G5 (`tests/claims.py`, `docs/reference/claims.md`). Every `physics` and `crosscheck` test has a `KNOWN_FROM`, and a missing one fails `selftest.claims` |
| Review requests on every new pull request | built 2026-11-02: on `opened`, `reopened` and `ready_for_review` the author gate asks every login in `.github/guard-reviewers.txt` and `.github/approved-developers.txt`, the author and bots excepted |
| Volunteers' channel | built 2026-11-02: the job `notify` posts new pull requests (and releases) to Discord when the secret `DISCORD_WEBHOOK_URL` exists; `Volunteer check:` comments from `.github/volunteers.txt` are counted and release nothing |

The earlier wording "79 checks" referred to `cmake/verify_build.sh`, the
autotools-parity script, which reports 79 passes in the default
configuration; it is a nightly parity job, not the per-change gate. The
per-change gate is the ctest suite (613 tests at r4493).
