# AGENTS.md — rules for AI-assisted work on this repository

This file governs any AI agent (a coding assistant editing this repo
directly, or the automated "CMake AI Guard" reviewing suggested
uploads) working on `PaulWatford/diagham-cmake-poc`. It is a first cut:
it encodes the rules from the Guard specification so they're versioned
with the code they govern, not just in a separate document. The
automated Guard pipeline itself (self-hosted runner, mailer, GitHub
status check) is not yet built — see "Build order" below for what
exists today versus what's planned.

## Locked facts

These are the two claims this repository stands on. **No agent may
weaken either of them**, in code, in tests, or in what it reports about
them:

- **The ctest suite passes in full** (`ctest --preset default`: 905 tests at
  DiagHam r4493, of which 250 physics goldens + 1 Python cross-check, 3-4
  cross-checks (one needs LAPACK), 11-12 regression (one needs GSL), 599 `--help` smokes, 1 install, 24 self-tests, 17 known-bug; the
  list of test names is the locked manifest `tests/manifest.txt`, checked by
  `selftest.manifest`). `cmake/verify_build.sh`, the autotools-parity script,
  reports 79 passes in the default configuration (FQHE and FTI on) and is a
  separate, nightly check.
- **The Hubbard ED demo matches the analytic ground state to machine
  precision** (see `docs/tutorials/hubbard-walkthrough.md` and the log `benchmarks/BENCHMARK.md`).

If a change affects either number, that is Damage (see below) even if
every check still reports green.

> **Resolved 2026-09-26:** the notes below described the proof-of-concept state; a ctest suite now exists (`docs/reference/tests.md`) and the 79 figure was reproduced against a real autotools build — in the default configuration, not core-only. Kept for the record.
>
> **Flagged, not corrected, 22/09 (a second independent audit):**
> two things about the "79/79" fact above do not currently hold and
> should be re-verified by whoever next has a real autotools toolchain
> available, before anyone relies on it as enforced:
> 1. **There is no `ctest`/`enable_testing()` wiring anywhere in this
>    build** (`grep -c "enable_testing\|add_test" CMakeLists.txt`
>    returns 0, confirmed directly). `ctest` reports 0 tests. The
>    "manifest of 79 test names" this file and the pipeline spec below
>    describe does not exist as `add_test` entries — `verify_build.sh`
>    is a hand-rolled shell script, run directly, not through `ctest`.
> 2. **The 79 figure itself is unconfirmed.** A second independent audit
>    reported getting 73 real passes (0 fails) running `verify_build.sh`
>    core-only against a real autotools build, and that 79 only appears
>    when 6 FTI libraries are also present — which would contradict the
>    stated `FTI=OFF` scope. This session could not independently
>    reproduce either number: `autoconf`/`automake`/`libtool` are not
>    installed in this sandbox and there is no root to add them, so no
>    real autotools build was available to compare against. Until
>    someone runs both builds side by side and gets a repeatable number,
>    treat "79" as **not currently enforced by anything** rather than as
>    a locked, machine-checked fact — the enforcement this file describes
>    (step 1: "Worktree + cmake + ctest + 79-name manifest") does not
>    exist yet, only `verify_build.sh` does.
> Also worth knowing: `verify_build.sh`'s checks are existence checks,
> per-library `nm`-based text-symbol *counts*, and one program's stdout
> comparison — not the byte-for-byte / binary comparison this file and
> `README.md` describe elsewhere.

## Purpose

Any AI review of a suggested change here is a **maintainer mailer and
holding queue**, not a merge robot and not a second DiagHam physics
engine. Its job is to make it possible for a human maintainer to see, in
one place: what the change does to the build system and tests, whether
it risks the 79-check suite or the Hubbard golden, whether it's useful
for the CMake port, and a clear hold — nothing reaches `main` until a
named maintainer replies `approve`.

## Scope

**In scope** for AI-assisted review or editing: `CMakeLists.txt`, CMake
presets, `cmake/Find*.cmake` modules, toolchain files, `ctest` wiring,
demo targets, goldens, CI workflows, and documentation of how to
configure and build.

**Out of scope**: computing new eigenvalues or other physics results,
changing DiagHam physics defaults, approving anything because "an LLM
likes it," and mixing review gates in from any other project.

## Authority

**Pass/fail of the suite is owned by CMake + ctest, not by any model.**
An AI agent may write a description of a diff and an opinion on its
usefulness. It may not set a pass bit, edit a golden file, or merge
anything.

Required human action on any suggested change: approve, reject, or ask
the author to restore tests. Silence is not approval — the hold lasts
until a named maintainer replies.

## Pipeline

1. **Intake.** A suggested upload (PR, patch, or dropped zip) is copied
   into an isolated worktree. `main` is never built in place.
2. **Machine gates (blocking).**
   - `cmake -S . -B build` (out-of-source)
   - `cmake --build build`
   - `ctest --test-dir build --output-on-failure`
   - Hubbard ED demo vs. the locked analytic golden, machine precision
   - Manifest check: the set of test names must match the locked
     list exactly, even if ctest still reports green
3. **Damage scan (blocking for the mail class, never an auto-fix).**
   Flag even when ctest is still green:
   - `add_test` / `enable_testing` lines removed
   - a test name disappears from the suite
   - the demo target is no longer registered with CTest
   - golden file bits changed
   - an option now defaults `OFF` so a DiagHam tree is silently unbuilt
   - an in-source build is reintroduced
4. **Code analysis (description).** Diff the incoming tree against the
   last approved `main` and produce a plain-language change list: files
   touched; new/removed/renamed CMake options and cache variables; new
   `find_package`/`FetchContent` pins; compiler or generator
   assumptions; install rules, export files, presets; any test added,
   skipped, or renamed. **The description must come from the diff plus
   the gate output — it must not invent files that aren't in the diff.**
5. **Hold.** The worktree stays quarantined. `main` and the protected
   branch are unchanged. Status reads `HELD_FOR_APPROVAL` until a
   maintainer sends the approve token or clicks the held Action.

## What an AI agent must never do here

- Auto-merge anything because a description sounded positive.
- Rewrite a golden file to make ctest green.
- Count the suite as passing if test names were dropped and new names
  filled the count back up.
- Run DiagHam itself as a replacement for ctest.
- Upload large object files or unpublished notes.
- Treat silence as approval.

## Usefulness, defined

Useful: easier to configure or build DiagHam with CMake, support for
more compilers, clearer options, extra tests that add to the 79. Not
useful: style-only changes, vendored junk, physics code that isn't
wired to ctest, or a green run that quietly dropped tests.

## Review-mail fields (for the future automated Guard)

Once the mailer pipeline exists, one email goes out per incoming
change (not per commit), to the people who run the repo, with these
fixed fields:

| Field | Content |
|---|---|
| Break | `PASS n/n` / `FAIL <first test>` / `BUILD BROKEN` |
| Damage | none / tests removed / checks skipped / golden edited / option default flipped |
| Demo | Hubbard ED vs golden OK / DRIFT / NOT RUN |
| Change description | What the upload will cause if merged, written from the diff. No physics claims. |
| Usefulness advice | Max eight lines, labelled `ADVICE` |
| Hold | `HELD_FOR_APPROVAL` until named maintainer replies `approve` \| `reject` \| `restore-tests` |
| Attachments | `LastTest.log` tail, `CMakeError.log` if configure died, test-name add/remove list, short diffstat |

Subject line format: `[diagham-cmake-poc] HOLD | <short title> | PASS n/n or FAIL k/n or BUILD BROKEN or TESTS DROPPED`

## System prompt for the reviewing model (use as written, once wired up)

> You explain diagham-cmake-poc CI and diffs. You describe what a patch
> will change in the CMake build and tests. You do not compute
> eigenvalues. You do not approve merges. If a number is not in ctest
> output or the locked golden, say it is not in the log. Prefer CMake
> cache variables, `find_package` failures, and missing `add_test`
> lines. Eight lines of usefulness advice, then stop.

## Approval

A change leaves hold only when a maintainer on the (separately kept,
non-public) maintainer list replies with one of:

- `approve` — merge allowed; quarantine discarded after merge
- `reject` — close; worktree deleted
- `restore-tests` — author must put the test names and golden back; the
  Guard re-runs

No AI agent is on that list, now or ever.

## Locked artefacts

- The manifest of test names (`ctest -N`, committed as `tests/manifest.txt`
  and enforced by `selftest.manifest`).
- Hubbard ED demo golden (analytic ground state, machine-precision
  tolerance).
- The regression references in `tests/data/regression/` (regenerated only
  for an understood, intended change, said so in the commit).
- Default CMake options that decide which DiagHam trees are built.

Changing any locked artefact is Damage even if the build stays green.
It still gets flagged and still holds.

## Build order (what exists vs. what's planned)

This repository today has the CMake build, the 79-check
`verify_build.sh`, and the Hubbard golden — i.e., step 1 below. The
automated mailer/hold pipeline (steps 2-5) is not yet implemented; this
file is the specification an implementation must follow when it is
built.

1. Worktree + `cmake` + `ctest` + test-name manifest on a self-hosted
   runner. *(exists — this is what `verify_build.sh` already checks
   locally; running it on a self-hosted runner with an isolated
   worktree is the remaining piece.)*
2. Damage scan and email template with an empty advice line.
3. A local LLM fills in the Change description and Usefulness advice.
4. Hold label + required GitHub status check.
5. Maintainer reply releases the hold.

If step 1 is not green and mailed on every failure, do not add the
model yet.

## Success looks like

A PR that breaks one of the 79, drifts the Hubbard demo, or deletes a
test name is held, and the maintainers get one message that states the
break, the damage, a plain description of the code change, and a
usefulness line that never pretends to be a merge decision.


The design rationale and implementation status of this pipeline: `docs/explanation/guard-design.md`.
