# CLAUDE.md: gates for any AI working on this repository

This file is for every AI coding tool that edits this repository (Claude Code reads it directly; other tools are
pointed at it from `AGENTS.md`). It is written so that the same file, with only the command table changed, can be
dropped into the ALPS repository when the ALPS team wants the same gates. The rules in `AGENTS.md` apply in full;
this file is the checklist that turns them into pass/fail.

@AGENTS.md

## The one rule

**Nothing an AI writes here is done until every gate below has been run and has passed, and the run is shown.**
"Shown" means the command and its last lines pasted into the pull request or the hand-off note, verbatim. A gate
that was not run is reported as NOT RUN, never as passed. A number that is not in a command's output or in a file
of this repository is reported as "not read".

## Gates

Run them in this order, from a clean, out-of-source build. Every one must exit 0.

| # | Gate | DiagHam command | What it proves |
|---|---|---|---|
| G1 | Configure and build, out of source | `cmake --preset default && cmake --build --preset default` | the tree still builds |
| G2 | Full suite with the manifest | `ctest --preset default --output-on-failure` | every test still passes, and `selftest.manifest` proves the set of test names equals `tests/manifest.txt` |
| G3 | Damage scan | `python3 scripts_cmake/damage_scan.py main HEAD` | no locked artefact changed (test names, goldens, `tests/data/`, `benchmarks/`, option defaults, `patches/`, in-source build). Exit 2 is allowed only if every finding is listed under "Damage, explained" in the pull request |
| G4 | Claims are sealed | `ctest --preset default -R selftest.claims` | every `physics` test states where its answer is known from, carries a tolerance, and uses only the labels in `tests/labels.txt` (the kinds are in the vocabulary below) |
| G5 | Generated pages are current | `python3 tests/coverage.py build/default --write docs/reference/test-coverage.md && python3 tests/claims.py build/default --write docs/reference/claims.md && python3 tests/numbers_of_record.py build/default --check` | the coverage page, the claim index and the canonical numbers page match this build; a stale page is a failure |
| G6 | Documentation builds | `python3 scripts_cmake/gen_mkdocs_nav.py --check && mkdocs build --strict` | no broken link, no page outside the nav |
| G7 | Disclosure | the pull request template's "AI assistance" section has exactly one box ticked | the Guard's author gate fails every pull request without it |

G4 fails the suite on a `physics` or `crosscheck` test without `KNOWN_FROM`; `python3 tests/claims.py
build/default --missing` lists any. A new test sets `KNOWN_FROM`, and `REASON` when its tolerance is looser than `rel:1e-12`, 4 ulp or
an absolute 1e-9 (`docs/how-to/develop/add-a-golden-test.md`).

## Before touching these, stop and ask the maintainer

- Anything under `src/`, `FQHE/`, `FTI/`, `Spin/`, `QuantumDots/`, `Base/` other than `CMakeLists.txt`,
  `*.cmake` and `Makefile.am`. These are DiagHam's physics sources; a change there is "for maintainer review"
  (`GOVERNANCE.md`) and needs a DiagHam author's agreement.
- An expected value, tolerance or `CHECK` line of any test; any file under `tests/data/` or `benchmarks/`;
  `tests/manifest*.txt`; the default of any `option()`; anything under `patches/`.
- `.github/approved-developers.txt`, `.github/guard-reviewers.txt`, `.github/volunteers.txt`, `CODEOWNERS`,
  the rulesets, the workflows.

Asking means: describe the change and why, in the pull request or hand-off note, and leave the files unchanged.

## Never

- Edit a golden, a reference or a tolerance to make a red test green.
- Remove, rename or skip a test to make a run green; put a name back if the diff shows it gone.
- Report a count that was not measured in this run, or carry a number from an earlier session.
- Claim a result is "validated", "proven" or "correct"; say which test passed and what it checks.
- Call anything a "bug" or "defect"; it is a finding, registered in `docs/reference/open-findings.md`.
- Merge, approve, release a hold, or add a label that releases one. No AI is on any approver list.
- Add a dependency, a download at build time, or a vendored copy without a line in `CHANGELOG.md` and the reason.

## Vocabulary

Test labels: `physics` (answer known independently of DiagHam), `crosscheck` (two programs or codes agree),
`regression` (a saved spectrum is reproduced; change detection only), `smoke` (`--help` exits 0), `selftest`
(the suite checks itself), `open-finding` (a `WILL_FAIL` reproducer that turns red when the finding is fixed).
A new test is one of these, and says so. Outputs of the migration are described as candidate, in review or
adopted; nothing is "published" until the authors adopt it.

## How to hand work back

One pull request per change, from a branch, never a push to `main`. The description has: what changed and why;
the G1 to G7 block; "Damage, explained" if G3 exited 2; "for maintainer review" in the title if a physics source
or a locked artefact had to change. The Guard (`docs/how-to/develop/guard.md`) then holds the pull request, asks
the reviewers and the approved developers to look, and posts it to the volunteers' channel. Only a listed
reviewer releases it.

## ALPS

The same file applies to the ALPS repository once its team adopts it. Only the command table changes: ALPS's
configure, build and test commands, its manifest or equivalent, and its generated pages. Those commands are not
read here and are left to the ALPS maintainers to fill in; the gates, the stop-and-ask list, the never list and
the vocabulary are the same. The shared Guard pieces (author gate, disclosure section, damage scan, Discord
post) are in this repository's `.github/` and `scripts_cmake/` and run on any repository with a ctest suite and a
manifest file.
