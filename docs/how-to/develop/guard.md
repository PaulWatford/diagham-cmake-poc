# Run, read and release the Guard (author gate and damage scan)

Purpose: how the Guard workflow treats a pull request, how to read its two comments, how a maintainer releases a hold, and how to run the scan by hand.
Source: `.github/workflows/guard.yml`, `.github/scripts/author_gate.py`, `.github/scripts/guard_report.py`, `scripts_cmake/damage_scan.py`; the rules are `AGENTS.md` and `GOVERNANCE.md`, the design `explanation/guard-design.md`. Written 2 November 2026.

## What runs on every pull request

Two jobs, both meant to be required status checks on `main`:

| Check | What it decides | Red means |
|---|---|---|
| `guard / author gate` | who must look at the change | the pull request is held (`held-for-approval`), or the "AI assistance" section of its text has no box ticked |
| `guard / damage scan` | whether a locked artefact changed | the scan flagged something and the hold has not been released |

Neither job builds or tests anything; the CI workflow does that. Neither job merges, approves, edits a file or rewrites a golden. They label, they request reviews, and each keeps one comment on the pull request up to date.

The workflow runs on `pull_request_target`, so it can label and comment on pull requests from forks. It therefore checks out the base branch only. The pull request's commits are fetched as git data and read by `scripts_cmake/damage_scan.py` from the base checkout; nothing from the pull request is executed. Keep it that way when editing the workflow.

## The author gate

Read from the base branch:

- `.github/approved-developers.txt`: logins whose pull requests are not held for being external.
- `.github/guard-reviewers.txt`: the named maintainers who are asked to review a held pull request and who alone can release it.

A pull request is **held** when its author is not an approved developer (or is a bot), or when AI assistance is declared in the pull request text or detected in the commits (AI co-author trailers, bot committers, "generated with" markers). Held pull requests get the labels `held-for-approval` and, as applicable, `external-contribution` and `ai-assisted`, a review request to every guard reviewer, and a comment with four fields: Author, AI assistance, DiagHam physics sources touched, Hold.

"DiagHam physics sources touched" is checked for AI-assisted changes only: a file under `src/`, `FQHE/`, `FTI/`, `Spin/`, `QuantumDots/` or `Base/` that is not a `CMakeLists.txt`, `.cmake` or `Makefile.am` is reported **for maintainer review**, the wording `GOVERNANCE.md` uses for physics-affecting changes, because `AGENTS.md` puts DiagHam's physics code out of bounds for automated changes.

The "AI assistance" section of the pull request template must have exactly one box ticked. With none ticked the gate fails for every author, approved or not, and says so; tick a box (editing the pull request text re-runs the gate). A pull request by an approved developer with "No AI tools" ticked passes with no label and no comment.

## The damage scan

`scripts_cmake/damage_scan.py BASE HEAD` compares the merge base with the pull request head and reports the classes `AGENTS.md` lists under "Damage scan":

| Class | Flagged when |
|---|---|
| tests removed | a name in `tests/manifest.txt` or `tests/manifest-core-programs.txt` is gone, even if the count is unchanged |
| demo unregistered | a `physics.hubbard.2x2.*` test (the analytic Hubbard golden and its companions) is gone |
| checks skipped | `enable_testing()` or `add_test()` calls removed; a `tests/*.cmake` file deleted; a locked name added to `tests/manifest-optional.txt`; `DIAGHAM_BUILD_TESTS` no longer defaults ON |
| golden edited | a file under `tests/data/` or a reference spectrum in `benchmarks/` modified, renamed or deleted; a line carrying a value or tolerance changed in the test registrations |
| option default flipped | an `option()` that defaulted ON now defaults OFF, or a `DIAGHAM_BUILD_*` option was removed |
| in-source build | the `CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR` guard is gone from `CMakeLists.txt` |
| patches changed | anything under `patches/` |

Additions are information, not damage: new test names, new options, options switched on, new `WILL_FAIL` registrations. The report ends with the test-name add/remove list and a diffstat, the attachments `AGENTS.md` asks for.

A flagged pull request gets the labels `damage` and `held-for-approval` and the check stays red. The scan is deliberately literal: it says what changed, not whether the change is right. A regenerated regression reference, an intended golden correction or a new carried fix will all be flagged, and the maintainer releases them after reading the diff. It never edits anything to make a run green.

Run it by hand from any checkout:

    python3 scripts_cmake/damage_scan.py origin/main HEAD
    python3 scripts_cmake/damage_scan.py v2026.9.0 main --markdown scan.md --json scan.json

Exit 0 is clean, 2 is damage, 1 means it could not run. `scripts_cmake/damage_scan_selftest.sh` clones the repository into a temporary directory, commits one synthetic change per class and checks that each is caught; it also runs the gate scripts on the fixtures in `.github/scripts/fixtures` without talking to GitHub.

## Releasing a hold

A login in `.github/guard-reviewers.txt` releases a hold in either way:

1. add the label `approved-by-maintainer` to the pull request (the label event re-runs the workflow, which turns both checks green and removes `held-for-approval`); or
2. approve the pull request through a review of its current head commit, then re-run the workflow (a new push also re-runs it, and a push after the approval needs a new approval).

The release label added by anyone not on the list is removed again and named in the comment. Silence, a comment, or a thumbs-up is not a release. The second approval mechanism, GitHub's own, is also what the branch ruleset requires: one approving review from a code owner.

## Rulesets

`.github/rulesets/main.json` and `.github/rulesets/tags.json` are the repository rulesets in GitHub's import format (Settings, Rules, Rulesets, New ruleset, Import a ruleset). The first requires pull requests with one code-owner approval and the eight status checks (five CI presets, the container job, the two Guard jobs) on `main`, blocks force pushes and deletion, and lets a repository administrator bypass only through a pull request. The second makes `v*` and `migration-r*` tags immutable. They are files so that the settings can be reviewed and restored; GitHub does not read them from the repository.

## When the repository moves to GitLab

The workflow and the rulesets are GitHub-specific. The scan itself is not: `.gitlab-ci.yml` can run `scripts_cmake/damage_scan.py $CI_MERGE_REQUEST_DIFF_BASE_SHA $CI_COMMIT_SHA` as a merge-request job, and a GitLab approval rule replaces the author gate's lists. See `move-to-gitlab.md`.
