# Governance

Purpose: say who decides what in this repository, and how a change travels between it and DiagHam itself.
Source: new (2026-09-28), following the JOSS review criteria (stated support and governance expectations) and the FAIR4RS checklist.

## Two codebases, two owners

- **DiagHam** — the exact-diagonalisation toolkit itself — is developed by Nicolas Regnault (ENS Paris), Gunnar Möller (University of Kent) and their collaborators in a Subversion repository at nick-ux.org. They own the physics: Hilbert spaces, Hamiltonians, algorithms, programs. Nothing in this repository changes that.
- **This repository** owns the build system (CMake), the test suite, the documentation set and the release process around DiagHam. Its `upstream` branch is a faithful conversion of the Subversion history and is never edited by hand; `main` carries the migration and, as ordinary commits with an `Upstream-Patch:` trailer, a small number of fixes to DiagHam sources.

## Roles

| Role | Who | Decides |
|---|---|---|
| Maintainer | Paul Watford (`@PaulWatford`, `CODEOWNERS`) | everything in `main` that is not DiagHam physics: build, tests, docs, releases, merging |
| DiagHam authors | Nicolas Regnault, Gunnar Möller and the DiagHam committers (`AUTHORS`) | any change that alters a computed number; whether a fix carried here is adopted upstream |
| Contributors | anyone | propose changes by pull request, report defects, add tests and documentation |

A change is **physics-affecting** when it changes a computed number rather than whether the code compiles, links or runs. Such a change is marked *for maintainer review* in its commit and pull request and is not merged without a DiagHam author's agreement; until then it lives on a branch or behind a `known-bug` test. Everything else the maintainer merges on the strength of the test suite.

## How a fix travels

1. A defect in DiagHam's sources is registered in `docs/reference/known-defects.md` with a reproducer, and where possible a `known-bug` test (`WILL_FAIL`) that turns red the day it is fixed.
2. A fix is one commit on `main` with the root cause and the verification in its message and the trailers `Upstream-Patch: <name>` / `Upstream-Base: SVN r<revision>`; `docs/explanation/upstream-fixes.md` indexes it.
3. The fix is offered upstream as a bug report with the reproducer (`docs/how-to/develop/send-a-fix-upstream.md`, when written). When upstream adopts it, the next synchronisation drops the local copy.
4. New Subversion revisions arrive through `upstream` only, by `scripts_cmake/sync_upstream.sh` run by a person and reviewed as a pull request (`docs/how-to/develop/sync-upstream.md`; `scripts_cmake/upstream-revision.txt` records the revision carried; the weekly CI job opens an `upstream-drift` issue when trunk moves on); `main` merges `upstream`, never the other way round.

## Decision rules

- The test suite is the arbiter of "does not break anything": a change that turns a `physics`, `crosscheck` or `regression` test red is not merged until the reason is understood and written down. Goldens are never edited to make a test pass.
- Regression references (`tests/data/regression/`) are regenerated only for an understood, intended change, said so in the commit.
- Documentation follows the code in the same pull request; the coverage page is generated and CI refuses a stale one.
- Disagreements about physics are settled by the DiagHam authors; about the build and infrastructure by the maintainer; either can be asked in a GitHub Discussion or issue.

## Repository settings the maintainer keeps

These are GitHub settings, not files, and are listed here so they can be checked:

- `main` protected: pull requests only, CI required (`default`, `full`, `hpc`, `core`), no force pushes (the only force pushes so far were the history rewrites of 2026-09-27, before any fork existed, and are recorded in `CHANGELOG.md`).
- Issues and Discussions on; wiki off (the knowledge base is `docs/`, versioned with the code and published as a site).
- Releases carry a Zenodo DOI once the release process (`docs/how-to/develop/release.md`, when written) is in place.

## Changing this document

By pull request, like everything else; the maintainer merges it.
