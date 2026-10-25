# Make a release

Purpose: the checklist that turns the state of `main` into a numbered, citable release.
Source: new (2026-09-28); `scripts_cmake/check_version.py`, `CITATION.cff`, `codemeta.json`, `.zenodo.json`, `CHANGELOG.md`. The first release made with it is 2026.9.0. Zenodo archiving added after the team's decision of 20 October 2026; `RELEASING.md` at the root is the short version of this page.

## Versioning

Calendar versions, `YEAR.MONTH.PATCH` (e.g. `2026.9.0`, `2026.9.1`, `2027.1.0`), git tag `v2026.9.0`. The DiagHam revision a release carries is in `scripts_cmake/upstream-revision.txt` and in the release notes; the `migration-r<rev>` tags mark synchronisations with upstream and are separate from releases. The one version number lives in `project(DiagHam VERSION …)` in `CMakeLists.txt`; `scripts_cmake/check_version.py` (run by CI) fails when `CITATION.cff`, `codemeta.json` or the newest `CHANGELOG.md` heading disagree.

## Checklist

1. `main` is green on every CI job (default, full, hpc, core, container, docs) and `scripts_cmake/sync_upstream.sh --dry-run` says what upstream revisions, if any, are *not* in this release.
2. Regenerate what is generated and commit if changed: program reference (`gen_program_reference.py build/lapack docs/reference/programs`), coverage page and manifest (`tests/coverage.py`), site nav (`gen_mkdocs_nav.py`).
3. `CHANGELOG.md`: rename `[Unreleased]` to `[x.y.z] — YYYY-MM-DD`, add an empty `[Unreleased]` above it, and write the release summary at the top of the section (revision carried, test numbers, what changed for users).
4. Set the version in `CMakeLists.txt`, `CITATION.cff` (`version`, `date-released`), `codemeta.json` (`version`, `dateModified`); `python3 scripts_cmake/check_version.py` passes.
5. Write `docs/history/release-notes/x.y.z.md` (the text of the GitHub Release: what it is, numbers, what changed, open findings, how to cite); regenerate the nav.
6. Commit as "Release x.y.z", tag: `git tag -a vx.y.z -m "DiagHam (CMake migration) x.y.z, DiagHam r<rev>"`, push `main` and the tag.
7. GitHub → Releases → *Draft a new release* → choose the tag, paste the release notes, publish.
8. **Stable releases only** (a version 1.0, or a version described in a paper): archive it on Zenodo. On GitHub, switch the repository on in Zenodo's GitHub integration *before* publishing the release; publishing it then creates the record from `.zenodo.json`. On GitLab, upload the tagged source archive to Zenodo by hand with the same metadata. After the first archived release, add the concept DOI to `CITATION.cff` (`identifiers:`) and to the README, and commit.
9. Announce where the users are (the DiagHam authors, the group).

## What a release is not

A release is **not a new author list**: the creators of a Zenodo record are the DiagHam authors as listed in `CITATION.cff` (`check_version.py` checks `.zenodo.json` against it). Not every release is archived, only stable ones, and each paper links the version it describes. A release does not change `upstream`, and it is not a promise that every program is correct: `docs/reference/test-coverage.md` says which programs have their physics checked, and `docs/reference/open-findings.md` lists what is known to be wrong.
