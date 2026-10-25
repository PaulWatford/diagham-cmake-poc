# Move the repository to GitLab

Purpose: the checklist for moving this repository from GitHub to its GitLab home with its history, tags and CI intact, and for checking afterwards that nothing was lost.
Source: written 2026-10-25. The team confirmed GitLab as the public host on 20 October 2026, partly for ethical reasons; the move waits for the approval of the GitLab structure (instance, group, project path, runners). The GitLab pipeline is `.gitlab-ci.yml`, checked against GitLab's CI schema and with `gitlab-ci-local` (see "CI on GitLab" below). Numbers are those of GitHub `main` at `2a40cb0`.

## What there is to move

| What | On GitHub today | Notes |
|---|---|---|
| `main` | 4,575 commits | the CMake build, tests and documentation on top of DiagHam r4493 |
| `upstream` | 4,477 commits | the faithful conversion of the DiagHam Subversion trunk, r1 to r4493 |
| Tags | `migration-r4493`, `v2026.9.0` | annotated tags |
| Pull request #1 | closed, not merged | the proof-of-concept audit; keep for the record |
| Issues, discussions | GitHub only | GitLab's GitHub importer brings issues and pull requests across |
| `svn-branch-testing_nr`, `svn-branch-testing_zp`, `svn-tag-first-release` | **not** on GitHub | old Subversion branches (last commits 2009, 2009, 2003) in the conversion repository; decide whether GitLab carries them |

Commit hashes do not change in a move. The tags, the release notes and any citation of a commit stay valid.

## Decisions to make first

1. **Where:** the GitLab instance, the group and the project path. The addresses in step 6 depend on them.
2. **Runners:** Docker executor; the images `ubuntu:24.04` and `python:3.12-slim`; at least 4 cores and a 2-hour job timeout for the preset jobs. A runner that allows Docker-in-Docker is needed only for the container job.
3. **The old Subversion branches:** carry them or not (table above).
4. **GitHub afterwards:** archive it read-only with a pointer to GitLab (recommended), or keep it as a push mirror of GitLab.
5. **Zenodo releases:** who maintains them after the move. GitLab has no Zenodo integration, so archiving a release becomes a manual upload (`RELEASING.md`).

## The move

1. **Freeze `main` on GitHub.** Announce it, and note the `main` and `upstream` hashes.
2. **Create an empty project on GitLab.** No README and no licence file, so that the first push does not meet an existing history.
3. **Copy branches and tags exactly.** Do not use `git push --mirror`: it would also push GitHub's pull-request refs (`refs/pull/*`), which GitLab refuses.

       git clone --mirror https://github.com/PaulWatford/diagham-cmake-poc.git
       cd diagham-cmake-poc.git
       git push <gitlab-url> 'refs/heads/*:refs/heads/*' 'refs/tags/*:refs/tags/*'

   To carry the old Subversion branches, push them from the conversion repository the same way. To bring issues and pull requests across as well, use GitLab's GitHub importer instead of these commands, then still do step 4.
4. **Check that nothing changed.** The two lists must be identical, and the counts must match the table above:

       git ls-remote https://github.com/PaulWatford/diagham-cmake-poc.git 'refs/heads/*' 'refs/tags/*' > github.refs
       git ls-remote <gitlab-url> 'refs/heads/*' 'refs/tags/*' > gitlab.refs
       diff github.refs gitlab.refs && echo identical
       git rev-list --count main upstream     # in a fresh clone from GitLab: 4575, then 4477 on its own

5. **Set the project up as GitHub is set up today** (`GOVERNANCE.md`, "Repository settings the maintainer keeps"):
   - protect `main`: merge through merge requests only, nobody pushes directly, force push off;
   - merge requests: "Pipelines must succeed" on; GitLab has no per-job required checks like GitHub, and this setting covers them all;
   - protect the tags `v*` and `migration-r*`;
   - issues on, wiki off (the knowledge base is `docs/`), Pages on.
6. **Update the addresses in one merge request.** These files name the GitHub repository or its Pages site:
   - `mkdocs.yml` (`site_url`, `repo_url`, `edit_uri`), `README.md`, `docs/README.md`, `docs/how-to/develop/docs-site.md`, `docs/explanation/checklists.md`;
   - `CITATION.cff` (`repository-code`, twice), `codemeta.json` (`codeRepository`, `url`, `issueTracker`, `readme`), `.zenodo.json` (the documentation link);
   - `container/Dockerfile` and `container/diagham.def` (the source labels);
   - `.github/ISSUE_TEMPLATE/`: GitLab reads issue templates from `.gitlab/issue_templates/`, so move them there;
   - `GOVERNANCE.md` ("Repository settings the maintainer keeps"), `docs/explanation/migration-history.md` (the host section), and the comments at the top of `.gitlab-ci.yml` and `scripts_cmake/ci.sh`.

   Leave `docs/history/release-notes/2026.9.0.md` as it is: it records the release as published.
7. **Start CI.** The first pipeline on `main` must be green: the five preset jobs (`default`, `full`, `hpc`, `core-programs` with GCC, `core` with Clang) and `pages`. Then:
   - create a weekly pipeline schedule (GitHub's runs on Mondays at 04:17 UTC) for `upstream-drift`, which fails as a warning when the DiagHam Subversion trunk has moved past the revision this repository carries;
   - set the project variable `DIAGHAM_CONTAINER_JOB` to `1` if a Docker-in-Docker runner exists, which switches on the container job.
8. **Close GitHub** as decided: archive it with a pointer to GitLab in its README, or set it up as a push mirror from GitLab.
9. **Announce the new address** to the DiagHam authors and the group, and update any links in the wiki.

## CI on GitLab

`.gitlab-ci.yml` has the same jobs as the GitHub workflows, and both call `scripts_cmake/ci.sh`:

| GitHub | GitLab | When |
|---|---|---|
| `ci.yml` build-and-test (5 presets) | `build-and-test` (5-job matrix, the same compilers) | merge requests, `main`, tags, manual runs |
| `ci.yml` container | `container` | when `DIAGHAM_CONTAINER_JOB` is `1` |
| `ci.yml` upstream-drift (opens an issue) | `upstream-drift` (fails as a warning instead) | the weekly schedule, manual runs |
| `pages.yml` build (pull requests) | `docs-check` | merge requests that touch the docs |
| `pages.yml` deploy | `pages` | `main` |

Checked on 25 October 2026: the file validates against GitLab's CI schema. `gitlab-ci-local` 4.56.2 expands each pipeline to the jobs in the table: a push to `main`, a merge request with and without docs changes, a schedule, a manual run, and the container switch. The `upstream-drift` script ran against the live Subversion trunk (r4493, no drift). The jobs themselves have not run on a GitLab runner yet; the first pipeline after step 7 is that test.
