# Migration history and the production plan

Purpose: how this repository relates to the canonical DiagHam, what was agreed for the production migration, and what has been done.
Source: `MIGRATION_ROADMAP.md` (2026-09-18, reconciling the proof of concept with the call with Gunnar Möller; corrected 2026-09-22), brought up to date 2026-09-26. Changed: the five agreed points are kept; each now states what has actually happened; the mirror procedure and the audited external contribution are added.

## The five points agreed with Gunnar Möller (University of Kent)

### 1. Git first, CMake on top — done

The agreed order was to move DiagHam's history into Git first and land the
CMake build on that, rather than finish CMake against autotools and migrate
version control afterwards.

**Where it stands.** The `upstream` branch is a `git svn` conversion of the
canonical Subversion repository at nick-ux.org: trunk, the branches
`testing_nr` and `testing_zp` and the tag `first-release`, 4,477 trunk
commits from 2003-05-03 (itself a CVS import) to 2026-09-18, SVN r4493.
Authors are mapped by name; the history stores no e-mail addresses, so none
are recorded. The converted tree at r4114 is byte-identical to the GitHub
mirror `guysoft/DiagHam` that the proof of concept was built against — which
is how it was established that the mirror stopped in October 2020 and that
the live SVN was 379 revisions (six years, ~200 fix commits) ahead. Every
CMake commit on `main` is now on top of r4493.

### 2. Host: University of Kent GitLab — pending

The production home is to be Kent's GitLab instance. This GitHub repository
is the development home in the meantime. Nothing here is GitHub-specific:
the CI pipeline is one script (`scripts_cmake/ci.sh`) that both the GitHub
workflow and the GitLab file call, so the move is a transfer, not a rework.
The Kent cluster recipe is a scaffold awaiting the site's values
(`docs/drafts/kent-cluster-build.md`).

### 3. Single repository — done

DiagHam stays one repository; no split by module. The generator walks one
tree and there are no submodule boundaries.

### 4. Full SVN history preserved — done

The conversion above carries the complete revision history. Earlier notes
here said the reference mirror had "one commit"; that was an artefact of
shallow clones and was corrected on 2026-09-22 — the mirror has 4,098 real
commits to r4114. Until the maintainers move to Git themselves, the SVN
remains the source of truth and `upstream` mirrors it (procedure below).

### 5. Static linking retained — done

`diagham_add_library` produces static archives only, as autotools'
`noinst_LIBRARIES` did.

## Branch layout

- `upstream` — the pure conversion. **Never edited by hand.**
- `main` — rooted in the original CMake proof-of-concept commits, then one
  merge bringing in the whole `upstream` history, then the build layer,
  the upstream fixes as individual commits, the generated files and the
  documentation. Synchronising with upstream is a merge into `main`, never a
  force-push.

Tags: `migration-r<revision>` at each synchronisation with the SVN.

## Keeping `upstream` current (the mirror procedure)

In the conversion repository (a `git svn` clone with the author map):

```
git svn fetch                      # appends new SVN revisions as commits
git branch -f upstream svn/trunk
git push origin upstream           # one-way: SVN -> Git
```

then in this repository `git merge upstream` into `main`, regenerate the
per-directory CMake files if any `Makefile.am` changed, build, run ctest,
tag. One direction only: nothing is ever pushed back into the SVN from Git.

## Where the proof of concept came from, and the audited contribution

The CMake build began as a personal proof of concept developed against the
GitHub mirror (see `CHANGELOG.md`, "Baseline"). On 2026-09-24 an automated
session opened a large draft pull request against it (Spin/QuantumDots, all
flags, ctest, install, CI) without instruction. It was audited in full; the
tests and the spinful-torus physics fix were verified correct, several
defects were found (dropped programs, a Windows-broken generator, install
collisions, unfixed regressions, no review gates), the PR was closed and its
verified parts were re-landed on r4493 as separate, reviewable commits.

## Still to do

- Kent GitLab hosting and the Kent cluster recipe.
- Maintainer review of the physics fix (`docs/explanation/torus-su2-coulomb-defect.md`)
  and of the Dice-lattice defaults; upstream submission of the fixes.
- The "best copy of each" review of duplicated and orphaned upstream code
  (`docs/explanation/deferred-code.md`).
- The review pipeline in CI (`docs/explanation/guard-design.md`).
