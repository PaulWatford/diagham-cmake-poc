# Migration roadmap: from this PoC to the production build

This document reconciles the CMake work in this repository with what was
actually agreed on the call with Gunnar Möller, so the PoC's scope claims
don't drift from the real plan. It covers five points from that
discussion that aren't build-system work per se, but constrain how this
PoC's output gets adopted.

## 1. Sequencing: Git migration before/alongside CMake, not after

The agreed order is to move DiagHam's history into Git first, and land
the CMake build on top of that Git repository, rather than finishing
CMake against autotools and only then migrating version control.

**Where this PoC currently stands relative to that:** everything here
was developed and verified against
[`guysoft/DiagHam`](https://github.com/guysoft/DiagHam), a **read-only
mirror** used purely as ground truth for current file contents. It is
not the production repository and should not be treated as one — see
point 3 below for why.

## 2. Host: University of Kent GitLab, not GitHub

The production repository's home is Kent's own GitLab instance, not
GitHub. This PoC repository (`PaulWatford/diagham-cmake-poc`, on GitHub)
is a personal scratch space for developing and demonstrating the CMake
port; it is explicitly not where the migrated project will live.

**Action for a later step:** once this CMake work is judged ready, it
needs to be re-hosted (or re-applied) against a Kent GitLab repository,
not merged into anything on GitHub. Nothing in this PoC assumes GitHub
specifically (no GitHub Actions, no GitHub-only tooling), so this is a
hosting/transfer step, not a rework.

## 3. Single repository, not split by module

DiagHam stays as one repository rather than being split into per-module
repos (e.g. separate repos for FQHE, FTI, Spin, QuantumDots). This
matches what the PoC already assumes: `extract_autotools.py` walks a
single source tree and there are no submodule or subtree boundaries
anywhere in `CMakeLists.txt` or `cmake/`. No change needed here — flagging
it so the "Beyond this iteration" section in the README doesn't imply a
multi-repo split when it lists Spin/QuantumDots as future scope; that's
future *build coverage*, not a future *repository*.

## 4. Full SVN history must be preserved

The production Git conversion needs to carry over DiagHam's complete SVN
revision history (24 years of commits), not a single squashed import.

**This is the one point where the PoC's reference material is explicitly
not representative of the production plan**, and it's worth stating
plainly since `DEFERRED.md` leans on this mirror's git log directly:
`guysoft/DiagHam` has exactly **one commit**
(`ed78a30`, 2020-10-23, an SVN-to-git import sweep), not real history.
`DEFERRED.md` already says this outright ("commit-level authorship
history for the deferred `.cc` files themselves is not available from
this snapshot") and points to the canonical SVN
(`nick-ux.org/diagham/websvn`) as the source that does have it — that
caveat was correct and should stay, it just needed to be tied explicitly
to this decision so it reads as "expected, and the production conversion
fixes it" rather than "a gap in this PoC."

**Action for a later step:** the real migration should use a
history-preserving SVN→Git conversion (`git svn`, `svn2git`, or
equivalent) against the canonical SVN, not a fresh import sweep like the
mirror used here.

## 5. Static linking is retained

The production build keeps static linking; this was not up for
reconsideration. The PoC already matches this: `diagham_add_library()` in
`cmake/DiagHamHelpers.cmake` calls `add_library(${target_name} STATIC
${ARG_SOURCES})` unconditionally, matching the autotools build's
`noinst_LIBRARIES` (`.a` archives). No change needed; recorded here as a
confirmed-consistent decision rather than an open question.

## Summary: what changes as a result of this reconciliation pass

| Decision | PoC status | Action |
|---|---|---|
| Git before/with CMake | CMake developed against a git mirror already | None — order is about the production repo, not this PoC |
| Kent GitLab host | PoC lives on GitHub (intentionally, as scratch space) | Re-host/transfer once this work is adopted |
| Single repo | Already assumed throughout | None |
| Full SVN history | Reference mirror has none (1 commit) — documented in `DEFERRED.md` | Production conversion must use a history-preserving SVN→Git tool |
| Static linking | Already the only mode `diagham_add_library` produces | None |

No CMake or patch content changes as a result of this pass — this was a
documentation/scope reconciliation, not a code change. The one
substantive addition is making explicit, in one place, which of this
PoC's reference material (the git mirror) is *not* representative of the
production plan and why, so a reviewer doesn't mistake the PoC's
single-commit mirror for the intended history-handling approach.
