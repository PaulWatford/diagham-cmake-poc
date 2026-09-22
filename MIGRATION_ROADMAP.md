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

> **Correction (22/09, audited by a second Claude session):** the
> paragraph originally here claimed `guysoft/DiagHam` "has exactly one
> commit" and used that to argue this PoC's reference material was not
> representative of production history needs. **That claim was false**,
> traced to this engagement's repeated use of `git clone --depth 1`
> (shallow clone) against the mirror, which truncates history to the tip
> commit regardless of what the mirror actually has. A full,
> non-shallow clone shows **4,098 real commits spanning 2003-05-03 to
> 2020-10-23** (top authors: regnault 2317, moller 714, repellin 375,
> sterdyni 306, papic 140), i.e. `guysoft/DiagHam` already carries real,
> rich SVN-derived history — it was never a single squashed import; only
> this PoC's *view* of it was truncated, by tooling choice, not by what
> the mirror contains. `DEFERRED.md` has the same correction and a
> worked example (`QHEFermionsTorusWithSpin.cc`'s real last-touch date,
> recovered with `git log --follow` against the full clone). This
> changes the framing below: the production conversion should still use
> a history-preserving SVN→Git conversion for the canonical authority
> (`nick-ux.org/diagham/websvn`, which has finer revision detail than
> even the full git mirror), but it is not correcting a single-commit
> mirror — it is carrying forward or re-deriving from a mirror that
> already has most of the real history, just not all of it, and not
> everything DEFERRED.md said about individual files' last-touch dates
> was actually checked against that fuller picture in this pass.

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
