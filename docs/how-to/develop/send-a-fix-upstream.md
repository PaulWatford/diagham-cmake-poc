# Offer a fix or a finding to the DiagHam authors

Purpose: turn a fix carried on `main`, or a defect found by the tests, into something the DiagHam authors can act on in their Subversion repository.
Source: new (2026-09-28); `scripts_cmake/upstream_pack.py`, `docs/upstream-reports/`, `patches/upstream/`, `GOVERNANCE.md`.

The physics code is theirs (`GOVERNANCE.md`): a fix is only finished when it is in Subversion, and until then this repository carries it as one commit with the trailers `Upstream-Patch:` / `Upstream-Base:`. Nothing in this repository sends anything; the maintainer does, by e-mail or by whatever route the authors prefer.

## What to send for a fix

`python3 scripts_cmake/upstream_pack.py` regenerates, from those commits:

- `patches/upstream/NN-name.diff` — the DiagHam source files only, no path prefixes: `svn patch NN-name.diff` in a trunk working copy at the base revision applies it;
- `docs/upstream-reports/NN-name.md` — the report: what, why, how it was verified (the commit message), files, how to apply.

Send the report and the diff together. If the fix changes a computed number (marked *PHYSICS FIX, FOR MAINTAINER REVIEW* in its subject), say so first and include the cross-check that supports it (`docs/reference/tests.md` names the tests).

## What to send for a defect without a fix

`docs/upstream-reports/open-*.md` hold the reproducers of registered defects that are not fixed here (U29 disk hang, U31 two-band checkerboard segfault, U34 wrong disk pseudopotential energies): the exact command, what happens, what was expected and why, the build it was seen in. The `known-bug` tests in ctest carry the same reproducers and turn red when upstream fixes them.

## When upstream has taken a fix

The next `scripts_cmake/sync_upstream.sh` run will either merge cleanly (upstream applied the same change) or conflict (upstream fixed it differently). In both cases: drop our commit's change if it is now redundant, update the row in `docs/explanation/upstream-fixes.md` to "fixed upstream at r…", and let the known-bug test (if any) be promoted to an ordinary test.

## Writing a new report by hand

Every report has the same four parts — what happens, how to reproduce it (one command that runs in seconds), what was expected and why it is known, what was changed and how that was verified — because those are the four things a maintainer needs to decide in one reading.
