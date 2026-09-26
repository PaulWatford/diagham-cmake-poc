# Contributing

This is the human workflow for the DiagHam CMake migration. The rules for
AI-assisted changes are in `AGENTS.md`; they apply on top of this.

## Where things go

- `upstream` is the pristine conversion of the DiagHam Subversion
  repository. **Never commit to it.** It is refreshed with
  `git svn fetch` and merged into `main`.
- `main` is where all work lands: build-system changes, tests, documentation,
  and accepted fixes to DiagHam sources.
- The canonical DiagHam source is still the Subversion repository at
  nick-ux.org. A fix to DiagHam code made here should also be reported
  upstream (see "Fixes to DiagHam code").

## Before you open a change

1. Build and test with the default preset:
   `cmake --preset default && cmake --build --preset default -j && ctest --preset default`.
   All tests must pass; `docs/reference/tests.md` explains what they check.
2. If you touched the build layer, also run the `lapack` preset (or `full`
   if you have the libraries).
3. If you added or removed a program or library, regenerate the
   per-directory CMake files (`python3 scripts_cmake/extract_autotools.py .`)
   and commit the result; do not edit generated `CMakeLists.txt` by hand.
4. Do not edit a golden value or a test tolerance to make a test pass. A
   changed golden is a physics change and needs a maintainer's review.

## Fixes to DiagHam code

- One fix per commit, with a message that states the root cause, what was
  compared, and how it was verified. Add the trailers
  `Upstream-Patch: <name>` and `Upstream-Base: SVN r<revision>`.
- Record the fix in `patches/PATCHES.md` (the audit trail).
- Physics-affecting changes (anything that changes a computed number rather
  than whether the code compiles) are marked **for maintainer review** in the
  commit and in the pull request, and are not merged without a DiagHam
  maintainer's sign-off.

## Documentation

- The knowledge base is organised as tutorials, how-to guides, reference and
  explanation under `docs/`; keep one purpose per page and cross-link.
- Every page adapted from the upstream wiki or an older document carries a
  provenance line (source, date, what was changed and why).
- Numbers (counts, revisions, results) are stated in one place and linked to,
  not re-typed.

## Pull requests

- Keep build-system changes, test changes and DiagHam-source fixes in
  separate commits so each can be reviewed and reverted on its own.
- CI must be green; a green CI is necessary, not sufficient — a named
  maintainer approves every merge.
