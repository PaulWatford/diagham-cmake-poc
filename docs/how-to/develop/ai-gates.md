# Run the gates an AI tool must pass before handing work back

Purpose: what each of the seven gates in `CLAUDE.md` checks, what a failure means, how to fix it without touching a locked artefact, and how to report the run in a pull request.
Source: `CLAUDE.md` (the gates and the stop-and-ask list), `AGENTS.md` (the rules they enforce), `GOVERNANCE.md`, `docs/how-to/develop/guard.md` (what happens to the pull request afterwards). Written 2 November 2026.

## The rule

Nothing an AI tool writes here is done until every gate has been run and has passed, and the run is shown: the command and its last lines pasted into the pull request or the hand-off note. A gate that was not run is reported as NOT RUN, with the reason. A number that is not in a command's output or in a file of this repository is reported as "not read".

## The gates

| # | Gate | Command |
|---|---|---|
| G1 | Configure and build, out of source | `cmake --preset default && cmake --build --preset default` |
| G2 | Full suite with the manifest | `ctest --preset default --output-on-failure` |
| G3 | Damage scan | `python3 scripts_cmake/damage_scan.py main HEAD` |
| G4 | Claims are sealed | `ctest --preset default -R selftest.claims` |
| G5 | Generated pages are current | `python3 tests/coverage.py build/default --write docs/reference/test-coverage.md && python3 tests/claims.py build/default --write docs/reference/claims.md && python3 tests/numbers_of_record.py build/default --check` |
| G6 | Documentation builds | `python3 scripts_cmake/gen_mkdocs_nav.py --check && mkdocs build --strict` |
| G7 | Disclosure | exactly one box ticked under "AI assistance" in the pull request text |

**G1, the build.** A failure is a compile or configure error. Fix the build files or the code you changed; an in-source build is refused by `CMakeLists.txt` on purpose, so configure into `build/`.

**G2, the suite.** Every test must pass and `selftest.manifest` must find exactly the names of `tests/manifest.txt`. A red physics test means the change altered a computed number; find out why, and do not touch its expected value or tolerance. A red `selftest.manifest` after adding a test means the manifest must be regenerated in the same change (`python3 tests/coverage.py build/default --manifest tests/manifest.txt`); after removing one it means a test is gone, which is damage (G3).

**G3, the damage scan.** It compares the branch with `main` and exits 2 when a locked artefact changed: a test name removed, a golden or reference file edited, a value or tolerance line changed, an option default flipped, the in-source guard removed, `patches/` touched. Exit 2 is acceptable only when every finding is intended and listed under "Damage, explained" in the pull request; the maintainer then decides. Never undo a finding by editing the scan.

**G4, sealed claims.** `selftest.claims` checks that every `physics` test says where its answer is known from, carries a tolerance (or the reason it has a looser one), uses only the labels of the vocabulary, and that open-finding tests are `WILL_FAIL` and registered. Until it exists the gate is reported NOT RUN (`docs/explanation/guard-design.md`, status table).

**G5, generated pages.** The coverage page, the claim index and the numbers page are written from the build. A stale page is a failure: regenerate it with the command and commit the result. Do not edit a generated page by hand. Until `tests/claims.py` exists that part of the gate is reported NOT RUN.

**G6, documentation.** `gen_mkdocs_nav.py --check` fails when a page is missing from the navigation (run it without `--check` to rewrite `mkdocs.yml`); `mkdocs build --strict` fails on a broken link.

**G7, disclosure.** The Guard's author gate fails every pull request whose "AI assistance" section has no box ticked. Work done by an AI tool ticks "AI-assisted" or "AI-generated", which holds the pull request for a named maintainer.

## Stop and ask

These are never changed on an AI tool's own decision; describe the change and why in the pull request or hand-off note and leave the files as they are:

- DiagHam's physics sources: anything under `src/`, `FQHE/`, `FTI/`, `Spin/`, `QuantumDots/`, `Base/` other than `CMakeLists.txt`, `*.cmake` and `Makefile.am`;
- an expected value, tolerance or `CHECK` line of a test; any file under `tests/data/` or `benchmarks/`; `tests/manifest*.txt`; the default of an `option()`; anything under `patches/`;
- `.github/approved-developers.txt`, `.github/guard-reviewers.txt`, `.github/volunteers.txt`, `CODEOWNERS`, the rulesets, the workflows.

## Writing the gate block in a pull request

One line per gate, the command's last line or lines after it, in a code block:

```
G1 build:          [100%] Built target ...
G2 suite:          100% tests passed, 0 tests failed out of N
G3 damage scan:    Damage: none
G4 claims:         NOT RUN (selftest.claims not built yet)
G5 pages:          docs/reference/numbers.md is current; claims.py NOT RUN (not built yet)
G6 docs:           mkdocs.yml nav up to date; mkdocs exit 0
G7 disclosure:     AI-generated ticked
```

Paste the real output; the N above stands for the number your run printed. If G3 exited 2, add a section "Damage, explained" with one line per finding. Put "for maintainer review" in the title when a physics source or a locked artefact had to change.

## ALPS

`CLAUDE.md` is written so that the ALPS repository can adopt it by changing only the command table (its configure, build and test commands, its manifest or equivalent, its generated pages). Those commands are not read here; the gates, the stop-and-ask list and the vocabulary are shared.
