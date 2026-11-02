<!-- Thank you. The checklist is what the maintainer will look for; ticking it plainly saves a round trip. -->

## What this changes and why

## AI assistance

<!-- Tick exactly one. The Guard workflow reads this section; a pull request with no box ticked fails its
     author gate until one is. AI-assisted and AI-generated changes are held for a named maintainer
     (AGENTS.md); changes to DiagHam's physics sources made with AI tools are marked "for maintainer review". -->

- [ ] No AI tools were used for this change.
- [ ] AI-assisted: a tool helped (name it and say what it did, e.g. "Claude Code drafted the CMake changes; I reviewed and tested them"):
- [ ] AI-generated: the change was produced by a tool and is submitted for maintainer review (name the tool and the prompt or task):

## Checklist

- [ ] `ctest --preset default` is green locally (say which other presets you ran).
- [ ] No generated file, build log or output file is committed (`build/`, `*.log`, `*.vec`).
- [ ] If a DiagHam source file changed: one fix per commit, root cause and verification in the message, trailers `Upstream-Patch:` and `Upstream-Base:`; a row in `docs/explanation/upstream-fixes.md`.
- [ ] If a computed number changes (physics-affecting): the commit and this PR say **for maintainer review**, and a DiagHam author has been asked (see `GOVERNANCE.md`).
- [ ] If a test was added: its expected value is known independently of DiagHam and the comment or the expected file says why; `tests/manifest.txt` and `docs/reference/test-coverage.md` regenerated (`python3 tests/coverage.py build/default --manifest tests/manifest.txt --write docs/reference/test-coverage.md`).
- [ ] If a regression reference was regenerated: the change is intended and understood, and the commit says so.
- [ ] Documentation follows the change in this PR (every page has `Purpose:` and `Source:` lines).
