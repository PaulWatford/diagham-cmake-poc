<!-- Thank you. The checklist is what the maintainer will look for; ticking it honestly saves a round trip. -->

## What this changes and why

## Checklist

- [ ] `ctest --preset default` is green locally (say which other presets you ran).
- [ ] No generated file, build log or output file is committed (`build/`, `*.log`, `*.vec`).
- [ ] If a DiagHam source file changed: one fix per commit, root cause and verification in the message, trailers `Upstream-Patch:` and `Upstream-Base:`; a row in `docs/explanation/upstream-fixes.md`.
- [ ] If a computed number changes (physics-affecting): the commit and this PR say **for maintainer review**, and a DiagHam author has been asked (see `GOVERNANCE.md`).
- [ ] If a test was added: its expected value is known independently of DiagHam and the comment or the expected file says why; `tests/manifest.txt` and `docs/reference/test-coverage.md` regenerated (`python3 tests/coverage.py build/default --manifest tests/manifest.txt --write docs/reference/test-coverage.md`).
- [ ] If a regression reference was regenerated: the change is intended and understood, and the commit says so.
- [ ] Documentation follows the change in this PR (every page has `Purpose:` and `Source:` lines).
