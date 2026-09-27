# Scored against the field's checklists

Purpose: say, item by item, where this repository stands against the guidelines a research-software reviewer would apply, with the evidence — and where it does not.
Source: the JOSS review checklist, the FAIR4RS-based checklist (TU Delft / Netherlands eScience Center), The Turing Way's testing and project-design checklists, and CERTIFY-ED (arXiv 2605.11787), scored on 2026-09-28 at release 2026.9.0. Scores: **done**, **partial**, **not done**, **n/a**.

## JOSS review checklist

| Item | Score | Evidence |
|---|---|---|
| Source in an open repository | done | github.com/PaulWatford/diagham-cmake-poc |
| OSI-approved licence file | done | `LICENSE` (GPL-2.0-or-later), `NOTICE` |
| Contribution, issue and support guidelines | done | `CONTRIBUTING.md`, `.github/ISSUE_TEMPLATE/`, `SUPPORT.md` |
| Installation documented; dependencies listed | done | `docs/how-to/build/`, `docs/reference/build-system.md`, `codemeta.json` |
| Automated tests | done | 738 tests; what they prove: `docs/explanation/verification.md` |
| Functional claims confirmed | partial | 34 of 603 programs physics-checked; the claim made is exactly that |
| Usage examples | done | three tutorials, 132 manuals, generated reference of every program |
| Evidence of sustained development, community engagement | partial | one maintainer so far; DiagHam itself has 23 years of history on `upstream`; Discussions to be enabled |
| Statement of need, state of the field, software paper | not done | out of scope for this cycle; material exists (README, `why-cmake.md`, `publications.md`) |
| AI-usage disclosure | done | `AGENTS.md`, `CONTRIBUTING.md` (AI-assisted contribution rules) |

## FAIR research software (TU Delft / eScience Center checklist)

| Item | Score | Evidence |
|---|---|---|
| Version control, public repository, meaningful commits | done | 4,500 upstream commits + 53 migration commits, each with its rationale |
| Issues, contribution guidelines, code of conduct | done | templates, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md` |
| README, licence, citation file | done | `README.md`, `LICENSE`, `CITATION.cff` (points at the DiagHam authors), `codemeta.json` |
| Installation instructions, user and developer docs, tutorials | done | `docs/` (Diátaxis) |
| Documentation hosted | done | GitHub Pages site (`.github/workflows/pages.yml`) — live after the first push with Pages enabled |
| API reference from docstrings | not done | DiagHam's C++ classes are not documented with Doxygen here; the program reference covers the command line |
| Installation/execution verified | done | `install.find_package_consumer`, smoke test of every program |
| Integration and unit tests | partial | physics, cross-check and regression tests at program level; no C++ unit tests of classes |
| Continuous integration | done | four presets, container, docs, coverage/manifest/version checks |
| Code coverage tracked | partial | per-program coverage page (`test-coverage.md`, CI-enforced); no line-coverage tool |
| Dependencies recorded | done | `codemeta.json`, presets, `environment.sh` records what was built with |
| Linters / formatters | not done | deliberately: upstream sources are not reformatted (`AGENTS.md`); no linter on the CMake/Python layer yet |
| Versioning, tagged releases, changelog | partial | calendar versioning, `CHANGELOG.md`, release checklist and metadata ready; first tag after the push (`docs/how-to/develop/release.md`) |
| DOI for releases | not done, by decision | DiagHam is the DiagHam authors' work; a DOI on this repository would present it as the maintainer's. Theirs to mint, if ever |
| Registry upload (PyPI/conda) | n/a | C++ programs; containers instead (`container/`) |

## The Turing Way — testing checklist

| Item | Score | Evidence |
|---|---|---|
| Smoke tests | done | 599 |
| Unit tests | not done | no C++ unit tests; the oracles unit-test *themselves* (`selftest.*_oracle`) |
| Integration tests | done | program-level goldens are integration tests of Hilbert space + Hamiltonian + diagonaliser |
| System tests on common paths | partial | the most used programs (sphere/torus two-body, Hubbard, spin chains) are covered; most of the 603 are not |
| Regression tests | done | 11 (+1), labelled honestly |
| Runtime/self tests | done | harness self-tests, manifest, oracle checks |
| Tests documented and runnable in one command | done | `ctest --preset default`; `docs/reference/tests.md` |
| Run in CI | done | |
| Coverage monitored | partial | per program, not per line |
| Code review | partial | PR template and CODEOWNERS; one maintainer |

## CERTIFY-ED layers (verification of exact-diagonalisation codes)

| Layer | Score | Evidence |
|---|---|---|
| Analytic limits / closed forms | done | zero-mode counting, free fermions, AKLT, Haldane–Shastry, Coulomb pseudopotentials |
| Independent reference solver | done | `tests/oracles/sphere_ed.py`, `spin_ed.py`, `hubbard_ed.py` |
| Alternative algorithms agree | done | Lanczos vs full diagonalisation; LAPACK vs internal |
| Algebraic invariants / conservation laws | partial | dimension counting per sector; no trace/sum-rule tests |
| Arbitrary-precision reference | not done | the Wigner-symbol oracle is exact-rational internally, the spectra are double precision |
| Error injection detected | done | `selftest.checker.*`, `selftest.runner.*`, and manual injections in every prompt's audit |
| Certificates / hashed results | not done | |

## What this scoring is not

Not a claim of completeness: the "partial" and "not done" rows are the work list. The numbers on this page are those of release 2026.9.0 and are not updated automatically; the coverage page is.

## Audit record

2026-09-28, fresh clone at commit e18bdd61: version and navigation checks, the sync dry run, the pack check and the strict site build pass; every preset built from scratch and tested — core 45/45, lapack 739/739, full 740/740, default 738/738 — and the coverage page reproduced. The hpc preset (MPI, ScaLAPACK) is exercised by CI only; this machine has no MPI.
