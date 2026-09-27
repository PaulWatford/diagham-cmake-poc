# Getting help

- **Something does not build or a test fails**: open an issue with the *Bug report* template; it asks for the preset, the compiler, the exact command and the output. The build guides are in [`docs/how-to/build/`](docs/how-to/build/README.md).
- **A program gives a result you do not expect**: first check [`docs/reference/known-defects.md`](docs/reference/known-defects.md) and [`docs/reference/test-coverage.md`](docs/reference/test-coverage.md) (has this program's physics been checked at all?). Then open a *Bug report* with the program, its options and what you expected and why.
- **How do I …?**: [`docs/README.md`](docs/README.md) is organised by task; the program reference (`docs/reference/programs/`) has the `--help` of every program and the manuals adapted from the DiagHam wiki. Questions that are not defects go to GitHub Discussions.
- **DiagHam physics questions** (what a program computes, which state a Hamiltonian targets): the DiagHam authors' wiki and publications list are the primary sources (`docs/reference/publications.md`); this repository's maintainer can say what the build and tests do, not always what the physics means.

Response times: this is a research project maintained alongside other work; issues are read, not always answered the same week. Pull requests with a test are the fastest way to get something changed (see `CONTRIBUTING.md`).
