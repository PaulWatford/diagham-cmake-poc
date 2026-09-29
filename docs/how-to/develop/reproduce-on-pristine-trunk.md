# Reproduce the registered defects on a pristine trunk

Purpose: show that a defect in `docs/reference/known-defects.md` is in DiagHam's Subversion trunk as DiagHam's own build produces it, not in this repository's conversion or its CMake build.
Source: `scripts_cmake/reproduce_on_pristine_trunk.sh`; first run 2026-09-29 against trunk r4493.

## What the script does

```
scripts_cmake/reproduce_on_pristine_trunk.sh [BUILD_DIR] [REV] [WORK]
```

1. `svn info` on `https://www.nick-ux.org/diagham/svn/DiagHam/trunk`: is the revision this repository carries (`scripts_cmake/upstream-revision.txt`) still the head?
2. `svn export -r REV` of trunk into `WORK` (default `~/diagham/pristine-r<REV>`), straight from the server.
3. A whole-tree diff of that export against the `upstream` branch. Expected differences: empty directories (git cannot store them) and Subversion keyword expansions such as `$Id$`; no source file.
4. DiagHam's own build: `bootstrap.sh`, `configure --enable-fqhe --enable-fti` (no GSL, the default of both build systems), `make -k`. Upstream's own build failures (the ungated LAPACK calls, U01 to U03) stop a plain `make`, so the script keeps going past them.
5. Every reproducer, run with those binaries, next to a positive control (a sphere Coulomb spectrum against `tests/oracles/sphere_ed.py`, which must pass) and next to what the register says. Only `tests/check_spectrum` and the reference files of this repository are involved.

It needs `svn`, autotools and a C++ compiler, and about half an hour for the build; a second run reuses the export and the build.

## The run of 2026-09-29 against r4493

- Trunk head: r4493, last changed 19 September 2026.
- Export versus `upstream`: four empty directories and one `$Id$` expansion in `mkinstalldirs`; no source file differs.
- The control passes to 10⁻⁹.
- U34: `FQHEDiskBosonsTwoBodyGeneric` gives 0, 1.5951, 1.9825, 2.0000, 2.1174, 2.8833, 5.9842 for three bosons at L_z = 6 with V₀ (exact: 0, 3/2, 15/8, 33/16, 9/4, 3, 6), the same digits as the CMake build; two particles give the single level 2.
- U29: N = 3 hangs after `start`; N = 2 at L_z = 5 segfaults.
- U32: both symmetrisation flags give dimension 0; one flag at a time gives 25 and 28.
- U30: −17.8677 and −18.4122 in two momentum sectors, below the exact −16 (the sources on `main` carry the fix; the pristine build does not).
- U33: segmentation fault after "Hilbert space dimension = 12" (same: fixed on `main`, present in the pristine build).

The registers (`docs/reference/known-defects.md`, `docs/explanation/verification.md`) say "present at r4493" for these; this is the evidence behind that column. Rerun the script when trunk moves.
