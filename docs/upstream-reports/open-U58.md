# Open finding U58: the Josephson-pairing program segfaults with a fixed particle number and off-diagonal terms

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U58), found 2026-10-26 while writing the program's first independent tests; r4493 sources, GCC 15. The same crash on a pristine autotools build of r4493.

## What happens

`FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing --nbr-particles N` (or `--all-fixednbrparticles`) with a confining file that has an off-diagonal line (m != n) takes the branch without momentum conservation (`FQHE/src/Programs/FQHEOnCylinder/FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing.cc`, from line 562). There the Hilbert space for a fixed particle number is commented out (line 581, `// Space = new QuasiholeOnSphereWithSpinAllMomenta (...)`), `Space` stays null, and line 588 calls `Space->GetLargeHilbertSpaceDimension()`. An address-sanitizer build reports the null read at line 588; a normal build ends with a segmentation fault (exit status 139). With only diagonal lines, a fixed particle number works, and without `--nbr-particles` off-diagonal lines work.

Example: 2S = 6, perimeter 6, quasihole data from `FQHECylinderQuasiholeMatrixElements -l 6 --cylinder-perimeter 6 --offdiagonal-density 2`, a confining file with the diagonal terms and the lines (0, 1) and (1, 3), `--nbr-particles 4`.

## Fix direction

Restore the fixed-N space, or stop with a message that `--nbr-particles` is not available with off-diagonal terms.

## Consequence for the tests

`openfinding.fqhe.cylinder.quasiholes_time_reversal_josephson.fixed_n_offdiagonal_crash.offdiagonal_fixed_n4_2s6_perimeter6_sz0` (WILL_FAIL) expects the spectrum of that case (`tests/oracles/cylinder_quasiholes.py`) and fails while the program crashes.
