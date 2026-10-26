# Open finding U57: the Josephson-pairing program reads a confining-file line with m > n as the wrong coupling

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U57), found 2026-10-26 while writing the program's first independent tests; r4493 sources, GCC 15. The same spectra on a pristine autotools build of r4493.

## What happens

The wiki page of `FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing` describes the confining file as encoding "any term such as" sum_{m,n,sigma} V_{m,n,sigma} c+_{m,sigma} c_{n,sigma}, one line (m, n, V_up, V_down) per term. The program (`FQHE/src/Programs/FQHEOnCylinder/FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing.cc`, lines 277 to 296) stores an off-diagonal line in `OffDiagonalOneBodyPotential...[CreationIndices[i]][|m - n| - 1]` and the Hamiltonian adds the Hermitian conjugate of each stored term. That slot is read as the coupling between orbital `CreationIndices[i]` and the orbital |m - n| above it. For m < n this is the pair (m, n), as intended. For m > n it is the pair (m, 2m - n): the line (1, 0) couples orbitals 1 and 2.

So a Hermitian V written out in full, as the formula invites (both (0, 1) and (1, 0)), gets a spurious coupling, and a line with m > n alone never means what it says. Measured against an independent construction (`tests/oracles/cylinder_quasiholes.py`; 2S = 6, perimeter 6, the two-layer problem with pairing and charging energy):
- a single line (0, 1, 0.2, 0) gives the spectrum of 0.2 (c+_0 c_1 + h.c.) to 5e-14;
- the lines (0, 1, 0.2, 0) and (1, 0, 0.2, 0) together give the spectrum of 0.2 (c+_0 c_1 + h.c.) + 0.2 (c+_1 c_2 + h.c.) to 5e-14, where the formula means a single 0.2 (c+_0 c_1 + h.c.); the spectra differ by up to 0.13.

## Fix direction

Store a line with m > n at `[AnnihilationIndices[i]][m - n - 1]` with the conjugate amplitude (it is the conjugate partner of (n, m)), and either skip a pair already given or say on the wiki page that only m < n lines should be written.

## Consequence for the tests

`openfinding.fqhe.cylinder.quasiholes_time_reversal_josephson.lower_triangle_line.both_triangles_2s6_perimeter6_sz0` (WILL_FAIL) gives both triangles of a Hermitian V. The physics tests write only m <= n lines.
