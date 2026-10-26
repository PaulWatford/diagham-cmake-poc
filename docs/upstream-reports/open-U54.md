# Open finding U54: FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity counts the central orbital in Q_R

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U54), found 2026-10-26 while writing the program's first independent tests; r4493 sources, GCC 15. The lines involved are unchanged in this repository (patch 03 only adds the built-in diagonaliser for builds without LAPACK); there was no pristine autotools build with LAPACK to run them on.

## What happens

The program's wiki page defines the charge imbalance as

    Q_L - Q_R = sum_{n < (N_phi + 1) / 2} sum_sigma (c+_{n,sigma} c_{n,sigma} - c+_{N_phi - n,sigma} c_{N_phi - n,sigma})

which pairs each orbital n with its mirror N_phi - n, so that a central orbital (N_phi even) drops out. The code (`FQHE/src/Programs/FQHEOnSphere/FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity.cc`, lines 577 and 635, the real and complex branches) instead adds the orbitals with `MomentumIndex < (LzMax + 1) / 2` in integer arithmetic and subtracts all the others. When N_phi is odd (an even number of orbitals, as in the wiki's N_phi = 21 example) the two agree. When N_phi is even, the central orbital N_phi / 2 is subtracted, so the printed "Charge imbalance between left and right (orbital cut)" is the wiki's value minus the central orbital's occupation in both layers. A state that is symmetric under the mirror then gets a non-zero imbalance.

Example: 2S = 6, the two-layer Hamiltonian of `tests/oracles/quasiholes_trs.py` case `pairing_2s6_sz0`, the three lowest states of 2Lz = 0. The program prints the imbalance eigenvalues -2.0039, -0.2284 and 1.9445; the wiki's definition gives -1.8769, 0.3544 and 2.0102. The printed values are, to 1e-14, those of the code's definition (the wiki's operator minus n_{3,up} + n_{3,down}), and every density matrix element and integrated charge printed by the same run matches the oracle to 5e-14.

## Fix direction

Skip the central orbital (`MomentumIndex == LzMax / 2` when `LzMax` is even), or say on the wiki page that it belongs to Q_R.

## Consequence for the tests

`openfinding.fqhe.sphere.quasiholes_time_reversal_density.imbalance_central_orbital.pairing_2s6_sz0_lz0` (WILL_FAIL) checks the wiki's imbalance for 2S = 6; the physics test for 2S = 5 (an even number of orbitals) checks it too, and passes.
