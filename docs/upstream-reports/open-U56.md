# Open finding U56: FQHECylinderConfiningPotentialCoefficients puts the left confining potential on the right

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U56), found 2026-10-26 while writing the program's first independent tests; r4493 sources, GCC 15, built with GSL. The lines involved are as in r4493; the pristine autotools build of r4493 has no GSL, so there it writes zeros (U50).

## What happens

The help and the wiki define the left potential as V_L(x) = V0,L (x0,L - x)^alpha_L for x < x0,L and zero above. Its matrix element for orbital k is V0 / sqrt(pi) * int_0^inf t^alpha exp(-(t + X_k - x0)^2) dt, largest for the orbitals on the left (small k).

`FQHECylinderComputePolynomialLeftPotentialCoefficients` (`FQHE/src/Programs/FQHEOnCylinder/FQHECylinderConfiningPotentialCoefficients.cc`, from line 606) is a line-for-line copy of the right-potential function, and the caller passes it the same orbital index and the left offset unchanged. So it computes the right potential V0,L (x - x0,L)^alpha_L for x > x0,L instead.

Example: `-s 24 --cylinder-perimeter 10.0 --confining-rightstrength 0 --confining-leftstrength 1 --confining-leftpower 2 --confining-leftoffset -4.5`. The program writes coefficients that grow to the right, from 3.3e-4 at k = 0 to 145.5 at k = 24. The left potential gives 9.74 at k = 0 and falls to below 1e-17 from k = 14 on.

## Fix direction

Mirror the geometry in the left function (or at the call): the left matrix element of orbital k is the right one with X_k and x0 replaced by -X_k and -x0, i.e. pass `-orbitalIndex` and `-LeftShift` to the right function.

## Consequence for the tests

In builds with GSL, `openfinding.fqhe.cylinder.confining_potential_coefficients.left_potential_on_the_right.power2_2s24_perimeter10_x-4.5` (WILL_FAIL) compares that example with the integral (`tests/oracles/confining.py`).
