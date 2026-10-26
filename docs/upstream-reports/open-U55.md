# Open finding U55: FQHECylinderConfiningPotentialCoefficients gets a sign wrong for orbitals before the offset

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U55), found 2026-10-26 while writing the program's first independent tests; r4493 sources, GCC 15, built with GSL. The lines involved are as in r4493 (this repository changes only the output precision of the file); the pristine autotools build of r4493 has no GSL, so there it writes zeros (U50).

## What happens

For the polynomial right potential V(x) = V0 (x - x0)^alpha (x > x0), the matrix element of orbital k, centred at X_k = 2 pi (k - N_phi / 2) / L, is

    V0 / sqrt(pi) * int_c^inf (u - c)^alpha exp(-u^2) du,   c = x0 - X_k,

and expanding (u - c)^alpha gives sum_i C(alpha, i) (-c)^(alpha - i) u^i. `FQHECylinderComputePolynomialRightPotentialCoefficients` (`FQHE/src/Programs/FQHEOnCylinder/FQHECylinderConfiningPotentialCoefficients.cc`, from line 554) has two branches. For c < 0 (the orbital centre beyond the offset) it uses `pow(-cutPosition, alpha - i)` and is right. For c >= 0 (the centre before the offset, line 563) it uses `pow(cutPosition, alpha - i)`, dropping the sign (-1)^(alpha - i), so every term with alpha - i odd enters with the wrong sign. The left potential has the same line (line 615; see also U56).

The wiki's own example, `-s 24 --cylinder-perimeter 10.0 --confining-rightpower 2 --spinful --error 1e-7 --confining-rightoffset 4.5`, reproduces exactly and shows the effect: the coefficients rise and fall (0.518 at k = 18, then 0.311 at k = 19) where a confining potential grows monotonically. The integral gives, for k = 13 to 19, 6.3e-10, 8.8e-8, 6.0e-6, 2.0e-4, 0.0036, 0.035 and 0.198; the program writes 1.4e-6, 9.9e-5, 0.0032, 0.043, 0.246, 0.518 and 0.311. From k = 20 (centres beyond 4.5) on, the two agree to every printed digit. Any alpha >= 1 is affected (alpha = 0 takes a separate, correct `erf` path).

## Fix direction

Use `pow(-cutPosition, alpha - i)` in the c >= 0 branch as well (that is, (-c)^(alpha - i)), in both the right and the left functions.

## Consequence for the tests

In builds with GSL, `openfinding.fqhe.cylinder.confining_potential_coefficients.sign_before_offset.wiki_example` (WILL_FAIL) compares the wiki's example with the integral (`tests/oracles/confining.py`). Physics tests cover the correct branch (every orbital centre beyond the offset, alpha = 1, 2, 3) and, in every build, the step potential (alpha = 0).
