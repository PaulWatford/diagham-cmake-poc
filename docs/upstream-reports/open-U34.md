# Open defect U34: the disk pseudopotential Hamiltonian gives wrong energies for N ≥ 3

Purpose: reproducer for the DiagHam authors; no fix is carried here (three `known-bug` tests, `knownbug.fqhe.disk.bosons_two_body_generic.*`, hold the exact spectra).
Source: `docs/reference/known-defects.md` (U34), found 2026-09-29 while writing the disk goldens; r4493 sources, `default` (no GSL) and `full` (GSL) presets, GCC 15.

## What happens

```
printf 'Pseudopotentials = 1 0 0 0 0 0 0\n' > v0.dat
FQHEDiskBosonsTwoBodyGeneric -p 3 --minimum-momentum 6 --maximum-momentum 6 \
    --interaction-file v0.dat --interaction-name v0 --full-diag 5000
```

writes, for three bosons at total angular momentum L_z = 6 with the V₀ pseudopotential:

```
0  1.5951314  1.9824860  2.0000000  2.1174040  2.8832800  5.9841990
```

## What was expected

```
0  1.5  1.875  2.0625  2.25  3  6
```

that is 0, 3/2, 15/8, 33/16, 9/4, 3, 6. Two independent computations give these numbers (`tests/oracles/geometry_ed.py`):

- second quantisation: the pair states (z₁−z₂)^m (z₁+z₂)^M expanded in the orbitals z^k e^{−|z|²/4}/√(2π 2^k k!) with binomials, H = 2 Σ_m V_m Σ_{i<j} P_m(ij) assembled on the Fock basis exactly as for the sphere and torus oracles that reproduce DiagHam's sphere and torus programs to 10⁻⁹;
- first quantisation with exact rational arithmetic: P_m(ij) keeps the (z_i−z_j)^m component of a symmetric polynomial as a function of (z_i+z_j, z_i−z_j) (the relative-angular-momentum decomposition is orthogonal, so no metric is needed); H acts on the monomial symmetric polynomials of degree 6; the 7×7 rational matrix has the eigenvalues above.

The two agree to 10⁻¹⁵. Other cases: four bosons at L_z = 12 with V₀ = 1, V₂ = ½ (34 states) differ by up to 0.42; three bosons at L_z = 9 with V₂ = 1 give 0.0586, 0.5685, 0.8436 … instead of 3/32, 3/4, 27/32 ….

## What is known

- The zero modes are right: the Laughlin state and its quasihole countings (this repository's `physics.fqhe.disk.bosons_laughlin_1_2.*` goldens) come out exactly, and every zero-mode count matches the oracle in every case tried. Only the nonzero energies are wrong.
- N = 2 is right at every L_z: with a single V_m = 1 the pair of relative angular momentum m sits at 2 for L_z = m, m+2, m+5 (m = 0, 2, 4). So the pair states (`EvaluateCGCoefficient`, formula B7 of PRB 95, 245117) have the right norm; the problem appears when several pairs share particles, i.e. in the many-body assembly (`EvaluateInteractionFactors`, the bosonic branch with `BosonicFactor` 4/2/1) or in the operator it is applied with.
- The same numbers come out of the `default` and `full` builds (the "in-built routine" branch is used in both; `UseGSL` is not switched on by the program).
- `FQHEDiskFermionsTwoBodyGeneric` uses the same Hamiltonian class but hangs (U29), so the fermionic side could not be tested.
- `FQHEDiskBosonsDelta` is affected in some other way: two bosons at L_z = 2 give two nonzero levels (2 and 2.159), which no contact interaction can produce; not understood.
