# Open finding U60: FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis tests degeneracies on unsorted energies

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U60), found 2026-10-26 while writing the program's first independent tests; r4493 sources, GCC 15. The same on a pristine autotools build of r4493.

## What happens

The wiki page says that `--nbr-states` sets the maximum number of states to keep, and that "the actual number of kept states can actually be lower to preserve the spectrum degeneracies". The program sorts the two-layer energies into `TwoLayerSortedEnergies`, with the permutation in `TwoLayerIndices` (`FQHE/src/Programs/FQHEOnCylinder/FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis.cc`, line 282). The cut that should avoid splitting a degenerate level (line 311) instead compares `TwoLayerEnergies[EffectiveSubspaceDimension - 1]` with `TwoLayerEnergies[EffectiveSubspaceDimension]`: two entries of the unsorted array, in the order the levels were generated. So the cut ignores the degeneracy it should protect, and it can also stop early when two unrelated neighbours in generation order happen to be equal (mirror Lz sectors often are).

Example (`tests/oracles/trs_generate_basis.py`, case `2s6_n3_ec0.3_n0_3_keep20`): made-up single-layer spectra for N = 0 to 3 particles on 2S = 6, charging energy 0.3 around N0 = 3, `--nbr-states 20`. The 20th and 21st lowest two-layer levels are both at E = 3.336 (a level of three states); the program keeps 20 states, one of the three, where the rule gives 19. With `--nbr-states 35` and `50` on two other cases it keeps 35 and 50 where the rule gives 32 and 49. The full spectrum (`--write-fullspectrum`) and the kept states themselves are right: they are the lowest levels, with the right energies and single-layer indices.

## Fix direction

Compare `TwoLayerSortedEnergies[EffectiveSubspaceDimension - 1]` with `TwoLayerSortedEnergies[EffectiveSubspaceDimension]` (or the energies through `TwoLayerIndices`).

## Consequence for the tests

`openfinding.fqhe.cylinder.quasiholes_time_reversal_generate_basis.degeneracy_cut.2s6_n3_ec0.3_n0_3_keep20` (WILL_FAIL) applies the wiki's rule. The physics tests check the full spectrum and that the kept states are the lowest levels with consistent energies and indices.
