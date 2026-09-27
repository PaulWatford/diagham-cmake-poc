# Open defect U31: FCICheckerboardLatticeModel segfaults without --single-band

Purpose: reproducer for the DiagHam authors; no fix is carried here (`knownbug.fci.checkerboard.two_band_model_segfault` in ctest holds it).
Source: `docs/reference/known-defects.md` (U31), found 2026-09-28; r4493 sources, `default` and `full` presets, GCC 15.

## What happens

```
FCICheckerboardLatticeModel -p 3 -x 3 -y 3 --flat-band --full-diag 5000
```

writes the tight-binding file, prints `nbr interaction = 0`, a line of `=`, `start`, then **segmentation fault** (exit 139). The same for `-p 2 -x 2 -y 3`, with or without `--flat-band`, and for `-p 4 -x 4 -y 3`.

## What was expected

The two-band checkerboard model for three fermions on 3×3 unit cells is a 816-state problem that should diagonalise in a fraction of a second. With `--single-band` (projection onto the lowest band) the same command runs normally and its 84-state spectrum is the regression case `regression.fci_checkerboard_singleband_n3_3x3`.

## What is known

`nbr interaction = 0` before the crash suggests the two-band interaction terms are not being generated for this Hamiltonian and a later loop reads past an empty array. Not investigated further; the FCI module's authors will know the intended path.

## Related, not investigated

`FQSHKagomeModel -p 2 -x 2 -y 3 --flat-band` also segfaults and `FQSH2DBHZModel -p 2 -x 2 -y 3 --flat-band` exits 255; the inputs may simply be invalid for those models.
