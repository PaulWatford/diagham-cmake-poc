# Open defect U37: 3D --export-onebodytext writes only the kz = 0 slice

Purpose: reproducer for the DiagHam authors; no fix is carried here (`knownbug.fti.simple_ti3d.export_onebodytext_kz0_slice_only` in ctest holds it).
Source: `docs/reference/known-defects.md` (U37), found 2026-09-29 while writing the band-structure goldens; r4493 sources, `default` preset, GCC 15.

## What happens

```
FTI3DSimpleTI -p 2 -x 3 -y 3 -z 2 --four-bands --mass 1.5 --singleparticle-spectrum --export-onebodytext
```

`..._tightbinding.dat` has the two-dimensional layout (a `kx ky E_0 ...` header) and 9 rows, one per (kx, ky) at kz = 0. The plain `.dat` written by the same run has all 18 momenta.

## What was expected

18 rows with a kz column, or an explicit statement that only kz = 0 is exported.

## What is known

`Abstract3DTightBindingModel` derives from `Abstract2DTightBindingModel` and does not override `WriteBandStructureASCII` (`FTI/src/Tools/FTITightBinding/Abstract2DTightBindingModel.cc`), so the 2D writer is used for a 3D model. The band tests read the plain `.dat` for this reason.
