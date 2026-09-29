# Open defect U36: FTI3DHopf segfaults after writing its band file

Purpose: reproducer for the DiagHam authors; no fix is carried here (`knownbug.fti.hopf3d.segfault_after_band_file` in ctest holds it).
Source: `docs/reference/known-defects.md` (U36), found 2026-09-29 while writing the band-structure goldens; r4493 sources, `default` preset, GCC 15, and a pristine autotools build of r4493.

## What happens

```
FTI3DHopf -p 2 -x 3 -y 3 -z 3 --singleparticle-spectrum
```

The program writes `fermions_singleband_hopf_p_2_x_3_y_3_z_3_l_1.dat` and then exits with status 139 (SIGSEGV). The same happens for 2x2x2, for other `--lambda` values, and with `--export-onebodytext` added. The crash is identical on the CMake build and on the autotools build of the untouched r4493 tree, so it is not caused by the migration.

## What was expected

A clean exit after the spectrum has been written, as `FTI3DSimpleTI` does.

## What is known

The band file is complete before the crash, so the failure is in the code that runs after `WriteAsciiSpectrum` in `FTI/src/Programs/FTI/FTI3DHopf.cc`. The cause has not been isolated here.
