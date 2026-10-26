# Open finding U50: without GSL, the sphere and disk real-space partition programs write all-zero weights and exit 0

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U50), found 2026-10-26 while writing the programs' first independent tests; r4493 sources, GCC 15. The same on a pristine autotools build of r4493 (configured without GSL).

## What happens

`FQHESphereRealSpacePartitionCoefficients` and `FQHEDiskRealSpacePartitionCoefficients` compute the orbital weights with GSL's incomplete beta and gamma functions. In a build without GSL they print `gsl is required` once per orbital, then still write the weight file, with every weight 0, and exit with status 0:

```
FQHESphereRealSpacePartitionCoefficients -s 15 --theta 0.5     # writes OrbitalSquareWeights = 0 0 ... 0, exit 0
FQHEDiskRealSpacePartitionCoefficients -s 10 --radius 2.5      # the same
```

A script, or a real-space entanglement run that reads the file, goes on with zero weights. With GSL both programs agree with the closed forms (incomplete beta for a sphere cap, incomplete gamma for a disk) to better than 1e-14 (`tests/oracles/realspace.py`).

## Fix direction

Stop with an error (and a non-zero exit status) when GSL is not available, or compute the weights without GSL: both have closed forms as finite sums of positive terms.

## Consequence for the tests

In builds without GSL, `openfinding.fqhe.sphere.realspace_partition_coefficients_zero_without_gsl.sphere_2s15_hemisphere` and `openfinding.fqhe.disk.realspace_partition_coefficients_zero_without_gsl.disk_2s10_radius2.5` (WILL_FAIL) fail while the zeros are written. In builds with GSL the same cases, and three more, are physics tests.

Found later, the same pattern in `FQHECylinderConfiningPotentialCoefficients`: without GSL every power >= 1 prints "error, GSL library is needed" once per term and writes all-zero coefficients, exit status 0 (the pristine autotools build of r4493, configured without GSL, does the same). The step potential (power 0) uses `erf` and is unaffected. `openfinding.fqhe.cylinder.confining_potential_coefficients.zero_without_gsl.power2_right_2s24_perimeter10_x-9` (WILL_FAIL, builds without GSL) documents it. With GSL the program has two further findings, U55 and U56.
