# FQHEDiskRealSpacePartitionCoefficients — manual

Source: DiagHam wiki page `FQHEDiskRealSpacePartitionCoefficients`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHEDiskRealSpacePartitionCoefficients](../FQHEOnDisk/FQHEDiskRealSpacePartitionCoefficients.md)
[FQHEDiskRealSpacePartitionCoefficients](FQHEDiskRealSpacePartitionCoefficients.md) is the analogue of [FQHECylinderRealSpacePartitionCoefficients](FQHECylinderRealSpacePartitionCoefficients.md) for the disk geometry. A typical usage is 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnDisk/FQHEDiskRealSpacePartitionCoefficients -s 27 --radius 5*

The number of flux qunanta is specified by -s and the radius of the circular region $A$ centred at the origin is set by --radius. The code generates a text file realspace_disk_radius_5.000000_2s_27.dat that should contain

```text
# real space coefficients for disk cut of radius 5 on a disk with N_phi=27
OrbitalSquareWeights = 0.99999627334683 0.99994969018218 0.99965854540311 0.99844544215699 0.99465449451286 0.98517712540256 0.96543260642275 0.93017453681595 0.87508380305533 0.79856889505446 0.70292526005338 0.59423931118852 0.48102478112095 0.37216465605598 0.2749681158194 0.19397099895558 0.13069200140573 0.084163326736715 0.051851747105456 0.030594128926997 0.017308117565459 0.0093997774693063 0.0049064024146738 0.0024643507545475 0.0011924488482317 0.0005564978950738 0.00025075224451712 0.00010920333222236
```

This file can be used with codes such as [FQHESphereFermionEntanglementEntropyParticlePartition](FQHESphereFermionEntanglementEntropyParticlePartition.md) to perform the real space entanglement spectrum.
