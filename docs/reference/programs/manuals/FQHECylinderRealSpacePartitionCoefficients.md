# FQHECylinderRealSpacePartitionCoefficients — manual

Source: DiagHam wiki page `FQHECylinderRealSpacePartitionCoefficients`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHECylinderRealSpacePartitionCoefficients](../FQHEOnCylinder/FQHECylinderRealSpacePartitionCoefficients.md)
[FQHECylinderRealSpacePartitionCoefficients](FQHECylinderRealSpacePartitionCoefficients.md) computes the geometrical coefficients required to evaluate a real space entanglement spectrum on the cylinder geometry. The cut splits the cylinder exactly into two equal parts. A typical usage is 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderRealSpacePartitionCoefficients -s 17 --cylinder-perimeter 10*

The cylinder is defined by its number of flux quanta (set by -s, here 17) and either its perimeter (using the --cylinder-perimeter option) or the aspect ratio (using the -r option). The previous example produces a text file realspace_cylinder_l_0.000000_perimeter_10.000000_2s_17.dat that should contain

```text
# real space coefficients for a sharp cut at x=0 on a cylinder with perimeter L=10 and N_phi=17
OrbitalSquareWeights = 0.99999999999998 0.99999999998671 0.9999999961694 0.9999994885244 0.99996814014954 0.99906462046345 0.98683946253913 0.90871192544881 0.67158291820171 0.32841708179829 0.091288074551193 0.013160537460871 0.00093537953654699 3.1859850463778e-05 5.1147559676501e-07 3.8306023353307e-09 1.3294254586071e-11 2.1260770921572e-14
```

Such a file should be used as the input file for real space entanglement spectrum codes such as [FQHESphereFermionEntanglementEntropyParticlePartition](FQHESphereFermionEntanglementEntropyParticlePartition.md).


[FQHECylinderRealSpacePartitionCoefficients](FQHECylinderRealSpacePartitionCoefficients.md) is also compatible with the MPS code, including for iMPS calculations. For that purpose, we just have to set the number of flux quanta to zero. Only the orbitals for the left hand side and having a non-zero or non equal to one geometry weights are computed. The error with which a weight is considered to be zero or one is set by the --error option. For example,

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderRealSpacePartitionCoefficients --cylinder-perimeter 14 -s 0 --error 1e-8*

will work for a cylinder of perimeter 14. It will produce a text realspace_cylinder_l_0.000000_perimeter_14.000000_2s_0.dat that should lokk like 

```text
# real space coefficients for a sharp cut at x=0 on an infinite cylinder with perimeter L=14
OrbitalSquareWeights = 0.99999980886422 0.99999556192682 0.99993000484254 0.99924690387519 0.99443817547374 0.97155094306578 0.89784991301619 0.7371871581773 0.5
```
