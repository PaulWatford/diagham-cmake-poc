# FQHESphereRealSpacePartitionCoefficients — manual

Source: DiagHam wiki page `FQHESphereRealSpacePartitionCoefficients`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereRealSpacePartitionCoefficients](../FQHEOnSphere/FQHESphereRealSpacePartitionCoefficients.md)
[FQHESphereRealSpacePartitionCoefficients](FQHESphereRealSpacePartitionCoefficients.md) is the analogue of [FQHECylinderRealSpacePartitionCoefficients](FQHECylinderRealSpacePartitionCoefficients.md) for the sphere geometry. It supports both patches preserving the rotational symmetry along z and patches breaking the full rotational symmetries.

## Sphere caps

For example, the weight file for a (north) hemisphere cut can be produced with the following command line

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereRealSpacePartitionCoefficients -s 15 --theta 0.5*

The output file has the default name (can be tuned with the --output-file option) realspace_sphere_theta_0.500000_2s_15.dat and should look like

```text
# real space coefficients for a cut at theta=0.5 on a sphere with N_phi=15
OrbitalSquareWeights = 0.99998474121094 0.99974060058594 0.99790954589844 0.98936462402344 0.96159362792969 0.89494323730469 0.77275085449219 0.59819030761719 0.40180969238281 0.22724914550781 0.10505676269531 0.038406372070313 0.010635375976563 0.0020904541015625 0.0002593994140625 1.52587890625e-05
```

Note that (square) weigths are sorted from the north orbital to the south one.


## Full rotational symmetry breaking patches

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereRealSpacePartitionCoefficients -s 15 --theta 0.5 --phi 0.5*
