# FQHESphereWithSU2SpinLLLProjectedLLLTimesManyLL — manual

Source: DiagHam wiki page `FQHESphereWithSU2SpinLLLProjectedLLLTimesManyLL`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereWithSU2SpinLLLProjectedLLLTimesManyLL](../FQHEOnSphere/FQHESphereWithSU2SpinLLLProjectedLLLTimesManyLL.md)

> **Options in this manual that the current program does not list**: `--use-alt`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
[FQHESphereWithSU2SpinLLLProjectedLLLTimesManyLL](FQHESphereWithSU2SpinLLLProjectedLLLTimesManyLL.md) generates spinful composite fermion states projected onto the lowest Landau level.

## Basic usage

To generate the bosonic spinless $\nu=2/3$ CF state on the sphere geometry, we just have to use

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereWithSU2SpinLLLProjectedLLLTimesManyLL -p 6 --nbr-ll 2 -s 6 --interaction-name cf_2ll --normalize --disable-lzsymmetry*

This will create a binary vector named bosons_sphere_su2_cf2ll_n_6_2s_6_sz_6_lz_0.0.vec. The middle label in the file name is set by the --interaction-name. The number of Lambda level is fixed by the --nbr-ll. We set twice the total spin projection along z via the -s option. Here it matches the number of particles set by the -p option. By default, the code assumes filled shells. So in the previous example, we have two filled shells of CF. Beware that for spinful bosons, the --use-alt basis is implicit. In many codes (such as [FQHESphereShowBasis](FQHESphereShowBasis.md)), this option has to be set manually.

## Advanced usage

The code is fully parallel but the parallelization is optimized for a number of cores/processors that is a divisor of the number of particles.
