# FQHESphereUnnormalizeState — manual

Source: DiagHam wiki page `FQHESphereUnnormalizeState`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereUnnormalizeState](../FQHEOnSphere/FQHESphereUnnormalizeState.md)

> **Options in this manual that the current program does not list**: `--boson`, `--get-index`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
if you start from an unnormalized state and want to go to the sphere geometry

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereUnnormalizeState -i bosons_unnormalized_haldane_laughlin2_n_8_2s_14_lz_0.0.vec -o bosons_conformal_haldane_laughlin2_n_8_2s_14_lz_0.0.vec --haldane --reference-file laughlin2_n_8_2s_14.dat --normalize*


----

if you start from an unnormalized state and want to go to the conformal limit

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereUnnormalizeState -i bosons_unnormalized_haldane_laughlin2_n_8_2s_14_lz_0.0.vec -o bosons_conformal_haldane_laughlin2_n_8_2s_14_lz_0.0.vec --haldane --reference-file laughlin2_n_8_2s_14.dat --symmetry-factor --normalize*

In order for the reduced density matrix to be normalized to one, you should normalize the vector to one (with respect to the basic scalar product)

*\$PATHTODIAGHAM/build/src/Programs/NormalizeVector bosons_conformal_haldane_laughlin2_n_8_2s_14_lz_0.0.vec*


----


if you start from an normalized state and want to go to the conformal limit


*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereUnnormalizeState -i bosons_laughlin2_n_8_2s_14_lz_0.0.vec -o bosons_conformal_laughlin2_n_8_2s_14_lz_0.0.vec --symmetry-factor --normalization 1635*

beware that unless the first component of your state is large (which is not the case in this example, since it is zero), you should pick a reference component using the --normalization option. A good component might be the one of the potential root configuration (here in this example 1635). The component index can be open using FQHESphereShowBasis

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereShowBasis -p 8 -l 14 -z 0 --boson --get-index "1 0 1 0 1 0 1 0 1 0 1 0 1 0 1"*  
memory requested for Hilbert space = 11Mo  
1635 : 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1  

----

if you start from a state in the conformal limit and want to go to the sphere geometry, just add the --conformal-limit compare to the the unnormalize to normalize case

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereUnnormalizeState -i bosons_unnormalized_haldane_laughlin2_n_8_2s_14_lz_0.0.vec -o bosons_conformal_haldane_laughlin2_n_8_2s_14_lz_0.0.vec --haldane --reference-file laughlin2_n_8_2s_14.dat --normalize --conformal-limit*

----

--symmetry-factor refers to the 1/n! for bosons.
