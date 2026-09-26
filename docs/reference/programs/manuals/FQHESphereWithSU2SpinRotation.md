# FQHESphereWithSU2SpinRotation — manual

Source: DiagHam wiki page `FQHESphereWithSU2SpinRotation`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereWithSU2SpinRotation](../FQHEOnSphere/FQHESphereWithSU2SpinRotation.md)

> **Options in this manual that the current program does not list**: `--all-sz`, `--eigenstate`, `--interaction-file`, `--interaction-name`, `--nbr-lz`, `--no-cross`, `--sum`, `--use-alt`, `--use-lapack`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
[FQHESphereWithSU2SpinRotation](FQHESphereWithSU2SpinRotation.md) can perform a spin rotation of any spinful bosonic state with a well defined total Sz. For fermions, you should use [FQHESphereFermionsWithSpinRotation](FQHESphereFermionsWithSpinRotation.md).

As an example, we can perform a pi/2 rotation of the ferromagnetic Laughlin 1/2 state computed in the largest Sz.

First we generate the fully polarized Laughlin 1/2 state

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereBosonsWithSpin -s 4 --use-alt --use-lapack -p 4 -l 6 --nbr-lz 1 --interaction-file pseudopotentials_v0.dat --interaction-name v0 --eigenstate -n 1*

Note that we used the --use-alt option of [FQHESphereBosonsWithSpin](FQHESphereBosonsWithSpin.md). [FQHESphereWithSU2SpinRotation](FQHESphereWithSU2SpinRotation.md) only works with data that has been generated with this option. We can now perform the rotation, providing the rotation angle through the -a option (in unit of pi)

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereWithSU2SpinRotation -i bosons_sphere_su2_v0_n_4_2s_6_sz_4_lz_0.0.vec -a 0.5*

This creates a file bosons_sphere_su2_v0_n_4_2s_6_syrotation_0.500000_lz_0.0.vec written is the spinful basis without Sz conservation. To test the code, we can generate all the spinful zero mode of the V0 interaction and check that the rotated state belong to that subspace. To generate the zero energy manifold

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereBosonsWithSpin --all-sz --use-alt --use-lapack -p 4 -l 6 --nbr-lz 1 --interaction-file pseudopotentials_v0.dat --interaction-name v0 --eigenstate -n 18*

and testing the overlap 

*\$PATHTODIAGHAM/build/src/Programs/GenericOverlap --sum --no-cross bosons_sphere_su2_v0_n_4_2s_6_syrotation_0.500000_lz_0.0.vec bosons_sphere_su2_v0_n_4_2s_6_lz_0.*.vec*

```text
File 0  bosons_sphere_su2_v0_n_4_2s_6_syrotation_0.500000_lz_0.0.vec
File 1  bosons_sphere_su2_v0_n_4_2s_6_lz_0.0.vec
...
Overlap |<0|17>|^2 = 0.02386494271652
Overlap |<0|18>|^2 = 0.0625
Total overlap = 1
```

If the rotated is known to have a Sz<->-Sz symmetry, this can be used to express the rotated state in the symmetrized basis using the option  --szsymmetrized-basis (and potentially the --minus-szparity option depending on which sector you are interested in). For example

*\$PATHTODIAGHAM/FQHE/src/Programs/FQHEOnSphere/FQHESphereWithSU2SpinRotation -i bosons_sphere_su2_v0_n_6_2s_10_sz_6_lz_0.0.vec -a 0.5 --szsymmetrized-basis*

will create a state bosons_sphere_su2_szsym_1_v0_n_6_2s_10_syrotation_0.500000_lz_0.0.vec expressed in the symmetrized Sz<->-Sz sector.
