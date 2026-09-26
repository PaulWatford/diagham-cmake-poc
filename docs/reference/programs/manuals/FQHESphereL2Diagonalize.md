# FQHESphereL2Diagonalize — manual

Source: DiagHam wiki page `FQHESphereL2Diagonalize`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereL2Diagonalize](../FQHEOnSphere/FQHESphereL2Diagonalize.md)
[FQHESphereL2Diagonalize](FQHESphereL2Diagonalize.md) diagonalizes the $L^2$ operator of the sphere. Its usage is similar to [QHEFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md) up to some minor differences. It can handle both bosons and fermions. The statistics is set through the *--statistics* option that can take either of the two values : bosons or fermions. The diagonalization can be done in a squeezed basis. In that case the root configuration has to be provided through a text file identical to the those used for [FQHESphereJackGenerator](FQHESphereJackGenerator.md). Here is an example for the bosonic Laughlin state with 4 particles

```text
NbrParticles=4
LzMax=6
ReferenceState=1 0 1 0 1 0 1
```


If this file is save under *laughlin_n_4_2s_6.dat*, one can generate the corresponding state using the property the Laughlin state is the unique densest $L=0$ having this root configuration

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereL2Diagonalize -p 4 -l 6 -z 0 --statistics bosons --haldane --reference-file laughlin_n_4_2s_6.dat --use-lapack --eigenstate -n 1*

*-z* allows to set the $2 L_z$ and  *--haldane* indicates that we are using a squeezed basis. Running this command line generates both the spectrum (*bosons_l2_n_4_2s_6_lz.dat*) and the groundstate (*bosons_l2_n_4_2s_6_lz_0.0.vec*) which is that precise case, the Laughlin state.
