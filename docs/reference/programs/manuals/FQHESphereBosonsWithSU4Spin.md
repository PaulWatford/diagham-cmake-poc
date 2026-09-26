# FQHESphereBosonsWithSU4Spin — manual

Source: DiagHam wiki page `FQHESphereBosonsWithSU4Spin`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereBosonsWithSU4Spin](../FQHEOnSphere/FQHESphereBosonsWithSU4Spin.md)
[FQHESphereBosonsWithSU4Spin](FQHESphereBosonsWithSU4Spin.md) is the analogue of [FQHETorusBosonsWithSU4SpinTwoBodyGeneric](FQHETorusBosonsWithSU4SpinTwoBodyGeneric.md) for bosons having a SU(4) internal degree of freedom on the sphere geometry. It shares the options of both [FQHETorusBosonsWithSU4SpinTwoBodyGeneric](FQHETorusBosonsWithSU4SpinTwoBodyGeneric.md) for the SU(4) options and [FQHESphereFermionsWithSpin](FQHESphereFermionsWithSpin.md) for the sphere geometry related options. The interaction file that describes the two body interaction is identical to the one used for [FQHETorusBosonsWithSU4SpinTwoBodyGeneric](FQHETorusBosonsWithSU4SpinTwoBodyGeneric.md).

A typical usage of [FQHESphereBosonsWithSU4Spin](FQHESphereBosonsWithSU4Spin.md) can be exemplified by generating the generalized SU(4) Halperin state [2,1]. Creating an interaction file pseudopotentials_su4_21.dat containing a single line


```text
Pseudopotentials = 1.0
```


is enough to define the interaction that produces the exact interaction related to the generalized Halperin state [2,1], irrespective to the number of particles or flux quanta. The command line 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereBosonsWithSU4Spin -p 8 -l 8 --interaction-name halperin_21 --interaction-file pseudopotentials_su4_21.dat --use-lapack --nbr-n1 2 --nbr-n2 2 --nbr-n3 2 --nbr-n4 2 --memory 1000 -n 1 --show-itertime  -S --processors 2 --eigenstate --nbr-lz 2*

generates the ground state in the angular momentum sectors $L_z=0$ and $L_z=1$. The energy spectrum bosons_sphere_su4_halperin_21_n_8_2s_8_sz_0_iz_0_pz_0.dat should look like


```text
# Lz E
0 -5.3290705182008e-15
1 0.96485914575916
```

As we can, the energy of the ground state in the  $L_z=0$ sector has expected for the model interaction of the generalized Halperin state [2,1], while the energy of the ground state in the  $L_z=1$ is non zero.
