# FQHESphereExtractRootPartitions — manual

Source: DiagHam wiki page `FQHESphereExtractRootPartitions`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereExtractRootPartitions](../FQHEOnSphere/FQHESphereExtractRootPartitions.md)

> **Options in this manual that the current program does not list**: `--eigenstate`, `--full-diag`, `--interaction-file`, `--interaction-name`, `--nbr-lz`, `--use-lapack`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
[FQHESphereExtractRootPartitions](FQHESphereExtractRootPartitions.md) extracts admissible configurations (i.e. root configurations) from a set of eigenstates. It can be used in two different ways. As an example, let's generate some quasihole states of the Halperin (332) phase using [FQHESphereFermionsWithSpin](FQHESphereFermionsWithSpin.md)

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereFermionsWithSpin -p 6 -l 13 --interaction-name halperin_332 --nbr-lz 2 --use-lapack --interaction-file pseudopotentials_332_2s_13.dat --full-diag 8000 -n 20 --eigenstate*

where we input the pseudo-potential file pseudopotentials_332_2s_13.dat

```text
Pseudopotentials = 1 1 0 0 0 0 0 0 0 0 0 0 0 0
```

We get 7 quasihole states at $S_z=0$, 4 in the $L_z=0$ sectors and 3 in the  $L_z=1$ sectors. Now providing the spectrum to [FQHESphereExtractRootPartitions](FQHESphereExtractRootPartitions.md)

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereExtractRootPartitions fermions_sphere_su2_halperin_332_n_6_2s_13_sz_0_lz.dat*

will produce a file fermions_sphere_su2_halperin_332_n_6_2s_13_sz_0_lz.root that contains the root configurations of the 7 states

```text
total nbr of root partitions : 10
nbr of root partitions per Lz sector : 
4 3 

-----------------------------------------
Lz = 0 : 
d 0 u 0 0 d 0 0 u 0 0 d 0 u 
d 0 u 0 0 u 0 0 d 0 0 d 0 u 
d 0 0 d 0 u 0 0 d 0 u 0 0 u 
u 0 0 d 0 u 0 0 d 0 u 0 0 d 
nbr states = 4
-----------------------------------------
Lz = 1 : 
d 0 u 0 0 0 d 0 u 0 0 d 0 u 
d 0 0 d 0 u 0 0 u 0 0 d 0 u 
u 0 0 d 0 u 0 0 d 0 0 d 0 u 
nbr states = 3
-----------------------------------------
```

The first line gives the total number of states, assuming the spectrum is $L_z \leftrightarrow -L_z$ symmetric (so here 4+3+3=10). The third line gives the number of quasihole states per positive $L_z$ sector. Then for each $L_z$ sector, we obtain the root configurations.

Instead of providing an energy spectrum, we can provide a list of states via the --eigenstate-list option. This option uses a single column text file that lists all the states that have to be considered.
