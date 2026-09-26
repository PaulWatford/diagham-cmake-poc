# FTIEntanglementEntropyParticlePartition — manual

Source: DiagHam wiki page `FTIEntanglementEntropyParticlePartition`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FTIEntanglementEntropyParticlePartition](../FTI/FTIEntanglementEntropyParticlePartition.md)

> **Options in this manual that the current program does not list**: `--export-densitymatrix`, `--import-densitymatrix`, `--selected-blocks`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
[FTIEntanglementEntropyParticlePartition](FTIEntanglementEntropyParticlePartition.md) computes the particle entanglement spectrum for systems in momentum space. Usage is similar to [FQHESphereBosonEntanglementEntropy](FQHESphereBosonEntanglementEntropy.md). 

## Basic usage

The code autodetects most of the system information such as the statistics or system size from the input vector names. The only exception if the number of bands (potentially including additional degrees of freedom such as spin or valley) that should be provided through the option --nbr-subbands (see e.g. [FTIDensity](FTIDensity.md)).
 
Compare to [FQHESphereBosonEntanglementEntropy](FQHESphereBosonEntanglementEntropy.md), the total density matrix will generally rely on more than one projector i.e.

$\rho=\frac{1}{N}\sum_{i=1}^{N} |\psi_i \rangle \langle \psi_i |$

In that case, you will need to use the --degenerated-groundstate option to give the list of states. Such a single column text file typically looks like

```text
fermions_singleband_checkerboardlattice_n_6_x_6_y_3_v_0.000000_t1_1.000000_t2_0.292893_gx_0.000000_gy_0.000000_kx_3_ky_0.0.vec
fermions_singleband_checkerboardlattice_n_6_x_6_y_3_v_0.000000_t1_1.000000_t2_0.292893_gx_0.000000_gy_0.000000_kx_3_ky_0.1.vec
fermions_singleband_checkerboardlattice_n_6_x_6_y_3_v_0.000000_t1_1.000000_t2_0.292893_gx_0.000000_gy_0.000000_kx_3_ky_0.2.vec
```


The command line has to be

*FTI/src/Programs/FTI/FTIEntanglementEntropyParticlePartition --degenerated-groundstate ground.dat --use-lapack --density-matrix fermions_singleband_checkerboardlattice_n_6_x_6_y_3_v_0.000000_t1_1.000000_t2_0.292893_gx_0.000000_gy_0.000000_kx_3_ky_0.0.full.parent --show-time*

where ground.dat is the column formatted text file described above. If you want to look at the entanglement spectrum, you should use the --density-matrix option to save the eigenvalues of the reduced density matrix and then use [FTIEntanglementSpectrum](FTIEntanglementSpectrum.md) to parse the output file. --show-time allows to show the amount of time to compute each block of the reduced density matrix. [FTIEntanglementEntropyParticlePartition](FTIEntanglementEntropyParticlePartition.md) supports the SMP mode. This can greatly decrease the amount of time to compute each block of the reduced density matrix.

fermions_singleband_checkerboardlattice_n_6_x_6_y_3_v_0.000000_t1_1.000000_t2_0.292893_gx_0.000000_gy_0.000000_kx_3_ky_0.0.full.parent conatins the eigenvalues of the reduced density matrix. The first column is the number of particles in the part A, the second column is the momentum along x for the part A, the third column is the momentum along x for the part A. Finally, the fourth column contains the eigenvalues themselves. This file typically looks like

```text
#  N    Kx    Ky    lambda
1 0 0 0.09590561524389
1 0 1 0.081030474931752
1 0 2 0.081030474931752
1 1 0 0.082982084898121
1 1 1 0.082012981462358
1 1 2 0.084055539822119
1 2 0 0.0796494013569
```

[FTIEntanglementEntropyParticlePartition](FTIEntanglementEntropyParticlePartition.md) also produces a text file containing the entanglement entropy for each cut. The file name is built from the first vector file name (in this example fermions_singleband_checkerboardlattice_n_6_x_6_y_3_v_0.000000_t1_1.000000_t2_0.292893_gx_0.000000_gy_0.000000_kx_3_ky_0.0.partent ).  The first column is the number of particles in the part A, the second column is the entanglement entropy, the third column is the trace of the reduced density matrix, the fourth column is $1-{\rm Tr}_A (\rho_A)$ .

## Advanced usage

[FTIEntanglementEntropyParticlePartition](FTIEntanglementEntropyParticlePartition.md) has a full support of MPI *(wiki page MPI)* and Scalapack *(wiki page Scalapack)* as opposed to other particle entanglement entropy codes. 

The default behavior of [FTIEntanglementEntropyParticlePartition](FTIEntanglementEntropyParticlePartition.md) is to consider all diagonal blocks of the reduced density. If we want to focus on very specific blocks (beyond choosing the particle cut via --min-na and --max-na), we can use the --selected-blocks option. This latest takes as an argument a column formatted ascii file that indicates which blocks of the reduced density matrix should be computed. Its format depends on the system symmetries. For example, the following file

```text
#Na Kxa Kya Sza SzSyma
2 0 1 0 -1
3 1 1 1 0
```

is for a system with the translation symmetry both in the x and y directions (respectively second and third columns), the $S_z$ symmetry (fourth column) and the $S_z \leftrightarrow -S_z$ symmetry (fifth column). The first column should always be the number of particles in the region $A$. If the translation symmetries are present, they should be the second and third column (if the the system is also invariant under the y translation). The next column is the $S_z$ symmetry (its position depends on the absence or presence of the translation symmetries). The last column (optional) is the $S_z \leftrightarrow -S_z$ symmetry.

[FTIEntanglementEntropyParticlePartition](FTIEntanglementEntropyParticlePartition.md) can also read and write diagonal blocks of the reduced density. This capability should be used with the --selected-blocks option to select the diagonal block. Then --export-densitymatrix compute the diagonal block without diagonalizing it and store it as a binary matrix. The output name is either deduced from --density-matrix option (is enabled) or a default name like densitymatrix_na_2_kxa_0_kya_1_sza_0_szsyma_1.mat (containing the diagonal block quantum numbers) is used. Then this file can be read from the drive to perform the diagonalization. We just have to use the --import-densitymatrix option. It takes the matrix file name as an argument. All the other options should be identical to those used during the generation phase.
