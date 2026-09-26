# FQHESphereFermionEntanglementEntropyParticlePartition — manual

Source: DiagHam wiki page `FQHESphereFermionEntanglementEntropyParticlePartition`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereFermionEntanglementEntropyParticlePartition](../FQHEOnSphere/FQHESphereFermionEntanglementEntropyParticlePartition.md)

> **Options in this manual that the current program does not list**: `--fermion`, `--normalize`, `--radius`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
FQHESphereFermionEntanglementEntropyParticlePartition allows to compute both the particle entanglement spectrum and the real space entanglement spectrum for a state on the sphere geometry. It is the analogous of [FQHESphereBosonEntanglementEntropy](FQHESphereBosonEntanglementEntropy.md) which handles the orbital entanglement spectrum. Note that most of this code valid for any genus zero geometry (i.e. the sphere, the disk and the cylinder geometries).


## Particle partition

As an example, let generate the Laughlin state $\nu=1/3$ on the sphere geometry using [FQHESphereJackGenerator](FQHESphereJackGenerator.md)

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator --reference-file laughlin3_n_8_2s_21.dat -o fermions_haldane_laughlin3_n_8_2s_21_lz_0.0.vec -a -2 --fermion --normalize*

where the root file laughlin3_n_8_2s_21.dat is

```text
NbrParticles = 8
LzMax = 21
ReferenceState = 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1
```


This produces our input vector fermions_haldane_laughlin3_n_8_2s_21_lz_0.0.vec . We can now generate the particle entanglement spectrum using

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereFermionEntanglementEntropyParticlePartition fermions_haldane_laughlin3_n_8_2s_21_lz_0.0.vec --haldane --reference-file laughlin3_n_8_2s_21.dat --use-lapack --use-svd --density-matrix fermions_haldane_laughlin3_n_8_2s_21_lz_0.0.full.parent --show-time*

Since the input vector given by [FQHESphereJackGenerator](FQHESphereJackGenerator.md) is in the squeezed basis, we need to provide the two options -haldane --reference-file. If we want the spectrum of the reduced density matrix, the option --density-matrix with the name of the output file (here fermions_haldane_laughlin3_n_8_2s_21_lz_0.0.full.parent ) has to be provided. The command produces two files fermions_haldane_laughlin3_n_8_2s_21_lz_0.0.partent and fermions_haldane_laughlin3_n_8_2s_21_lz_0.0.full.parent. The former contains the entanglement entropies and should look like

```text
1 3.0910424533583 1 4.5519144009631e-15
2 5.2181973385854 1 4.7739590058882e-15
3 6.5870098686786 1 3.663735981263e-15
4 7.1171383310127 0.99999999999999 1.0658141036402e-14
```

The fist column is the number of kept particles (i.e. $N_A$) which defines how the system is cut into two, the second column is the entanglement entropy for that specific cut, the third column is the trace of the reduced density matrix (should be 1) and the fourth column is the error on the trace.

The entanglement spectrum is contained in the file fermions_haldane_laughlin3_n_8_2s_21_lz_0.0.full.parent which looks like

```text
#  N    Lz    lambda
1 -21 0.045454545454545
1 -19 0.045454545454545
1 -17 0.045454545454545
1 -15 0.045454545454545
1 -13 0.045454545454545
```

The first column is the number of kept particles (i.e. $N_A$) , the second is the angular momentum of the part A (i.e. $L_{z,A}$) and the third columns gives the eigenvalues of the reduced density matrix.

While the --use-svd option is not mandatory, it greatly improves the calculation accuracy. It can be used when considering a single state (i.e. it cannot be used at the same time than the --degenerated-groundstate option).

## Real space partition

The real space partition works very similarly to the particle partition up to the additional --realspace-cut option

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereFermionEntanglementEntropyParticlePartition fermions_haldane_laughlin3_n_8_2s_21_lz_0.0.vec --haldane --reference-file laughlin3_n_8_2s_21.dat --use-lapack --use-svd --density-matrix fermions_haldane_laughlin3_n_8_2s_21_lz_0.0.full.rsent --show-time --realspace-cut*

The entanglement spectrum will be stored in fermions_haldane_laughlin3_n_8_2s_21_lz_0.0.full.rsent and its contain is identical to the one described  in the particle partition section. The entropy file is a bit different.The output name is fermions_haldane_laughlin3_n_8_2s_21_lz_0.0._thetabot_0.000000_thetabot_90.000000.realent  but the the columns have a different meaning.

```text
0 1.2716482141353e-08 5.9881528831023e-10 0.99999999940118
1 0.00047384283241268 4.4786234605411e-05 0.99995521376539
2 0.078728943275143 0.013572509551709 0.98642749044829
3 0.87012003614535 0.22346904767986 0.77653095232014
4 1.7779680830002 0.52582731186886 0.47417268813114
5 0.87012003614534 0.22346904767986 0.77653095232014
6 0.078728943275139 0.013572509551708 0.98642749044829
7 0.00047384283241246 4.4786234605392e-05 0.99995521376539
8 1.2716482140749e-08 5.9881528828039e-10 0.99999999940118
# Entanglement entropy = 3.676613752939
# Total trace = 0.99999999999883
```


The first column is still the number of particles in A, the second (resp. third) column is the contribution to the entropy (resp. to the trace) from the block of the reduced density matrix with $N_A$ particles. Thus to get the entanglement entropy (resp. trace) for the given real space cut, one has to sum up all the entries of the second (resp. third) column. These two sums are actually given by the last two commented lines.


By default the cut is set to the hemisphere. It can be changed by using the --realspace-theta-bot option. One can thus define the cut by the azimuthal angle ranging from 0 to --realspace-theta-bot  (set to 90 degrees by default). Notice that the starting angle can also be changed by using the --realspace-theta-top option.

## Other geometries preserving one symmetry (rotation like)

As already mentioned, most of this code valid for any genus zero geometry (i.e. the sphere, the disk and the cylinder geometries). In particular, the particle partition can be used with any input file describing a state in any o these geometries. For the real space partition, we need to provide some external information about the geometry. The option --realspace-generic allows to provide the geometrical weight of each orbital on the region $A$. More precisely, let $\phi_m(\vec{r})$ be the single particle wavefunction associated to the $m$-th orbital, its geometrical weight $\alpha_m$ over $A$ is defined as 

$\alpha_m=\int_A d^2 r \left|\phi_m(\vec{r})\right|^2$

All the $\alpha_m$ are positive numbers lower or equal to one. They should be provided in a text file that should look like

```text
OrbitalSquareWeights = 1.0 0.8 0.5 0.2 0.0
```

for a system of five orbitals (4 flux quanta). If the number of weights is lower than the number of orbitals, an equal number of weights one to the left and weights 0 to the right will be added. 

For the cylinder geometry, these coefficients can be computed using [FQHECylinderRealSpacePartitionCoefficients](FQHECylinderRealSpacePartitionCoefficients.md). For the disk geometry, we can use [FQHEDiskRealSpacePartitionCoefficients](FQHEDiskRealSpacePartitionCoefficients.md) and [FQHESphereRealSpacePartitionCoefficients](FQHESphereRealSpacePartitionCoefficients.md) for the sphere.

Let's look at a simple example. Assuming that we have generated the Laughlin $\nu=1/3$ state on the disk geometry for 8 fermions using both [FQHESphereJackGenerator](FQHESphereJackGenerator.md) and [FQHEDiskUnnormalizeState](FQHEDiskUnnormalizeState.md), we can compute the geometrical weights with 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnDisk/FQHEDiskRealSpacePartitionCoefficients -s 21 --radius 4*

To compute the corresponding real space entanglement spectrum, we just have to use

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereFermionEntanglementEntropyParticlePartition fermions_disk_normalized_n_8_2s_21_lz_0.0.vec --haldane --reference-file laughlin3_n_8_2s_21.dat --realspace-cut --realspace-generic realspace_disk_radius_4.000000_2s_21.dat --use-svd -o fermions_disk_normalized_n_21_2s_8_lz_0.0.fgenericcut.full.realent*

## Patches breaking all rotation/translation symmetries
