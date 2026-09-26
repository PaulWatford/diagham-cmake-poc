# FCIKagomeLatticeModel — manual

Source: DiagHam wiki page `FCIKagomeLatticeModel`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FCIKagomeLatticeModel](../FCI/FCIKagomeLatticeModel.md)

> **Options in this manual that the current program does not list**: `--eigenstates`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
FCIKagomeLatticeModel handles interacting fermions or bosons on the Kagome lattice. This exact diagonalization code offers similar calculation features to those of [FCICheckerboardLatticeModel](FCICheckerboardLatticeModel.md). 

## How to produce a bosonic Laughlin like state

We will show how to produce a Laughlin bosonic state at filling $\nu=1/2$. We will use only the on-site interaction and turn to zero the next nearest neighbor hopping term. To get the energy spectrum for 2 bosons on a 2x2 lattice, we just have to run

\$PATHTODIAGHAM/build/FTI/src/Programs/FCI/FCIKagomeLatticeModel -p 2 -x 2 -y 2 --flat-band --use-lapack --boson --t1 1.0 --l1 1.0 --t2 0.0 --l2 0.0 

This produces an ASCII text file *bosons_singleband_kagomelattice_n_2_x_2_y_2_t1_1_t2_0_l1_1_l2_0_gx_0_gy_0.dat* which contains

```text
# eigenvalues
# kx ky  E
0 0 -2.7755575615629e-17
0 0 0.03349364905389
0 0 0.41666666666667
0 0 0.46650635094611
0 1 0.25
0 1 0.33333333333333
1 0 0.25
1 0 0.33333333333333
1 1 0.25
1 1 0.33333333333333
```


For this particular choice of the lattice aspect ratio (a 2x2 lattice), the almost twofold degenerate groundstate is located in the momentum sector $(k_x,k_y)=(0,0)$. To get the two eigenstates corresponding to these two low energies, we can run the code as followed

\$PATHTODIAGHAM/build/FTI/src/Programs/FCI/FCIKagomeLatticeModel -p 2 -x 2 -y 2 --flat-band --use-lapack --boson --t1 1.0 --l1 1.0 --t2 0.0 --l2 0.0 --eigenstate -n 2 --only-kx 0 --only-ky 0

The --eigenstate option requires to compute the eigenstates and -n 2 indicates that we want only the two first eigenstates. Since we want to focus on the (0,0) sector, we have added the two options --only-kx 0 --eigenstates only-ky 0 . Running this command line will generate a new spectrum file and two binary vectors *bosons_singleband_kagomelattice_n_2_x_2_y_2_t1_1_t2_0_l1_1_l2_0_gx_0_gy_0_kx_0_ky_0.0.vec* and *bosons_singleband_kagomelattice_n_2_x_2_y_2_t1_1_t2_0_l1_1_l2_0_gx_0_gy_0_kx_0_ky_0.1.vec*. The extension *.x.vec* indicates this the the x-th eigenstate (zero being the groundstate in the given momentum sector).

The eigenstates are binary vectors that can be managed with the corresponding tools *(wiki page Binary vectors)*. Their decomposition onto the n-body basis can be read using [FTIShowBasis](FTIShowBasis.md).

If more than one eigenvalue or eigenstate per momentum are needed and the Hilbert space requires to use the Lanczos algorithm, it is safer to use the options --block-lanczos --block-size x --force-reorthogonalize . Here x should match the number of eigenvalues set by the -n option i.e.

\$PATHTODIAGHAM/build/FTI/src/Programs/FCI/FCIKagomeLatticeModel -p 8 -x 4 -y 4 --flat-band --use-lapack --boson --t1 1.0 --l1 1.0 --t2 0.0 --l2 0.0 -n 3 --block-lanczos --block-size 3 --force-reorthogonalize

## How to produce a fermionic Laughlin like state

We will show how to produce a Laughlin fermionic state at filling $\nu=1/3$. We will use only the nearest neighbor interaction and turn to zero the next nearest neighbor hopping term. To get the energy spectrum for 4 bosons on a 4x3 lattice, we just have to run

\$PATHTODIAGHAM/build/FTI/src/Programs/FCI/FCIKagomeLatticeModel -p 4 -x 4 -y 3 --t1 1.0 --l1 1.0 --t2 0.0 --l2 0.0 --use-lapack --flat-band -n 1 --eigenstate

According to the folding principle, the three almost degenerate states corresponding to the Laughlin groundstate manifold are located in the momentum sectors (kx,ky)=(2,0), (2,1) and (2,2). The three corresponding eigenstates have the following names :

```text
fermions_singleband_kagomelattice_n_4_x_4_y_3_t1_1_t2_0_l1_1_l2_0_gx_0_gy_0_kx_2_ky_0.0.vec
fermions_singleband_kagomelattice_n_4_x_4_y_3_t1_1_t2_0_l1_1_l2_0_gx_0_gy_0_kx_2_ky_1.0.vec
fermions_singleband_kagomelattice_n_4_x_4_y_3_t1_1_t2_0_l1_1_l2_0_gx_0_gy_0_kx_2_ky_2.0.vec
```

Their decomposition onto the n-body basis can be read using [FTIShowBasis](FTIShowBasis.md) and the particle entanglement spectrum can be computed using [FTIEntanglementEntropyParticlePartition](FTIEntanglementEntropyParticlePartition.md).
