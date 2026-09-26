# FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis — manual

Source: DiagHam wiki page `FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis](../FQHEOnCylinder/FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis.md)

> **Options in this manual that the current program does not list**: `--cddcd-matrices`, `--cducu-matrices`, `--column`, `--cylinder-perimeter`, `--eigenstate`, `--force-negativelz`, `--interaction-file`, `--interaction-name`, `--nbr-lz`, `--nbr-particles`, `--use-cylinder`, `--use-hilbert`, `--use-lapack`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
[FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis](FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis.md) builds an effective Hilbert space for the two layer system with time reversal symmetry starting from the single layer spectra and eigenstates. It focuses on a single momentum sector and it also allows to truncate the spectrum. 


## Generating the single layer data

The single layer data (both the spectra and the eigenstate) should be generated using [FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing](FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing.md). The command should look like

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing -l 21 --force-negativelz --use-lapack --use-cylinder --cylinder-perimeter 8.0 --directory nphi_21_cylinder_perimeter_8.000000 --interaction-name symlinearmomentum --interaction-file confining_momentum_cylinder_perimeter_8.000000_2s_21_alphar_1_x0r_0.000000_v0r_1.000000_alphal_1_x0l_0.000000_v0l_1.000000.dat --nbr-particles 6 -s 6 -n 30 --eigenstate*

We just have to set the number of particles to match the  $S_z$ value, and repeat this operation for all admissible values of the number of particles.

## Building the effective basis

The effective basis is built as a product of two decoupled single layer states i.e. 

$|\Psi_{N_{\uparrow},N_{\Phi},L_{z,\uparrow},\alpha_\uparrow}\rangle\otimes |\Psi_{N_{\downarrow},N_{\Phi},L_{z,\downarrow},\alpha_\downarrow}\rangle$

where 
A typical usage og [FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis](FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis.md) to generate this basis is 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis -s spectra.dat --charging-energy 0.24966151390041 --average-nbrparticles 30.0  --nbr-states 200*

where spectra.dat contains the list of single layer spectra

```text
fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_0_2s_21_sz_0_lz.dat
fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_1_2s_21_sz_1_lz.dat
fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_2_2s_21_sz_2_lz.dat
fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_3_2s_21_sz_3_lz.dat
fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_4_2s_21_sz_4_lz.dat
fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_5_2s_21_sz_5_lz.dat
fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_6_2s_21_sz_6_lz.dat
fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_7_2s_21_sz_7_lz.dat
fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_8_2s_21_sz_8_lz.dat
```

The  charging energy is set by the two options  --charging-energy and --average-nbrparticles (similar to [FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing](FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing.md)). The option --nbr-states set the maximum number of states to keep. The actual number of kept states can actually be lower to preserve the spectrum degeneracies.

The standard output provides some useful information

```text
total number of levels in the Sz=0 Ky=0 : 320622
require 1 eigenstates for N=5 and Ky=-5
require 1 eigenstates for N=5 and Ky=5
require 1 eigenstates for N=6 and Ky=-24
require 1 eigenstates for N=6 and Ky=-14
require 2 eigenstates for N=6 and Ky=-12
require 1 eigenstates for N=6 and Ky=-10
require 2 eigenstates for N=6 and Ky=-4
require 2 eigenstates for N=6 and Ky=-2
...
```

It gives the number of single layer eigenstates that will be needed per number of particles and per momentum. The code will generate two ASCII files

- fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_0_2s_21_sz_0_lz_0_effective_198.vectors contains all the single layer eigenvectors that are required to construct the effective basis.
- fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_0_2s_21_sz_0_lz_0_effective_198.basis contains the information about the effective basis, including its spectrum.

This latest file should look like

```text
# N Sz Lz N_u N_d Lz_u Lz_d E E-E_c file_index_up file_index_down
12 0 0 6 6 0 0 133.74203360766 52.851703103928 50 50
12 0 0 6 6 12 12 134.27436375013 53.384033246394 81 81
12 0 0 6 6 -12 -12 134.27436375013 53.384033246394 23 23
12 0 0 6 6 0 0 135.82486382332 54.934533319586 51 50
12 0 0 6 6 0 0 135.82486382332 54.934533319586 50 51
12 0 0 6 6 2 2 135.82679499461 54.936464490878 60 60
12 0 0 6 6 -2 -2 135.82679499461 54.936464490878 44 44
12 0 0 6 6 -10 -10 136.30154339982 55.411212896088 29 29
12 0 0 6 6 10 10 136.30154339982 55.411212896088 78 78
10 0 0 5 5 5 5 136.30648840747 36.441882847308 8 8
10 0 0 5 5 -5 -5 136.30648840747 36.441882847308 4 4
14 0 0 7 7 -7 -7 136.30648840747 72.393140848968 98 98
14 0 0 7 7 7 7 136.30648840747 72.393140848968 102 102
12 0 0 6 6 12 12 136.31109083121 55.420760327478 81 82
12 0 0 6 6 -12 -12 136.31109083121 55.420760327478 23 24
12 0 0 6 6 -12 -12 136.31109083121 55.420760327478 24 23
12 0 0 6 6 12 12 136.31109083121 55.420760327478 82 81
....
```

Each column has the following meaning:
- 1st column : total number of particles $N$.
- 2nd column : twice the $S_z$ value.
- 3rd column : twice the $L_z$ value.
- 4th column : number of particles in the up layer $N_\uparrow$.
- 5th column : number of particles in the down layer $N_\downarrow$.
- 6th column : total momentum for the up layer $L_{z,\uparrow}$.
- 7th column : total momentum for the down layer $L_{z,\downarrow}$.
- 8th column : total energy $E$.
- 9th column : confining energy i.e. $E-E_c$. where $E_c$ is the charging energy
- 10th column : index of the single layer eigenstate used for the up layer (as described in fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_0_2s_21_sz_0_lz_0_effective_198.vectors).
- 11th column : index of the single layer eigenstate used for the down layer (as described in fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_0_2s_21_sz_0_lz_0_effective_198.vectors).


Note that the full spectrum of the two decoupled layer problem (including the charging energy) can be stored by adding the --write-fullspectrum option. This creates a file (in this example named fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_0_2s_21_sz_0_lz_0_fullspectrum.dat) that is similar to effective basis description (i.e. fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_0_2s_21_sz_0_lz_0_effective_198.basis) but without truncation and without the three last columns.

The charging energy can also be automatically optimize using the --optimize-chargingenergy option. This option takes an integer as a parameter and will adjust the factor in front of the charging energy such that the ground state is in the total number of particles sector set by --optimize-chargingenergy while the  the ground state in the other particle number sectors that differ by two units from this one are identical.

## Building the eigenstates

[FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis](FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis.md) can also export the eigenstates associated to each requested eigenvalues in the two layer system. In that case, we just have to add the --build-eigenstates option i.e.

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis -s spectra.dat --charging-energy 0.24966151390041 --average-nbrparticles 30.0  --nbr-states 40 --directory nphi_21_cylinder_perimeter_8.000000 --build-eigenstates*

Note that the  --directory option (already described in [FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing](FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing.md)) should be provided. The code will generate the eigenstates whose file names are deduced from the single layer energy spectrum files. In this example, they will look like fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_effective_39_fixedn_12_n_0_2s_21_sz_0_lz_0.0.vec .

The list of generated eigenstates can be exported in a format that is compatible with the --use-hilbert option of [FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing](FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing.md). You just have to include the --write-eigenstatebasis i.e.

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis -s spectra.dat --charging-energy 0.24966151390041 --average-nbrparticles 30.0  --nbr-states 40 --directory ../../nphi_21_cylinder_perimeter_8.000000 --build-eigenstates --write-eigenstatebasis*

As a sanity check, we can use this basis without the superconducting order parameter using the sams confining potential and charging energy that were used to generate these eigenstates

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing -l 21  --directory nphi_21_cylinder_perimeter_8.000000 --interaction-name test --use-lapack --use-cylinder --cylinder-perimeter 8.0  --interaction-file confining_momentum_cylinder_perimeter_11.757129_2s_21_alphar_1_x0r_0.000000_v0r_1.000000_alphal_1_x0l_0.000000_v0l_1.000000.dat --use-hilbert fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_effective_39_n_0_2s_21_sz_0_lz_0.basis  --charging-energy 0.24966151390041 --average-nbrparticles 30.0 --nbr-lz 1*

The computed spectrum fermions_cylinder_perimeter_8.000000_su2_quasiholes_test_cenergy_0.249662_n0_30.000000_pairing_n_0_2s_21_sz_0_lz.dat should match the one contained in fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_0_2s_21_sz_0_lz_0_effective_39.basis .

## Building the effective Hamiltonian

[FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis](FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis.md) can directly compute all the building blocks to write the Hamiltonian discussed in [FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing](FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing.md) without using this code. We just need to add the --build-effectivehamiltonian option

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis -s spectra.dat --charging-energy 0.24966151390041 --average-nbrparticles 30.0  --nbr-states 200 --directory nphi_21_cylinder_perimeter_8.000000 --build-effectivehamiltonian*


This will generate many binary matrix file with three different suffixes

- cducu_m.mat for the matrix elements of $c^\dagger_{m,\uparrow} c_{m,\uparrow}$ in the effective basis.
- cddcd_m.mat for the matrix elements of $c^\dagger_{m,\downarrow} c_{m,\downarrow}$ in the effective basis.
- cddcdu_m.mat for the matrix elements of $c^\dagger_{m,\downarrow} c^\dagger_{-m,\uparrow} + {\rm h.c.}$ in the effective basis.

Three ASCII files with a .list extension are also produced. They contain the file names of all the generated binary matrix files in the first column and some 0.0 numerical coefficients in the second column. There is one file for each type of matrix (cducu, cddcd, cddcdu).

This way, the .cddcdu.list file can be used as a starting point for the input of [GenericHamiltonianDiagonalization](GenericHamiltonianDiagonalization.md) (in this case, the diagonal matrix coefficients can be set using the energy spectrum fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_n_0_2s_21_sz_0_lz_0_effective_198.basis, and --column 7 to specify the index of the column where the energy is stored). The .cducu.list and .cddcd.list files can be used to compute the charge imbalance of a vector generated using this method, using [FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity](FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity.md) and its options --cducu-matrices and --cddcd-matrices.
