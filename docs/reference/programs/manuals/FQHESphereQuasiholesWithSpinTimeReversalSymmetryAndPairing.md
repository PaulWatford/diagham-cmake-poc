# FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing — manual

Source: DiagHam wiki page `FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing](../FQHEOnSphere/FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing.md)
[FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing](FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing.md) deals with two copies of single Landau level with opposite chiralities and superconducting regions at the poles. The Hamiltonian is projected into the subspace of Laughlin 1/3 quasiholes and reads

$H=\sum_{m,\sigma} V_{m,n,\sigma} c^\dagger_{m,\sigma} c_{m,\sigma} + \sum_{m} \Delta_{m} c^\dagger_{m,\uparrow} c^\dagger_{m,\downarrow}\;+\; h.c.$

The options are similar to [FQHESphereFermionsWithSpinTimeReversalSymmetryAndPairing](FQHESphereFermionsWithSpinTimeReversalSymmetryAndPairing.md) (more details can also be found in the documentation of [FQHESphereFermionsWithSpin](FQHESphereFermionsWithSpin.md)). Note that there is no option to fix the number of particles (for compatibility reasons, any output file still contains a _n_0_ string). An interaction file has to be provided using the "--interaction-file" option. This file specifies the form of the pairing term (OneBodyPotentialPairing) as well as that of the confining potential (OneBodyPotentialUpUp and OneBodyPotentialDownDown). The number of coefficients in each of these lines has to match the number of orbitals. Unlike [FQHESphereFermionsWithSpinTimeReversalSymmetryAndPairing](FQHESphereFermionsWithSpinTimeReversalSymmetryAndPairing.md), the two-body pseudopotentials should be left undefined, since these terms are only implemented through the projection onto the quasihole subspace. Finally, the average number of particles N_0 can be set through the option "--average-nbrparticles", with a charge energy E_c (N - N_0)^2 whose amplitude E_c is fixed by the option "--charging-energy".

[FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing](FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing.md) does not generate the Hilbert space or matrix elements from scratch. Instead, one has to run [FQHESphereQuasiholeMatrixElements](FQHESphereQuasiholeMatrixElements.md) with the same number of flux quanta first. The files containing the numbers of quasiholes for each spin (fermions_qh_states_k_1_r_2_nphi_*.dat), the pairing matrix elements in binary format (fermions_qh_k_1_r_2_n_*_nphi_*_lz_*_c_*.mat), the one-body potential matrix elements (fermions_qh_k_1_r_2_n_*_nphi_*_lz_*_cdc_*.mat) should be present either in the directory where the program is run or in any directory that should be indicated using the --directory option.

A typical command line would look like this:

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing -l 6 --charging-energy 1 --average-nbrparticles 4 --interaction-file onebodypotentials_2s_6.dat --interaction-name test --use-lapack --full-diag 5000*

If all the required files are located in the directory nphi_6, you should use

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing -l 6 --charging-energy 1 --average-nbrparticles 4 --interaction-file onebodypotentials_2s_6.dat --interaction-name test --use-lapack --full-diag 5000 --directory nphi_6*

As an example, here is the file content of onebodypotentials_2s_6.dat

```text
  OneBodyPotentialUpUp = 1.0 0 0 0 0 0 1.0
  OneBodyPotentialDownDown = 1.0 0 0 0 0 0 1.0
  OneBodyPotentialPairing = 0 1.0 0 0 0 1.0 0
```

[FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing](FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing.md) has a full support of both lapack and scalapack *(wiki page scalapack)*. When the full spectrum is needed with moderate Hilbert space dimension (between 5000 and 60000), scalapack *(wiki page scalapack)* provides a clear speed boost over lapack.

Note that by default, [FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing](FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing.md) only computes the positive angular momentum sectors. If the potentials are not symmetric, you should include the --force-negativelz to compute all the angular momentum sectors.

## Cylinder geometry

Most of the [FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing](FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing.md) are geometry agnostic. If the default geometry is the sphere geometry, ee can switch to the cylinder geometry thanks to a few option. First, we need to generate the data on such a geometry with [FQHECylinderQuasiholeMatrixElements](FQHECylinderQuasiholeMatrixElements.md). For example, let start with a cylinder where the perimeter is set to 6 magnetic length and with $N_\Phi=6$ flux quanta

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHECylinderQuasiholeMatrixElements -l 6 --cylinder-perimeter 6.0*

All the data will automatically be stored in nphi_6_cylinder_perimeter_6.000000 . Now we just have to call

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing -l 6  --charging-energy 1 --average-nbrparticles 4 --cylinder-perimeter 6 --use-cylinder --directory nphi_6_cylinder_perimeter_6.000000 --interaction-file onebodypotentials_2s_6.dat*

to perform the same calculation that described for the sphere but this time on a cylinder.

## Fixed number of particles

It is also possible to fix the total number of particles via the --nbr-particles option. For example we can simulate a single layer by setting both the total spin and the number of particles to the same value

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing -l 26 -s 7 --nbr-particles 7 --interaction-file confining_momentum_cylinder_perimeter_8.000000_2s_26_alphar_1_x0r_0.000000_v0r_1.000000_alphal_1_x0l_0.000000_v0l_0.000000.dat --use-cylinder --cylinder-perimeter 8.0 --directory nphi_26_n_7_cylinder_perimeter_8.000000 --use-lapack --force-negativelz --full-diag 8000 --interaction-name righlinearmomentum*

Here we effectively have a single layer with 7 particles and 26 flux quanta on a cylinder with a perimeter of $L_y=8$, the linear confining can be generated by [FQHECylinderConfiningPotentialCoefficients](FQHECylinderConfiningPotentialCoefficients.md).

If we want to look at all the particle number sectors in a single run, we can use both --nbr-particles and --all-fixednbrparticles. In that case, [FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing](FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing.md) will consider all the particle number sectors from 0 (or 1 is $2S_z$ is odd) to the integer set by --nbr-particles with an increment of two units. For example, if we want to compute the 30 lowest energy states of the problem of a symmetric confining potential without superconducting coupling, with  $2S_z=0$ and for all the number of particles going from 0 to 16, we just have to use

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing -l 21 --directory nphi_21_cylinder_perimeter_8.000000 --interaction-file confining_momentum_cylinder_perimeter_8.000000_2s_21_alphar_1_x0r_0.000000_v0r_1.000000_alphal_1_x0l_0.000000_v0l_1.000000.dat --interaction-name symlinearmomentum --use-lapack --full-diag 12000  -s 0  --use-cylinder --cylinder-perimeter 8.0  --memory 5000 -n 30  --block-lanczos --block-size 30  --force-reorthogonalize -S --processors 4 --eigenstate --nbr-lz 1 --nbr-particles 16 --all-fixednbrparticles*

Most of the options are similar to the previous example with the exception of --all-fixednbrparticles.


## Using an effective subspace

As many DiagHam codes, the diagonalization can be performed in an effective subspace described by a series of binary vectors (as long they belong to the same subspace defined by the various quantum numbers). For example, we can first generate the eigenstates for the problem with only the confining potential, keep only the lowest energy states and then diagonalize the superconducting coupling within that effective Hilbert space. Nicely, getting the eigenstates for the problem with only the confining potential allows to use the particle number conservation. We can generate these states using the example mentioned previously. We then need to convert them to the full basis without the particle number conservation using [FQHESphereQuasiholesWithSpinTimeReversalSymmetryConvertStates](FQHESphereQuasiholesWithSpinTimeReversalSymmetryConvertStates.md) (let say in the $K_y=0$ sector). Then we create a text file (like basis.dat) that contains the list of states that we will use for the effective Hilbert space. This file is a single line that should look like

```text
Basis= fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_fixedn_0_n_0_2s_21_sz_0_lz_0.0.vec fermions_cylinder_perimeter_8.000000_su2_quasiholes_symlinearmomentum_cenergy_0.000000_n0_0.000000_pairing_fixedn_2_n_0_2s_21_sz_0_lz_0.0.vec ...
```


The way the effective subspace basis is sorted is the order in which the vectors appears in this file. Thus, it is highly recommended that the effective subspace basis should be sorted per particle number sector and within each sector from the smallest to the largest energy.

If we want to diagonalize the superconducting coupling acting only on the leftmost orbital, we create a simple interaction file confining_supra_0_2s_21.dat

```text
OneBodyPotentialPairing  =  1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
```

Then we just run the following command 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing -l 21 --directory nphi_21_cylinder_perimeter_8.000000 --interaction-file confining_supra_0_2s_21.dat --interaction-name supra_0_effective_30 --use-lapack --full-diag 12000  -s 0  --use-cylinder --cylinder-perimeter 8.0  --memory 16000 -S --processors 4 --nbr-lz 1 --use-hilbert basis.dat --export-binhamiltonian fermions_cylinder_perimeter_8.000000_su2_quasiholes_supra_0_effective_30_cenergy_0.000000_n0_0.000000_pairing_n_0_2s_21_sz_0_lz_0.mat*


Note taht we have asked to export the hamiltonian restricted to the effective space in the binary file fermions_cylinder_perimeter_8.000000_su2_quasiholes_supra_0_effective_30_cenergy_0.000000_n0_0.000000_pairing_n_0_2s_21_sz_0_lz_0.mat . This file can be processed and diagonalized using the generic code [GenericHamiltonianDiagonalization](GenericHamiltonianDiagonalization.md).
