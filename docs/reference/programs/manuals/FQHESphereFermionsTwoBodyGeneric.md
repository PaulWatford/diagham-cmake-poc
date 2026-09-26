# FQHESphereFermionsTwoBodyGeneric — manual

Source: DiagHam wiki page `FQHESphereFermionsTwoBodyGeneric`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereFermionsTwoBodyGeneric](../FQHEOnSphere/FQHESphereFermionsTwoBodyGeneric.md)
FQHESphereFermionsTwoBodyGeneric (formerly [QHEFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md)) provides exact diagonalization of any two-body hamiltonian for fermions on the sphere geometry. The interaction is provided through a file similar to those generated with the [CoulombPseudopotentials](CoulombPseudopotentials.md) code, that gives the Haldane pseudo-potentials. 


## Basic usage

As a first example, let's generate the complete spectrum for N=4 fermions with 9 flux quanta for the Coulomb interaction in the lowest Landau level. We assume the pseudo-potentials have been generated like described in [CoulombPseudopotentials](CoulombPseudopotentials.md). Then we run the command 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereFermionsTwoBodyGeneric -p 4 -l 9 --interaction-name coulomb_l_0 --interaction-file pseudopotential_coulomb_l_0_2s_9.dat --use-lapack*

This will diagonalize the coulomb interaction in each Lz sectors. Option -p sets the number of particles and -l the one-body maximum angular momentum (which matches the number of flux quanta in the lowest Landau level). The --use-lapack is not mandatory but improves both the speed and accuracy of the calculation. The name of the file that contains the pseudo-potentials is provided using the --interaction-file option. Running this command produces one file named fermions_coulomb_l_0_n_4_2s_9_lz.dat that looks like :

```text
# Lz E
0 3.7422768242912
0 3.9293432402432
0 3.9745134584252
0 4.0059903341452
0 4.1522333057227
```

This a two-column ASCII file. The first column is the Lz sector (rounded to the closest lower integer) and the second column is the energies sorted from the lowest to the highest. The name of the interaction appearing in the output file is set using the --interaction-name option. The name of the interaction can also be defined in the pseudo-potential file adding a line such as 

```text
Name = coulomb_l_0
```

In virtually all DiagHam's routines, the value of Lz is encoded as 2*Lz because the latter is always an integer. One can diagonalise many Lz sectors at once by using the options "--initial-lz <Lz0> --nbr-lz <n>". This will diagonalize Lz0 sector, then Lz0 + 2, then Lz0 + 4, ..., Lz0 + 2*n. It will store all the eigenvalues in a single file and the eigenvectors will be labelled by "...*lz_<value of Lz>.<index of vector>.vec" (see below for further information regarding eigenvectors). Quite often, one only needs the ground state, which means the calculation can be restricted to a single sector Lz=0 that we can specify using "--nbr-lz 1" (we can drop "--initial-lz 0" option as this is assumed by default). 

## Convention regarding 1/2 factor in energies

Note that in second-quantized notation, the Hamiltonian is usually defined with a prefactor 1/2, i.e., H = (1/2) \sum_{m_1,m_2,m_3,m_4} V_{m1,..,m4} c_m1^+ c_m2^+ c_m3 c_m4, but this prefactor is absent in DiagHam for historical reasons. Therefore, if the Hamiltonian only contains two-body interactions, the energies should be divided by 2 to get the correct values at the end of the calculation. However, if the Hamiltonian involves a mix of 2-body and 1-body terms, then one should take extra care (e.g., divide the Haldane pseudopotentials by 2, but not modify the 1-body terms). As a reminder, a one-body potential can be added by including a Onebodypotentials line in the interaction file (similar to what can be produced with [CoulombPseudopotentials](CoulombPseudopotentials.md)).



## Eigenstates

By default, [FQHESphereFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md) only computes the energy spectrum. One can extract the eigenstates turning on the --eigenstate option. The number of eigenvectors that will be computed is set thruugh the -n (or --nbr-eigen) option. In full diagonalization mode, this option sets the number of eigenvectors to evalute, while in Lanczos mode it  sets both the number of eigenvalues and eigenvectors. Vectors are saved in binary vector *(wiki page binary vectors)* format with names such as **fermions_coulomb_l_0_n_4_2s_9_lz_x.y.dat** where **x** stands for twice the Lz value and **y** is the eigenstate index in this Lz sector (0 being the lower energy state in this sector).

## Partial or full diagonalization

[FQHESphereFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md) completely diagonalizes the Hamiltonian until the Hilbert space dimension reaches the limit set through the --full-diag option (default value is set to 500). Above this limit, the code will switch to Lanczos algorithm and only compute the lowest part of the spectrum.

## Diagonalizing in a given subspace

In principle, [FQHESphereFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md)  uses the full Hilbert space constraint by the number of particles, flux quanta and total $L_z$. If you want to work in a given subspace that is defined by a set of binary vectors, you can use the --use-hilbert option. This option requires a text file that looks like

```text
Basis = vector1.vec vector2.vec vector3.vec
```

where *vector1.vec vector2.vec vector3.vec* is the list of binary vectors that define the subspace. Notice that this option also provide an easy way to compute the energy expectation value of a single state. In that case, one needs to only have a single state in the text file provided to the --use-hilbert option.
