# FQHESphereBosonEntanglementEntropy — manual

Source: DiagHam wiki page `FQHESphereBosonEntanglementEntropy`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereBosonEntanglementEntropy](../FQHEOnSphere/FQHESphereBosonEntanglementEntropy.md)
FQHESphereBosonEntanglementEntropy (for bosons) and FQHESphereFermionEntanglementEntropy (for fermions) generate the reduced density matrix, including its eigenvalues and eigenstates. They also compute the entanglement entropy

## Basic usage

Let's say you have generate the Gaffnian state using the Jack generator *(wiki page Jack generator)*. So you have a bosons_haldane_gaffnian_n_14_2s_18_lz_0.0.vec and the description of the root configuration in gaffnian_n_14_2s_18.dat

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereBosonEntanglementEntropy bosons_haldane_gaffnian_n_14_2s_18_lz_0.0.vec  --haldane --reference-file gaffnian_n_14_2s_18.dat --use-lapack*

This will compute the orbital entanglement entropy of the Gaffnian state for a number of orbitals going from 1 to 10 (i.e. up to the hemisphere cut). The result is stored in a text file named bosons_partial_unnormalized_haldane_gaffnian_n_14_2s_18_lz_0.0.ent

If you want the eigenvalues of the reduced density matrix, just add the --density-matrix option with the name of file where you want to store them (here bosons_partial_unnormalized_haldane_gaffnian_n_14_2s_18_lz_0.0.full.ent)

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereBosonEntanglementEntropy bosons_haldane_gaffnian_n_14_2s_18_lz_0.0.vec  --haldane --reference-file gaffnian_n_14_2s_18.dat --density-matrix bosons_partial_unnormalized_haldane_gaffnian_n_14_2s_18_lz_0.0.full.ent --use-lapack*

The data in bosons_partial_unnormalized_haldane_gaffnian_n_14_2s_18_lz_0.0.full.ent will look like

```text
# l_a    N    Lz    lambda
1 0 0 0.38469657408456
1 1 0 0.49376474656819
1 2 0 0.12153867934773
```

The file is column formatted (using a single space as separator) :
- the first column is the number of orbitals (l_a)
- the second column the number of particles (N)
- the third column is the projection of the angular momentum (Lz)
- the forth column is the eigenvalue (lambda)

Lines in the file are sorted by increasing l_a, N, Lz and decreasing lambda. 

The bosons_partial_unnormalized_haldane_gaffnian_n_14_2s_18_lz_0.0.full.ent file can be processed using [FQHESphereEntanglementSpectrum](FQHESphereEntanglementSpectrum.md) to extract the entanglement spectrum.



You can restrict the calculation to a given range of l_a values using the --min-la and --max-la options. In the following example

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereBosonEntanglementEntropy bosons_haldane_gaffnian_n_14_2s_18_lz_0.0.vec  --haldane --reference-file gaffnian_n_14_2s_18.dat --density-matrix bosons_partial_unnormalized_haldane_gaffnian_n_14_2s_18_lz_0.0.full.ent --use-lapack --min-la 9 --max-la 9*

the calculation will be done only for l_a=9.

To generate the spectrum when A has 2 edges, add the option --stripe-subsystem while calling FQHESphereBosonEntanglementEntropy. The option --shift-la then indicates the position of the edge of A closest to the north pole. For example,

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereBosonEntanglementEntropy bosons_haldane_gaffnian_n_14_2s_18_lz_0.0.vec  --haldane --reference-file gaffnian_n_14_2s_18.dat --density-matrix bosons_partial_unnormalized_haldane_gaffnian_n_14_2s_18_lz_0.0.full.ent --use-lapack --stripe-subsystem --shift-la 5*

will generate the eigenvalue spectrum for all possible sizes l_a, for A being the subsystem between orbitals 5 and 5+l_a.

## Advanced usage

Eigenstates of the reduced can be computed using the --density-eigenstate flag. Since this can generate a huge number of vectors, you restrict to  a given Na and Lza sector with the --na-eigenstate and --lza-eigenstate options. You may also restrict the calculation to the x first largest eigenvalues using --nbr-eigenstates x . For example, 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereBosonEntanglementEntropy bosons_haldane_gaffnian_n_14_2s_18_lz_0.0.vec  --haldane --reference-file gaffnian_n_14_2s_18.dat --density-matrix bosons_partial_unnormalized_haldane_gaffnian_n_14_2s_18_lz_0.0.full.ent --use-lapack --min-la 9 --max-la 9  --na-eigenstate 7 --lza-eigenstate 0 --nbr-eigenstates 2 --density-eigenstate*

Two vectors are generated : 
bosons_sphere_density_n_14_2s_18_lz_0_la_9_na_7_lza_0.0.vec  and bosons_sphere_density_n_14_2s_18_lz_0_la_9_na_7_lza_0.1.vec

Beware that the --lza-eigenstate uses the lz counting convention from the density matrix eigenvalue file and not the convention from the entanglement spectrum. As usual, --lza-eigenstate is twice the lz value to avoid half integer number. 

You can store the full density matrix using the --full-densitymatrix option

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereBosonEntanglementEntropy bosons_haldane_gaffnian_n_14_2s_18_lz_0.0.vec  --haldane --reference-file gaffnian_n_14_2s_18.dat --density-matrix bosons_partial_unnormalized_haldane_gaffnian_n_14_2s_18_lz_0.0.full.ent --use-lapack --min-la 9 --max-la 9  --na-eigenstate 7 --lza-eigenstate 0 --nbr-eigenstates 2 --full-densitymatrix matrix.dat*

This will store the reduced density matrix in the ASCII file matrix.dat for the sector la=9, Na=7 Lza=0. Indices of the matrix are associated to n-body using the standard DiagHam convention that can be obtained using [FQHESphereShowBasis](FQHESphereShowBasis.md) with 7 bosons, 8 flux quanta (i.e. la-1) and total lz=0.
