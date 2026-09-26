# FQHESphereCP2GenericInteraction — manual

Source: DiagHam wiki page `FQHESphereCP2GenericInteraction`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereCP2GenericInteraction](../FQHEOnSphere/FQHESphereCP2GenericInteraction.md)
## Basic usage

[FQHESphereCP2GenericInteraction](FQHESphereCP2GenericInteraction.md) allows one to diagonalize any two-body interaction in the lowest Landau level of the CP2 geometry. 

The form of the interaction is defined in an ASCII file provided in the --interaction-file option. It specifies the value of all two-body pseudopotentials. It should contain one line of the form:
```text
Pseudopotentials = 0 1 
```

One also has to provide the list of all Clebsch-Gordan coefficients necessary for the composition of SU(3) momenta. A particle in the lowest Landau level of the CP2 geometry is in a <MATH>(N_{\phi},0)</MATH> <MATH>SU(3)</MATH> representation, where <MATH>N_{\phi}</MATH> is the number of flux quanta specified by the -l option. The composition of two such representations writes:

<center>
<MATH>(N_{\phi},0) \otimes (N_{\phi},0) = \oplus_{p' = 0}^{N_{\phi}} \left( 2(N_{\phi} - p'), p' \right)</MATH>
</center>
One thus has to provide <MATH>N_{\phi} + 1</MATH> files containing the Clebsch-Gordan coefficients for each <MATH>( 2(N_{\phi} - p'), p'), p' = 0, ..., N_{\phi}</MATH> resulting representation. The name of these files specifies the original single-particle SU(3) representations <MATH>(p_1, q_1) = (N_{\phi}, 0)</MATH>, <MATH>(p_2, q_2) = (N_{\phi}, 0)</MATH>, and the resulting representation <MATH>(p, q)</MATH>. For example, for <MATH>N_{\phi} = 3</MATH>, the following files should be provided:
```text
CGCSU3_p1_3_q1_0_p2_3_q2_0_p_0_q_3.dat
CGCSU3_p1_3_q1_0_p2_3_q2_0_p_2_q_2.dat
CGCSU3_p1_3_q1_0_p2_3_q2_0_p_4_q_1.dat
CGCSU3_p1_3_q1_0_p2_3_q2_0_p_6_q_0.dat
```

Each of these files is a four-column file, taking the form
```text
#Q(M) Q(M1) Q(M2) CGC
1 1 7 -0.500000000000000
1 2 4 0.500000000000000
1 4 2 -0.500000000000000
1 7 1 0.500000000000000
2 1 8 -0.500000000000000
2 2 5 0.408248290463863
```

All Clebsch-Gordan files necessary for calculations up to <MATH>N_{\phi} = 10</MATH> are given [here](http://www.nick-ux.org/~repellin/cp2_data/ClebschGordanCoefficients.tgz).

The generic use of the code is:


*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereCP2GenericInteraction -p 2 -l 3 --interaction-name V1 --interaction-file pseudopotential_v1.dat --use-lapack*

This will diagonalize the interaction in each Tz and Y sector. Option -p sets the number of particles and -l the one-body maximum angular momentum (which matches the number of flux quanta in the lowest Landau level). Running this command produces one file named fermions_cp2_V1_n_2_2s_3.dat that looks like :

```text
# tz y E
1 -9 1
0 -6 1
0 -6 1
2 -6 1
1 -3 2.2204460492503e-16
1 -3 1
1 -3 1
0 0 -6.1309083055668e-31
0 0 1
0 0 1
3 -3 1.1102230246252e-16
3 -3 1
...
```

This is a three-column ASCII file. The first column is the tz sector, the second column is the y sector and the third column is the energies sorted from the lowest to the highest. The name of the interaction appearing in the output file is set using the --interaction-name option.

By default, the statistics is fermionic. Bosonic statistics can be activated thanks to the option --boson.


## Discrete symmetries

Two discrete symmetries are implemented for bosons.

The discrete symmetry <MATH>T_z \rightarrow -T_z </MATH> can be activated by using the option --tzsymmetrized-basis. The parity sector can be chosen by using the option --minus-tzparity

There is an additional Z3 symmetry, which commutes with the <MATH>T_z \rightarrow -T_z </MATH> symmetry. Activating both of these symmetries requires to use the option --tzZ3symmetrized-basis

If the eigenvectors are obtained using these symmetries, one can convert them to the full basis by using [FQHESphereCP2ConvertSymmetrizedState](FQHESphereCP2ConvertSymmetrizedState.md).
