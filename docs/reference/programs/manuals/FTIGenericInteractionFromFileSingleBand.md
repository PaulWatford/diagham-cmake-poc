# FTIGenericInteractionFromFileSingleBand — manual

Source: DiagHam wiki page `FTIGenericInteractionFromFileSingleBand`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FTIGenericInteractionFromFileSingleBand](../FTI/FTIGenericInteractionFromFileSingleBand.md)

> **Options in this manual that the current program does not list**: `--conserve-bandoccuption`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
[FTIGenericInteractionFromFileSingleBand](FTIGenericInteractionFromFileSingleBand.md) is similar to [FTIGenericInteractionFromFileTwoBands](FTIGenericInteractionFromFileTwoBands.md) (or [FTIGenericInteractionFromFileThreeBands](FTIGenericInteractionFromFileThreeBands.md)) for a single band. These two programs share most of the same option. While [FTIGenericInteractionFromFileTwoBands](FTIGenericInteractionFromFileTwoBands.md) can technically handle a single band problem by using the option --conserve-bandoccuption and choosing the proper quantum number sectors, [FTIGenericInteractionFromFileSingleBand](FTIGenericInteractionFromFileSingleBand.md) supports more momentum sectors per band and is usually faster for the hilbert space and hamiltonian generations. 

The typical usage is for the simplest case (no additional valley or spin quantum numbers) is

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIGenericInteractionFromFileSingleBand --use-lapack --interaction-file Ne12_matrix_element.csv --singleparticle-file Ne12_one_body_term.csv --force-reorthogonalize --block-lanczos --block-size 5 -n 5 -x 6 -y 6 -p 12 --memory 1500000 --interaction-name pentalayer_without_q0 -S --processors 32 --eigenstate --lanczos-precision 1e-12 --show-itertime --selected-sectors sectors.dat*

The options are similar to [FTIGenericInteractionFromFileTwoBands](FTIGenericInteractionFromFileTwoBands.md). The only difference consists in the format for the different input files.

## Format of the one-body and two-body terms:

At its core, [FTIGenericInteractionFromFileSingleBand](FTIGenericInteractionFromFileSingleBand.md) diagonalizes the following hamiltonian

$\sum_k M_k c^\dagger_k c_k + \sum{k_1, k_2, k_3, k_4} U_{k_1, k_2, k_3, k_4} c^\dagger_{k_1} c^\dagger_{k_2} c_{k_3} c_{k_4}$

where $M$ and $U$ are provided through ascii text files and the options --interaction-file (for $U$) and --singleparticle-file (for $M$). Here we did not include any additional degrees of freedom. In the above command line example, the file content for Ne12_matrix_element.csv look like 

```text
# kx_1 ky_1 kx_2 ky_2 kx_3 ky_3 kx_4 ky_4 matrix_element
0 0 0 0 0 0 0 0 (0.0000000000000000,0.0000000000000000)
0 1 0 5 0 0 0 0 (-0.0024110134775047,-0.0006607366248230)
0 2 0 4 0 0 0 0 (0.0004708271743521,-0.0009749655343664)
0 3 0 3 0 0 0 0 (0.0000401621379381,0.0005855673161462)
0 4 0 2 0 0 0 0 (0.0004708271743583,-0.0009749655343856)
...
```


Complex entries are given as (Re,Im) and should include no space (space is used a separator). Note that apart from momentum conservation and hermiticity, no other assumptions should be made about the two-body matrix elements. For example, no ordering should be assumed at this stage, the code will take care automatically of such optimization.

For the one-body matrix elements, the file of Ne12_one_body_term.csv is simpler 

```text
# kx ky one_body_term
0 0 0.0465304246569367
0 1 0.0435119280864994
0 2 0.0335189166933237
0 3 0.0272147438234846
0 4 0.0338461361172866
0 5 0.0428872838337418
1 0 0.0435102234072057
...
```

Note that due to hermiticity, entries are real and thus have a simpler syntax.


## Sector selection:

By default, [FTIGenericInteractionFromFileSingleBand](FTIGenericInteractionFromFileSingleBand.md) will compute all the momentum sectors. While some options can be used to specify some sectors through the command line, the most efficient method is to use an ascii input file and the option --selected-sectors. For a single band without any additional degree of freedom, the sector selection input file (*sectors.dat* in the command line example) should look like

```text
# kx ky
5 0
5 1
5 2
5 3
5 4
5 5
```

It is possible to change the number of eigenvalues/eigenvectors computed/requested per quantum number sector. To do so, add a column containing the number of requested eigenvalues in the and include the option *--variable-nbreigenvalues* to the command line. Note that when using block Lanczos, you should set the block size to the maximal number of eigenvalues to compute. Still, the code will make sure the block size does not exceed the number of eigenvalues to compute. For example, in the case of the spinless/valleyless situation, the sector selection input file should include a third column

```text
# kx ky deg
5 0 5
5 1 1
5 2 1
5 3 1
5 4 1
5 5 3
```
