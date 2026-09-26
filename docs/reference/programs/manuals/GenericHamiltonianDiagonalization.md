# GenericHamiltonianDiagonalization — manual

Source: DiagHam wiki page `GenericHamiltonianDiagonalization`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [GenericHamiltonianDiagonalization](../core/GenericHamiltonianDiagonalization.md)
[GenericHamiltonianDiagonalization](GenericHamiltonianDiagonalization.md) is a simple program that diagonalizes a matrix (or a sum of matrice) given as an input file.

## Text file matrices

## Binary  matrices

The basic usage for an hamiltonian stored as a binary file is

*\$PATHTODIAGHAM/build/src/Programs/GenericHamiltonianDiagonalization --use-lapack -o spectrum.dat --bin-hamiltonian hamiltonian.mat*

where hamiltonian.mat is the matrix binary file and spectrum.dat is the output file for the spectrum. The hamiltonian can be rescaled using the --rescaling-factor option. The diagonal matrix elements can be modified by providing an external column based ASCII text file:

*\$PATHTODIAGHAM/build/src/Programs/GenericHamiltonianDiagonalization --use-lapack -o spectrum.dat --bin-hamiltonian hamiltonian.mat --rescaling-factor 2.0 --add-diagonal diagonal.dat --column-diagonal 2 *

In this example, diagonal.dat is the text file containing the contribution that will be added to the diagonal entries of hamiltonian.mat. Note that the rescaling is performed *before* changing the diagonal entries. The column of diagonal.dat that should be considered is set by --column-diagonal (0 being the first column, the default value of --column-diagonal).

A linear combination of matrices can be used as the input hamiltonian. For example

*\$PATHTODIAGHAM/build/src/Programs/GenericHamiltonianDiagonalization --use-lapack -o test2 --multiple-binhamiltonians hamiltonians.dat --column-diagonal 3 --add-diagonal diagonal.dat  --rescaling-factor 2.0*

will first read a list of matrices from the column formatted text file hamiltonians.dat, add them (with the optional weight), rescale by --rescaling-factor and then add diagonal.dat to the diagonal entries. hamiltonians.dat should look like

```text
hamiltonian1.mat
hamiltonian2.mat
hamiltonian3.mat
```


If the second column is present, it will be used as weights i.e.

```text
hamiltonian1.mat 0.25
hamiltonian2.mat 0.75
```

will correspond to a hamiltonian (hamiltonian1.mat * 0.25) + (hamiltonian2.mat * 0.75).

If the third column is present, it will be used as phases (in unit of $2 \pi$) to transform the real symmetric matrices in the first column into hermitian matrices. For example:

```text
hamiltonian1.mat 0.25 0.0
hamiltonian2.mat 0.75 0.1
```

will correspond to a hamiltonian (hamiltonian1.mat * 0.25) + (hamiltonian2.mat.herm * 0.75), where hamiltonian2.mat.herm is the matrix formed by multiplying all of the upper triangular elements of hamiltonian2.mat by a factor $e^{2 \pi i * 0.1}$ (and the lower triangular elements by a factor $e^{-2 \pi i * 0.1}$).


The eigenstates can be obtained using the standard options of DiagHam (such as --eigenstate) except that a prefix for the eigenstate file name should be provided using the --eigenstate-file option.
