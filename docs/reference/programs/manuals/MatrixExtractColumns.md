# MatrixExtractColumns — manual

Source: DiagHam wiki page `MatrixExtractColumns`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [MatrixExtractColumns](../core/MatrixExtractColumns.md)
MatrixExtractColumns allows to extract a series of columns from a binary matrix and store them as binary vectors. A typical usage is

*\$PATHTODIAGHAM/build/src/Programs/MatrixExtractColumns --first-column 1 --nbr-columns 2 -i spin_1_periodicaklt_n_6_sz_0_invsym_1_szsym_1_k_0.eigenvec.mat -o vector*

In this example, MatrixExtractColumns extracts the second and the third column of and store them in vector.1.vec and vector.2.vec. The index of the first column is set by the --first-column option (the first column having an index of zero) and the number of column to extract is fixed by --nbr-columns. If the matrix is complex, you should add the -c option.
