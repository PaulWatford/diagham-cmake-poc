# BuildSuperPosition — manual

Source: DiagHam wiki page `BuildSuperPosition`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [BuildSuperPosition](../core/BuildSuperPosition.md)
BuildSuperPosition builds a linear superposition of binary vectors. Typical usage looks like

*\$PATHTODIAGHAM/build/src/Programs/BuildSuperPosition -f list.dat -o output.vec*

list.dat is an ASCII column formatted text, the first column being the binary file names and the second column being the multiplicative coefficients. The -o option sets th output file name. By default, the output vector is normalized to one (this can be turned off using the --no-normalize option). BuildSuperPosition can work with complex vector if the --complex option is set. In that case, if any third column in list.dat will be used to define the imaginary part of the multiplicative coefficients.
