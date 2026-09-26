# VectorBinary2Ascii — manual

Source: DiagHam wiki page `VectorBinary2Ascii`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [VectorBinary2Ascii](../core/VectorBinary2Ascii.md)
VectorBinary2Ascii converts a binary vector into a text file


*\$PATHTODIAGHAM/build/src/Programs/VectorBinary2Ascii -i vector.vec -o vector.vec.txt*

If you just want to see all components of vector.vec without saving them into a text file, just type

*\$PATHTODIAGHAM/build/src/Programs/VectorBinary2Ascii -i vector.vec -s*

The --add-index option adds an index in front of each component

*\$PATHTODIAGHAM/build/src/Programs/VectorBinary2Ascii -i vector.vec --add-index -s*

The indices use the C convention and thus are zero based. 

It is possible to only convert a part of the binary vector. The two options --min-range and --max-range  allows to define the range of components that have to be converted

*\$PATHTODIAGHAM/build/src/Programs/VectorBinary2Ascii -i vector.vec -s --min-range 12 --max-range 20*

This will display components 12 to 20 of vector.vec.
