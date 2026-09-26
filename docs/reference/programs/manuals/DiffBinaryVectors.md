# DiffBinaryVectors — manual

Source: DiagHam wiki page `DiffBinaryVectors`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [DiffBinaryVectors](../core/DiffBinaryVectors.md)
DiffBinaryVectors compares two binary vectors. The ypical usage is

*\$PATHTODIAGHAM/build/src/Program/DiffBinaryVectors vector1.vec vector2.vec*

If you don't want to compare components that are strictly zero in the vector1.vec to those of vector2.vec, just add the --discard-zero option. The error threshold is set through the --error option. If it is set to 0, it means that strict equality between two floating point numbers is required.
