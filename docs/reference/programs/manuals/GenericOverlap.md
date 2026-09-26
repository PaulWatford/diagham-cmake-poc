# GenericOverlap — manual

Source: DiagHam wiki page `GenericOverlap`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [GenericOverlap](../core/GenericOverlap.md)
GenericOverlap allows to compute overlaps (i.e. scalar products) between binary vectors. The typical usage is

*\$PATHTODIAGHAM/build/src/Programs/GenericOverlap vector1.vec vector2.vec*

If the vectors are complex, one should add the -c option. One can also compute overlaps for a group of vectors

*\$PATHTODIAGHAM/build/src/Programs/GenericOverlap vector1.vec vector2.vec vector3.vec vector4.vec*

or equivalently, if only vector1.vec vector2.vec vector3.vec vector4.vec are in the current directory

*\$PATHTODIAGHAM/build/src/Programs/GenericOverlap vector*.vec*
