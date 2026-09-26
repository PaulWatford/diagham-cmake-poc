# VectorPhaseMultiply — manual

Source: DiagHam wiki page `VectorPhaseMultiply`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [VectorPhaseMultiply](../core/VectorPhaseMultiply.md)
VectorPhaseMultiply multiplies a binary vector with a complex phase factor. The input vector can be real or complex while the output vector is always a complex vector. Typical usage is 

*\$PATHTODIAGHAM/build/src/Programs/VectorPhaseMultiply -i input.vec -o output.vec -p 3.14*

The phase is set through the -p option and should be expressed in radian or in unit of pi if the --pi-unit option is set.

VectorPhaseMultiply can be used to convert real vector into complex phase by setting the phase to zero (default value).
