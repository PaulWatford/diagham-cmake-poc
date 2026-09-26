# FQHEMPSEMatrix — manual

Source: DiagHam wiki page `FQHEMPSEMatrix`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHEMPSEMatrix](../FQHEOnSphere/FQHEMPSEMatrix.md)
FQHEMPSEMatrix diagonalizes the MPS tranfer matrix of various models states.


As a simple example, here is how to compute both the left and right eigenstates (for the largest eigenvalue) for the Laughlin state 1/3

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHEMPSEMatrix --p-truncation 14 --trim-qsector --ematrix-memory 0 -S --use-nonrational --fixed-qsector --laughlin-index 3 --show-itertime --right-eigenstates --left-eigenstates --diagonal-block --normalize-cylinder --processors 4 --cylinder-perimeter 20.0 --boson -S*

Beware that some code changes occurred around 2021. Even for a fermionic Laughlin state, the --boson option should be used to benefit from a more optimal MPS construction.
