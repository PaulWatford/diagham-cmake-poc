# FQHESphereGetDimension — manual

Source: DiagHam wiki page `FQHESphereGetDimension`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereGetDimension](../FQHEOnSphere/FQHESphereGetDimension.md)
FQHESphereGetDimension computes the Hilbert space dimension for particles on the sphere geometry. For example,

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereGetDimension -n 12 -s 21  --ground-only *

will display the size of the Hilbert space for 12 fermions and 21 flux quanta in the Lz=0 sector. The --boson option has to be used to switch to bosons

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereGetDimension -n 12 -s 10  --ground-only --boson*

FQHESphereGetDimension can also compute Hilbert space dimensions with an internal degree of freedom (like spin) or with multiple Landau levels. This can be done using the following options :

- --su2-spin : consider particles with SU(2) spin
- --su3-spin : consider particles with SU(3) spin
- --su2su2-spin : consider particles with SU(2)xSU(2) spin
- --su4-spin : consider particles with SU(4) spin
- --2-ll : consider particles within two Landau levels
- --3-ll : consider particles within three Landau levels
