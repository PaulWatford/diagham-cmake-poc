# FQHECylinderQuasiholeMatrixElements — manual

Source: DiagHam wiki page `FQHECylinderQuasiholeMatrixElements`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHECylinderQuasiholeMatrixElements](../FQHEOnCylinder/FQHECylinderQuasiholeMatrixElements.md)

> **Options in this manual that the current program does not list**: `--ratio`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
[FQHECylinderQuasiholeMatrixElements](FQHECylinderQuasiholeMatrixElements.md) is the analogue of [FQHESphereQuasiholeMatrixElements](FQHESphereQuasiholeMatrixElements.md) for the cylinder geometry. Its usage and outputs are similar to its sphere counterpart except that it supports ths SMP/parallelization mode. 

There are a few specific options. The cylinder geometry can be specified either by its  perimeter $L_y$ (using the --cylinder-perimeter option) or its aspect ratio $L_y/L_x$ (using the --ratio option). Unless --cylinder-perimeter is set to a positive number, the aspect ratio is used as the default way to set the geometry.

There is also the possibility to focus on a specific number of particles using the --fix-nbrparticles option. In that case, only the matrix elements involving this specific number of particles and the vacuum are computed. For example,

``\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderQuasiholeMatrixElements -l 26 --fix-nbrparticles 7 --cylinder-perimeter 8.0''

generates all the matrix elements involving 7 particles on cylinder with 26 flux quanta and a perimeter of  $L_y=8$. The data are stored in a folder named nphi_26_n_7_cylinder_perimeter_8.000000.
