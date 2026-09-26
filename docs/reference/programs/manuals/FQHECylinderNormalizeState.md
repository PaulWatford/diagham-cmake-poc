# FQHECylinderNormalizeState — manual

Source: DiagHam wiki page `FQHECylinderNormalizeState`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHECylinderNormalizeState](../FQHEOnCylinder/FQHECylinderNormalizeState.md)
[FQHECylinderNormalizeState](FQHECylinderNormalizeState.md) converts a state from the unnormalized basis to the cylinder geometry. It is the analogue of [FQHESphereUnnormalizeState](FQHESphereUnnormalizeState.md) for the cylinder geometry. Its typical usage is 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderNormalizeState -i fermions_unnormalized_haldane_n_8_2s_21_lz_0.0.vec --haldane --reference-file laughlin3_n_8_2s_21.dat --cylinder-perimeter 9.0 --normalize --output-file fermions_haldane_cylinder_perimeter_9.000000_n_8_2s_21_lz_0.0.vec*
