# FQHESphereBosonsWithSpin — manual

Source: DiagHam wiki page `FQHESphereBosonsWithSpin`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereBosonsWithSpin](../FQHEOnSphere/FQHESphereBosonsWithSpin.md)
FQHESphereBosonsWithSpin is the analogue of FQHETorusBosonsWithSU2SpinTwoBodyGeneric for bosons having a SU(2) internal degree of freedom on the sphere geometry. It also has similarities with the fermionic version [FQHESphereFermionsWithSpin](FQHESphereFermionsWithSpin.md) or the SU(3)/SU(4) codes [FQHESphereBosonsWithSU3Spin](FQHESphereBosonsWithSU3Spin.md)/[FQHESphereBosonsWithSU4Spin](FQHESphereBosonsWithSU4Spin.md).

## Basis choice and compatibility

A big warning about the data generated with this code. Unless you know what you are doin, you should always include the --use-alt option. It has internally two ways to organize the many-body Hilbert space resulting in different eigenstate components. This is for historical reason. Most codes using the eigenstate either allow to switch from one convention to another (using the --use-alt) or will assume that the data is generated with the --use-alt (see, e.g., [FQHESphereWithSU2SpinRotation](FQHESphereWithSU2SpinRotation.md)).
