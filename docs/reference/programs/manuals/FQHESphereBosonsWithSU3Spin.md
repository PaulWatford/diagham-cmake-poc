# FQHESphereBosonsWithSU3Spin — manual

Source: DiagHam wiki page `FQHESphereBosonsWithSU3Spin`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereBosonsWithSU3Spin](../FQHEOnSphere/FQHESphereBosonsWithSU3Spin.md)
[FQHESphereBosonsWithSU3Spin](FQHESphereBosonsWithSU3Spin.md) is the analogue of [FQHESphereBosonsWithSU4Spin](FQHESphereBosonsWithSU4Spin.md) for bosons having a SU(4) internal degree of freedom on the sphere geometry. It shares the options of [FQHESphereBosonsWithSU4Spin](FQHESphereBosonsWithSU4Spin.md) except for the SU(3) specific options. Similar to [FQHETorusBosonsWithSU3SpinTwoBodyGeneric](FQHETorusBosonsWithSU3SpinTwoBodyGeneric.md), the SU(3) sectors are defined by the two quantum numbers $T_z$ and $Y$, they are related to the number of particles per component $N_1, N_2$ and $N_3$ through the relation 

$T_z=\frac{N_1-N_2}{2}$ and $Y=\frac{N_1+N_2 -2 N_3}{3}$

The SU(3) sector can be set through the --total-tz and --total-y options (for $T_z$ and $Y$) or through --nbr-n1, --nbr-n2, --nbr-n3 (for $N_1, N_2$ and $N_3$ ). Unless the sum of the values given by --nbr-n1, --nbr-n2, --nbr-n3 matches the total number of particles, --total-tz and --total-y are used by default.

The interaction file that describes the two body interaction is identical to the one used for [FQHETorusBosonsWithSU3SpinTwoBodyGeneric](FQHETorusBosonsWithSU3SpinTwoBodyGeneric.md).
