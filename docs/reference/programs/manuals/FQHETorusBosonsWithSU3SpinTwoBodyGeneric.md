# FQHETorusBosonsWithSU3SpinTwoBodyGeneric — manual

Source: DiagHam wiki page `FQHETorusBosonsWithSU3SpinTwoBodyGeneric`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHETorusBosonsWithSU3SpinTwoBodyGeneric](../FQHEOnTorus/FQHETorusBosonsWithSU3SpinTwoBodyGeneric.md)
[FQHETorusBosonsWithSU3SpinTwoBodyGeneric](FQHETorusBosonsWithSU3SpinTwoBodyGeneric.md) allows to consider bosons with a SU(3) degree of freedom on the torus geometry. Only the translation along the y direction is taken into account. Many options are similar to the ones of [QHEFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md). The SU(3) sectors are defined by the two quantum numbers $T_z$ and $Y$, they are related to the number of particles per component $N_1, N_2$ and $N_3$ through the relation 

$T_z=\frac{N_1-N_2}{2}$ and $Y=\frac{N_1+N_2 -2 N_3}{3}$

The SU(3) sector can be set through the --total-tz and --total-y options (for $T_z$ and $Y$) or through --nbr-n1, --nbr-n2, --nbr-n3 (for $N_1, N_2$ and $N_3$ ). Unless the sum of the values given by --nbr-n1, --nbr-n2, --nbr-n3 matches the total number of particles, --total-tz and --total-y are used by default.

The interaction file that describes the two body interaction is similar to the one used for [FQHESphereFermionsWithSpin](FQHESphereFermionsWithSpin.md). Interactions between each species are defined by Pseudopotentials11, Pseudopotentials12, Pseudopotentials13, Pseudopotentials22, Pseudopotentials23, Pseudopotentials33. If an additional generic Pseudopotentials is provided, it will replace any missing  Pseudopotentialsij term. Note that the number of pseudo-potentials per interaction does not need not match the number of flux quanta. For example, the following file


```text
Pseudopotentials = 1.0
```


is enough to define the interaction that produces the exact interaction related to the generalized Halperin state [222,111], irrespective to the number of particles or flux quanta.
