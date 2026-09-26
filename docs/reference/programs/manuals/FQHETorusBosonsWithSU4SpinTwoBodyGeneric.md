# FQHETorusBosonsWithSU4SpinTwoBodyGeneric — manual

Source: DiagHam wiki page `FQHETorusBosonsWithSU4SpinTwoBodyGeneric`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHETorusBosonsWithSU4SpinTwoBodyGeneric](../FQHEOnTorus/FQHETorusBosonsWithSU4SpinTwoBodyGeneric.md)
FQHETorusBosonsWithSU4SpinTwoBodyGeneric is the analogue of [FQHETorusBosonsWithSU3SpinTwoBodyGeneric](FQHETorusBosonsWithSU3SpinTwoBodyGeneric.md) for bosons having a SU(4) internal degree of freedom on the torus geometry. Only the translation along the y direction is taken into account. Many options are thus similar to the ones of [FQHETorusBosonsWithSU3SpinTwoBodyGeneric](FQHETorusBosonsWithSU3SpinTwoBodyGeneric.md). The SU(4) sectors are defined by the three quantum numbers $S_z$ (spin), $I_z$ (isospin) and $P_z$, they are related to the number of particles per component $N_1, N_2, N_3$ and $N_4$ through the relation 

$S_z=\frac{N_1+N_2-N_3-N_4}{2}$, $I_z=\frac{N_1+N_3 -N_2-N_4}{2}$ and $P_z=\frac{N_1+N_4 -N_2-N_3}{2}$ 

The SU(4) sector can be set through the --total-sz (for $S_z$), --total-isosz (for $I_z$) and --total-entanglement options (for $P_z$) or through --nbr-n1, --nbr-n2, --nbr-n3, --nbr-n4 (for $N_1, N_2, N_3$ and $N_4$ ). Unless the sum of the values given by --nbr-n1, --nbr-n2, --nbr-n3, -nbr-n4  matches the total number of particles, --total-sz, --total-isosz and --total-entanglement are used by default.

The interaction file that describes the two body interaction is similar to the one used for [FQHESphereFermionsWithSpin](FQHESphereFermionsWithSpin.md) or [FQHETorusBosonsWithSU3SpinTwoBodyGeneric](FQHETorusBosonsWithSU3SpinTwoBodyGeneric.md). Interactions between each species are defined by ten pseudo-potentials :

- PseudopotentialsUpPlusUpPlus
- PseudopotentialsUpPlusUpMinus
- PseudopotentialsUpPlusDownPlus
- PseudopotentialsUpPlusDownMinus
- PseudopotentialsUpMinusUpMinus
- PseudopotentialsUpMinusDownPlus
- PseudopotentialsUpMinusDownMinus
- PseudopotentialsDownPlusDownPlus
- PseudopotentialsDownPlusDownMinus
- PseudopotentialsDownMinusDownMinus

Here UpPlus stands for the type 1 particles, UpMinus stands for the type 2 particles, DownPlus stands for the type 3 particles and DownMinus stands for the type 4 particles. If an additional generic Pseudopotentials is provided, it will replace any missing pseudo-potential term. Note that the number of pseudo-potentials per interaction does not need not match the number of flux quanta. For example, the following file


```text
Pseudopotentials = 1.0
```


is enough to define the interaction that produces the exact interaction related to the generalized Halperin state [2,1], irrespective to the number of particles or flux quanta.
