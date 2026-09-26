# FQHETorusFermionsWithSpinTwoBodyGeneric — manual

Source: DiagHam wiki page `FQHETorusFermionsWithSpinTwoBodyGeneric`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHETorusFermionsWithSpinTwoBodyGeneric](../FQHEOnTorus/FQHETorusFermionsWithSpinTwoBodyGeneric.md)
FQHETorusFermionsWithSpinTwoBodyGeneric allows to consider fermions with a SU(2) degree of freedom on the torus geometry. Only the translation along the y direction is taken into account. Many options are similar to the ones of [QHEFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md) or [FQHESphereFermionsWithSpin](FQHESphereFermionsWithSpin.md). The aspect ratio of the torus is by default set to 1 and can be tuned using --ratio optio. The two vectors that define the torus are orthogonal. Unless the --redundant-kymomenta option is set, the program  only computes the spectrum in momentum sector not related by the center of mass symmetry.

Similar to [FQHESphereFermionsWithSpin](FQHESphereFermionsWithSpin.md), a file describing the interaction should be provided. Contrary to [FQHESphereFermionsWithSpin](FQHESphereFermionsWithSpin.md), the number of pseudo-potentials does not have to match the number of flux quanta. For example, the following file pseudopotentials_torus_331.dat

```text
PseudopotentialsUpUp = 0 1
PseudopotentialsDownDown = 0 1
PseudopotentialsUpDown = 1
```

is a perfectly valid file to describe the model interaction that gives rise to the Halperin (331) state.

A typical usage of FQHETorusFermionsWithSpinTwoBodyGeneric is 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusFermionsWithSpinTwoBodyGeneric -p 4 -l 8 --interaction-file pseudopotentials_torus_331.dat --interaction-name 331 --use-lapack*

This will generate the spectrum for the model interaction described in pseudopotentials_torus_331.dat for 4 fermions at filling factor nu=1/2 .

Additionally, it is possible to specify different combinations of one-body potentials, using "OneBodyPotentialUpUp", "OneBodyPotentialDownDown" (useful for adding a Zeeman term), while the Sx term can be defined using "OneBodyPotentialUpDown" (this is also useful for implementing tunnelling term in a bilayer simulation). Note that each line of the one-body potential should have as many numbers as the number of flux quanta. For example, a typical Zeeman contribution on the first and last orbitals in a system with 8 flux quanta should look like

```text
 OneBodyPotentialUpUp = 1.0 0 0 0 0 0 0 1.0
 OneBodyPotentialDownDown = -1.0 0 0 0 0 0 0 -1.0
```
