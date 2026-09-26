# CoulombPseudopotentials — manual

Source: DiagHam wiki page `CoulombPseudopotentials`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [CoulombPseudopotentials](../FQHEOnSphere/CoulombPseudopotentials.md)
CoulombPseudopotentials generates the pseudopotentials on the sphere for the Coulomb interaction. For example,

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/CoulombPseudopotentials -s 9 -l 0*

will generate the Coulomb pseudopotentials in the lowest Landau level for 9 flux quanta. The result is stored in pseudopotential_coulomb_l_0_2s_9.dat and looks like 


```text
# pseudopotentials on the sphere for coulomb interaction
# in the Landau level N=0 for 2S=9 flux quanta
#
# Pseudopotentials = V_0 V_1 ...

Pseudopotentials = 0.97624284500207 0.50131389337944 0.386727860607 0.33203907223833 0.29990625879591 0.27922306853413 0.26543427502627 0.25633367131108 0.25076120019562 0.2481076425216
```


the Landau level can be set through the -l option, 0 being the lowest Landau level. The number of flux quanta is set using the -s option. The files that are produced with CoulombPseudopotentials can be directly provided to programs such as [FQHESphereFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md), [FQHESphereBosonsTwoBodyGeneric](FQHESphereBosonsTwoBodyGeneric.md) ...

The Coulomb interaction between two layers separated from a distance d (or the ZDS interaction) can be computed using the -d or --layer-separation options i.e.

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/CoulombPseudopotentials -s 9-l 0 -d 1.4*

## One-body potentials

In addition to two-body pseudopotentials, there is an option to include one-body potentials in the same pseudopotential file (e.g., basic codes such as FQHEFermionsOnSphereTwoBodyGeneric support the one-body potential of the form \sum_m U_m c_m^+ c_m). A typical example when this is needed is if we wish to include an impurity, usually positioned at one of the poles of the sphere to preserve Lz symmetry. CoulombPseudopotentials has the following options to include a delta potential impurity at one or both poles ("--north-potential" or "--south-potential"). This will produce a line in the pseudopotential file: "Onebodypotentials = U_0 U_1 U_2 ... U_{2Q}".
