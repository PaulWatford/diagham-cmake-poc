# FQHESphereFermionsWithSpin — manual

Source: DiagHam wiki page `FQHESphereFermionsWithSpin`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereFermionsWithSpin](../FQHEOnSphere/FQHESphereFermionsWithSpin.md)
FQHESphereFermionsWithSpin provides exact diagonalization of any two-body hamiltonian for fermions with SU(2) on the sphere geometry. The way to use this code close to the one of [FQHESphereFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md). FQHESphereFermionsWithSpin works in a given sector of total spin projection along z (Sz) that can be set using the -s or --total-sz option. The code also supports the option to work in the space containing all possible Sz components, this can be activated with the option --all-sz.


The interaction is specified through its Haldane pseudo-potentials. Such a file can be generated using [CoulombPseudopotentials](CoulombPseudopotentials.md) for SU(2) invariant interaction. In the most generic case, it looks like 

```text
PseudopotentialsUpUp = 1.0 1.0 1.0 0.0 0.0 0.0
PseudopotentialsDownDown = 1.0 1.0 1.0 0.0 0.0 0.0
PseudopotentialsUpDown = 1.0 0.0 0.0 0.0 0.0 0.0
```

This example, for the Halperin (331) state, shows the naming convention to provide the pseudo-potentials between two spin up, two spin down and spin up / spin down.

Additionally, it is possible to specify different combinations of one-body potentials, using "OneBodyPotentialUpUp", "OneBodyPotentialDownDown" (useful for adding a Zeeman term), while the Sx term can be defined using "OneBodyPotentialUpDown" (this is also useful for implementing tunnelling term in a bilayer simulation).
