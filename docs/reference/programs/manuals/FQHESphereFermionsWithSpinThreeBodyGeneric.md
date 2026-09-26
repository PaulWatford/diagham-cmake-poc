# FQHESphereFermionsWithSpinThreeBodyGeneric — manual

Source: DiagHam wiki page `FQHESphereFermionsWithSpinThreeBodyGeneric`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereFermionsWithSpinThreeBodyGeneric](../FQHEOnSphere/FQHESphereFermionsWithSpinThreeBodyGeneric.md)
FQHESphereFermionsWithSpinThreeBodyGeneric handles SU(2) symmetric three-body interactions on the sphere geometry. A typical pseudo-potential file should look like

```text
ThreebodyPseudopotentials32 = 0 0 0 1
ThreebodyPseudopotentials12 = 0 1 1
```

where ThreebodyPseudopotentials32 gives the pseudo-potential in the spin 3/2 sector (starting from relative angular momentum 0) and ThreebodyPseudopotentials12 gives the pseudo-potential in the spin 1/2 sector. For example, this  pseudo-potential file reproduces the model interaction for the Spin-Charge separated state of [arXiv:1107.2232](http://arxiv.org/abs/1107.2232) 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereFermionsWithSpinThreeBodyGeneric -p 8 -l 9 --use-lapack --interaction-name spinchargeseparated --interaction-file pseudopotentials_spinchargeseparated_fermions.dat --full-diag 4000*

The code can handle Hilbert spaces with a fixed value of S_z (*-s*) or without Sz conservation (*--all-sz* option). 

Note that, similar to [FQHESphereFermionsThreeBodyGeneric](FQHESphereFermionsThreeBodyGeneric.md), this code does *not* include an overall normalization which would guarantee that 3 particles only have energies 0 and 1 in the appropriate momentum sector. 

## Including 2-body interactions

In addition to the 3-body pseudopotentials described above, the code can also include 2-body interaction. This can be defined in two ways. The first option is the standard definition of spinful pseudopotentials (given by PseudopotentialsUpUp, PseudopotentialsUpDown, PseudopotentialsDownDown arrays). The second option, relevant for cases without Sz conservation, is to use arbitrary 2-body interactions similar to FQHESphereFermionsWithSpinFull *(wiki page FQHESphereFermionsWithSpinFull)*. In this case, we parametrize the interactions using PseudopotentialsUpUpUpUp, PseudopotentialsDownDownDownDown,PseudopotentialsUpUpDownDown, PseudopotentialsDownDownUpUp, PseudopotentialsUpDownUpDown (this covers a wide range of interactions that still preserve inversion symmetry). If inversion is broken, one could also add terms such as PseudopotentialsDownDownDownUp (but note the code has not been tested in such cases). 

**Note that if you wish to use the second type of 2-body interaction, you must include a command line option *--twobody-full*.**

To illustrate a possible use of these pseudopotentials, a bilayer interaction given in terms of pseudopotentials V^{inter} and V^{intra} in the usual z-basis can be rotated to the x-basis according to the following linear transformation (note the important 1/2 factors in certain terms):

PseudopotentialsUpUpUpUp, PseudopotentialsDownDownDownDown = V^{intra}/2 + V^{inter}/2, 
PseudopotentialsUpUpDownDown, PseudopotentialsDownDownUpUp = V^{intra}/2 - V^{inter}/2, 
PseudopotentialsUpDownUpDown= V^{intra} + V^{inter}.

This will produce the same spectrum as the usual bilayer interaction (V^{intra}, V^{inter}), but the eigenstates will be rotated from z- to x-basis. 

Finally, a 1-body potential can also be included in the pseudopotential file, using OneBodyPotentialUpUp, OneBodyPotentialDownDown, OneBodyPotentialUpDown (these are defined in the standard z-basis).
