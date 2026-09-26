# FQHESphereFermionsThreeBodyGeneric — manual

Source: DiagHam wiki page `FQHESphereFermionsThreeBodyGeneric`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereFermionsThreeBodyGeneric](../FQHEOnSphere/FQHESphereFermionsThreeBodyGeneric.md)
This code diagonalizes a general 3-body FQHE Hamiltonian (with the optional 2-body and 1-body terms). Most of the options are identical to [FQHESphereFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md). The interaction file ("--interaction-file") should contain a line such as:

ThreebodyPseudopotentials= 0 0 0 1

(and optionally one can also provide the 2-body pseudopotentials using "Pseudopotentials=..." and 1-body terms using "Onebodypotentials=..."). 

The numbers following "ThreebodyPseudopotentials=" are the 3-body pseudopotentials giving energy penalty to three fermions in a state with relative angular momentum L=0, 1, 2, 3, etc. For more information, see: https://journals.aps.org/prb/pdf/10.1103/PhysRevB.75.195306 The example above is the interaction that produces the fermionic Moore-Read state at flux 2Q=2N-3, since three fermions are closest together when L=3. 

## Degeneracy

Note that starting with L=9, the space of states for three fermions can be degenerate (e.g., at L=9 there is a two-fold degeneracy). This code may not work properly in such cases.

## Normalization of 3-body pseudopotentials

Unlike the usual two-body (Haldane) pseudopotentials, the present code for fermionic three-body pseudopotentials is not normalized in the standard way (the energies of 3 fermions should be either 0 or 1, regardless of the size of the sphere). 

In the absence of a convenient analytical formula for the normalization, it is recommended to set the normalization via a simple numerical procedure: fix the flux 2Q one wishes to study and diagonalize the three-body Hamiltonian with only N=3 particles (we can set Lz=0 or, alternatively, Lz can be chosen depending on what relative angular momentum we wish to penalize by the 3-body interaction). For the Moore-Read example above, the resulting spectrum will consist of a bunch of exact zero energy states and a single non-zero energy (possibly degenerate depending on Lz). The non-zero energy is the norm we seek; if we divide our 3-body pseudopotential with this energy, we will obtain a correctly normalized Moore-Read parent Hamiltonian. Below is a convenient list of these normalization factors for the Moore-Read state:

   #2Q   3-particle Energy
```text
5 0.269022290028 
7 0.242789392187
9 0.231142245955
11 0.224566229889
13 0.220342378021
15 0.217400402154
17 0.215233737353
19 0.213571626874
21 0.212256227478
23 0.211189311495
25 0.210306552053 
27 0.209564060678
29 0.208930867795
31 0.208384495334
33 0.207908233218
35 0.207489400848
37 0.207118201623
39 0.206786947174
41 0.206489519221
43 0.206220988258
45 0.20597733821
```

Thus one should define the input pseudopotential by diving with the energy at the given value of 2Q to have the properly normalized Moore-Read interaction. For example, if we want the normalized Moore-Read interaction at 2Q=17, we should use "ThreebodyPseudopotentials= 0 0 0 4.646" since 1/0.215233737353 = 4.646.
