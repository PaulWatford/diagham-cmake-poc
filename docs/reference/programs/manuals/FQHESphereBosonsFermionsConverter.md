# FQHESphereBosonsFermionsConverter — manual

Source: DiagHam wiki page `FQHESphereBosonsFermionsConverter`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereBosonsFermionsConverter](../FQHEOnSphere/FQHESphereBosonsFermionsConverter.md)
FQHESphereBosonsFermionsConverter computes the bosonic (resp. fermionic) counterpart of a given fermionic (resp. bosonic) expressed in <em>the unnormalized </em> basis.

More precisely, if the input states are fermionic (resp. bosonic), they will be divided (multiplied) by a Jastrow factor.

Typical usage is :

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereBosonsFermionsConverter fermions_laughlin3_n_6_2s_15_lz_0.0.vec -o fermions_laughlin2_n_6_2s_10_lz_0.0.vec*

More than one input state can be given but they must all be either bosonic or fermionic.

The squeezed basis can be used but, in this case, the reference-file to give as a parameter after the option --reference-file must correspond to the bosonic basis. For instance, the same result of the previous command will be obtained with:

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereBosonsFermionsConverter fermions_laughlin3_n_6_2s_15_lz_0.0.vec -o fermions_laughlin2_n_6_2s_10_lz_0.0.vec --haldane --reference-file laughlin2_n_6_2s_10.dat*

where laughlin2_n_6_2s_10.dat is a text file with:
```   NbrParticles = 6
```text
LzMax = 10
ReferenceState = 1 0 1 0 1 0 1 0 1 0 1
```
```

This program can be used in a parallelized way via the options -S and --processors (followed by the number of thread to use). When the input states are fermionic, the code is completely parallelized over the Hilbert Space as it consists in computing a vector for each Fock states and multiply the input states by the matrix whose columns are those vectors whereas  when the input states are bosonic, this matrix must be inverted which is not parallelized.
