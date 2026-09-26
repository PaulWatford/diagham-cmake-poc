# FQHESphereLLLStatesProduct — manual

Source: DiagHam wiki page `FQHESphereLLLStatesProduct`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereLLLStatesProduct](../FQHEOnSphere/FQHESphereLLLStatesProduct.md)
FQHESphereLLLStatesProduct allows to compute the product of two states (either fermionic or bosonic) that belongs to the lowest Landau level.

## Basic usage

The first way to use this program is to have only one input state. In this case, the square of this state will be computed.

\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereLLLStatesProduct bosons_unnormalized_laughlin2_n_6_2s_10_lz_0.0.vec -o bosons_unnormalized_laughlin4_n_6_2s_20_lz_0.0.vec

Note that the input vector is given in the unnormalized basis. In this case, one can use Haldane basis by adding the option --haldane and --reference-file followed by the reference file corresponding to the input state.

The second way is to provide two input states. 

\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereLLLStatesProduct bosons_unnormalized_laughlin2_n_6_2s_10_lz_0.0.vec fermions_unnormalized_laughlin3_n_6_2s_15_lz_0.0.vec  -o fermions_unnormalized_laughlin5_n_6_2s_25_lz_0.0.vec

The option --normalize indicates that the product state will be normalized on the sphere. 

## Important remarks

The input states must be given in the unnormalized basis.

In the case two input states are provided, the haldane basis option are not yet supported.

Up to now, this program only supports FermionOnSphere and BosonOnSphereShort Hilbert spaces.
