# FQHESphereFermionsTruncatedSchmidtDecomposition — manual

Source: DiagHam wiki page `FQHESphereFermionsTruncatedSchmidtDecomposition`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereFermionsTruncatedSchmidtDecomposition](../FQHEOnSphere/FQHESphereFermionsTruncatedSchmidtDecomposition.md)
FQHESphereFermionsTruncatedScmidtDecomposition (for fermions) generates a state with its Schmidt decomposition truncated from a given fermionic state.

## Basic usage

Let's say that we've generated the 1/3 Laughlin state using the Jack generator *(wiki page Jack generator)*. We therefore have  fermions_haldane_laughlin_n_11_2s_30_lz_0.0.vec and the description of the root configuration in laughlin_n_11_2s_30.dat

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereFermionsTruncatedSchmidtDecomposition fermions_haldane_laughlin_n_11_2s_30_lz_0.0.vec  --haldane --reference-file laughlin_n_11_2s_30.dat --use-lapack -a 15 -c 13*

This will compute the reduced density matrix for subsystem A with 15 (as specified by option -a) orbitals, throw away those eigenvectors with Schmidt number greater than exp(-13) (given by the option -c) for all possible values of the number of particles in A and z-angular momentum of A and return the vector with the truncated spectrum in

*fermions_haldane_laughlin_n_11_2s_30_lz_0.0.trunc.vec*

by default.
