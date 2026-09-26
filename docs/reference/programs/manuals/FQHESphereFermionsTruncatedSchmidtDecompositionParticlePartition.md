# FQHESphereFermionsTruncatedSchmidtDecompositionParticlePartition — manual

Source: DiagHam wiki page `FQHESphereFermionsTruncatedSchmidtDecompositionParticlePartition`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereFermionsTruncatedSchmidtDecompositionParticlePartition](../FQHEOnSphere/FQHESphereFermionsTruncatedSchmidtDecompositionParticlePartition.md)
FQHESphereFermionsTruncatedSchmidtDecompositionParticlePartition (for fermions) generates a state using its Schmidt decomposition from a given fermionic state. The Schmidt decomposition is performed in the particle space, contrary to [FQHESphereFermionsTruncatedSchmidtDecomposition](FQHESphereFermionsTruncatedSchmidtDecomposition.md) which works in the orbital space.

## Basic usage

Let's say that we've generated the 1/3 Laughlin state using the Jack generator *(wiki page Jack generator)*. We therefore have  fermions_haldane_laughlin_n_11_2s_30_lz_0.0.vec and the description of the root configuration in laughlin_n_11_2s_30.dat

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereFermionsTruncatedSchmidtDecompositionParticlePartition fermions_haldane_laughlin_n_11_2s_30_lz_0.0.vec  --haldane --reference-file laughlin_n_11_2s_30.dat --use-lapack --use-svd --na 5 -c 13*

This will compute the reduced density matrix for subsystem A with 5 (as specified by option --na)  particles, throw away those eigenvectors with Schmidt number greater than exp(-13) (given by the option -c) for all possible values z-angular momenta of A. Notice that both the --use-lapack and --use-svd are mandatory for the time being. The vector built from the truncated Schmidt decomposition is stored by default in *fermions_haldane_laughlin_n_11_2s_30_lz_0.0.trunc.vec*. The code will also compute for free the eigenvalues of the full (i.e. not truncated) reduced density matrix. The corresponding data will be stored in the ASCII file *fermions_haldane_laughlin_n_11_2s_30_lz_0.0.full.parent*.
