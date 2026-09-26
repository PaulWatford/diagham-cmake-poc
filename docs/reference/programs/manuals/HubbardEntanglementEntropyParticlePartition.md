# HubbardEntanglementEntropyParticlePartition — manual

Source: DiagHam wiki page `HubbardEntanglementEntropyParticlePartition`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [HubbardEntanglementEntropyParticlePartition](../HubbardModels/HubbardEntanglementEntropyParticlePartition.md)
HubbardEntanglementEntropyParticlePartition computes the particle entanglement spectrum for Hubbard models. Its usage is almost identical to the one of [FTIEntanglementEntropyParticlePartition](FTIEntanglementEntropyParticlePartition.md). Eigenstates of the reduced density matrix can be exported using the --density-eigenstate . By default, all eigenstates corresponding to an eigenvalue greater than $10^{-14}$ are exported. This behaviour can be changed by setting the option --nbr-eigenstates to the required number of eigenstates.
