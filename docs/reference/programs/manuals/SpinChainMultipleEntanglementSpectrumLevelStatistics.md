# SpinChainMultipleEntanglementSpectrumLevelStatistics — manual

Source: DiagHam wiki page `SpinChainMultipleEntanglementSpectrumLevelStatistics`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [SpinChainMultipleEntanglementSpectrumLevelStatistics](../Spin/SpinChainMultipleEntanglementSpectrumLevelStatistics.md)
[SpinChainMultipleEntanglementSpectrumLevelStatistics](SpinChainMultipleEntanglementSpectrumLevelStatistics.md) allows to manipulate entanglement spectrum files produced by [SpinChainMultipleEntanglementSpectra](SpinChainMultipleEntanglementSpectra.md).

## Extracting a single entanglement spectrum

Let say we have produce the entanglement spectrum file spin_1_periodicaklt_n_12_sz_0_invsym_1_szsym_1_k_0.eigenvec_la_6.full.ent and we want to extract the entanglement spectrum for the eigenstate with index 4 (the lowest energy state having index 0). We just need to run 

*\$PATHTODIAGHAM/build/Spin/src/Programs/SpinChainMultipleEntanglementSpectrumLevelStatistics -s spin_1_periodicaklt_n_12_sz_0_invsym_1_szsym_1_k_0.eigenvec_la_6.full.ent --extract-singlespectrum 4*

A text file  spin_1_periodicaklt_n_12_sz_0_invsym_1_szsym_1_k_0.eigenvec_la_6.full.spec_4.ent will be produced that should look like


```text
-12 7.79406248455e-22 48.603509821978
-10 1.4498534754712e-14 31.864728801989
-10 1.2350354885308e-17 38.932846875579
-10 6.4075538433171e-19 41.891639184495
-10 3.0916308287941e-19 42.620418038893
-10 9.2953922072206e-22 48.427353230081
-10 7.7940615982093e-22 48.603509935698
-8 1.1538772156394e-09 20.580138073447
-8 1.4406580497704e-10 22.660750941728
...
```

The first column is the $S_z$ for the region A, the second column contains the reduced density matrix eigenvalue and the third column contains the entanglement energies.
