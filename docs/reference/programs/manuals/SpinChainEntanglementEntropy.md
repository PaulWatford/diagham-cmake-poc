# SpinChainEntanglementEntropy — manual

Source: DiagHam wiki page `SpinChainEntanglementEntropy`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [SpinChainEntanglementEntropy](../Spin/SpinChainEntanglementEntropy.md)
SpinChainEntanglementEntropy computes the entanglement spectrum and entanglement entropy of a spin chain. Its typical usage is

*\$PATHTODIAGHAM/build/Spin/src/Programs/SpinChainEntanglementEntropy spin_2_periodicaklt_p3_1.000000_p4_1.000000_n_8_sz_0_s_0_invsym_1_szsym_1_k_4.21.vec --use-lapack --use-svd*

The --use-svd greatly improves the accuracy but only works when computing the entanglement spectrum of a single state. The code produces two output files whose name is deduced from the input vector, here spin_2_periodicaklt_p3_1.000000_p4_1.000000_n_8_sz_0_s_0_invsym_1_szsym_1_k_4.21.ent and  spin_2_periodicaklt_p3_1.000000_p4_1.000000_n_8_sz_0_s_0_invsym_1_szsym_1_k_4.21.full.ent. The first one contains the entanglement entropies for all cuts up to the half the system size (unless the --min-la and/or --max-la options are used) and should look like

```text
1 1.6094379124341 1
2 2.8249257205672 1
3 3.7717695164438 1
4 4.0513979917824 1
```

- The first column is number of sites in the region A.
- The second column is the Von Neumann entanglement entropy.
- The third column is the trace of the reduced matrix and should be equal to one.

The second file spin_2_periodicaklt_p3_1.000000_p4_1.000000_n_8_sz_0_s_0_invsym_1_szsym_1_k_4.21.full.ent stores the entanglement spectrum

```text
# l_a    Sz    lambda
1 -4 0.2
1 -2 0.2
1 0 0.2
1 2 0.2
1 4 0.2
2 -8 0.0060679611650486
2 -6 0.027912621359223
2 -6 0.0060679611650486
2 -4 0.095873786407767
...
```

- The first column is number of sites in the region A.
- The second column is twice the $S_z$ of the region A.
- The third column contains the eigenvalues of the reduced density matrix.
