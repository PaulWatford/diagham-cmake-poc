# FQHETorusBosonsEntanglementEntropy — manual

Source: DiagHam wiki page `FQHETorusBosonsEntanglementEntropy`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHETorusBosonsEntanglementEntropy](../FQHEOnTorus/FQHETorusBosonsEntanglementEntropy.md)
FQHETorusBosonsEntanglementEntropy allows to compute the orbital entanglement entropy and the corresponding entanglement spectrum for bosons on the torus geometry. Its usage is similar to [FQHESphereBosonEntanglementEntropy](FQHESphereBosonEntanglementEntropy.md).  Notice that the states have to be provided without the relative momentum $K_x$ (conversion can be done via [FQHETorusComputeKx](FQHETorusComputeKx.md)) but the vectors can either real or complex (-c to indicate complex vectors).

As an example, we will show how to compute the orbital entanglement entropy for a state defined through a linear combination of two vectors. we consider the two zero energy Moore-Read ground states in the momentum sectors $(0,0)$ and $(\pi ,0)$ for 4 bosons bosons_torus_kysym_3body_hardcore_n_4_2s_4_ratio_1.000000_kx_0_ky_0.0.vec and bosons_torus_kysym_3body_hardcore_n_4_2s_4_ratio_1.000000_kx_2_ky_0.0.vec . These states have been generated using [FQHETorusBosonsWithTranslationsGenericNBody](FQHETorusBosonsWithTranslationsGenericNBody.md) and then converting to the $K_y$ basis using [FQHETorusComputeKx](FQHETorusComputeKx.md). If we want to look at the superposition defined with the same weight, we shall use the following command line

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusBosonsEntanglementEntropy --use-lapack --use-svd -c --degenerated-groundstate basis.dat --min-la 2 --max-la 2 -o oes_pfaffian_n_4_theta_pi_4.dat*

where basis.dat contains the description of the complex linear superposition

```text
bosons_torus_kysym_3body_hardcore_n_4_2s_4_ratio_1.000000_kx_0_ky_0.0.vec (0.707106781186547573,0)
bosons_torus_kysym_3body_hardcore_n_4_2s_4_ratio_1.000000_kx_2_ky_0.0.vec (0.707106781186547573,0)
```

Here we override the output file name where the entropies are stored, setting it to oes_pfaffian_n_4_theta_pi_4.dat. Moreover, we exactly focus on the half cut, keeping only two orbitals, thanks to the options --min-la 2 --max-la 2. Notice that as long as we compute teh entanglement entropy for a single state (be it defined as a linear superposition of severals wave functions) it is safe and more accurate to use the --use-svd option. The contain of the output file should be

```text
2 0.34387962854592 1 -2.2204460492503e-16
```

The first column is the number of orbitals that you keep (i.e. $l_a$), the second is the Von Neumann entropy, the third one is the trace of the reduced density matrix (should be one), the fourth is the error on the reduced density matrix trace (should be 0 up to the machine precision).
