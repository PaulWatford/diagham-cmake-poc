# FQHETorusComputeKx — manual

Source: DiagHam wiki page `FQHETorusComputeKx`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHETorusComputeKx](../FQHEOnTorus/FQHETorusComputeKx.md)
FQHETorusComputeKx allows convert any state to and from the $(k_x,k_y)$.

## Convert from  the $(k_x,k_y)$ basis to the $(k_y)$ basis

Assuming we want to convert the state fermions_torus_4body_hollowcore_n_9_2s_15_ratio_1.000000_kx_0_ky_12.0.vec, we just have to use

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusComputeKx -i fermions_torus_4body_hollowcore_n_9_2s_15_ratio_1.000000_kx_0_ky_12.0.vec  --invert-real --invert --interaction-name readrezayi3*

This will generate an output vector fermions_torus_kysym_readrezayi3_n_9_2s_15_ratio_1.000000_kx_0_ky_12.0.vec . Notice that the option --invert-real should be used only when the input vector is real (up to a global phase factor). This is true when $k_x=0,\pi$ and the torus is rectangular.
