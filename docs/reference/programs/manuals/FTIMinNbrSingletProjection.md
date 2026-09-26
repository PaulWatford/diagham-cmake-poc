# FTIMinNbrSingletProjection — manual

Source: DiagHam wiki page `FTIMinNbrSingletProjection`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FTIMinNbrSingletProjection](../FTI/FTIMinNbrSingletProjection.md)
[FTIMinNbrSingletProjection](FTIMinNbrSingletProjection.md) projects a fermionic spinful state onto a space with a minimum number of on-site singlets. Its usage is either

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIMinNbrSingletProjection --input-state fermions_su2_realspace_hofstadter_X_4_Y_2_q_1_n_6_ns_24_x_1_y_3_u_-20_gx_0_gy_0_sz_0_szsym_-1_kx_0_ky_0.0.vec --min-nbrsinglets 3*

for a single state or 

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIMinNbrSingletProjection --degenerate-states basis.dat --min-nbrsinglets 4*

when several states have to be processed (in this example, basis.dat is a single column ASCII file containing the vector list). The minimum number of on-site singlets is set by the --min-nbrsinglets option. The typical output will look like

```text
intermediate Hilbert space dimension = 1293121600
Hilbert space dimension = 161651036
memory requested for Hilbert space = 1849Mo
memory requested for lookup table = 8Mo
intermediate Hilbert space dimension = 35960
Hilbert space dimension = 9024
memory requested for Hilbert space = 105ko
memory requested for lookup table = 8Mo
fermions_su2_realspace_hofstadter_X_4_Y_2_q_1_ns_32_n_8_x_1_y_4_u_-100_gx_0_gy_0_sz_0_szsym_1_kx_0_ky_0.0.vec : 0.99462566426357
fermions_su2_realspace_hofstadter_X_4_Y_2_q_1_ns_32_n_8_x_1_y_4_u_-100_gx_0_gy_0_sz_0_szsym_1_kx_0_ky_0.1.vec : 0.99463120531123
```

The two last lines of this example give the weights of each state in the projected space. The projected states are stored in binary vectors with an additional string indicating the minimum number of on-site singlets (in this example, the stored vectors will have names such as fermions_su2_realspace_hofstadter_X_4_Y_2_q_1_ns_32_minnbrsinglet_4_n_8_x_1_y_4_u_-100_gx_0_gy_0_sz_0_szsym_1_kx_0_ky_0.0.vec). The storage of the projected states can be bypassed by adding the --only-weights option.
