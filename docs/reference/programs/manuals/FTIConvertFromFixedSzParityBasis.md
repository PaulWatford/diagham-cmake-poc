# FTIConvertFromFixedSzParityBasis — manual

Source: DiagHam wiki page `FTIConvertFromFixedSzParityBasis`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FTIConvertFromFixedSzParityBasis](../FTI/FTIConvertFromFixedSzParityBasis.md)
FTIConvertFromFixedSzParityBasis converts spinful states written in real space from a basis with a defined $S_z$ parity to a basis without such a constraint. Its usage is 

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIConvertFromFixedSzParityBasis -i fermions_su2_realspace_hofstadter_X_4_Y_2_q_1_n_6_ns_24_x_1_y_3_u_-20_gx_0_gy_0_sz_0_szsym_-1_kx_0_ky_0.0.vec*

The converted state has the same output name than the input file without the string szsym_xx_ .

FTIConvertFromFixedSzParityBasis can be used in the opposite direction i.e. from a basis whithout a defined $S_z$ parity to a basis with this constraint. This is automatically detected from the input file name. The $S_z$ parity is also computed automatically but this is time consuming (but SMP mode can be turned on). Otherwise the sector can be manually set with the option . The reverse mode can be used as

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIConvertFromFixedSzParityBasis -i fermions_su2_realspace_hofstadter_X_4_Y_2_q_1_n_6_ns_24_x_1_y_3_u_-20_gx_0_gy_0_sz_0_kx_0_ky_0.0.vec -S  --processors 4*

This will create a vector with the correct naming convention fermions_su2_realspace_hofstadter_X_4_Y_2_q_1_n_6_ns_24_x_1_y_3_u_-20_gx_0_gy_0_sz_0_szsym_-1_kx_0_ky_0.0.vec. The standard output should look like

```text
Convert from fixed Sz parity, kx, ky basis to full basis
Nbr particles=6 Nbr sites=24 kx = 0 ky = 0 Sz = 0
intermediate Hilbert space dimension = 4096576
Hilbert space dimension = 1365568
memory requested for Hilbert space = 15Mo
memory requested for lookup table = 6Mo
Sz parity of fermions_su2_realspace_hofstadter_X_4_Y_2_q_1_n_6_ns_24_x_1_y_3_u_-20_gx_0_gy_0_sz_0_kx_0_ky_0.0.vec = (-1,0)
intermediate Hilbert space dimension = 4096576
Hilbert space dimension = 683124
memory requested for Hilbert space = 7Mo
memory requested for lookup table = 6Mo
dimension = 683124 683124
```

Notice that the $S_z$ parity computed from the input vector is displayed. If it is far from either +1 or -1, you should double check that your input state is not part of a degenerate subspace.
