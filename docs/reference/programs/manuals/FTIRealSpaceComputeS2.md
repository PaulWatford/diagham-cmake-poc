# FTIRealSpaceComputeS2 — manual

Source: DiagHam wiki page `FTIRealSpaceComputeS2`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FTIRealSpaceComputeS2](../FTI/FTIRealSpaceComputeS2.md)
FTIRealSpaceComputeS2 computes the expectation value of $S^2$ for lattice models written in real space. A typical usage is

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIRealSpaceComputeS2 -i fermions_realspace_checkerboardlatticewithspin_n_12_ns_12_x_3_y_2_kx_2_ky_1_sz_0.0.vec*

The output should look like 

```text
intermediate Hilbert space dimension = 853776
Hilbert space dimension = 142224
memory requested for Hilbert space = 1Mo
memory requested for lookup table = 6Mo
fermions_realspace_checkerboardlatticewithspin_n_12_ns_12_x_3_y_2_kx_2_ky_1_sz_0.0.vec.ref :  < S^2 >=(2,1.49205e-16) < S >=1
```
