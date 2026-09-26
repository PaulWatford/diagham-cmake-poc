# FTIGenerateMagnon — manual

Source: DiagHam wiki page `FTIGenerateMagnon`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FTIGenerateMagnon](../FTI/FTIGenerateMagnon.md)
FTIGenerateMagnon creates magnon-like states on top of any many-body state. Depending on the the initial Sz (or Bz/Pz), it constructs all bilinears \$\gamma^\dagger_{q1+q',\downarrow} \gamma_{q1,\uparrow}|\Psi\rangle\$ (if Sz>=0) or \$\gamma^\dagger_{q1+q',\uparrow} \gamma_{q1,\downarrow}|\Psi\rangle\$ (if Sz<0) in the momentum sector \$q+k\$ (modulo the system size) where \$k\$ is the ground state momentum. A typical usage is


*FTIGenerateMagnon --eigenstate-file fermions_twoband_u_1.000_test_n_14_ns_42_x_21_y_1_kx_0_ky_0_bz_-14.0.vec*


where *fermions_twoband_u_1.000_test_n_14_ns_42_x_21_y_1_kx_0_ky_0_bz_-14.0.vec* is the state on top of which bilinears should be built. By default, bilinears will be generated for all momentum sectors, a behavior that can be overridden by the options --only-kx and --only-ky.

All generated states are normalized to 1 (if a bilinear operator annihilates the state, no output state is recorded). The output file name for the generated states should look like *fermions_twoband_u_1.000_test_n_14_ns_42_x_21_y_1_kx_0_ky_0_bz_-14.0_magnon_qx_20_qy_0_q1x_0_q1y_0_bz_-12.vec*. In particular the total momentum of the state (corresponding to the sum of \$q\$ and the momentum of the ground state modulo the \$N_x\$ and \$N_y\$ is given in qx and qy. What is denoted by q1x and q1y correspond to the formulas given above. Note that except in some specific case, generated states in a given momentum sector do not have to be orthogonal. You can use [ExtractLinearlyIndependentVectors](ExtractLinearlyIndependentVectors.md) to extract an orthogonal basis.
