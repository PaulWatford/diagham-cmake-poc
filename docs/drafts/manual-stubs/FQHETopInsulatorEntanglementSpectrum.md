# FQHETopInsulatorEntanglementSpectrum (superseded wiki page)

Source: DiagHam wiki page `FQHETopInsulatorEntanglementSpectrum`, as of 2026-09-24. Reason: the program is now `FTIEntanglementSpectrum` and the wiki has a page under that name, which is the attached manual; this older page is kept for comparison.

[FTIEntanglementSpectrum](FTIEntanglementSpectrum.md) provides an easy way to parse the entanglement entropy and the eigenvalues of the reduced density matrices produced by [FTIEntanglementEntropyParticlePartition](FTIEntanglementEntropyParticlePartition.md). Usage is similar to [FQHESphereEntanglementSpectrum](FQHESphereEntanglementSpectrum.md). The typical usage is 

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIEntanglementSpectrum bosons_singleband_kagomelattice_n_6_x_4_y_3_t1_1_t2_0_l1_1_l2_0_gx_0_gy_0_kx_0_ky_0.0.full.parent -n 3 --particle-entanglement*

*bosons_singleband_kagomelattice_n_6_x_4_y_3_t1_1_t2_0_l1_1_l2_0_gx_0_gy_0_kx_0_ky_0.0.full.parent* is the file that contains the reduced density matrix eigenvalues obtained through [FTIEntanglementEntropyParticlePartition](FTIEntanglementEntropyParticlePartition.md). [FTIEntanglementSpectrum](FTIEntanglementSpectrum.md) will create a file named *bosons_singleband_kagomelattice_n_6_x_4_y_3_t1_1_t2_0_l1_1_l2_0_gx_0_gy_0_kx_0_ky_0.0.na_3.parentspec* that will contain the particle entanglement spectrum for $N_A=3$ .
