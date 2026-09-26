# FQHETopInsulatorEntanglementEntropyParticlePartition (superseded wiki page)

Source: DiagHam wiki page `FQHETopInsulatorEntanglementEntropyParticlePartition`, as of 2026-09-24. Reason: the program is now `FTIEntanglementEntropyParticlePartition` and the wiki has a page under that name, which is the attached manual; this older page is kept for comparison.

[FTIEntanglementEntropyParticlePartition](../../reference/programs/manuals/FTIEntanglementEntropyParticlePartition.md) computes the particle entanglement spectrum for systems with a Chern insulator in the single band approximation. Usage is similar to [FQHESphereBosonEntanglementEntropy](../../reference/programs/manuals/FQHESphereBosonEntanglementEntropy.md). The code autodetects most of the system information such as the statistics or system size from the input vector names.
 
Compare to [FQHESphereBosonEntanglementEntropy](../../reference/programs/manuals/FQHESphereBosonEntanglementEntropy.md), the total density matrix will generally rely on more than one projector i.e.

$\rho=\frac{1}{N}\sum_{i=1}^{N} |\psi_i  \psi_i |$

In that case, you will need to use the --degenerated-groundstate option to give the list of states. Such a single column text file typically looks like

```text
fermions_singleband_checkerboardlattice_n_6_x_6_y_3_v_0.000000_t1_1.000000_t2_0.292893_gx_0.000000_gy_0.000000_kx_3_ky_0.0.vec
fermions_singleband_checkerboardlattice_n_6_x_6_y_3_v_0.000000_t1_1.000000_t2_0.292893_gx_0.000000_gy_0.000000_kx_3_ky_0.1.vec
fermions_singleband_checkerboardlattice_n_6_x_6_y_3_v_0.000000_t1_1.000000_t2_0.292893_gx_0.000000_gy_0.000000_kx_3_ky_0.2.vec
```


The command line has to be

*FTI/src/Programs/FTI/FTIEntanglementEntropyParticlePartition --degenerated-groundstate ground.dat --use-lapack --density-matrix fermions_singleband_checkerboardlattice_n_6_x_6_y_3_v_0.000000_t1_1.000000_t2_0.292893_gx_0.000000_gy_0.000000_kx_3_ky_0.0.full.parent --show-time*

ground.dat . If you want to look at the entanglement spectrum, you should use the --density-matrix option to save the eigenvalues of the reduced density matrix and then use [FTIEntanglementSpectrum](../../reference/programs/manuals/FTIEntanglementSpectrum.md) to parse the output file. --show-time allows to show the amount of time to compute each block of the reduced density matrix. [FTIEntanglementEntropyParticlePartition](../../reference/programs/manuals/FTIEntanglementEntropyParticlePartition.md) supports the SMP mode. This can greatly decrease the amount of time to compute each block of the reduced density matrix.
