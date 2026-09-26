# FQHESphereMPSEntanglementSpectrumParticlePartition — manual

Source: DiagHam wiki page `FQHESphereMPSEntanglementSpectrumParticlePartition`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereMPSEntanglementSpectrumParticlePartition](../FQHEOnSphere/FQHESphereMPSEntanglementSpectrumParticlePartition.md)

> **Options in this manual that the current program does not list**: `--left-eigenstates`, `--right-eigenstates`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
[FQHESphereMPSEntanglementSpectrumParticlePartition](FQHESphereMPSEntanglementSpectrumParticlePartition.md) allows to compute various type of entanglement spectra for model states that have an exact MPS description. The selection of the model is identical to one of [FQHESphereMPSCreateState](FQHESphereMPSCreateState.md).


## Orbital entanglement spectrum

The orbital entanglement spectrum (OES) can be computed for either a finite or an infinite cylinder. 

### Finite cylinder

Let's start with the Laughlin state $\nu=\frac{1}{3}$. The calculation can be done with or without the padding (i.e. the extra empty orbitals such that the total number of orbitals is commensurate with the number of particles) using the --use-padding option. 

*FQHESphereMPSEntanglementSpectrumParticlePartition --p-truncation 10 --laughlin-index 3 --orbital-es  --normalize-cylinder --cylinder-perimeter  10  --nbr-fluxquanta 11 --use-padding --trim-qsector*

The number of flux quanta is fixed through the --nbr-fluxquanta. Since we use the padding in this example, the number of particles is 4. The perimeter of cylinder in magnetic length unit is set via the --cylinder-perimeter option. The filling factor is fixed via the --laughlin-index option. Parallelization can be turned on setting the -S option and the number of processors via --processors. Beware that the parallel efficiency is far from perfect. So the optimal number of CPUs should be tuned by hand. 

The code produces two files fermions_cylinder_laughlin3_perimeter_10.000000_plevel_10_n_4_2s_11_lz_0.0.ent and fermions_cylinder_laughlin3_perimeter_10.000000_plevel_10_n_4_2s_11_lz_0.0.full.ent . The first file contains the entanglement entropy. Unless the option is used, only the entropy for the half system cut is computed. Here fermions_cylinder_laughlin3_perimeter_10.000000_plevel_10_n_4_2s_11_lz_0.0.full.ent contains a single line

```text
6 0.91994333151389 1
```

where the first column is the number of kept orbitals and the second column is the entanglement entropy (following the conventions of [FQHESphereBosonEntanglementEntropy](FQHESphereBosonEntanglementEntropy.md)). fermions_cylinder_laughlin3_perimeter_10.000000_plevel_10_n_4_2s_11_lz_0.0.full.ent contains the entanglement spectrum and should look like

```text
# x    Q    P    lambda    -ln(lambda)
0 8 0 0.032198003340378 3.4358508364172
0 11 0 0.64597356129258 0.43699670282873
0 11 1 0.2776840852705 1.2812711953843
0 11 2 0.0095343028465402 4.6528591578057
0 11 2 0.0017675431377849 6.3381647543896
0 11 3 0.0005233955469789 7.5551730757781
0 11 4 0.00012110522485752 9.0188507633523
0 14 0 0.032198003340378 3.4358508364172
```


Instead of the usual way to label the entanglement energies, as used in [FQHESphereBosonEntanglementEntropy](FQHESphereBosonEntanglementEntropy.md), [FQHESphereMPSEntanglementSpectrumParticlePartition](FQHESphereMPSEntanglementSpectrumParticlePartition.md) relies on the quantum numbers $(x,Q,P)$ that label the auxiliary states. $Q$ measures the charge in unit of the quasihole charge, up to a global shift. For the $\nu=\frac{1}{m}$ Laughlin state, the relation between the truncation in the CFT level $P_{max}$, the number of orbitals $N_{orb,A}$ in the part A and the number of particles $N_A$ in the part A is $Q=N_{orb,A}+P_{max}-m N_A$. $x$ labels the field of the neutral sector. For tha Laughlin case, it is always equal to zero (i.e. the identity field). The relation between the CFT level and the total momentum of the part A is a bit more involved. $P$ measures the total momentum with respect to the densest root configuration for $N_{A}$ particles.


We now consider another example, the Gaffnian state with 14 fermions (and thus 31 flux quanta). In that case, one should use

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereMPSEntanglementSpectrumParticlePartition --p-truncation 8 --k-2 --r-index 3 --orbital-es --use-nonrational -S --normalize-cylinder --matrices-cft cft_data_gaffnian_num --nbr-fluxquanta 31 --cylinder-perimeter 15 --trim-qsector*

Here the CFT data are provided through the --matrices-cft option. CFT data are available here *(wiki page FQHESphereMPSCreateState#CFT data)*. The --use-nonrational flag is used to indicates that the CFT data located in the cft_data_gaffnian_num directory use real numbers and not rational numbers.

### Infinite cylinder

To compute the OES of an infinite cylinder, one has to provide the transfer matrix left and right eigenstates corresponding to the largest eigenvalue. These can be generated using [FQHEMPSEMatrix](FQHEMPSEMatrix.md). For example 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHEMPSEMatrix --p-truncation 10 --trim-qsector  --ematrix-memory 2000 -S  --use-nonrational --fixed-qsector --laughlin-index 3 --show-itertime --right-eigenstates --left-eigenstates --diagonal-block --normalize-cylinder --processors 4 --cylinder-perimeter 15*

generates the left and right eigenstates ematrix_diagblock_cylinder_laughlin3_fixedq_0_3_perimeter_15.000000_plevel_10_left.0.vec and ematrix_diagblock_cylinder_laughlin3_fixedq_0_3_perimeter_15.000000_plevel_10.dat. With these vectors, we now use [FQHESphereMPSEntanglementSpectrumParticlePartition](FQHESphereMPSEntanglementSpectrumParticlePartition.md) to get the infinite OES

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereMPSEntanglementSpectrumParticlePartition --p-truncation 10 --right-eigenstate ematrix_diagblock_cylinder_laughlin3_fixedq_0_3_perimeter_15.000000_plevel_10_right.0.vec --left-eigenstate ematrix_diagblock_cylinder_laughlin3_fixedq_0_3_perimeter_15.000000_plevel_10_left.0.vec --orbital-es --laughlin-index 3 --infinite-cylinder --fixed-qsector --trim-qsector --use-nonrational --diagonal-block -S --normalize-cylinder --cylinder-perimeter 15*

Beware than the same parameters (such as --cylinder-perimeter, --laughlin-index 3, --p-truncation, --fixed-qsector, --trim-qsector and --diagonal-block),  have to be used for both [FQHEMPSEMatrix](FQHEMPSEMatrix.md) and [FQHESphereMPSEntanglementSpectrumParticlePartition](FQHESphereMPSEntanglementSpectrumParticlePartition.md). 

Two files are produced : fermions_infinite_cylinder_laughlin3_fixedq_0_3_perimeter_15.000000_plevel_10_n_0_2s_0_lz_0.0.full.ent that contains the entanglement spectrum and fermions_infinite_cylinder_laughlin3_fixedq_0_3_perimeter_15.000000_plevel_10_n_0_2s_0_lz_0.0.ent that gives the entanglement entropy. The content of these files is similar to the finite cylinder case.

## Particle entanglement spectrum

### Finite cylinder

### Infinite cylinder

## Real space entanglement spectrum

### Finite cylinder

### Infinite cylinder
