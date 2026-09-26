# FQHESphereMPSCreateState — manual

Source: DiagHam wiki page `FQHESphereMPSCreateState`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereMPSCreateState](../FQHEOnSphere/FQHESphereMPSCreateState.md)
FQHESphereMPSCreateState generates the decomposition of several model wave functions using their matrix product state description. The state can be obtained in the unnormalized basis (i.e. the Jack normalization) or on the cylinder geometry.

## Laughlin states

We first consider the Laughlin state $\nu=\frac{1}{3}$. The calculation can be done with or without the padding (i.e. the extra empty orbitals such that the total number of orbitals is commensurate with the number of particles) using the --use-padding option. 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereMPSCreateState --laughlin-index 3 --reference-file laughlin3_n_4_2s_11.dat --use-padding --p-truncation 1 -t fermions_unnormalized_haldane_plevel_1_laughlin3_n_4_2s_11_lz_0.0.txt*

The root configuration is specified via the --reference-file option which uses the same kind of root description file than [FQHESphereJackGenerator](FQHESphereJackGenerator.md).  Since we use the padding in this example, the file laughlin3_n_4_2s_11.dat should look like

```text
NbrParticles = 4
LzMax = 11
ReferenceState = 0 1 0 0 1 0 0 1 0 0 1 0
```

Here we have set the CFT truncation to 1  (--p-truncation 1) and we have chosen to export the resulting state decompostion in an ASCII text file using the -t option. The output file fermions_unnormalized_haldane_plevel_1_laughlin3_n_4_2s_11_lz_0.0.txt contains

```text
1 [10,7,4,1]
-3 [10,7,3,2]
-3 [10,6,5,1]
6 [10,6,4,2]
-3 [9,8,4,1]
9 [9,8,3,2]
6 [9,7,5,1]
-12 [9,7,4,2]
```


A binary output can be obtained by using the -o option instead of the -t option. The filling factor is fixed via the --laughlin-index option and thus we can build the $\nu=\frac{1}{m}$ in a similar way. 

If we now want to do the same calculation on the cylinder geometry, we just have to set the two options --normalize-cylinder and --cylinder-perimeter 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereMPSCreateState --laughlin-index 3 --reference-file laughlin3_n_4_2s_11.dat --use-padding --p-truncation 1 -t fermions_cylinder_perimeter_5.000000_plevel_1_laughlin3_n_4_2s_11_lz_0.0.txt --normalize-cylinder --cylinder-perimeter 5.0*

We work on a torus with a perimeter of 5 magnetic lengths (--cylinder-perimeter 5.0). The size of the cylinder can also be set through the aspect ratio  (--aspect-ratio) instead of the perimeter. fermions_cylinder_perimeter_5.000000_plevel_1_laughlin3_n_4_2s_11_lz_0.0.txt should look like

```text
0.97634612782278 [10,7,4,1]
-0.12448136712101 [10,7,3,2]
-0.12448136712101 [10,6,5,1]
0.0021812391027076 [10,6,4,2]
-0.12448136712101 [9,8,4,1]
0.015871021883264 [9,8,3,2]
0.0021812391027076 [9,7,5,1]
-3.8221013579932e-05 [9,7,4,2]
```

Note that the truncated state is normalized to 1. So the coefficients that can be reached with this CFT truncation, are identical to the full state coefficients up to a global factor.

The bosonic Laughlin states can be generated in a similar way. We have to specify the --boson option and to set a maximum occupation for any orbital via the --boson-truncation option.

## (k=2,r) clustered states

These states can be generated using the --k-2 option. They include the Moore-Read state (r=2) and the Gaffnian state (r=3). The r index is set through the --r-index option.

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereMPSCreateState --k-2 --r-index 2 --reference-file pfaffian_n_6_2s_9.dat --p-truncation 2 -t fermions_pfaffian_n_6_2s_9_lz_0.0.txt*


The root configuration is provided through the pfaffian_n_6_2s_9.dat file. In that case, it reads

```text
NbrParticles = 6
LzMax = 9
ReferenceState = 1 1 0 0 1 1 0 0 1 1
```


The output file fermions_pfaffian_n_6_2s_9_lz_0.0.txt will contain 

```text
1 [9,8,5,4,1,0]
-2 [9,8,5,3,2,0]
10 [9,8,4,3,2,1]
-2 [9,7,6,4,1,0]
4 [9,7,6,3,2,0]
2 [9,7,5,4,2,0]
-16 [9,7,5,3,2,1]
-14 [9,6,5,4,3,0]
28 [9,6,5,4,2,1]
10 [8,7,6,5,1,0]
-16 [8,7,6,4,2,0]
28 [8,7,6,3,2,1]
28 [8,7,5,4,3,0]
-6 [8,7,5,4,2,1]
```


The CFT calculations that are involved in the creation of the B matrices are time consuming. You can provide the name of a directory where all these calculations are stored. This can be done using the --matrices-cft option. If the code does not find the required file in this directory, it computes the CFT calculations and store them in this directory. The previous example using this option will look like


*PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereMPSCreateState --k-2 --r-index 2 --reference-file pfaffian_n_6_2s_9.dat --p-truncation 2 -t fermions_pfaffian_n_6_2s_9_lz_0.0.txt  --matrices-cft cft_data_pfaffian*

where cft_data_pfaffian is the directory name. Sets of data are available below. By default, the CFT calculations are performed using arbitrary precision numbers. This can be turned off using the --use-nonrational option. The calculations are then faster but they might be less reliable. Beware that if you use the --matrices-cft option on top of --use-nonrational, the code will create (or use) another series of files in the directory (inlcuding a _num_ tag in the file names).

## Read-Rezayi k=3 states

## N=1 Superconformal

## CFT data

| State | arb. accuracy | double accuracy |
|---|---|---|
| Moore-Read (vacuum sector) | [http://www.nick-ux.org/~regnault/mpsdata/cft/cft_data_pfaffian.tgz Pmax=13] | [http://www.nick-ux.org/~regnault/mpsdata/cft/cft_data_pfaffian_num.tgz Pmax=15] |
| Moore-Read (quasihole sector) | [http://www.nick-ux.org/~regnault/mpsdata/cft/cft_data_pfaffian_qh.tgz Pmax=13] | [http://www.nick-ux.org/~regnault/mpsdata/cft/cft_data_pfaffian_qh_num.tgz Pmax=15] |
| Gaffnian (vacuum sector) | [http://www.nick-ux.org/~regnault/mpsdata/cft/cft_data_gaffnian_num.tgz Pmax=13] | [http://www.nick-ux.org/~regnault/mpsdata/cft/cft_data_gaffnian_num.tgz Pmax=15] |
| Gaffnian (quasihole sector) | [http://www.nick-ux.org/~regnault/mpsdata/cft/cft_data_gaffnian_qh.tgz Pmax=13] | [http://www.nick-ux.org/~regnault/mpsdata/cft/cft_data_gaffnian_qh_num.tgz Pmax=14] |
| Read-Rezayi k=3 (vacuum sector) | [http://www.nick-ux.org/~regnault/mpsdata/cft/cft_data_readrezayi3.tgz Pmax=13] | [http://www.nick-ux.org/~regnault/mpsdata/cft/cft_data_readrezayi3_num.tgz Pmax=14] |
| Read-Rezayi k=3 (quasihole sector) | [http://www.nick-ux.org/~regnault/mpsdata/cft/cft_data_readrezayi3_qh.tgz Pmax=13] | [http://www.nick-ux.org/~regnault/mpsdata/cft/cft_data_readrezayi3_qh_num.tgz Pmax=13] |
| Jack (k=2, r=6) (vacuum sector) | [http://www.nick-ux.org/~regnault/mpsdata/cft/cft_data_jack_k_2_r_6.tgz Pmax=12] | [http://www.nick-ux.org/~regnault/mpsdata/cft/cft_data_jack_k_2_r_6_num.tgz Pmax=14] |
| Ising tricritical (vacuum sector) |  | [http://www.nick-ux.org/~regnault/mpsdata/cft/cft_data_clustered_k_2_r_6_ising_tricritical_num.tgz Pmax=11] |
