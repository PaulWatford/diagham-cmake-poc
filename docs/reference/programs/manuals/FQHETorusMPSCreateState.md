# FQHETorusMPSCreateState — manual

Source: DiagHam wiki page `FQHETorusMPSCreateState`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHETorusMPSCreateState](../FQHEOnTorus/FQHETorusMPSCreateState.md)

> **Options in this manual that the current program does not list**: `--use-gmp`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
[FQHETorusMPSCreateState](FQHETorusMPSCreateState.md) generates the decomposition of several model wave functions using their matrix product state description on the torus geometry. It is thus similar to [FQHESphereMPSCreateState](FQHESphereMPSCreateState.md) and handles the same model states.

## Laughlin states

We first consider the Laughlin state $\nu=\frac{1}{3}$. To generate the corresponding two particle Laughlin state on a torus with aspect ratio 2, we just have to use

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusMPSCreateState --p-truncation 3 --nbr-fluxquanta 6 --trim-qsector --aspect-ratio 2.0 --topological-sector 0*

In this case, the program generates on of the three Laughlin state in a binary file fermions_torus_kysym_mps_plevel_3_laughlin3_n_2_2s_6_ratio_2.000000_ky_1.0.vec . If you want to generate a text file instead of a binary file, just add the -t option. in the present example, this will generate a file fermions_torus_kysym_mps_plevel_3_laughlin3_n_2_2s_6_ratio_2.000000_ky_1.0.vec.txt that will contain

```text
-0.99793678381685 0 0 1 0 0 1 
0.045399204317316 0 0 0 1 1 0 
-0.045399204317316 1 1 0 0 0 0 
```

Each line contains the coefficient and its corresponding configuration in the occupation number basis.
 
The number of flux quanta is set via the --nbr-fluxquanta option and the number of particles is automatically evaluated. The torus aspect ratio is tuned through the --aspect-ratio option. The CFT is set via the --p-truncation option.

The two other Laughlin states can be generate by setting other values of --topological-sector. For the Laughlin state $\nu=\frac{1}{m}$, the value of --topological-sector is between 0 and $m-1$. The angular momentum along the y axis is 1 in the above mentioned example. It depends on the number of flux quanta, the truncation in the CFT levels and the topological sector and is automatically evaluated by the code. Looking at other Laughlin states $\nu=\frac{1}{m}$ can be done by providing $m$ through the --laughlin-index option.

The bosonic Laughlin states can be generated in a similar way. We have to specify the --boson option and to set a maximum occupation for any orbital via the --boson-truncation option. 

The program is fully parallel. So you should turn on the -S option and set the number of processors via the --processors options. If you have enough memory, you can also tweak the --precalculation-blocksize option with a larger value than one to precalculate products of B matrices (the value should be of course lower than the number of orbitals/flux quanta, usually a fraction of it). For example, a full command line for the Laughlin state $\nu=\frac{1}{5}$ looks like

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusMPSCreateState --p-truncation 10 --nbr-fluxquanta 30  --trim-qsector --aspect-ratio 1.0 --topological-sector 1  -S --processors 4 --precalculation-blocksize 5 --laughlin-index 5*

(This example takes 4m42 on an intel i7-3517U with 10Gb of RAM and gcc version 4.8.3).

By default, the code uses a rectangular torus. This can be changed by setting the --angle option to a non-zero value. This latest option allows to set the angle (in pi units) between the two vectors 1 and $\tau$ that span the torus. Notice that the length of $\tau$ (or more precisely the ratio of lengths between $\tau$ and 1 is set via the --aspect-ratio option. A typical usage would look like 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusMPSCreateState --p-truncation 9 --nbr-fluxquanta 6 --topological-sector 0  --trim-qsector --angle 0.4 --aspect-ratio 0.75*

This command line produces a complex vector fermions_torus_kysym_mps_plevel_9_laughlin3_n_2_2s_6_ratio_0.750000_angle_1.256637_ky_1.0.vec where the angle is encoded (in radian unit) in the file name (here 1.256637). Beware that codes such as [FQHEFermionsTorusWithTranslation](FQHETorusFermionsWithTranslations.md) set the angle in radian units. Enough digits have to be provided to make an exact comparison between the MPS and the ED code.

## (k=2,r) clustered states and Read-Rezayi k=3 states

The (k=2,r) clustered  states can be generated using the --k-2 option. They include the Moore-Read state (r=2) and the Gaffnian state (r=3). The r index is set through the --r-index option. The Read-Rezayi k=3 states can be obtained by using the --rr-3 option. In both cases, the CFT calculations that are involved in the creation of the B matrices are time consuming. You can provide the name of a directory where all these calculations are stored. This can be done using the --matrices-cft option. If the code does not find the required file in this directory, it computes the CFT calculations and store them in this directory. Notice that archives with the CFT calculations for several states are already provided in FQHESphereMPSCreateState#CFT_data *(wiki page FQHESphereMPSCreateState#CFT_data)* (the data are identical for the genus 0 and genus 1 surfaces).


A typical usage to get the Moore-Read state in the vacuum sector would look like

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusMPSCreateState --p-truncation 6 --nbr-fluxquanta 8 --topological-sector 0 --trim-qsector --k-2 --matrices-cft ../cft_data_pfaffian*

This produces a binary vector fermions_torus_kysym_mps_plevel_6_clustered_k_2_r_2_n_4_2s_8_ratio_1.000000_ky_2.0.vec. Here we assume that use arbitrary accuracy CFT data located in the ../cft_data_pfaffian directory. Remember to add the --use-nonrational option if you use double accuracy data (otherwise don't forget to compile diagham with the --use-gmp option). For the fermionic Moore-Read in the vacuum sector, they are 4 possible values for the topological sector. Beware that two of them will give the same $k_y$ value and thus the same output file name.

To get the quasihole sector, just add the --quasihole-sector option. The output vector file name should look like fermions_torus_kysym_mps_plevel_6_clustered_k_2_qh_r_2_n_4_2s_8_ratio_1.000000_ky_4.0.vec (notice the extra _qh_). Only two topological sector values are possible in that case.

Similarly to the Laughlin case, we can move away from the rectangular torus using the --angle option to a non-zero value.
