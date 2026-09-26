# HubbardEntanglementEntropyRealSpacePartition — manual

Source: DiagHam wiki page `HubbardEntanglementEntropyRealSpacePartition`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [HubbardEntanglementEntropyRealSpacePartition](../HubbardModels/HubbardEntanglementEntropyRealSpacePartition.md)

> **Options in this manual that the current program does not list**: `--density-matrix`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
HubbardEntanglementEntropyRealSpacePartition computes the entanglement entropy and entanglement spectrum with a real space partition for the Hubbard models. For example, let's compute the entanglement spectrum of the state fermions_kitaev_heisenberg_gutzwiller_stripe_n_6_x_6_t_1_tK_1_j1_0.49_j2_0.51.0.vec produced by [HubbardKitaevHeisenbergModel](HubbardKitaevHeisenbergModel.md) on a single hexagon (6 sites) where we want to keep only the sites 0 and 2 and to trace out the remaining sites.
The usage is

*\$PATHTODIAGHAM/build/FTI/src/Programs/HubbardModels/HubbardEntanglementEntropyRealSpacePartition fermions_kitaev_heisenberg_gutzwiller_stripe_n_6_x_6_t_1_tK_1_j1_0.49_j2_0.51.0.vec --use-lapack --density-matrix fermions_kitaev_heisenberg_gutzwiller_stripe_n_6_x_6_t_1_tK_1_j1_0.49_j2_0.51.0.full.ent --kept-sites sites_0_2.dat*

where sites_0_2.dat is a column-based ASCII file containing

```text
0
2
```

It lists all the sites that have to be kept. System parameters such as the number of particles or Gutzwiller projection are detected from the input file name. HubbardEntanglementEntropyRealSpacePartition produces two files fermions_kitaev_heisenberg_gutzwiller_stripe_n_6_x_6_t_1_tK_1_j1_1_j2_0.0.ent and fermions_kitaev_heisenberg_gutzwiller_stripe_n_4_x_6_t_1_tK_1_j1_0.49_j2_0.51.0.full.ent. The first file contains the entanglement entropy and the second file contains the entanglement spectrum. This latest looks like

```text
#  N    lambda
2 0.71712927295533
2 0.094290242348223
2 0.094290242348223
2 0.094290242348223
```

The first column corresponds to the number of particles that are kept the second are the reduced density eigenvalues.
