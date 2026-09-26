# SpinChainMultipleComputeS2 — manual

Source: DiagHam wiki page `SpinChainMultipleComputeS2`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [SpinChainMultipleComputeS2](../Spin/SpinChainMultipleComputeS2.md)

> **Options in this manual that the current program does not list**: `--all-eigenstates`, `--full-diag`, `--initial-sz`, `--momentum`, `--nbr-sz`, `--use-lapack`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
[SpinChainMultipleComputeS2](SpinChainMultipleComputeS2.md) computes the total $S^2$ of a series of eigenstates stored in a binary matrix. For example, if we generate all the eigenstates of the periodic AKLT spin chain with 12 sites at momentum $k=0$ and total spin projection $S_z=0$ using [PeriodicSpinChainAKLT](PeriodicSpinChainAKLT.md) the following way

*\$PATHTODIAGHAM/build/Spin/src/Programs/PeriodicSpinChainAKLT -p 12 --all-eigenstates --full-diag 10000 --use-lapack  --momentum 1 --initial-sz 0 --nbr-sz 1*

This will create two binary matrices spin_1_periodicaklt_n_12_sz_0_szsym_-1_k_1.eigenvec.mat and spin_1_periodicaklt_n_12_sz_0_szsym_1_k_1.eigenvec.mat (one for each sector of the symmetry $S_z\leftrightarrow -S_z$) and one spectrum file spin_1_periodicaklt_n_12.dat. If we want to compute all the $S^2$ in the sector $k=0$  $S_z=0$ and $S_{z,{\rm sym}}=-1$, we need the command line

*\$PATHTODIAGHAM/build/Spin/src/Programs/SpinChainMultipleComputeS2 -S --processors 2 --spectrum spin_1_periodicaklt_n_12.dat -c --multiple-states spin_1_periodicaklt_n_12_sz_0_szsym_-1_k_1.eigenvec.mat*

Note that we use the -c option since the matrix  spin_1_periodicaklt_n_12_sz_0_szsym_-1_k_1.eigenvec.mat is complex. Such an option is not needed when the momentum is 0 or $\pi$ because the eigenstate matrix is real.

[SpinChainMultipleComputeS2](SpinChainMultipleComputeS2.md) can also export the eigenstates sorted by their S value. In that case, we just need to add the --export-eigenstate option

*\$PATHTODIAGHAM/build/Spin/src/Programs/SpinChainMultipleComputeS2 -S --processors 2 --spectrum spin_1_periodicaklt_n_12.dat -c --export-eigenstate --multiple-states spin_1_periodicaklt_n_12_sz_0_szsym_-1_k_1.eigenvec.mat*

This will create a binary file another eigenstate matrix spin_1_periodicaklt_n_12_sz_0_szsym_-1_k_1.eigenvec.s2sorted.mat with eigenstates sorted by their S value.
