# FCICheckerboardLatticeModel — manual

Source: DiagHam wiki page `FCICheckerboardLatticeModel`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FCICheckerboardLatticeModel](../FCI/FCICheckerboardLatticeModel.md)
[FCICheckerboardLatticeModel](FCICheckerboardLatticeModel.md) handles interacting fermions or bosons on the checkerboard lattice. This exact diagonalization code offers similar calculation features to those of  [QHEFermionsTwoBodyGeneric](FQHESphereFermionsTwoBodyGeneric.md). The typical usage is 

*\$PATHTODIAGHAM/build/FTI/src/Programs/FCI/FCICheckerboardLatticeModel -p 4 -x 4 -y 3 --flat-band --single-band --use-lapack*

The -p option sets the number of fermions, -x (resp. -y) the number of sites in the s (resp. y) direction. The --flat-band allows to work in the flat band limit and --single-band projects in the Hamiltonian on the lowest energy subband (greatly improved the calculation speed). The code will produce an output file with a name such as **fermions_singleband_checkerboardlattice_n_4_x_4_y_3_v_0.000000_t1_1.000000_t2_0.292893_tpp_0.207107_gx_0.000000_gy_0.000000.dat**. It will look like

```text
# eigenvalues
# kx ky E
0 0 0.22722637070671
0 0 0.24524464039421
0 0 0.2498205060043
0 0 0.27409716067925
0 0 0.44636044726849
0 0 0.47911382030328
0 0 0.51799401735293
0 0 0.52053705316647
0 0 0.57405393558818
0 0 0.59094271653235
0 0 0.62156107864437
0 0 0.63528364243047
0 0 0.68277333254545
0 0 0.70081134893839
```

The output is an ASCII column formatted text file. The first column is the total kx momentum, the second column is the total ky momentum and the third column contains the sorted energy eigenvalues.

The default behavior is to compute the spectrum in each momentum sector. One can restrict to a fixed kx and/or ky momentum sector using the --only-kx and --only-ky options.

## Interaction

Both two body and three body interactions are available in the single band approximation. File:checkerboardthreebody.png *(wiki page File:checkerboardthreebody.png)*


## One-body hamiltonian

Most of the programs related to fractional Chern insulators or fractional topological insulators provide a way to retrieve information about the band stucture or the one-body mode hamiltonian. First, one has to turn on the --singleparticle-spectrum option. This will store the band stucture in an ASCII text file. The system sizes that are doable in this mode, are much larger than those reachable in the interacting case. The full information about the one body model, i.e. the energies and the eigenstates, can be dumped in a binary file (using the --export-onebody option) or in a text file (using the --export-onebodytext option). The default output file name is similar to the one of the interacting spectrum with a *_tightbinding.dat* extension. It can be replaced by a user define file name through the --export-onebodyname .

## Memory usage

FQHECheckerboardLatticeModel uses as much memory as possible to store the matrix elements. The amount of memory (in MBytes) can be tune using the --memory option. If this option is set to 0, then no matrix elements are stored which will greatly increase the computation time. You should never provide more memory than the physical amount of RAM on your computer. To know if you have provided enough memory, just look at a typical run 

*\$PATHTODIAGHAM/build/FTI/src/Programs/FCI/FCICheckerboardLatticeModel -p 8 -x 6 -y 4 --single-band --flat-band -n 1 --show-itertime --memory 100*


---
```text
start
nbr interaction = 3225676
reduction factor=1
------------------------------------------------------------------

time = 0.426645
fast = 62Mb start
------------------------------------------------------------------

time = 0.526953
----------------------------------------------------------------
 Hilbert space dimension = 30709
Run Lanczos Algorithm
0.93090256321753 0.93090256321753 1 (0.069106 s for step 5)
0.78123005742037 0.78123005742037 0.16078213951828 (0.069303 s for step 6)
0.66740307853967 0.66740307853967 0.14570225223611 (0.069505 s for step 7)
```


A reduction factor of 1 means that the program has enough memory to store the complete hamiltonian (here it requires 62 MBytes). If we now set --memory 50, we get 

---
```text
start
nbr interaction = 3225676
reduction factor=2
------------------------------------------------------------------

time = 0.427688
fast = 31Mb start
------------------------------------------------------------------

time = 0.269822
----------------------------------------------------------------
 Hilbert space dimension = 30709
Run Lanczos Algorithm
0.93090256321753 0.93090256321753 1 (0.318085 s for step 5)
0.78123005742037 0.78123005742037 0.16078213951828 (0.322451 s for step 6)
0.66740307853968 0.66740307853968 0.14570225223611 (0.317772 s for step 7)
```

The reduction factor is now 2, which implies that only half of the hamiltonian can be stored. One can see that the amount of time per iteration is almost multiply by 5. As a thumb rule, there is almost no gain between a reduction factor of 4 and disabling the hamiltonian storage.

Notice that the total memory usage during a run depends on several factors beyond the (partial) storage of matrix elements. In particular, the options that are set for the Lanczos algorithm may greatly increase the memory usage.


## Running the code in parallel mode

The code is fully parallel. On a machine with several cores or processors, you can turn on parallelization using the -S option and setting --processors to the total number of cores you can use.

If you want to run the code across multiple nodes, DiagHam has to configure and compile with the MPI *(wiki page MPI)* option enabled.

## Real space

FCICheckerboardLatticeModel can also work in real space, keeping the two. It supports both spinless fermions and boson (with or without the hardcore constraint), and with or without translations. For example, 

*\$PATHTODIAGHAM/build/FTI/src/Programs/FCI/FCICheckerboardLatticeModel -p 6 -x 4 -y 3 --use-lapack --boson --u-potential 5.0 --real-space -n 3 --block-lanczos --block-size 3 --force-reorthogonalize --eigenstate --memory 4000 -S --processors 4 --gutzwiller*

will consider a checkerboard lattice with 4x3 unit cells (24 sites) with 6 bosons using the hardcore constraint (--gutzwiller option). The spectrum bosons_realspace_gutzwiller_checkerboardlattice_n_6_ns_24_x_4_y_3_u_5.000000_v_0.000000_t1_1.000000_t2_0.292893_tpp_0.207107_gx_0.000000_gy_0.000000.dat should look like

```text
# eigenvalues
# kx ky  E
0 0 -12.223408402675
0 0 -11.821248147999
0 0 -11.801980073472
0 1 -11.85234570983
0 1 -11.757051023166
0 1 -11.713935150874
0 2 -11.85234570983
0 2 -11.757051023166
0 2 -11.713935150874
1 0 -11.884877664428
1 0 -11.735772253883
1 0 -11.723202200313
1 1 -11.884289796641
1 1 -11.759313765328
1 1 -11.693554456905
1 2 -11.884289796641
1 2 -11.759313765328
1 2 -11.693554456905
2 0 -12.226313403825
2 0 -11.864887534923
2 0 -11.745326608006
...
```

The (almost) twofold degeneracy is expected for the Laughlin $\nu=\frac{1}{2}$ state. Adding the --no-translation disables the translations and will give an energy spectrum bosons_realspace_gutzwiller_notranslation_checkerboardlattice_n_6_ns_24_x_4_y_3_u_5.000000_v_0.000000_t1_1.000000_t2_0.292893_tpp_0.207107_gx_0.000000_gy_0.000000.dat

```text
# eigenvalues
# kx ky  E
0 0 -12.226313403825
0 0 -12.223408402675
0 0 -11.884877664428
```
