# FTIGenericInteractionFromFileTwoBands — manual

Source: DiagHam wiki page `FTIGenericInteractionFromFileTwoBands`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FTIGenericInteractionFromFileTwoBands](../FTI/FTIGenericInteractionFromFileTwoBands.md)
FTIGenericInteractionFromFileTwoBands is a generic program for strongly interacting systems written in momentum space with two bands and potentially with a spin and/or a valley degrees of freedom. It requires external files providing the interaction matrix elements and optional one-body terms. If spin is added, the interaction is assumed to be SU(2) symmetric, namely the interaction is spin independent. Most of the running options are similar to those of [FQHECheckerboardLatticeModel](FCICheckerboardLatticeModel.md). 

 
By default, the code does not include any spin or valley degrees of freedom. Each can be included by adding the running options --add-spin and --add-valley

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIGenericInteractionFromFileTwoBands --use-lapack -x 4 -y 2 -p 8 --interaction-file interaction.txt --interaction-name singled_gated_coulomb --full-diag 1000 -u 1.0 --memory 1400000 --show-itertime -S --processors 8 --add-valley --add-spin --use-valleyspin --sz1-value 0 --sz2-value 0 --pz-value 0 --real-interaction --disable-pzsymmetry --disable-szsymmetry  --lanczos-precision 1e-12 --singleparticle-file onebody.txt --only-kx 0 --only-ky 0 --eigenstate -n 3 --block-lanczos --block-size 3 --force-reorthogonalize*

The formats of the interaction file interaction.txt and one-body term file onebody.txt are described below. The name of the interaction, as it should appear in the output file name is set by the option --interaction-name (here singled_gated_coulomb). In this example, we compute 3 eigenstates using a safe version of the block Lanczos algorithm. The system has 8 particles (set by -p) on a 4x2 discretized Brillouin zone (-x 4 -y 2), and we solely focus on the kx=ky=0 momentum sectors using the options --only-kx 0 --only-ky 0. If one of these two or both options are absent, then all the corresponding momentum sectors are computed. 

The number of particles per valley is set by --pz-value, Pz being the twice valley polarization, it should be understood as the difference between the number of particles in valley eta=+1 and the number of particles in the valley eta=-1.

The spin projection along z (aka Sz) can be fixed in two different manners. One can use the --sz-value which fixes twice the total Sz and --ez-value which fixes twice Ez defined as Ez=(N_{up, eta=+1}+N_{down, eta=-1}-N_{down, eta=+1}-N_{up, eta=-1})/2. Sometimes, it is more convenient to provide the Sz value per valley. In that case, we can use the --use-valleyspin  and provide the the twice Sz for the two valleys through --sz1-value (for eta=-1) and --sz2-value (for eta=+1).

Add the moment, spin flip and valley flip discrete symmetries are not available. As such, the two options --disable-pzsymmetry --disable-szsymmetry should always be added.

## Format of the two-body term:

Here we use the notation m_i for the band index (either 0 or 1), eta_i for the valley index (either -1 or +1) and kx_i, ky_i for the momentum index (from 0 to Nx-1 or Ny-1). When taking into account the valley degree of freedom for a real hamiltonian, the ascii file providing the interaction matrix elements when valley is added should look like  


```text
# m_1 eta_1 kx_1 ky_1 m_2 eta_2 kx_2 ky_2 m_3 eta_3 kx_3 ky_3 m_4 eta_4 kx_4 ky_4 matrix_element
0 -1 2 1 0 -1 2 1 0 -1 0 0 0 -1 0 0 -1.2029024074361752
0 -1 2 1 0  1 2 1 0  1 0 0 0 -1 0 0 -1.2136487227435184
```


while the complex version should read

```text
# m_1 eta_1 kx_1 ky_1 m_2 eta_2 kx_2 ky_2 m_3 eta_3 kx_3 ky_3 m_4 eta_4 kx_4 ky_4 matrix_element
1  1 1 1 0  1 3 1 1  1 0 0 1  1 0 0 (-1.3725397342421886e-13,-0.23684972098099083)
1 -1 1 1 1 -1 3 1 1 -1 0 0 1 -1 0 0 (0.056951985858466,-4.8304189142295314e-14)
```

Note that the code checks if the U(1) valley symmetry is conserved and fails in case it is broken. By default, the interaction is assumed to be complex. Use the --real-interaction option to force it to be real (greatly improve speed and memory usage). If the valley degree of freedom is not used, the interaction file should not include the four eta columns.


## Format of the one-body term:

By default the one-body term assumes conservation of all quantum numbers (i.e. a fully diagonal term). If the option --full-singleparticle is selected, then only the momentum and valley degree of freedom should be conserved. In that case, each matrix element is defined through 5 indices, with m and n standing for the band indices. The ascii file for the one-body term should look like

```text
# kx    ky      m       n       eta     one_body_term
1 1 0 0 1 (-209.8414774464715,1.4387741762747747e-15)
1 1 0 1 -1 (-1.7174948001364406e-13,-1.1811546809422595)
...
```

Matrix elements can be either complex or real, but the choice should be consistent over the whole file. If complex entries are provided, the running option --complex-singlebody should be added to the command line. 

If the valley degree of freedom is not required, the file is simplified by dropping the eta column

```text
# kx    ky      m       n       one_body_term
1 1 0 0 (-209.8414774464715,1.4387741762747747e-15)
1 1 1 0 (1.7174948001364406e-13,1.1811546809422595)
1 1 0 1 (1.7200543150922606e-13,-1.1811546809422595)
1 1 1 1 (-209.8414774464714,-1.912571380449884e-16)
...
```


## Hilbert space generation

Generating the Hilbert space might take a substantial amount of time, especially when involving options such as --max-band1, --min-band0, --allowed-orbitals, ... When computing a phase diagram, Hilbert spaces are usually identical for the different hamiltonian parameters. In that case, it is suggested to use the --hilbert-directory. You should provide a directory as a parameter to the option where the Hilbert spaces will be stored. During the first run, all Hilbert spaces will be generated and saved there, while any other run will directly read the Hilbert spaces from this folder. As an example on a laptop computing 3 momentum sectors, we get

```text
''time $PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIGenericInteractionFromFileTwoBands --full-diag 2000 --lanczos-precision 1e-11 --use-lapack -m 1500000 --show-itertime -S --processors 32 --interaction-file EDInt_twoband.txt  -x 4 -y 6 -p 16 -n 2 --block-lanczos --block-size 2 --force-reorthogonalize -u 1.0 --full-singleparticle --complex-singlebody --interaction-name test --max-band0 -1 --max-band1 0 --eigenstate --selected-sectors ./sectorseclection.by._n24 --3band-convention --singleparticle-file EDSP_twoband.txt --hilbert-directory .''
One-body diagonalization done in 0.00121 s
(kx=0,ky=0) : 
...
time = 3.042403
----------------------------------------------------------------
------------------------------------
Total dim=91957

real	8m39.690s
user	9m26.203s
sys	0m8.092s
```


Three files are written in the current directory that look like   *fermions_klattice_su2_n_16_nx_4_ny_6_kx_0_ky_0_minband1_0_minband1_0_maxband0_24_maxband1_0.hil*.  A second run of the same code is now much faster

```text
time = 3.318556
----------------------------------------------------------------
------------------------------------
Total dim=91957

real	0m12.292s
user	1m0.289s
sys	0m6.068s
```


There are a few caveats to keep in mind:
- If you run a phase diagram using a script potentially on a cluster, be sure that the folder is shared by all nodes. Also to avoid potential conflict when writing the Hilbert space on disk, make sure a single code is actually performing that operation (e.g. by running a single point first).
- If you use the option --allowed-orbitals, the Hilbert space file names will have the same name irrespective of which orbitals are conserved. So any time  you use a new orbital filtering either wipe out the previous files or make sure that you do not use the same folder.

## Options:


- To test the matrix elements have been properly implemented (both two-body and one-body), the --test-hermitian option checks that the hamiltonian is indeed hermitian. Note that the accuracy of the test can be tuned through the --testhermitian-error option.

- If the band occupation is conserved, the --conserve-bandoccuption can be added to speed up the calculations and decrease the memory usage.

- If only some quantum number sectors shall be computed, these sectors can be provided through a text file whose name should be set by the option --selected-sectors. This ascii column-formatted text format depends on the options that have been selected. 

   - Without spin and valley degrees of freedom, only the total momentum per direction should be provided

```text
# kx ky 
0 1
1 0
1 1
2 0
2 1
```

Any additional column will be discarded. if both spin and valley are used, then kx, ky, 2Sz, 2Pz and 2Ez (in that specific order) should be provided as a 5 column text.


   - With valley but without spin, you can use three columns for kx, ky and Pz respectively, i.e.,

```text
# kx ky Pz
0 1 6
1 0 6
1 1 6
2 0 6
2 1 6
```


   - When --conserve-bandoccuption is used, the structure of the selected sector file is slightly different. For example in the 2 band model without valley nor spin degrees of freedom, 10 particles on a 5x6 momentum lattice, we consider the following file 

```text
# kx ky n1 n2
0 0 10 0
0 1 10 0
0 2 10 0
0 3 10 0
1 0 10 0
1 1 10 0
1 2 10 0
1 3 10 0
1 4 10 0
1 5 10 0
2 0 10 0
2 1 10 0
2 2 10 0
2 3 10 0
2 4 10 0
2 5 10 0
3 0 10 0
3 1 10 0
3 2 10 0
3 3 10 0
```


There we only focus on the momentum sectors that are not related by the inversion symmetry and on the case where all 10 particles are in band 1

   - It is possible to change the number of eigenvalues/eigenvectors computed per quantum number sector. To do so, add a column containing the number of requested eigenvalues in the and include the option *--variable-nbreigenvalues* to the command line. Note that when using block Lanczos, you should set the block size to the maximal number of eigenvalues to compute. Still, the code will make sure the block size does not exceed the number of eigenvalues to compute.
