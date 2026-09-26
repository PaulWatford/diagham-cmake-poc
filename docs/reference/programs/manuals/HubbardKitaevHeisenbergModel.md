# HubbardKitaevHeisenbergModel — manual

Source: DiagHam wiki page `HubbardKitaevHeisenbergModel`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [HubbardKitaevHeisenbergModel](../HubbardModels/HubbardKitaevHeisenbergModel.md)
$$
H_{SH}= -t \sum_{ <i,j>} [c_{i\sigma}^{\dagger}d_{j\sigma}+d_{j\sigma}^{\dagger}c_{i\sigma}] 
-t'\sum_{ <i,j> } [c_{i\sigma}^{\dagger}d_{j\sigma'}\sigma_{\sigma\sigma'}^{\alpha}+d_{i\sigma}^{\dagger}c_{j\sigma'}\sigma_{\sigma\sigma'}^{\alpha}] 
+J_1 \sum_{ <i,j>} \vec{S}_i\cdot\vec{S}_j 
+J_2 \sum_{ <i,j>} [S_i^{\alpha}S_j^{\alpha}-S_i^{\beta}S_j^{\beta}-S_i^{\gamma}S_j^{\gamma}] 
+U \sum_{i} n_{i\uparrow} n_{i\downarrow}
$$


*\$PATHTODIAGHAM/build/FTI/src/Programs/HubbardModels/HubbardKitaevHeisenbergModel -p 6 -x 6 --gutzwiller --j1 0.49 --j2 0.51 --use-lapack --stripe --eigenstate -n 1*

Exact diagonalisation program for a system of spin 1/2 fermions on the Honeycomb lattice Kitaev-Heisenberg model.
The number of particles and the number of sites are set through the -p option and -x respectively. The Gutzwiller projection is turned on via the --gutzwiller option. In that case, a reduced Hilbert space forbidding the occupation of a site by two fermions is used. The parameters of the Hamiltonian are set using :

- --isotropic-t for the spin isotropic nearest neighbor hopping amplitude i.e. $t$
- --anisotropic-t for the spin anisotropic nearest neighbor hopping amplitude i.e. $t'$
- --j1 for the neareast neighbor Heisenberg interaction i.e. $J_1$
- --j2 for the neareast neighbor anisotropic spin-spin interaction i.e. $J_2$
- --u-potential for the on-site Hubbard potential i.e. $U$. This option is useless if the Gutzwiller projection is turned on.

For the parallelization options, they are similiar to those of many diagonalization codes such as FCICheckerboardLatticeModel *(wiki page FCICheckerboardLatticeModel#Running_the_code_in_parallel_mode)*.

## Open parallelogram configuration
The *--open* option sets the geometry of the system to be a parallelogram with open boundary conditions (zigzag alignment in both x and y directions). The options *--nbrsites-x* and *--nbrsites-y* set the number of hexagons in the x and y directions.

The figure below gives an example of the bond types and sites indices for an open system of length 3 and of width 2. The bond type is given by the following color code: red ($x$), green ($y$), blue ($z$)

300px *(wiki page File:honeycombopen.png)*

## Stripe configuration
The --stripe option sets the geometry of the system to be a stripe: the boundary conditions are periodic in the x direction, and open in the y direction.
- The option *--nbrsites-x* sets the number of unit cells (number or hexagons in the direction, or length of the stripe).
- The option *--nbrsites-y* sets the number of hexagons in the y direction (width of the stripe).

The figure below gives an example of the bond types and sites indices for a stripe of length 3 and of width 2. The bond type is given by the following color code: red ($x$), green ($y$), blue ($z$)


300px *(wiki page File:Honeycombstripe.png)*

## Torus configuration
The --torus option sets the geometry of the system to a torus (periodic boundary conditions in both directions).
- The number of unit cells in the x direction is specified using the option --nbrsites-x
- The number of unit cells in the y direction is specified using the option --nbrsites-y
The total number of sites (set using the option --nbr-sites) should be equal to twice the product of the two. 

For instance:
*\$PATHTODIAGHAM/build/FTI/src/Programs/HubbardModels/HubbardKitaevHeisenbergModel -p 12 -x 12 --use-lapack  --gutzwiller --j1 0.3 --j2 0.25 --isotropic-t 1.0 --anisotropic-t 0.5 --full-diag 10000 --torus --nbrsites-x 3 --nbrsites-y 2*

The sites are labeled in a way that is consistent with translational invariance. For instance, for 3 unit cells in the x directions and 2 unit cells in the y direction:

250px *(wiki page File:honeycombtorus.png)*

## Generic geometry
When the --stripe and --torus options are not activated, one can specify the geometry of the system by providing a three-column text file through the option --geometry-file. The two first columns give two site indices, and the third column gives the nature of the link between them (0 if link $x$, 1 if link $y$, 2 if link $z$). For instance, the following input geometry file corresponds to a single hexagon:

0 1 2  
0 2 1  
1 3 0  
2 4 0  
3 5 1  
4 5 2  

The name of the output files (spectrum and vectors) is specified using the option --geometry-name. For instance 

*\$PATHTODIAGHAM/build/FTI/src/Programs/HubbardModels/HubbardKitaevHeisenbergModel -p 4 -x 6 --j1 1 --j2 1 --isotropic-t 0 --anisotropic-t 0 --use-lapack  --eigenstate -n 1 --geometry-file hexagon.dat --geometry-name hexagon*

will generate a file fermions_kitaev_heisenberg_hexagon_n_4_x_6_t_0_tK_0_j1_1_j2_1.dat

## Periodic boundary conditions
The current version can resolve momentum in one or two directions for a system with periodic boundary conditions. We define the real space vectors **a**<sub>1</sub>, **a**<sub>2</sub> and reciprocal space vectors **a' ** <sub>1</sub>, **a' **<sub>2</sub> in the following figure. When the periodic boundary conditions are implemented, a translation by a vector N<sub>x</sub> **a**<sub>1</sub> or N<sub>y</sub> **a**<sub>2</sub> leaves the system invariant (where N<sub>x</sub> and N<sub>y</sub> are defined by the *--nbrsites-x* and *--nbrsites-y* options).

300px *(wiki page File:honeycombPeriodicBoundaryConditions.png)*

For one direction, this is enabled with the --xperiodic-boundary option. The periodicity is expressed using the site index convention. For example, in the stripe geometry with boundary conditions, the periodicity has to be set to 4 since the site $n$ and $n+4$ are related by translation. This can be done via the --max-xmomentum option which sets the number of possible momentum values in the x direction. For the stripe geometry, this should be the numbers of sites divided by 4. A typical usage would be 

*\$PATHTODIAGHAM/build/FTI/src/Programs/HubbardModels/HubbardKitaevHeisenbergModel -p 10 -x 12 --gutzwiller --stripe --isotropic-t 1.0 --anisotropic-t 1.0 --j1 1.0 --j2 1.0 --use-lapack --memory 500 -n 3 --block-lanczos --block-size 3 --force-reorthogonalize --xperiodic-boundary  --max-xmomentum 3 --nbrsites-x 3 --nbrsites-y 1*

A specific momentum sector can be computed by setting to the --x-momentum option to the corresponding momentum sector. Note that these options are not restricted to the stripe geometry. Any geometry where the periodity can be written as a relation between site $n$ and $n+T_x$ where $T_x$ is the periodicity (that has to be provided through the --max-xmomentum option) can use these options.

For two (decoupled) directions, one needs to enable the --2dperiodic-boundaries. In that case both --max-xmomentum and --max-ymomentum options have to be set. Note that the site numbering should satisfy the above mentioned  condition for the x direction. For the y direction the site numbering should be such that two sites related by a single translation in the y direction satisfy $n$ and $n+T_Y$, with $N_Y T_Y=T_X$, $N_Y$ being the number of unit cells in the y direction. Here is an example for a torus having 3x3 hexagons and periodic boundary conditions

*\$PATHTODIAGHAM/build/FTI/src/Programs/HubbardModels/HubbardKitaevHeisenbergModel -p 18 -x 18  --use-lapack --geometry-file torus_nx_3_ny_3.dat --geometry-name torus --gutzwiller --j1 0.3 --j2 0.25 --isotropic-t 1.0 --anisotropic-t 0.5   --full-diag 7000 --max-ymomentum 3 --max-xmomentum 3 --2dperiodic-boundaries  --memory 6000*

If using the --torus option instead of specifying the geometry by an external file, the values of --nbrsites-x and --nbrsites-y are automatically set to those of --max-xmomentum and --max-ymomentum.
