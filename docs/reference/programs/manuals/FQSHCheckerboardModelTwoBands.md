# FQSHCheckerboardModelTwoBands — manual

Source: DiagHam wiki page `FQHEQuantumSpinHallCheckerboardModelTwoBands`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build). The wiki page was named `FQHEQuantumSpinHallCheckerboardModelTwoBands`; the program is now `FQSHCheckerboardModelTwoBands`.

Generated option reference: [FQSHCheckerboardModelTwoBands](../FTI/FQSHCheckerboardModelTwoBands.md)
FQHEQuantumSpinHallCheckerboardModelTwoBands provides an exact diagonalization code for FQSH based on two copies of the checkerboard lattice model. Most of the running options are similar to those of [FQHECheckerboardLatticeModel](FCICheckerboardLatticeModel.md). 

## Decoupled layers

In its simplest version, [FQHEQuantumSpinHallCheckerboardModelTwoBands](FQSHCheckerboardModelTwoBands.md) involved two decoupled copies of the checkerboard lattice model (spin up and down). The typical usage is

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FQSHCheckerboardModelTwoBands -p 4 -x 3 -y 2 --use-lapack --flat-band --project-fourbands --decoupled*

The eigenvalue output file looks like 

```text
# eigenvalues
# kx ky Sz E
0 0 -4 1.7063718059659
0 0 -4 2.1234062994611
0 0 -4 2.1434774759191
0 0 -2 0.80880078458526
0 0 -2 0.84834788180414
0 0 -2 0.85652103153087
0 0 -2 0.86641822256772
0 0 -2 0.915790768175
0 0 -2 1.0302577225926
0 0 -2 1.0356483557114
```

Compared to [FCICheckerboardLatticeModel](FCICheckerboardLatticeModel.md), there is an additional column (the third one) which is two times the Sz value (Sz being a good quantum number in the case of two decoupled copies).

The interaction can be tuned through three paramaters : --u-potential sets the repulsive nearest neighbor potential strength between two identical spins (note that this option is inactive in flat-band mode, as U sets the energy scale of the system), --v-potential is the repulsive on-site potential strength between opposite spins and --w-potential allows to change the repulsive nearest neighbor potential strength between opposite spins.

## Tilted periodic boundary conditions
Tilted periodic boundary conditions are available for the spin conserved (decoupled) mode. The typical usage is

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHETopInsulator/FQHEQuantumSpinHallCheckerboardModelTwoBands -p 6 -x 9 -y 1 --use-lapack --flat-band --project-fourbands --decoupled  --nx1 3 --ny1 0 --nx2 1 --ny2 3 --offset 3*


*nx1* and *ny1* (respectively *nx2* and *ny2*) set the coordinates of **T1** (respectively **T2**), the periodicity vectors (a translation by either of these vectors leaves the system invariant). They must verify the condition nx1*ny2 - nx2*ny1 = Ns where Ns is the total number of unit cells. Nx and Ny must be set (using the *-x* and *-y* options) to 

$Nx = Ns/GCD(Ns, nx2, ny2)$ 

$Ny = GCD(Ns, nx2, ny2)$


*--offset* is an integer that must verify the relations

$(offset*ny2 - ny1)  mod  Nx = 0$

$(nx1 - offset*nx2)  mod   Nx = 0$



For nx1 = Nx, ny2 = Ny, ny1 = nx2 = 0, and offset = 0, one recovers the non-tilted periodic boundary conditions.


In tilted mode, one can also look at a time-reversal breaking system of two Checkerboard Chern insulator models. The *--break-timereversal* option activates the bilayer mode.

## Coupled layers

The coupled layer model is defined through the following hamiltonian


$$
H=\left( {\begin{array}{cc}
 H_0(k) & \Delta \; \sigma_y^*  \\
 \Delta^* \; \sigma_y &  H_0^*(-k) \\
 \end{array} } \right)
$$


where $H_0(k)$ is the one body hamiltonian of the checkerboard lattice and 
$$
\sigma_y=\left( {\begin{array}{cc}
 0 & i  \\
 -i &  0 \\
 \end{array} } \right)
$$
.

The code can be run using only the two lowest energy bands or the full four bands using the --four-bands option. The coupling parameter $\Delta$ is set through the two options --mixing-norm $\left(\left|\Delta\right|\right)$ and --mixing-arg $\left(arg(\Delta) / 2 \pi\right)$. The typical usage is 


*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHETopInsulator/FQHEQuantumSpinHallCheckerboardModelTwoBands -p 4 -x 3 -y 2 --flat-band --use-lapack*

The eigenvalue output file looks like 

```text
# eigenvalues
# kx ky  E
0 0 0.27058366707565
0 0 0.27271072155586
0 0 0.38540586407147
0 0 0.53925025231446
0 0 0.58915542718289
0 0 0.61372286112607
0 0 0.70253212436423
0 0 0.74470356767789
```


The --four-bands is much more demanding than the default two band mode. As a safety check, one can use the combination of --four-bands and --project-fourbands options to compare the result of the four band calculation projected onto the two lower band model.
