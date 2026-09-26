# PeriodicSpinChainAKLT — manual

Source: DiagHam wiki page `PeriodicSpinChainAKLT`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [PeriodicSpinChainAKLT](../Spin/PeriodicSpinChainAKLT.md)
PeriodicSpinChainAKLT diagonalize the AKLT Hamiltonian using periodic boundary conditions. The AKLT Hamiltonian is defined by

$H=\sum_{i=0}^{N-1} \left( \vec{S}_{i}\cdot\vec{S}_{i+1} + \frac{1}{3} \left(\vec{S}_{i}\cdot\vec{S}_{i+1}\right)^2 \right)$

## Basic usage

A typical usage is

*\$PATHTODIAGHAM/build/Spin/src/Programs/PeriodicSpinChainAKLT -p 6 --full-diag 10000 --use-lapack*

PeriodicSpinChainAKLT computes all momentum sectors and all $S_z$ sectors. To focus on a single momentum sector, we just have to include the option --momentum x where x is the momentum sector (in $2\pi/N$ units and should be between 0 and  $N-1$). Similarly, we can restrict to one or several $S_z$ sectors with the options --initial-sz (twice the initial $S_z$ sector that should be computed) and --nbr-sz (the actual number of $S_z$ sectors that should be computed).

By default the code uses both the parity and inversion symmetry when possible. To disable these symmetries, we can use the options --disable-szsymmetry and --disable-inversionsymmetry. Note that for the momentum sector $k=0$ and  $k=\pi$, PeriodicSpinChainAKLT uses a real hamiltonian to improve the performance. If you want to produce complex eigenstates, you can disable this feature by adding the --disable-realhamiltonian option. For example, 

*\$PATHTODIAGHAM/build/Spin/src/Programs/PeriodicSpinChainAKLT -p 6 --all-eigenstates --full-diag 10000 --use-lapack  --momentum 0 --initial-sz 8 --nbr-sz 1  --disable-szsymmetry  --disable-inversionsymmetry*

will diagonalize the AKLT Hamiltonian in the momentum sector $k=0$ without using the $S_z$ parity or the inversion symmetry. This example also computes all the eigenstates. To extract an individual eigenstates out this binary matrix, you can use [MatrixExtractColumns](MatrixExtractColumns.md).

## Advanced  options

[PeriodicSpinChainAKLT](PeriodicSpinChainAKLT.md) allows to set the coefficients freely in front of the bilinear / biquadratic terms. Thus we can manage the following Hamiltonian

$H=\sum_{i=0}^{N-1}  \left( \lambda_1 \vec{S}_{i}\cdot\vec{S}_{i+1} + \lambda_2 \left(\vec{S}_{i}\cdot\vec{S}_{i+1}\right)^2 \right)$

$\lambda_1$ is set by the --linear-factor option while $\lambda_2$ is set by the --quadratic-factor option. For example, running 

*\$PATHTODIAGHAM/build/Spin/src/Programs/PeriodicSpinChainAKLT -p 8 --use-lapack --linear-factor 3.0  --quadratic-factor 1.0*

will generate a spectrum file spin_1_periodicaklt_linear_3.000000_quadratic_1.000000_n_8.dat. Note that this mode is only activated when $\lambda_1$ is different from 1. If you want to only tune $\lambda_2$, you should rely on the other option --additional-quadratic. In that case, the hamiltonian is 


$H=\sum_{i=0}^{N-1}  \left(\vec{S}_{i}\cdot\vec{S}_{i+1} + \left(\lambda + \frac{1}{3}\right) \left(\vec{S}_{i}\cdot\vec{S}_{i+1}\right)^2 \right)$

where  $\lambda$ is set by the --additional-quadratic option and the spectrum output file name should look like spin_1_periodicaklt_quadratic_-0.300000_n_8.dat if we run a command such as

*\$PATHTODIAGHAM/build/Spin/src/Programs/PeriodicSpinChainAKLT -p 8 --use-lapack  --additional-quadratic -0.3*
