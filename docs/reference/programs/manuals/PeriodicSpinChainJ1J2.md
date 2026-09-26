# PeriodicSpinChainJ1J2 — manual

Source: DiagHam wiki page `PeriodicSpinChainJ1J2`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [PeriodicSpinChainJ1J2](../Spin/PeriodicSpinChainJ1J2.md)
PeriodicSpinChainJ1J2 diagonalizes the $J_1-J_2$ Hamiltonian using periodic boundary conditions. The $J_1-J_2$ Hamiltonian is defined by

$H=J_1\sum_{i=0}^{N-1} \vec{S}_{i}\cdot\vec{S}_{i+1} \; + \; J_2\sum_{i=0}^{N-1} \vec{S}_{i}\cdot\vec{S}_{i+2}$

A typical usage is

*\$PATHTODIAGHAM/build/Spin/src/Programs/PeriodicSpinChainJ1J2 -p 6 --full-diag 10000 --use-lapack*

By default, $J_1=1$ and $J_2=\frac{1}{2}$. This can be tuned using the options --j1 and --j2. Most of the options are identical to [PeriodicSpinChainAKLT](PeriodicSpinChainAKLT.md).
