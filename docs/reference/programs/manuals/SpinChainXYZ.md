# SpinChainXYZ — manual

Source: DiagHam wiki page `SpinChainXYZ`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [SpinChainXYZ](../Spin/SpinChainXYZ.md)
SpinChainXYZ diagonalizes the XYZ hamiltonian  :

$H_{\rm XYZ}=-\sum_{i=1}^{N-1} \left( J_x S_{x,i}S_{x,i+1} + J_y S_{y,i}S_{y,i+1}  + J_z S_{z,i}S_{z,i+1} \right) - J_x B \left(\prod_i \sigma_{z,i}\right) S_{x,N}S_{x,1} - J_y B \left(\prod_i \sigma_{z,i}\right) S_{y,N}S_{y,1}  - J_z |B| S_{z,N}S_{z,1}$

where $B$ is the boundary condition (can be either 0,1 or -1).
