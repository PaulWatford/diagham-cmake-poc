# Open finding U47: HubbardKitaevHeisenbergModel's J1 and J2 multiply Pauli matrices, four times the S.S of its wiki page

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U47), found 2026-10-26 while writing the program's first independent tests; r4493 sources, GCC 15. The program's output is byte-identical on a pristine autotools build of r4493.

## What happens

The program's wiki page writes the spin terms as J_1 sum S_i.S_j + J_2 sum [S_i^a S_j^a - S_i^b S_j^b - S_i^c S_j^c] (a the type of the bond). The program's `--j1` and `--j2` multiply Pauli matrices instead, J_1 sum sigma_i.sigma_j + J_2 sum [sigma^a sigma^a - sigma^b sigma^b - sigma^c sigma^c], which is four times the wiki's Hamiltonian. The help text ("strength of the nearest neighbor Heisenberg interaction") does not say which.

Example, the single hexagon of the wiki page (geometry file `0 1 2 / 0 2 1 / 1 3 0 / 2 4 0 / 3 5 1 / 4 5 2`) at half filling, where the Gutzwiller projection leaves a pure spin model:

```
HubbardKitaevHeisenbergModel -p 6 --gutzwiller --geometry-file hexagon.dat --geometry-name hexg --isotropic-t 0 --anisotropic-t 0 --j1 0.3 --j2 0.5 --full-diag 5000 --use-lapack
```

The 64 eigenvalues agree to 1e-13 with an independent exact diagonalisation of the Pauli-matrix form (`tests/oracles/core_tools.py`), and differ by up to 3.4 from the S.S form of the wiki. The bond-dependent signs and the bond types (0 = x, 1 = y, 2 = z) are as the wiki describes.

## Fix direction

Either is enough: write the wiki formula with Pauli matrices (or say that J_1 and J_2 are four times the S.S couplings), or say so in the help text of `--j1` and `--j2`.

## Consequence for the tests

`physics.hubbard.kitaev_heisenberg.hexagon_spin_j1_0.3_j2_0.5` takes the factor of four from the program, as stated in `tests/oracles/core_tools.py`; everything else in it is independent. No `open-finding` test: the code is consistent, only its description differs.
