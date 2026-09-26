# PeriodicSpinChainGeneralizedAKLT — manual

Source: DiagHam wiki page `PeriodicSpinChainGeneralizedAKLT`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [PeriodicSpinChainGeneralizedAKLT](../Spin/PeriodicSpinChainGeneralizedAKLT.md)
[PeriodicSpinChainGeneralizedAKLT](PeriodicSpinChainGeneralizedAKLT.md) performs the exact diagonalization of the spin 2 analogue of the AKLT Hamiltonian. We denote $P^{(s)}_{i,j}$ the projector on the total spin $s$ sector for the two spins located on sites $i$ and $j$. The generalized AKLT hamiltonian for a chain of \$N\$ spins with periodic boundary conditions reads

$H=\sum_{i=0}^{N-1} \lambda^{(3)} P^{(3)}_{i,i+1} + \lambda^{(4)}  P^{(4)}_{i,i+1}$

[PeriodicSpinChainGeneralizedAKLT](PeriodicSpinChainGeneralizedAKLT.md) supports the same features than [PeriodicSpinChainAKLT](PeriodicSpinChainAKLT.md), including all discrete symmetries, real hamiltonians / vectors when available, ... The two coefficients $\lambda^{(3)}$ and $\lambda^{(4)}$ are set by the options --p3-factor and --p4-factor . Their values appear in the output file names. For example 

*\$PATHTODIAGHAM/build/Spin/src/Programs/PeriodicSpinChainGeneralizedAKLT -p 2 --p3-factor 1.0 --p4-factor 1.0 --use-lapack*

will generate a file spin_2_periodicaklt_p3_1.000000_p4_1.000000_n_2.dat containing the spectrum for these parameters. The file content should look like

```text
#  periodic spin 2 chain with 2 sites 
# 2Sz K SzSym InvSym  E
0 0 1 1  0
0 0 1 1  2.0816681711722e-17
0 0 1 1  2
0 1 -1 1  0
0 1 -1 1  2
2 0 0 1 1.1102230246252e-16
2 0 0 1 2
2 1 0 1 -1.6653345369377e-16
2 1 0 1 2
4 0 0 1 -1.1102230246252e-16
4 0 0 1 2
4 1 0 1 2
6 0 0 1 2
6 1 0 1 2
8 0 0 1 2
```
