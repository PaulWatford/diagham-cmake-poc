# GenericPeriodicSpinChain — manual

Source: DiagHam wiki page `GenericPeriodicSpinChain`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [GenericPeriodicSpinChain](../Spin/GenericPeriodicSpinChain.md)
GenericPeriodicSpinChain diagonalizes the Heisenberg spin chain with periodic boundary conditions. For example

\$PATHTODIAGHAM/build/Spin/src/Programs/GenericOpenSpinChain -s 1 -p 4 --use-lapack

diagonalizes a 4 sites (-p option) spin 1/2 (the -s option should be set two twice the spin per site) with coupling constant J=1 . The code only looks at positive Sz sector and all momentum sectors. The spectrum is stored in a file whose name looks like spin_1_2_periodicchain_n_4.dat with the following content


```text
#  periodic spin 1/2 chain with 4 sites
# 2Sz K  E
0 0 -2
0 0 1
0 1 0
0 2 -1
0 2 7.4983036091107e-33
0 3 0
2 0 1
2 1 -6.1230317691119e-17
2 2 -1
2 3 1.8369095307336e-16
4 0 1
```
