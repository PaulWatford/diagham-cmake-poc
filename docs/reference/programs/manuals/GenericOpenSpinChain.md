# GenericOpenSpinChain — manual

Source: DiagHam wiki page `GenericOpenSpinChain`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [GenericOpenSpinChain](../Spin/GenericOpenSpinChain.md)
[GenericOpenSpinChain](GenericOpenSpinChain.md) diagonalizes the Hamiltonian

$H=\sum_{i=1}^{N-1} \left( J S_{x,i}S_{x,i+1} + J S_{y,i}S_{y,i+1}  + (J + \Delta J_z)S_{z,i}S_{z,i+1} \right) + \sum_{i=1}^{N} h_{z,i} S_{z,i}$

Note that this Hamiltonian has the total spin projection along z as a good quantum number. It shares many options with [FullGenericOpenSpinChain](FullGenericOpenSpinChain.md).

At typical usage is 
*\$PATHTODIAGHAM/build/Spin/src/Programs/GenericOpenSpinChain -s 1 -p 4 --use-lapack*

diagonalizes a 4 sites (-p option) spin 1/2 (the -s option should be set two twice the spin per site) with coupling constant J=1 . The code only looks at positive Sz sector. The spectrum is stored in a file whose name looks like spin_1_2_openchain_n_4.dat with the following content 


```text
#  open spin 1/2 chain with 4 sites
# 2Sz E
0 -1.6160254037844
0 -0.95710678118655
0 -0.25
0 0.11602540378444
0 0.45710678118655
0 0.75
2 -0.95710678118655
2 -0.25
2 0.45710678118655
2 0.75
4 0.75
```

$\Delta J_z$ can be set by using the --djz-value option.
