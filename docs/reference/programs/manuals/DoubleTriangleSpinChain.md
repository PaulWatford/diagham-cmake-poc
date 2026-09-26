# DoubleTriangleSpinChain — manual

Source: DiagHam wiki page `DoubleTriangleSpinChain`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [DoubleTriangleSpinChain](../Spin/DoubleTriangleSpinChain.md)
DoubleTriangleSpinChain diagonalizes the Heisenberg spin chain with periodic boundary conditions and second nearest neighbor coupling. In the following, J1 denotes the nearest neighbor coupling constant and J2 the second nearest neighbor coupling constant. For example, 

*\$PATHTODIAGHAM/build/Spin/src/Programs/DoubleTriangleSpinChain -s 1 -p 4 --use-lapack*

diagonalizes a 4 sites (-p option) spin 1/2 (the -s option should be set two twice the spin per site) with coupling constant J1=1 and J2=0.5. By default, the coupling constants are set at the Majumdar-Gosh point. They can be tuned using the --j1-value and --j2-value. The code only looks at positive Sz sector and all momentum sectors. The spectrum is stored in a file whose name looks like spin_1_2_translations_doubletrianglechain_j1_1.000000_j2_0.500000_n_4.dat with the following content

```text
#  double triangle spin 1 / 2 chain with tanslations and 4 sites j1=1.000000 j2=0.500000
# 2Sz K  E
0 0 -1.5
0 0 1.5
0 1 -0.5
0 2 -1.5
0 2 -0.5
0 3 -0.5
2 0 1.5
2 1 -0.5
2 2 -0.5
2 3 -0.5
4 0 1.5
```
