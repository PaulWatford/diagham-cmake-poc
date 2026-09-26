# GenericSpinChainShowBasis — manual

Source: DiagHam wiki page `GenericSpinChainShowBasis`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [GenericSpinChainShowBasis](../Spin/GenericSpinChainShowBasis.md)
GenericSpinChainShowBasis displays the many-body basis. For example, if we want to look at the basis for the spin 1 chain with 6 sites, periodic boundary conditions, $2S_z=8$ and with zero momentum, we just have to use

*\$PATHTODIAGHAM/build/Spin/src/Programs/GenericSpinChainShowBasis -p 6 -s 2 --periodic-chain --momentum 0 --sz-value 8*

The output should look like

```text
memory requested for Hilbert space = 48
memory requested for lookup table = 96ko
1 1 0 1 1 0 
1 1 1 0 1 0 
1 1 1 1 0 0 
1 1 1 1 1 -1 
```

We can also output the components of an eigenstate behind each configuration. For example, if we consider a complex binary vector vector.2.vec, we can type

*\$PATHTODIAGHAM/build/Spin/src/Programs/GenericSpinChainShowBasis -p 6 -s 2 --periodic-chain --momentum 0 --complex-vector --state vector.2.vec*

and we should see something like

```text
memory requested for Hilbert space = 48
memory requested for lookup table = 96ko
1 1 0 1 1 0  : (-0.5,-0)
1 1 1 0 1 0  : (-0.35355339059327,-0)
1 1 1 1 0 0  : (0.35355339059327,-0)
1 1 1 1 1 -1  : (0.70710678118655,-0)
```
