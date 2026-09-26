# SpinChainComputeS2 — manual

Source: DiagHam wiki page `SpinChainComputeS2`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [SpinChainComputeS2](../Spin/SpinChainComputeS2.md)
[SpinChainComputeS2](SpinChainComputeS2.md) computes the S value of an eigenstate or a group of eigenstates. For example, if we generate an eigenstate spin_2_periodicaklt_p3_1.000000_p4_1.000000_n_8_sz_0_invsym_1_szsym_1_k_4.22.vec using [PeriodicSpinChainGeneralizedAKLT](PeriodicSpinChainGeneralizedAKLT.md)

*\$PATHTODIAGHAM/build/Spin/src/Programs/SpinChainComputeS2 -i spin_2_periodicaklt_p3_1.000000_p4_1.000000_n_8_sz_0_invsym_1_szsym_1_k_4.22.vec --force-real*

The output will look like

```text
can't guess x-momentum sector from file name spin_2_periodicaklt_p3_1.000000_p4_1.000000_n_8_sz_0_invsym_1_szsym_1_k_4.22.vec
can't guess x-momentum sector from file name spin_2_periodicaklt_p3_1.000000_p4_1.000000_n_8_sz_0_invsym_1_szsym_1_k_4.22.vec
N=8 Sz=0 2s=4 kx=4
memory requested for Hilbert space = 13ko
memory requested for lookup table = 16Mo
&lt;S^2&gt;=20 &lt;S&gt;=4
round(&lt;2S&gt;)=8
```


The first few lines are just information about testing symmetries (detected from the file name). The two last lines gives the expectation of $S^2$ and thus the S value. The last line is the rounded value of 2S. Note that we have used the --force-real to indicate that the eigenstate is real despite it has been computed using momentum conservation (since the state here is at momentum $\pi$).

If we have a set of eigenstates (for example corresponding to some degenerate energy), we can diagonalize them in one shot. Remember that diagonalization codes usually never produced $S^2$ eigenstates for degenerate energies, so this step is mandatory in that case to get the correct $S^2$ values. First ew need to provide a list of text through a text file (here basis.dat)


```text
spin_2_periodicaklt_p3_1.000000_p4_1.000000_n_8_sz_0_invsym_1_szsym_1_k_4.21.vec
spin_2_periodicaklt_p3_1.000000_p4_1.000000_n_8_sz_0_invsym_1_szsym_1_k_4.22.vec
```


Note that all the states should have the same quantum numbers. Then we just have to run

*\$PATHTODIAGHAM/build/Spin/src/Programs/SpinChainComputeS2 --degenerate-states basis.dat --force-real*

and we will get as an output

```text
can't guess x-momentum sector from file name spin_2_periodicaklt_p3_1.000000_p4_1.000000_n_8_sz_0_invsym_1_szsym_1_k_4.21.vec
can't guess x-momentum sector from file name spin_2_periodicaklt_p3_1.000000_p4_1.000000_n_8_sz_0_invsym_1_szsym_1_k_4.21.vec
N=8 Sz=0 2s=4 kx=4
memory requested for Hilbert space = 13ko
memory requested for lookup table = 16Mo
&lt;S^2&gt;=3.5687597400896e-15 &lt;S&gt;=3.5527136788005e-15
round(&lt;2S&gt;)=0
&lt;S^2&gt;=20 &lt;S&gt;=4
round(&lt;2S&gt;)=8
```

We now have the two  $S^2$ eigenvalues and thus the corresponding S values. If we want the corresponding eigenstates, we just have to include the --compute-eigenstates option

*\$PATHTODIAGHAM/build/Spin/src/Programs/SpinChainComputeS2 --degenerate-states basis.dat --force-real --compute-eigenstates*

This will create two new vectors spin_2_periodicaklt_p3_1.000000_p4_1.000000_n_8_sz_0_s_0_invsym_1_szsym_1_k_4.21.vec and spin_2_periodicaklt_p3_1.000000_p4_1.000000_n_8_sz_0_s_8_invsym_1_szsym_1_k_4.22.vec. Note that each file name includes the s value (here _s_0_ and _s_8_).
