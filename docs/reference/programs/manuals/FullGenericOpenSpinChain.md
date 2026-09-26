# FullGenericOpenSpinChain — manual

Source: DiagHam wiki page `FullGenericOpenSpinChain`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FullGenericOpenSpinChain](../Spin/FullGenericOpenSpinChain.md)
[FullGenericOpenSpinChain](FullGenericOpenSpinChain.md) performs the exact diagonalization of the Hamiltonian

$H=\sum_{i=1}^{N-1} \left( J_x S_{x,i}S_{x,i+1} + J_y S_{y,i}S_{y,i+1}  + J_z S_{z,i}S_{z,i+1} \right) + \sum_{i=1}^{N} \sum_{\alpha=x,y,z} h_{\alpha,i} S_{\alpha,i}$

The simplest usage is 

*\$PATHTODIAGHAM/build/Spin/src/Programs/FullGenericOpenSpinChain -p 10 --use-lapack --full-diag 2000*

This will generate the full spectrum for a spin 1/2 open chain with $J_x=J_y=J_z=1$ and $h_{\alpha,i}=0$. The spectrum is stored in the ASCII file spin_1_2_openchain_n_10_jx_1.000000_jy_1.000000_jz_1.000000.dat  that should contain 

```text
#  open spin 1/2 chain with 10 sites 
# E
-4.2580352072829
-3.9306735895016
...
```

We can jet the values $J_x,J_y$ and $J_z$ using the options --jx-value, --jy-value and --jz-value which are set to one by default. Periodic boundary conditions can be added by adding the --use-periodic option. To set the magnetic field, several possibilities are provided. For a uniform magnetic field $h_{\alpha,i}=h_{\alpha}$, we can use the --hx-value, --hy-value and --hz-value. The magnetic field can also be specified through a column formatted ASCII file. For example, if we use the following h filed defined by hfield.dat


```text
0.074582052417099 -0.028589533409104 0.18689403915778
-0.078413534443825 -0.044323729118332 -0.2201355881989
0.035255062393844 0.043060069903731 0.1460151784122
0.12923212815076 0.18053271784447 0.069194715237245
0.23684391356073 0.073561968049034 -0.060641295043752
0.22684709331952 -0.21963313990273 -0.030818197177723
0.24988580634817 0.15365192270838 -0.12620972911827
-0.2404253995046 -0.00093650259077549 0.086510408204049
-0.24202130222693 -0.22648968733847 0.13108505285345
0.21121106413193 0.17433513072319 -0.15865674638189
```


we just run for a system with periodic boundary conditions

*\$PATHTODIAGHAM/build/Spin/src/Programs/FullGenericOpenSpinChain -p 10 --use-lapack --full-diag 2000 --use-periodic --fullh-values hfield.dat*

We can also generate a random on-site field. Its magnitude in each direction is set by --random-hxvalue, --random-hyvalue and --random-hzvalue . It uses a uniform distribution (a gaussian distribution can be used instead is we use --random-gaussianhxvalue, --random-gaussianhyvalue and --random-gaussianhzvalue). For example 

*\$PATHTODIAGHAM/build/Spin/src/Programs/FullGenericOpenSpinChain -p 10 --use-lapack --full-diag 2000 --use-periodic --random-hxvalue 0.25 --random-hzvalue 0.25 --all-eigenstates*

The spectrum is stored in spin_1_2_closedchain_n_10_jx_1.000000_jy_1.000000_jz_1.000000_randomhx_0.250000_randomhy_0.000000_randomhz_0.250000_runid_0.dat. The random magnetic that is generated is available in an ASCII file spin_1_2_closedchain_n_10_jx_1.000000_jy_1.000000_jz_1.000000_randomhx_0.250000_randomhy_0.000000_randomhz_0.250000_runid_0.hvalues that can used as the input of the --fullh-values if we need to re-run the code with the same random field. Finally in this example, we have required to compute all the eigenstates by adding the --all-eigenstates options. So a file spin_1_2_closedchain_n_10_jx_1.000000_jy_1.000000_jz_1.000000_randomhx_0.250000_randomhy_0.000000_randomhz_0.250000_runid_0.eigenvec.mat is written, storing the eigenstates in a binary matrix. This matrix can be processed either with programs such as [MatrixExtractColumns](MatrixExtractColumns.md) to extract individual eigenstates more directly to compute properties duch as [SpinChainMultipleEntanglementSpectra](SpinChainMultipleEntanglementSpectra.md).

Note that as long as $h_{y,0}=0$, the hamiltonian is real and so are the eigenstates. This greatly speeds-up the calculations. Moreover, [FullGenericOpenSpinChain](FullGenericOpenSpinChain.md) has a full support of scalapack *(wiki page scalapack)*.
