# FTIDensity — manual

Source: DiagHam wiki page `FTIDensity`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FTIDensity](../FTI/FTIDensity.md)
FTIDensity is a generic code to compute density expectation values for FTI eigenstates. A typical usage is 

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIDensity fermions_twoband_u_0.600_twoband_eps_10_xi_20_th_3.7_m_0.6_V_20.8_w_-23.8_phi_-107.7_n_12_ns_18_x_3_y_3_kx_1_ky_2_pz_4.0.vec --nbr-subbands 4 --decoupled*

Most of the system parameters (wuth the exception of the number of subbands) are inferred from the input file (here *fermions_twoband_u_0.600_twoband_eps_10_xi_20_th_3.7_m_0.6_V_20.8_w_-23.8_phi_-107.7_n_12_ns_18_x_3_y_3_kx_1_ky_2_pz_4.0.vec*).

The code creates an ASCII output file, replacing the .vec extension of the input file with .rho.dat

```text
# kx ky spin sigma sigma' <c^+ c>
0 0 0 0 0 (0.016040169509158,0)
0 0 1 0 0 (0.012407150138962,0)
0 0 0 0 1 (-0.0072963145226506,-0.00030110869553328)
0 0 1 0 1 (-6.5829645858117e-05,2.7384256575587e-06)
0 0 0 1 0 (-0.0072963145226506,0.00030110869553328)
...
2 2 1 1 0 (-0.0069063507889334,-0.015948744985578)
2 2 0 1 1 (0.78655889170032,0)
2 2 1 1 1 (0.365714667446,0)
# partial density spin=0 sigma=0 = (0.98055461971139,0)
# partial density spin=0 sigma=1 = (7.0194453802891,0)
# partial density spin=1 sigma=0 = (0.42831734360924,0)
# partial density spin=1 sigma=1 = (3.5716826563911,0)
# total density = (12.000000000001,0)
```

Beyond the expectation values of each c^\dagger_{k,spin,sigma}c_{k,spin,sigma'} operator, the last lines (starting with a #) give the summation of the expectation values over the momentum index.

Note that --nbr-subbands refers to the total number of orbitals (including e.g. spin/valley degree of freedom, number of bands) while --decoupled indicates that one degree of freedom has a U(1) conservation rule. The quantum number (Sz or Pz) is inferred from the file name (like _pz_0 or _sz_1). Be sure that it only appears once in the file name to avoid potential issues.

FTIDensity is also able to evaluate density-density (normal ordered) expectation values by adding the --rhorho option. It will create an additional ASCII output file, replacing the .vec extension of the input file with .rhorho.dat that will contain something like

```text
# kx1 ky1 spin1 sigma1 kx2 ky2 spin2 sigma2 kx3 ky3 spin3 sigma3 kx4 ky4 spin4 sigma4 <c^+ c^+ c c>
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 (0,0)
0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 (0.026354001819583,0)
0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 (0.04637497979938,0)
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 (0,0)
0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 1 (-0.00010170659213893,0.00053068238161885)
0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 1 (0.00034267944226811,0.00059127337598877)
```


If you are only interested in the intra band terms (namely setting sigma1=sigma2=sigma3=sigma4), then add the option --intraband-only. This will greatly speed up the calculations in multiband systems.
