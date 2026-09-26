# FQHESphereFermionsCorrelation — manual

Source: DiagHam wiki page `FQHESphereFermionsCorrelation`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereFermionsCorrelation](../FQHEOnSphere/FQHESphereFermionsCorrelation.md)
FQHESphereFermionsCorrelation allows to compute the density-density correlation function and the density correlation. Typical usage is 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereFermionsCorrelation -e fermions_laughlin3_n_6_2s_15_lz_0.0.vec*.

The input vector set through the -e option, can be either real or complex. The density-density correlation function is recorder within a ascii text file whose default name is the same than the input vector with the extension .vec replaced by .rhorho. It looks like :

```text
# density-density correlation coefficients for fermions_haldane_laughlin3_phase_0.25_n_6_2s_15_lz_0.0.vec
#
# (l+S)    n_l
# 0 (-0.1853253351536,0)
# 1 (-0.17749894714969,0)
# 2 (-0.17022771010782,0)
# 3 (-0.16473643925751,0)
# 4 (-0.16441377280697,0)
# 5 (-0.17224229844679,0)
# 6 (-0.18854502150451,0)
# 7 (-0.21001956558218,0)
# 8 (-0.23006010567721,0)
# 9 (-0.23907930130314,0)
# 10 (-0.22450829648961,0)
# 11 (-0.17377823526627,0)
# 12 (-0.086889117633136,0)
# 13 (0,0)
# 14 (0,0)
# 15 (0,0)
0 0
0.0086036058143182 3.3167740786046e-15
0.017207211628636 2.1226254229208e-13
0.025810817442955 2.4175942316431e-12
0.034414423257273 1.3581987409501e-11
0.043018029071591 5.1803117124306e-11
```


The first column is the distance on the sphere from the north pole (at $\theta = 0$) and the second column is the correlation function at this point. First column can be switched to the chord distance (--chord option) or the theta angle (--radians option). The output file name can be set through the -o option. The function sampling can be set using the --nbr-points option (the default is set to 1000 steps).

The --density option allows to compute the density correlation instead of the density-density correlation. In that case, the default output file name replaces the .vec extension of the input vector with .rho. Thus a command line such as 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereFermionsCorrelation --density -e fermions_laughlin3_n_6_2s_15_lz_0.0.vec*

would give a output file fermions_laughlin3_n_6_2s_15_lz_0.0.rho which looks like

```text
# density correlation coefficients for fermions_laughlin3_n_6_2s_15_lz_0.0.vec
#
# (l+S)    n_l
# 0 0.375
# 1 0.375
# 2 0.375
# 3 0.375
# 4 0.375
# 5 0.375
# 6 0.375
# 7 0.375
# 8 0.375
# 9 0.375
# 10 0.375
# 11 0.375
# 12 0.375
# 13 0.375
# 14 0.375
# 15 0.375
0 0.47746482927569
0.0086036058143182 0.47746482927569
0.017207211628636 0.47746482927569
0.025810817442955 0.47746482927569
0.034414423257273 0.47746482927569
0.043018029071591 0.47746482927569
```

The first column is the distance on the sphere from the north pole (at $\theta = 0$) and the second column is the density at this point. Note that the density is normalized such that is equal to the number of particles once we integrate over the whole sphere, meaning $\int_0^{2\pi} d \phi \int_0^{\pi} \sin\theta d\theta \rho = N$. Note that the commented first lines gives all the matrix elements $\langle \Psi | c^{\dagger}_{S+l} c_{S+l} | \Psi \rangle$.
