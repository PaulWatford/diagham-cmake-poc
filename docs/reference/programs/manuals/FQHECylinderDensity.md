# FQHECylinderDensity — manual

Source: DiagHam wiki page `FQHECylinderDensity`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHECylinderDensity](../FQHEOnCylinder/FQHECylinderDensity.md)

> **Options in this manual that the current program does not list**: `--cylinder-perimeter`, `--fermion`, `--normalize`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
[FQHECylinderDensity](FQHECylinderDensity.md) computes the density profile and the integrated charge for states on the cylinder geometry. 

First we generate the Laughlin state $\nu=1/3$ using [FQHESphereJackGenerator](FQHESphereJackGenerator.md)
*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator --reference-file laughlin3_n_8_2s_21.dat --fermion -o fermions_unnormalized_haldane_n_8_2s_21_lz_0.0.vec*

where we have used the reference file

```text
NbrParticles=8
LzMax=21
ReferenceState=1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1
```


The unnormalized state can be converted to the cylinder geometry using [FQHECylinderNormalizeState](FQHECylinderNormalizeState.md)

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderNormalizeState -i fermions_unnormalized_haldane_n_8_2s_21_lz_0.0.vec --haldane --reference-file laughlin3_n_8_2s_21.dat --cylinder-perimeter 9.0 --normalize --output-file fermions_haldane_cylinder_perimeter_9.000000_n_8_2s_21_lz_0.0.vec*

To obtain the density profile, we just have to run

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderDensity --nbr-points 400 --haldane --reference-file laughlin3_n_8_2s_21.dat --input-state fermions_haldane_cylinder_perimeter_9.000000_n_8_2s_21_lz_0.0.vec*

The code will create a file fermions_haldane_cylinder_perimeter_9.000000_n_8_2s_21_lz_0.0.rho.dat (replacing the .vec extension with .rho.dat). This text file should look like

```text
# Length = 15.35889741755, Perimeter = 9, Aspect ratio = 0.58597956320198
# 0 0 0 (0.54756193124569,0)
# 1 0 0 (0.43662316354643,0)
# 2 0 0 (0.40006895675818,0)
# 3 0 0 (0.36705856641123,0)
# 4 0 0 (0.25514337026317,0)
# 5 0 0 (0.26513773727377,0)
# 6 0 0 (0.41770106572259,0)
# 7 0 0 (0.30891888141043,0)
# 8 0 0 (0.30625304866107,0)
# 9 0 0 (0.40465377729257,0)
# 10 0 0 (0.29087950141488,0)
# 11 0 0 (0.29087950141488,0)
# 12 0 0 (0.40465377729257,0)
# 13 0 0 (0.30625304866107,0)
# 14 0 0 (0.30891888141043,0)
# 15 0 0 (0.41770106572259,0)
# 16 0 0 (0.26513773727377,0)
# 17 0 0 (0.25514337026317,0)
# 18 0 0 (0.36705856641123,0)
# 19 0 0 (0.40006895675818,0)
# 20 0 0 (0.43662316354643,0)
# 21 0 0 (0.54756193124569,0)
-12.679448708775 1.1580477168475e-13 1.0638527610851e-14
-12.616209563095 2.2689176354567e-13 2.1094680462658e-14
-12.552970417416 4.4099979670168e-13 4.1484099343606e-14
-12.489731271736 8.5032613162762e-13 8.0937838883275e-14
-12.426492126056 1.6265231331468e-12 1.5664760134065e-13
...
```


The commented region provides the occupation numbers in orbital space. The first column is the orbital index and the fourth is the occupation number, while the second and the third column are the left and right state  index (both 0 here since we have a single state). Below the commented region is the density (integrated over the perimeter) in real space (second column) as a function of the position along the cylinder axis (first column). The third column is the integrated charge from minus infinity to the position given by the first column. The region along the cylinder axis where the density is computed is set by two parameters : 
- The length $L_x$ deduced from the number of flux quanta $N_\Phi$ and the cylinder perimeter $L_y$ (i.e. $L_x = 2 \pi (N_\Phi + 1) l_B^2/ L_y$.
- An offset set by the option --offset (with a default value of $5 l_B$)

The code evaluates the density with a regular grid of --nbr-points over the interval  $[-L_x-{\rm offset},L_x+{\rm offset}]$.
