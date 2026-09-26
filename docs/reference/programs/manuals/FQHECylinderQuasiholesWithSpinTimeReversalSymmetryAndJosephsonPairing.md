# FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing — manual

Source: DiagHam wiki page `FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing](../FQHEOnCylinder/FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing.md)
[FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing](FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing.md) is the analogue of [FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing](FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing.md) for the cylinder geometry in absence of the momentum conservation and Josephson effect. Indeed, it can handle any type of one body terms. For this reqason, its usage is slightly different than [FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing](FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing.md). 

## Translation invariant and Josephson effect

First we consider the case where the confining potential and the superconducting regions do not break the translation symmetry along the cylinder perimeter.

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing -l 10 --cylinder-perimeter 8.0 --directory ../nphi_10_cylinder_perimeter_8.000000 --use-lapack --charging-energy 0.195751 --average-nbrparticles 18.0 --interaction-name test --confining-file confining_momentum_nphi_10.dat --superconducting-file pairing_nphi_10.dat*

As explained previously, the confining potential and the superconducting region descriptions are now provided as separated files. The confining potential is provided by the --confining-file option and allows to encode any term such as $\sum_{m,n,\sigma} V_{m,n,\sigma} c^\dagger_{m,\sigma} c_{n,\sigma}$. The file should look like

```text
0 0 5.0 5.0
1 1 4.0 4.0
2 2 3.0 3.0
3 3 2.0 2.0
4 4 1.0 1.0
5 5 0.0 0.0
6 6 1.0 1.0
7 7 2.0 2.0
8 8 3.0 3.0
9 9 4.0 4.0
10 10 5.0 5.0
```

- The first column is the creation operator index i.e. $m$ (0 being the leftmost orbital).
- The second column is the annihilation operator index i.e. $n$ (0 being the leftmost orbital).
- The third column is the amplitude of $V_{m,n,\uparrow}$.
- The fourth column is the amplitude of $V_{m,n,\downarrow}$.

The mean field superconducting coupling is set by the --superconducting-file option. The text file allows to defined any superconducting term $\sum_{m} \Delta_{m} c^\dagger_{m,\uparrow} c^\dagger_{m,\sigma}\;+\; h.c.$. A typical example is 

```text
0 0 1.8 0.2
1 1 1.8 0.2
2 2 1.8 0.2
3 3 1.8 0.2
4 4 0.0 0.0
5 5 0.0 0.0
6 6 0.0 0.0
7 7 1.8 0.0
8 8 1.8 0.0
9 9 1.8 0.0
10 10 1.8 0.0
```

- The first column is the leftmost creation operator index i.e. $m$ (0 being the leftmost orbital).
- The second column is the rightmost creation operator index i.e. $m$ (0 being the leftmost orbital).
- The third column is the amplitude of $\Delta_{m}$.
- The fourth column is the phase of $\Delta_{m}$ in $\pi$ units.

## Translation breaking confining potential and/or superconducting term

For example, the command

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnCylinder/FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing -l 18 --directory ../nphi_18_cylinder_perimeter_8.000000 --cylinder-perimeter 8.0 --use-lapack --average-nbrparticles 20 --charging-energy 0.37246249700455 --full-diag 3000   --memory 8000 -s 4 --force-negativeky -n 10 --block-lanczos --block-size 10 --force-reorthogonalize --confining-file  confining_momentum_cylinder_perimeter_8.000000_2s_18_alphar_1_x0r_0.000000_v0r_1.000000_alphal_1_x0l_0.000000_v0l_1.000000_l_7.000000.dat  --superconducting-file confining_momentum_cylinder_perimeter_8.000000_2s_18_l_7.000000_alphar_0_x0r_2.600000_v0r_2.000000_phir_0.000000_alphal_0_x0l_-2.600000_v0l_2.000000_phil_0.000000.dat*

shows that the code requires two files to provide the description of the confining potential (the --confining-file option) and the description of the mean field superconducting term (the --superconducting-file option).
