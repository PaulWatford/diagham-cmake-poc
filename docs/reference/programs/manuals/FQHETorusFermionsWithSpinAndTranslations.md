# FQHETorusFermionsWithSpinAndTranslations — manual

Source: DiagHam wiki page `FQHETorusFermionsWithSpinAndTranslations`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHETorusFermionsWithSpinAndTranslations](../FQHEOnTorus/FQHETorusFermionsWithSpinAndTranslations.md)
FQHETorusFermionsWithSpinAndTranslations is the analogue of [FQHETorusBosonsWithSpinAndTranslations](FQHETorusBosonsWithSpinAndTranslations.md) for fermionic systems.

## Different types of interaction

This code supports double gated Coulomb interaction, i.e., V(q) = tanh(q * xi / 2)/q where xi is the screening length (in magnetic length). It can be selected by the following parameters to pseudopotential file

```text
DoubleGatedCoulombLandauLevel = 0
DoubleGatedCoulombScreeningLength = 1.0
```

Note that only the lowest Landau level is available at the moment. By default, the interaction name is doublegatedcoulomb_l_0_xi_x.xxxxx and can be overrideby the Name option.
This interaction can be added to any pseudopotential interactions and the amplitude of the double gated Coulomb interaction can be controlled in each channel by the three parameters ScalingFactorUpUp, ScalingFactorDownDown and ScalingFactorUpDown. Thus the most generic pseudopotential file should look like

```text
Name = doublegatedcoulomb_l_0_xi_1.0000_v0ud_1.0000
DoubleGatedCoulombLandauLevel = 0
DoubleGatedCoulombScreeningLength = 1.0
PseudopotentialsUpDown = 1.0
ScalingFactorUpUp = 1.0
ScalingFactorDownDown = 1.0
ScalingFactorUpDown = 0.5
```


The code also supports the Rytova-Keldysh interaction (with the double gate screening) V(q) = tanh(q * xi / 2)/(q * (1+l_RK q)) where the interlayer distance is fixed through the option RytovaKeldyshInterlayerDistance. Note that the interaction is automatically selected when RytovaKeldyshInterlayerDistance is present in the interaction file.

## Time-Reversal symmetric model

FQHETorusFermionsWithSpinAndTranslations can also be used to simulate a system where the two copies labelled by the spin degree of freedom are related by time reversal symmetry, i.e., to simulate a fractional topological insulator using Landau levels. The typical usage is

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusFermionsWithSpinAndTranslations build/FQHE/src/Programs/FQHEOnTorus/FQHETorusFermionsWithSpinAndTranslations -p 10 -l 15 -s 0 --time-reversal --interaction-file pseudopotentials.dat --use-lapack -S --processors 4 --memory 5000 -n 6 --block-lanczos --block-size 6 --force-reorthogonalize --show-itertime*

where the interaction file looks like 

```text
Name = fti
Pseudopotentials = 1.0 1.0
```


The spectrum is written in the file *fermions_torus_timereversal_fti_n_10_2s_15_sz_0_ratio_1.000000.dat* which should contain

```text
# Kx Ky E
0 0 5.457321578659
0 0 6.203434371726
0 0 6.2210106197006
0 0 6.3871588023526
0 0 6.4347390174024
0 0 6.5576585066811
0 1 6.2023429331221
0 1 6.2647565469827
0 1 6.3344409405428
...
```
