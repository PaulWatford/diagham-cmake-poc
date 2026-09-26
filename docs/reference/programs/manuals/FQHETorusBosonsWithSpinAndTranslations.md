# FQHETorusBosonsWithSpinAndTranslations — manual

Source: DiagHam wiki page `FQHETorusBosonsWithSpinAndTranslations`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHETorusBosonsWithSpinAndTranslations](../FQHEOnTorus/FQHETorusBosonsWithSpinAndTranslations.md)
FQHETorusBosonsWithSpinAndTranslations provides exact diagonalization for bosons with a spin 1/2 degree of freedom on a torus geometry and any generic two body interaction. It fully supports magnetic translations in both directions. It is the spinful counterpart of [FQHETorusBosonsWithTranslations](FQHETorusBosonsWithTranslations.md) (see this code for additional details). 

For SU(2)-symmetric interaction, the pseudopotential file is similar to [FQHETorusBosonsWithTranslations](FQHETorusBosonsWithTranslations.md). As an example, let's consider the Halperin (221) state. The exact interaction can be enconded with the following file

```text
Name = halperin_221
Pseudopotentials = 1.0
```

Then the command 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusBosonsWithSpinAndTranslations -p 4 -l 6 --use-lapack --interaction-file pseudopotentials_221.dat*

produces a file *bosons_torus_su2_halperin_221_n_4_2s_6_sz_0_ratio_1.000000.dat*. Its content  should look like

```text
# Kx Ky E
0 0 -1.3322676295502e-15
0 0 0.73294529995177
0 0 0.83998531343515
0 0 0.89853855304335
0 0 0.9399234125861
0 0 1.0150775727123
```

Kx is the relative momentum and Ky the centre of mass momentum. By default, only the Sz=0 sector is computed. Other sectors can be accessed using the -s option. This options takes twice Sz as an integer parameter.

For a generic (i.e.e non SU(2)-symmetric) interaction, the pseudopotentials can be set between per spin-component, i.e.  up-up, down-down or up-down. For example the Halperin (220) state can be obtained with the following pseudopotential file

```text
Name = halperin_220
PseudopotentialsUpUp = 1.0
PseudopotentialsDownDown = 1.0
PseudopotentialsUpDown = 0.0
```

Many options such as fixing the aspect ratio of the torus are identical to [FQHEFermionsTorusWithTranslation](FQHETorusFermionsWithTranslations.md). In particular, if you intend to use the Lanczos algorithm, we strongly advise to look at the explanations there.
