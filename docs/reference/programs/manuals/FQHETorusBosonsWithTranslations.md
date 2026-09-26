# FQHETorusBosonsWithTranslations — manual

Source: DiagHam wiki page `FQHETorusBosonsWithTranslations`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHETorusBosonsWithTranslations](../FQHEOnTorus/FQHETorusBosonsWithTranslations.md)
FQHETorusBosonsWithTranslations provides exact diagonalization for bosons on a torus geometry and any generic two body interaction. It fully supports magnetic translations in both directions. It is the bosonic counterpart of [FQHEFermionsTorusWithTranslation](FQHETorusFermionsWithTranslations.md) (see this code for additional details). Typical usage is 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusBosonsWithTranslations -p 4 -l 8 --use-lapack --interaction-file pseudopotentials_delta.dat --all-points *

pseudopotentials_delta.dat is the two body pseudopotential file that looks like

```text
Name = laughlin
Pseudopotentials = 1.0
```

The spectrum will be stored in a file *bosons_torus_laughlin_n_4_2s_8_ratio_1.000000.dat*. It should look like

```text
# Kx Ky E
0 0 -1.0113437864945e-15
0 0 1.168247084162
0 0 1.3436250395106
0 0 1.4347042589888
0 0 1.9575631506679
0 0 2.1809984955673
0 0 2.3431986010129
0 0 2.9925162798495
0 0 3.1006312362024
```
