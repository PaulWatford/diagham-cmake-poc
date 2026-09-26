# FTIGenericInteractionFromFileThreeBands — manual

Source: DiagHam wiki page `FTIGenericInteractionFromFileThreeBands`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FTIGenericInteractionFromFileThreeBands](../FTI/FTIGenericInteractionFromFileThreeBands.md)
[FTIGenericInteractionFromFileThreeBands](FTIGenericInteractionFromFileThreeBands.md) is the analog of [FTIGenericInteractionFromFileTwoBands](FTIGenericInteractionFromFileTwoBands.md) for models with three bands.


## With valley and spin


The quantum number sector selection can be done using a text file whose name should be set by the option --selected-sectors. It should contain four columns providing kx, ky, 2Sz (twice the spin projection) and 2Pz (twice the valley projection) respectively.

```text
# kx ky 2Sz 2Pz
0 0 4 4
0 1 -4 4
1 1 2 -2
...
```
