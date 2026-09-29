# Open defect U41: 2DToricCodeEntanglementEntropy loses its spectrum when --kept-sites has a path

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/known-defects.md` (U41), found 2026-09-29 while writing the entanglement goldens; r4493 sources, GCC 15.

## What happens

```
2DToricCodeEntanglementEntropy -x 3 -y 3 --kept-sites /abs/path/kept.txt --export-entspectrum
```

prints the entropy but writes no `.dat` file. The output name is built in `Spin/src/Programs/2DToricCodeEntanglementEntropy.cc` (lines 978-981) as `2dtoriccode_entspectrum_x_%d_y_%d_%s.dat` with the `--kept-sites` argument pasted in as given, so the name contains slashes (`2dtoriccode_entspectrum_x_3_y_3_/abs/path/kept.txt.dat`); the directory does not exist, `open` fails and the code neither checks nor reports it. With a bare file name in the working directory it works.

## Fix direction

Use the base name of the kept-sites file, and check that the output file opened.

## Consequence for the tests

`tests/entanglement.cmake` exercises the rectangle regions (`--nbra-sitex/--nbra-sitey`), whose output name is fixed; the kept-sites path is not tested.
