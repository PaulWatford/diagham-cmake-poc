# FQHESphereASCIIToState — manual

Source: DiagHam wiki page `FQHESphereASCIIToState`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereASCIIToState](../FQHEOnSphere/FQHESphereASCIIToState.md)
FQHESphereASCIIToState allows to write a binary vector corresponding to state whose description is given by an text file.  For example

*\$PATHTODIAGHAM/build/src/Programs/FQHEOnSphere/FQHESphereASCIIToState -p 3 -l 4 -z 0 --boson --ascii-state dummy.txt -o bosons_dummy_n_3_s_4_lz_0.0.vec*

will create a binary vector bosons_dummy_n_3_s_4_lz_0.0.vec for a state of 3 bosons with 4 flux quanta at Lz=0 which is described by dummy.txt. A typical text file looks like

```text
1 0 1 0 1 : 1
0 2 0 0 1 : -2
```

The text file should have the same syntax than [FQHESphereShowBasis](FQHESphereShowBasis.md). Only non-zero components are mandatory. By default, FQHESphereASCIIToState normalize the state to 1, bosons_dummy_n_3_s_4_lz_0.0.vec should look like (you can use [VectorBinary2Ascii](VectorBinary2Ascii.md))

```text
0.44721359549996
-0.89442719099992
0
0
0
```

If you don't want the state to be normalize to 1, just use the --no-normalization option.
