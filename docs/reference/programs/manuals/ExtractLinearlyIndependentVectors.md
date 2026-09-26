# ExtractLinearlyIndependentVectors — manual

Source: DiagHam wiki page `ExtractLinearlyIndependentVectors`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [ExtractLinearlyIndependentVectors](../core/ExtractLinearlyIndependentVectors.md)
ExtractLinearlyIndependentVectors extract a set of linearly independent vectors from a given set of vectors. The typical usage

*\$PATHTODIAGHAM/build/src/Programs/ExtractLinearlyIndependentVectors -b basis.dat --use-lapack*

where basis.dat is text file that gives the description of the initial set of vectors. A typical input file (for four vectors) is 

```text
Basis = intialvector0.vec intialvector1.vec intialvector2.vec intialvector3.vec
```

By default, linearly independent vectors are saved under names like vectorX.vec. The prefix vector in the file name might be changed using the --vector-prefix option. To only check the number of linearly independent vectors without computing them, use the --check-only option.
