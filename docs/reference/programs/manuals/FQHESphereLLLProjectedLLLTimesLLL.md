# FQHESphereLLLProjectedLLLTimesLLL — manual

Source: DiagHam wiki page `FQHESphereLLLProjectedLLLTimesLLL`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereLLLProjectedLLLTimesLLL](../FQHEOnSphere/FQHESphereLLLProjectedLLLTimesLLL.md)
[FQHESphereLLLProjectedLLLTimesLLL](FQHESphereLLLProjectedLLLTimesLLL.md) computes the product of wavefunctions in the lowest Landau level with or without reverse flux attachment.

```text
time  $PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereLLLProjectedLLLTimesLLL -1 bosons_haldane_pfaffian_n_8_2s_6_lz_0.0.vec -2 fermions_haldane_laughlin3_n_8_2s_21_lz_0.0.vec --reference-file1 pfaffian_n_8_2s_6.dat --reference-file2 laughlin3_n_8_2s_21.dat --reverse-flux --normalize --interaction-name tpfaffian3 --disable-lzsymmetry
memory requested for Hilbert space = 1Mo
memory requested for Hilbert space = 11Mo
memory requested for Hilbert space = 4Mo
generating state fermions_sphere_tpfaffian3_n_8_2s_15_lz_0.0.vec

real    6m35.087s
user    6m35.212s
sys     0m0.008s
```

Note that this code is fully parallel

```text
time  ~$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereLLLProjectedLLLTimesLLL -1 bosons_haldane_pfaffian_n_10_2s_8_lz_0.0.vec -2 fermions_haldane_laughlin3_n_10_2s_27_lz_0.0.vec --reference-file1 pfaffian_n_10_2s_8.dat --reference-file2 laughlin3_n_10_2s_27.dat --reverse-flux --normalize --interaction-name test --disable-lzsymmetry -S --processors 24
....
processing 1058 4
processing 1062 4
processing 1066 4

real    7388m35.517s
user    174896m23.236s
sys     0m1.616s
```
