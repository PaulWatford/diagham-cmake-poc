# FQHETorusShowBasis — manual

Source: DiagHam wiki page `FQHETorusShowBasis`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHETorusShowBasis](../FQHEOnTorus/FQHETorusShowBasis.md)
FQHESphereShowBasis displays the various n-body basis that are used in DiagHam on the torus geometry. For example,

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusShowBasis -p 3 -l 6 -y 0 --boson --no-translation*

shows the n-body basis for 4 bosons with 6 flux quanta and a total Ky momentum of zero. The --no-translation indicated that only the translation along one direction is taken into account (i.e. not the full magnetic translations). The display style uses occupation numbers, leftmost occupation being the corresponding to the one body state with the smallest Ky value, i.e. zero). Typical output will look like :

```text
memory requested for Hilbert space = 14ko
 (k_y = 0) :
0 0 1 0 0 2
0 0 0 1 1 1
1 1 0 0 0 1
0 0 0 0 3 0
1 0 1 0 1 0
0 2 0 0 1 0
1 0 0 2 0 0
0 1 1 1 0 0
0 0 3 0 0 0
3 0 0 0 0 0
```


If for a given eigenstate, you want to display each component with the corresponding occupation number, just you the --state option. For example, if you generate one eigenstate bosons_torus_kysym_delta_n_3_2s_6_ratio_1.000000_ky_0.0.vec using FQHETorusBosonsDelta *(wiki page FQHETorusBosonsDelta)*, just type 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnTorus/FQHETorusShowBasis -p 3 -l 6 -y 0 --boson --no-translation --state bosons_torus_kysym_delta_n_3_2s_6_ratio_1.000000_ky_0.0.vec*


```text
memory requested for Hilbert space = 14ko
 (k_y = 0) :
0 0 1 0 0 2  : -0.37430491565847
0 0 0 1 1 1  : 0.0650645634067
1 1 0 0 0 1  : 0.065064563406703
0 0 0 0 3 0  : -0.027911892041863
1 0 1 0 1 0  : 0.75143201752102
0 2 0 0 1 0  : -0.37430491565849
1 0 0 2 0 0  : -0.37430491565848
0 1 1 1 0 0  : 0.065064563406708
0 0 3 0 0 0  : -0.02791189204187
3 0 0 0 0 0  : -0.027911892041867
```
