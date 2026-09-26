# FQHESphereEntanglementSpectrum — manual

Source: DiagHam wiki page `FQHESphereEntanglementSpectrum`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereEntanglementSpectrum](../FQHEOnSphere/FQHESphereEntanglementSpectrum.md)
FQHESphereEntanglementSpectrum provides an easy way to parse the eigenvalues of the reduced density matrices produced by [FQHESphereBosonEntanglementEntropy](FQHESphereBosonEntanglementEntropy.md), FQHESphereFermionEntanglementEntropy *(wiki page FQHESphereFermionEntanglementEntropy)*,  or


## Basic usage

If you want to get the entanglement spectrum from the example bosons_haldane_gaffnian_n_14_2s_18_lz_0.0.full.ent that you have computed using [FQHESphereBosonEntanglementEntropy](FQHESphereBosonEntanglementEntropy.md), and look at the sector of 9 orbitals and 7 partilces, you just have to type

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereEntanglementSpectrum bosons_haldane_gaffnian_n_14_2s_18_lz_0.0.full.ent  -n 7 -l 9*

should do what this [website](http://www.nick-ux.org/~regnault/entropy) do for you by looking at Na=7 la=9. Still this is not 
a plot that is produced but a column formated file bosons_haldane_gaffnian_n_14_2s_18_lz_0.0.la_9_na_7.entspec

To plot the data you can use gnuplot with the following script

```text
set xlabel "{L_z}^A" font "Helvetica,34"
set ylabel "{/Symbols x}" font "Helvetica,34"
set size 1.3, 1.1
set xtics nomirror
set ytics nomirror
set yrange [:24]
set mxtics 5
show mxtics
set mytics 5
show mytics
set terminal postscript color eps enhanced "Helvetica" 28
set bars small
set output "bosons_haldane_gaffnian_n_14_2s_18_lz_0.0.la_9_na_7.entspec.eps"
plot "bosons_haldane_gaffnian_n_14_2s_18_lz_0.0.la_9_na_7.entspec" using 4:6:(0.35) notitle with xerrorbars lc 1 lt 1 lw 3 pt 0 ps 0
```
