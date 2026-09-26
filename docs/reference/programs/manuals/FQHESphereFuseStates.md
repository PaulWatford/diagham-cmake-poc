# FQHESphereFuseStates — manual

Source: DiagHam wiki page `FQHESphereFuseStates`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereFuseStates](../FQHEOnSphere/FQHESphereFuseStates.md)
FQHESphereFuseStates creates states according to the product rule symmetry. Typical usage is

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereFuseStates -i fuse.dat --haldane --reference-file laughlin3_n_4_2s_9.dat -o fermions_haldane_fused_n_4_2s_9_lz_0.0.vec*

The options --haldane --reference-file laughlin3_n_4_2s_9.dat are needed if the output state has to be computed in a squeezed basis. The -o option sets the output file name. fuse.dat is a text file that describes which has to be used for the product rule symmetry. A typical line is this column formatted text file should have the following form

```text
fermions_unnormalized_haldane_laughlin3_qh_0_n_1_2s_1_lz_1.0.vec fermions_unnormalized_haldane_laughlin3_qh_0_n_3_2s_7_lz_3.0.vec none laughlin3_qh_0_n_3_2s_7.dat none laughlin3_qh_0_n_3.hil 0 1.0
```

The first two columns are the vector binary states (in the unnormalized basis) use for the left part (first column) and right part (second column) of the product rule symmetry. This is the minimal information that has to be provided. If the states are stored in a squeezed basis, column 3 to 6 should be defined. Columns 3 and 4 give the names of file that describe the squeezed basis ("none" if it is a regular basis). Columns 5 and 6  give the names of file where the squeezed Hilbert is stored (it should be set to **none** if no such file is needed). Column 7 defines the padding i.e. the number of empty orbitals between the two fused states. It is optional unless column 8 exists but supersede the --padding option if it exists. Column 8 allows to define a global multiplicative factor to the resulting fused state. Notice that the final output state is the linear combination of the fused states described by each line of the input file.
