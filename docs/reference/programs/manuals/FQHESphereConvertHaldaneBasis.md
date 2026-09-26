# FQHESphereConvertHaldaneBasis — manual

Source: DiagHam wiki page `FQHESphereConvertHaldaneBasis`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereConvertHaldaneBasis](../FQHEOnSphere/FQHESphereConvertHaldaneBasis.md)
FQHESphereConvertHaldaneBasis allows to convert a state from a squeezed basis to the full Fock basis with a fixed total Lz value. Typical usage looks like

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereConvertHaldaneBasis bosons_haldane_laughlin_n_8_2s_14_lz_0.0.vec --reference-file laughlin_n_8_2s_14.dat -o bosons_laughlin_n_8_2s_14_lz_0.0.vec*

where laughlin_n_8_2s_14.dat is the squeezed basis description (as explained in [FQHESphereJackGenerator](FQHESphereJackGenerator.md)), bosons_haldane_laughlin_n_8_2s_14_lz_0.0.vec is the input state in binary format. The converted will be stored in bosons_laughlin_n_8_2s_14_lz_0.0.vec, using the -o option.

If you want to convert from a squeezed basis to another basis, you can use the --target-haldane option. In that case, the reference file of the target squeezed basis has to be provided through the --target-referencefile option.
