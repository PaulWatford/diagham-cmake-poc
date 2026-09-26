# FTIConvertConstrainedBasis — manual

Source: DiagHam wiki page `FTIConvertConstrainedBasis`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FTIConvertConstrainedBasis](../FTI/FTIConvertConstrainedBasis.md)

> **Options in this manual that the current program does not list**: `--target-max-band`, `--target-min-band`. They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.
FTIConvertConstrainedBasis converts states defined 2 band systems or more with restrictions on, e.g., the number of particles per band to another Hilbert space with different restrictions. For example, coonsider the state *fermions_threeband_u_1.000_test_maxband1_12_maxband2_1_n_12_ns_24_x_2_y_6_kx_0_ky_0.0.vec* defined in a Hilbert space with at most 12 particles in band 1 and 1 particle in band 2 that we want to project on a Hilbert space with now at most 5 particles in band 1. The command 

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIConvertConstrainedBasis fermions_threeband_u_1.000_test_maxband1_12_maxband2_1_n_12_ns_24_x_2_y_6_kx_0_ky_0.0.vec -s 3 --target-max-band1 5 --target-max-band2 1* 

will generate the following output

```text
fermions_threeband_u_1.000_test_maxband1_12_maxband2_1_n_12_ns_24_x_2_y_6_kx_0_ky_0.0.vec N=12 Nx=2 Ny=6 kx=0 ky=0
Temporary Hilbert space dimension=2721736
Hilbert space dimension 2721736
memory requested for Hilbert space = 31Mo
memory requested for lookup table = 3Mo
Temporary Hilbert space dimension=1325240
Hilbert space dimension 1325240
memory requested for Hilbert space = 15Mo
memory requested for lookup table = 3Mo
Weight of projected fermions_threeband_u_1.000_test_maxband1_12_maxband2_1_n_12_ns_24_x_2_y_6_kx_0_ky_0.0.vec : 0.99738429562518
```

and a projected vector fermions_threeband_projected_maxband1_5_maxband2_1_n_12_ns_36_x_2_y_6_kx_0_ky_0.0.vec . Note that the weight of the projected state is provided in the output and the vector is *NOT* normalized unless explicitly requested bu adding the option --normalize. By default, the output name contains the "projected", the rest of the file name is obtained form the target Hilbert space parameters. This default behavior can be superseded using the option --output-file-name.

Most of the parameters of the input vectors are inferred from the file names with the notable exception of the number of bands (set by -s or --nbr-subbands) and the filtered orbitals list that should be provided with the option --source-allowed-orbitals. For the target Hilbert space, options --target-max-bandX, --target-min-bandY and/or --target-allowed-orbitals. Note that these parameters are independent from the source Hilbert space parameters. Thus by default, there are no restrictions on the target Hilbert space. if more than one vectors should be converted (all defined in the same Hilbert space), you can provide a list of vectors using the option --degenerated-groundstate.
