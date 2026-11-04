# Open finding U63: FTIGenericInteractionFromFileSingleBand ignores --eigenstate-file

Purpose: reproducer for the DiagHam authors; no fix is carried here.
Source: `docs/reference/open-findings.md` (U63), found 2026-11-03 while turning Nicolas Regnault's example into a test; r4493 sources, GCC 15.

## What happens

The option is declared at line 111 of `FTI/src/Programs/FTI/FTIGenericInteractionFromFileSingleBand.cc` ("filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec"), but nothing reads it: the eigenstate names are always built from the eigenvalue file name (line 1543). With

```
FTIGenericInteractionFromFileSingleBand -x 3 -y 2 -p 2 --interaction-file x3_y2_interaction.dat --interaction-name oracle \
    --only-kx 0 --only-ky 0 --full-diag 5000 --eigenstate --eigenstate-file named
```

the vectors are `fermions_oneband_u_1.000_oracle_n_2_ns_12_x_3_y_2_kx_0_ky_0.0.vec` and `.1.vec`. No `named_kx_0_ky_0.0.vec` is written, and nothing says the option was ignored.

## Fix direction

Use the option when it is given, or remove it from the help.

## Consequence for the tests

`openfinding.fti.generic_interaction_single_band.eigenstate_file_ignored` (WILL_FAIL) asks GenericOverlap for the squared overlap of `named_kx_0_ky_0.0.vec` with itself. The eigenvector tests of Nicolas Regnault's example give their two runs different `--interaction-name` values, which name the files.
