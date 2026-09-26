# FTIProjectOntoBands — manual

Source: DiagHam wiki page `FTIProjectOntoBands`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FTIProjectOntoBands](../FTI/FTIProjectOntoBands.md)
[FTIProjectOntoBands](FTIProjectOntoBands.md) projects a state defined in 2 or more bands into a single band (the lowest band). Its usage a similar to, e.g., [FTIDensity](FTIDensity.md). A typical example to project a state defined in three bands is

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIProjectOntoBands fermions_threeband_u_1.000_EScomp_maxband1_1_maxband2_0_n_12_ns_36_x_9_y_2_kx_0_ky_0.0.vec -s 3*

The code gives the weight onto the projected (defined as the square norm of the projected state) before exiting

```text
memory requested for Hilbert space = 12ko
memory requested for lookup table = 19Mo
Temporary Hilbert space dimension=404868
Hilbert space dimension 404868
memory requested for Hilbert space = 4Mo
memory requested for lookup table = 4Mo
weight onto band 0 : 0.48635639298124
```


The resulting state, whose name is build by replacing the string threeband with singleband_proj_band0 is NOT normalized. As a consequence, computing the PES with [FTIEntanglementEntropyParticlePartition](FTIEntanglementEntropyParticlePartition.md) or the density with [FTIDensity](FTIDensity.md) will result in an incorrect result. If needed, we can add the option --normalize when running [FTIProjectOntoBands](FTIProjectOntoBands.md) to enforce the normalization of the projected state.
