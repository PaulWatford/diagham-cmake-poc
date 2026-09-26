# FTIParticleHoleConjugation — manual

Source: DiagHam wiki page `FTIParticleHoleConjugation`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FTIParticleHoleConjugation](../FTI/FTIParticleHoleConjugation.md)
[FTIParticleHoleConjugation](FTIParticleHoleConjugation.md) applies a particle hole conjugation to any FTI state. At the moment, only single band spinless cases are supported. Its usage is straightforward

*\$PATHTODIAGHAM/build/FTI/src/Programs/FTI/FTIParticleHoleConjugation fermions_singleband_u_1.000_EScomp_maxband1_0_maxband2_0_n_12_ns_36_x_9_y_2_kx_3_ky_0.0.vec*

It produces a state fermions_ph_singleband_u_1.000_EScomp_maxband1_0_maxband2_0_n_6_ns_36_x_9_y_2_kx_6_ky_1.0.vec using a "_ph" after femions and updating automatically the number of particles and momentum sector.
