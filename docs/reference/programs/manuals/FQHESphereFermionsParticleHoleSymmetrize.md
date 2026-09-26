# FQHESphereFermionsParticleHoleSymmetrize — manual

Source: DiagHam wiki page `FQHESphereFermionsParticleHoleSymmetrize`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [FQHESphereFermionsParticleHoleSymmetrize](../FQHEOnSphere/FQHESphereFermionsParticleHoleSymmetrize.md)
FQHESphereFermionsParticleHoleSymmetrize computes the particle-hole conjugate state of a given fermionic state. Typical usage is 

*\$PATHTODIAGHAM/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereFermionsParticleHoleSymmetrize fermions_laughlin3_n_4_2s_9_lz_0.0.vec *

This create the particle-hole conjugate state and save it under fermions_laughlin3_n_6_2s_9_lz_0.0.vec. The input vector has to be expressed in the normalized sphere basis. The default output file name  replace fermions_ with fermions_holes and set the correct value for n. It can be changed using the -o option. Squeezed Hilbert space or Lz symmetric space can be set using the options --haldane, --reference-file, --symmetrized-basis or --load-hilbert.
