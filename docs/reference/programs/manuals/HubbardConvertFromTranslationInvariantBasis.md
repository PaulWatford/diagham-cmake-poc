# HubbardConvertFromTranslationInvariantBasis — manual

Source: DiagHam wiki page `HubbardConvertFromTranslationInvariantBasis`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [HubbardConvertFromTranslationInvariantBasis](../HubbardModels/HubbardConvertFromTranslationInvariantBasis.md)
This code allows to unwarp a state expressed in a basis resolved in momentum to a basis where there is no constraint on the momentum. The typical usage is

*\$PATHTODIAGHAM/build/FTI/src/Programs/HubbardModels/HubbardConvertFromTranslationInvariantBasis -i fermions_kitaev_heisenberg_gutzwiller_xmomentum_4_stripe_n_10_x_12_t_0.4_tK_0.8_j1_0_j2_0_kx_0.0.vec*

where fermions_kitaev_heisenberg_gutzwiller_xmomentum_4_stripe_n_10_x_12_t_0.4_tK_0.8_j1_0_j2_0_kx_0.0.vec is the binary vector of the state to convert. The output name is automatically deduced from the input file name. Multiple convertions can be done at the same by providing an ascii file to the --degenerate-states. This file should containg one vector file name per line.
