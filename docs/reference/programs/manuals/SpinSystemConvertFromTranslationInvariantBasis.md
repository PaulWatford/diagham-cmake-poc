# SpinSystemConvertFromTranslationInvariantBasis — manual

Source: DiagHam wiki page `SpinSystemConvertFromTranslationInvariantBasis`, as of 2026-09-24 (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build).

Generated option reference: [SpinSystemConvertFromTranslationInvariantBasis](../Spin/SpinSystemConvertFromTranslationInvariantBasis.md)
SpinSystemConvertFromTranslationInvariantBasis converts a many-body spin chain state written in a basis with a fixed momentum to the full many-body basis. Its usage is 

*\$PATHTODIAGHAM/build/Spin/src/Programs/SpinSystemConvertFromTranslationInvariantBasis -i spin_1_periodicaklt_n_8_sz_0_k_1.0.vec*

This will create a binary vector spin_1_periodicaklt_n_8_sz_0.0.0.vec where the string pattern indicating the momentum has been removed (in this example  _k_1).

This code only works when the input vector is complex. Beware that this code is highly untested for one dimensional systems.
