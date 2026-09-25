# Auto-generated CMakeLists.txt for src/MathTools/NumericalAnalysis
# Regenerate with: python3 scripts_cmake/extract_autotools.py

diagham_add_library(NumericalAnalysis
    SOURCES
        AbstractNumericalInterval.cc
        FiniteRegularSubdivisionNumericalInterval.cc
        RNumericalInterval.cc
        Abstract1DRealFunction.cc
        Abstract1DComplexFunction.cc
        Abstract1DComplexTrialFunction.cc
        Abstract1DComplexFunctionOnSphere.cc
        Abstract1DComplexTrialFunctionOnSphere.cc
        AntisymmetrizedComplexFunction.cc
        SymmetrizedComplexFunction.cc
        Tabulated1DRealFunction.cc
        Constant1DRealFunction.cc
        Linear1DRealFunction.cc
)
