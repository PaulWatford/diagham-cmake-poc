# Auto-generated CMakeLists.txt for src/MathTools
# Regenerate with: python3 scripts_cmake/extract_autotools.py

add_subdirectory(NumericalAnalysis)
add_subdirectory(RandomNumber)

diagham_add_library(MathTools
    SOURCES
        FactorialCoefficient.cc
        ClebschGordanCoefficients.cc
        IntegerAlgebraTools.cc
        BinomialCoefficients.cc
        GroupedPermutations.cc
        PermutationElement.cc
        SimplePermutations.cc
        ClebschGordanDiskCoefficients.cc
        ArithmeticGeometricMean.cc
        JacobiThetaFunction.cc
        WignerSmallDMatrix.cc
        ThreeJSymbol.cc
        PartialThreeJSymbol.cc
        KCombinations.cc
        JacobiPolynomials.cc
        MonopoleHarmonics.cc
        NSphereParameters.cc
        Rational.cc
        LongRational.cc
        JacobiThetaDerivativeFunction.cc
        SU3ClebschGordanCoefficients.cc
        SU3IrreducibleRepresentations.cc
        ClebschGordanThinAnnulusCoefficients.cc
)
