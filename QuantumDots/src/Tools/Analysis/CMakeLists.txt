# Auto-generated CMakeLists.txt for QuantumDots/src/Tools/Analysis
# Regenerate with: python3 scripts_cmake/extract_autotools.py

# Programs listed in the upstream Makefile.am bin_PROGRAMS.
diagham_add_programs(PREFIX Analysis
    PROGRAMS
        DOS DOS.cc
        Spectrum Spectrum.cc
        SumSpectrum SumSpectrum.cc
        SumDOS SumDOS.cc
        SumDOS2 SumDOS2.cc
        XYReflexionSymmetryPolarizationDegree XYReflexionSymmetryPolarizationDegree.cc
        XYReflexionSymmetryOscillatorForce XYReflexionSymmetryOscillatorForce.cc
        DefaultXYReflexionSymmetryOscillatorForce DefaultXYReflexionSymmetryOscillatorForce.cc
        DefaultCylinderSymmetryOscillatorForce DefaultCylinderSymmetryOscillatorForce.cc
        DefaultCylinderSymmetryInMagneticFieldOscillatorForce DefaultCylinderSymmetryInMagneticFieldOscillatorForce.cc
        DotProbability DotProbability.cc
        OverlapComplex OverlapComplex.cc
        TetrapodDotAnalysis TetrapodDotAnalysis.cc
        PeriodicWaveFunctionValue PeriodicWaveFunctionValue.cc
        XYReflexionSymmetryWaveFunctionValue XYReflexionSymmetryWaveFunctionValue.cc
        PeriodicProbabilityDensity PeriodicProbabilityDensity.cc
        AbsorptionInQuantumWellBField AbsorptionInQuantumWellBField.cc
        EvaluateBroadening EvaluateBroadening.cc
        QuantumWellEscapeTime QuantumWellEscapeTime.cc
)