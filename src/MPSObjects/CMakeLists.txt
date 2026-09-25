# Auto-generated CMakeLists.txt for src/MPSObjects
# Regenerate with: python3 scripts_cmake/extract_autotools.py

diagham_add_library(MPSObjects
    SOURCES
        AbstractMPOperatorOBC.cc
        RealMPOperatorOBC.cc
        ComplexMPOperatorOBC.cc
        AbstractMPSSite.cc
        RealMPSSite.cc
        ComplexMPSSite.cc
        DMRGFiniteSizeRealOBCMainTask.cc
        DMRGFiniteSizeComplexOBCMainTask.cc
        DMRGInfiniteSizeComplexMainTask.cc
        CompletelyPositiveMap.cc
        LeftOrthonormalizeAlgorithm.cc
)
