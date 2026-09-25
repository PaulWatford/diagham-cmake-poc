# Auto-generated CMakeLists.txt for src/MathTools/RandomNumber
# Regenerate with: python3 scripts_cmake/extract_autotools.py

diagham_add_library(RandomNumber
    SOURCES
        AbstractRandomNumberGenerator.cc
        StdlibRandomNumberGenerator.cc
        RanluxRandomNumberGenerator.cc
        NumRecRandomGenerator.cc
        FileRandomNumberGenerator.cc
)
