# Auto-generated CMakeLists.txt for Spin/src/Tools
# Regenerate with: python3 scripts_cmake/extract_autotools.py

add_subdirectory(SpinFiles)

diagham_add_library(SpinTools
    SOURCES
        LatticeConnections.cc
        SUNSpinFileTools.cc
)
