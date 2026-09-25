# Auto-generated CMakeLists.txt for src/HilbertSpace
# Regenerate with: python3 scripts_cmake/extract_autotools.py

add_subdirectory(DMRGHilbertSpace)
add_subdirectory(ManyBodyHilbertSpace)

diagham_add_library(HilbertSpace
    SOURCES
        AbstractHilbertSpace.cc
        SpaceDecomposition.cc
        SubspaceSpaceConverter.cc
        UndescribedHilbertSpace.cc
)
