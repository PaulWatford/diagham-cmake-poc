# Auto-generated CMakeLists.txt for src/Interaction
# Regenerate with: python3 scripts_cmake/extract_autotools.py

add_subdirectory(InternalInteraction)

diagham_add_library(Interaction
    SOURCES
        AbstractInteraction.cc
        BasicInteraction.cc
        DiamondInteraction.cc
)
