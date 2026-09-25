# Auto-generated CMakeLists.txt for src/Architecture
# Regenerate with: python3 scripts_cmake/extract_autotools.py

add_subdirectory(ArchitectureOperation)
add_subdirectory(ClusterArchitecture)

diagham_add_library(Architecture
    SOURCES
        AbstractArchitecture.cc
        MonoProcessorArchitecture.cc
        SimpleMPIArchitecture.cc
        SMPArchitecture.cc
        MixedMPISMPArchitecture.cc
        ArchitectureManager.cc
)
