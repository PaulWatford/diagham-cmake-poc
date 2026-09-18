# ============================================================================
#  cmake/CompilerDefaults.cmake
# ============================================================================
#
#  Build-type validation and per-compiler warning flags.
#
#  Adapted (not copied verbatim) from Gunnar Möller & Lars Schonenberg's
#  BDMC_UFL project (cmake/DefaultHeader.cmake, cmake/DefaultFooter.cmake),
#  supplied by Gunnar 18/09 as a template. Two deliberate differences from
#  the BDMC_UFL original, both to avoid silently changing DiagHam-PoC
#  decisions already made elsewhere in this repo:
#    - Does NOT set a default CMAKE_BUILD_TYPE. The top-level CMakeLists.txt
#      already defaults to Release (BDMC_UFL defaults to Debug); this file
#      only validates whatever build type is in effect, it never picks one.
#      Include this file AFTER that default is set.
#    - Does NOT set CMAKE_CXX_FLAGS for the language standard. The top-level
#      CMakeLists.txt already pins CMAKE_CXX_STANDARD 11 (BDMC_UFL sets
#      -std=c++14 directly); duplicating that here risked a silent conflict.
#    - Does NOT add its own configuration-summary footer. CMakeLists.txt
#      already prints one (see "DiagHam CMake configuration:"); adding a
#      second one would just be noise.
#
#  What this file does add, which the PoC had nothing for previously:
#    - A validated list of build types (Debug/Release/RelWithDebInfo/
#      MinSizeRel/RelWithChecks), with a FATAL_ERROR on anything else,
#      matching BDMC_UFL's cmake-gui-friendly pattern.
#    - Per-compiler (GCC/Clang/Intel) warning flags, so build warnings are
#      visible instead of silently suppressed -- useful groundwork given
#      this migration's whole approach is "surface config bugs like the
#      LAPACK/MPI one before they happen" (see patches/PATCHES.md).
# ============================================================================

# ---- Build-type validation --------------------------------------------------
set(DiagHam_ValidBuildTypes "Debug" "Release" "MinSizeRel" "RelWithDebInfo" "RelWithChecks")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${DiagHam_ValidBuildTypes})

string(TOUPPER "${CMAKE_BUILD_TYPE}" DiagHam_UpperBuildType)
set(DiagHam_UpperValidBuildTypes "")
foreach(type IN LISTS DiagHam_ValidBuildTypes)
    string(TOUPPER ${type} upperType)
    list(APPEND DiagHam_UpperValidBuildTypes ${upperType})
endforeach()

list(FIND DiagHam_UpperValidBuildTypes "${DiagHam_UpperBuildType}" DiagHam_BuildTypeFound)
if(DiagHam_BuildTypeFound LESS 0)
    message(FATAL_ERROR
        "Unknown CMAKE_BUILD_TYPE '${CMAKE_BUILD_TYPE}'. Use one of: ${DiagHam_ValidBuildTypes}")
endif()

# ---- Per-compiler warning flags --------------------------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    message(STATUS "Using warning flags for Intel compiler")
    # -w2: warning level 2 (matches Gunnar's BDMC_UFL default)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w2")

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(STATUS "Using warning flags for Clang compiler")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wconditional-uninitialized")

elseif(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "Using warning flags for g++-compatible compilers")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

else()
    message(STATUS "Unrecognised compiler (${CMAKE_CXX_COMPILER_ID}) -- no warning flags set")
endif()
