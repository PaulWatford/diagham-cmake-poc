# ============================================================================
# DiagHamInstall.cmake: install rules and the find_package(DiagHam) export
# ============================================================================
#
# Included from the top-level CMakeLists.txt after every module has been
# added, so the DIAGHAM_ALL_LIBS registry is complete.
#
# Layout under CMAKE_INSTALL_PREFIX:
#
#   bin/                          every program (component Runtime)
#   lib/diagham/lib<NAME>.a       every static library (component Development)
#   include/diagham/...           every public header, one merged tree
#                                 (component Development)
#   lib/cmake/DiagHam/            DiagHamConfig.cmake + exported targets
#
# The merged header tree is safe: no header path relative to its module
# root (src/, Base/src/, FQHE/src/, ...) is duplicated across modules
# (checked against upstream), so "Vector/RealVector.h" and
# "HilbertSpace/FermionOnSphere.h" resolve exactly as they do in-tree.
#
# Downstream use:
#
#   find_package(DiagHam REQUIRED)
#   target_link_libraries(myprog PRIVATE DiagHam::DiagHam)
#
# DiagHam::DiagHam carries the include path, the compile definitions DiagHam
# headers need (HAVE_CONFIG_H, MACHINE_PRECISION, ...), every static library
# (wrapped in --start-group/--end-group on GNU-style linkers, because the
# libraries are mutually recursive), and the external dependencies.
# ============================================================================

include(CMakePackageConfigHelpers)

set(DIAGHAM_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR}/diagham)
set(DIAGHAM_INSTALL_CMAKEDIR   ${CMAKE_INSTALL_LIBDIR}/cmake/DiagHam)

# Module roots that are part of this configuration.
set(DIAGHAM_MODULE_ROOTS src Base/src)
if(DIAGHAM_BUILD_QUANTUMDOTS)
    list(APPEND DIAGHAM_MODULE_ROOTS QuantumDots/src)
endif()
if(DIAGHAM_BUILD_SPIN)
    list(APPEND DIAGHAM_MODULE_ROOTS Spin/src)
endif()
if(DIAGHAM_BUILD_FQHE)
    list(APPEND DIAGHAM_MODULE_ROOTS FQHE/src)
endif()
if(DIAGHAM_BUILD_FTI)
    list(APPEND DIAGHAM_MODULE_ROOTS FTI/src)
endif()

# ----------------------------------------------------------------------------
# Headers
# ----------------------------------------------------------------------------
foreach(root ${DIAGHAM_MODULE_ROOTS})
    install(DIRECTORY ${CMAKE_SOURCE_DIR}/${root}/
        DESTINATION ${DIAGHAM_INSTALL_INCLUDEDIR}
        COMPONENT Development
        FILES_MATCHING PATTERN "*.h"
        PATTERN "Programs" EXCLUDE
    )
endforeach()
install(FILES ${CMAKE_BINARY_DIR}/src/config_ac.h
    DESTINATION ${DIAGHAM_INSTALL_INCLUDEDIR}
    COMPONENT Development
)

# ----------------------------------------------------------------------------
# The umbrella target
# ----------------------------------------------------------------------------
get_property(_diagham_libs GLOBAL PROPERTY DIAGHAM_ALL_LIBS)

add_library(diagham INTERFACE)
add_library(DiagHam::DiagHam ALIAS diagham)
set_target_properties(diagham PROPERTIES EXPORT_NAME DiagHam)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_link_libraries(diagham INTERFACE
        -Wl,--start-group ${_diagham_libs} -Wl,--end-group)
else()
    # Linkers without group support (e.g. Apple ld) rescan archives anyway.
    target_link_libraries(diagham INTERFACE ${_diagham_libs})
endif()
target_link_libraries(diagham INTERFACE diagham_external_deps)

set(_diagham_build_includes ${CMAKE_BINARY_DIR}/src)
foreach(root ${DIAGHAM_MODULE_ROOTS})
    list(APPEND _diagham_build_includes ${CMAKE_SOURCE_DIR}/${root})
endforeach()
target_include_directories(diagham INTERFACE
    "$<BUILD_INTERFACE:${_diagham_build_includes}>"
    "$<INSTALL_INTERFACE:${DIAGHAM_INSTALL_INCLUDEDIR}>"
)
target_compile_definitions(diagham INTERFACE HAVE_CONFIG_H ${DIAGHAM_COMPILE_DEFS})
target_compile_features(diagham INTERFACE cxx_std_11)

install(TARGETS diagham diagham_external_deps
    EXPORT DiagHamTargets
    COMPONENT Development
)
install(EXPORT DiagHamTargets
    NAMESPACE DiagHam::
    DESTINATION ${DIAGHAM_INSTALL_CMAKEDIR}
    COMPONENT Development
)

# ----------------------------------------------------------------------------
# Package config + version files
# ----------------------------------------------------------------------------
configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/cmake/DiagHamConfig.cmake.in
    ${CMAKE_BINARY_DIR}/DiagHamConfig.cmake
    INSTALL_DESTINATION ${DIAGHAM_INSTALL_CMAKEDIR}
)
write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/DiagHamConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(FILES
    ${CMAKE_BINARY_DIR}/DiagHamConfig.cmake
    ${CMAKE_BINARY_DIR}/DiagHamConfigVersion.cmake
    DESTINATION ${DIAGHAM_INSTALL_CMAKEDIR}
    COMPONENT Development
)
# Find modules the config may need to re-run for its dependencies.
install(FILES
    ${CMAKE_SOURCE_DIR}/cmake/FindGMP.cmake
    ${CMAKE_SOURCE_DIR}/cmake/FindFFTW3.cmake
    ${CMAKE_SOURCE_DIR}/cmake/FindMKL.cmake
    DESTINATION ${DIAGHAM_INSTALL_CMAKEDIR}
    COMPONENT Development
)
