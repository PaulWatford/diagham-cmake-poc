# ============================================================================
# DiagHamHelpers.cmake: helpers for the DiagHam CMake build
# ============================================================================
#
# These macros let each src/<subdir>/CMakeLists.txt be a one-liner, e.g.:
#
#    diagham_add_library(Matrix
#        SOURCES Matrix.cc RealSymmetricMatrix.cc ...
#    )
#
# They produce the same static libraries as the autotools build:
# Base/src/<subdir>/lib<NAME>.a or src/<subdir>/lib<NAME>.a, with names
# matching the autotools targets one-for-one.
#
# This is deliberately a thin layer over add_library / add_executable,
# the modernisation should *not* hide CMake from contributors, only
# make repetitive boilerplate go away.
# ============================================================================

# A registry of all DiagHam libraries built so far, so executables can
# easily link against "everything" without enumerating per-binary.
set_property(GLOBAL PROPERTY DIAGHAM_ALL_LIBS "")

# A registry of every program target, used by tests/ (the --help smoke
# test) and by the install rules.
set_property(GLOBAL PROPERTY DIAGHAM_ALL_PROGRAMS "")

# Every DiagHam target installs under one export set so downstream projects
# can find_package(DiagHam); see cmake/DiagHamInstall.cmake.
include(GNUInstallDirs)

function(diagham_add_library target_name)
    cmake_parse_arguments(ARG "" "" "SOURCES" ${ARGN})

    if(NOT ARG_SOURCES)
        message(FATAL_ERROR "diagham_add_library(${target_name}): no SOURCES given")
    endif()

    add_library(${target_name} STATIC ${ARG_SOURCES})

    # Make the static archive output as lib<target>.a (which matches autotools).
    # CMake does this by default on Linux, but be explicit so the matching is
    # unambiguous when comparing binaries to the autotools build.
    set_target_properties(${target_name} PROPERTIES
        ARCHIVE_OUTPUT_NAME ${target_name}
    )

    # Track this library globally
    get_property(libs GLOBAL PROPERTY DIAGHAM_ALL_LIBS)
    list(APPEND libs ${target_name})
    set_property(GLOBAL PROPERTY DIAGHAM_ALL_LIBS "${libs}")

    install(TARGETS ${target_name}
        EXPORT DiagHamTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/diagham
        COMPONENT Development
    )
endfunction()


# A program is a single-file executable that links against (essentially)
# every static library. We rely on the linker's dead-stripping rather
# than enumerating per-binary LDADD lists, the autotools build already
# does this implicitly via the duplicated -l flags.
function(diagham_add_program target_name source)
    if(NOT DIAGHAM_BUILD_PROGRAMS)
        return()
    endif()

    add_executable(${target_name} ${source})

    # Link against all libraries built so far. The order matters because of
    # mutual recursion between Architecture and ArchitectureOperation, so we
    # link with --start-group/--end-group on GNU ld.
    get_property(all_libs GLOBAL PROPERTY DIAGHAM_ALL_LIBS)

    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        target_link_libraries(${target_name} PRIVATE
            -Wl,--start-group
            ${all_libs}
            -Wl,--end-group
        )
    else()
        target_link_libraries(${target_name} PRIVATE ${all_libs})
    endif()

    # External dependencies (libm, pthreads, LAPACK/MKL, MPI, GSL, GMP,
    # FFTW, bz2, ScaLAPACK, ...) are collected on one interface target in
    # the top-level CMakeLists.txt, so adding an optional library there is
    # the only change needed to link it into every program.
    target_link_libraries(${target_name} PRIVATE diagham_external_deps)

    get_property(programs GLOBAL PROPERTY DIAGHAM_ALL_PROGRAMS)
    list(APPEND programs ${target_name})
    set_property(GLOBAL PROPERTY DIAGHAM_ALL_PROGRAMS "${programs}")

    install(TARGETS ${target_name}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT Runtime
    )
endfunction()


# Whether a program is on the upstream-exclusion list (see
# docs/explanation/deferred-code.md).
function(_diagham_program_excluded prog_name out_var)
    set(${out_var} FALSE PARENT_SCOPE)
    if(DEFINED DIAGHAM_UPSTREAM_EXCLUDED_PROGRAMS)
        list(FIND DIAGHAM_UPSTREAM_EXCLUDED_PROGRAMS "${prog_name}" excluded_idx)
        if(NOT excluded_idx EQUAL -1)
            message(STATUS "DiagHam: skipping ${prog_name} (upstream issue, see patches/PATCHES.md)")
            set(${out_var} TRUE PARENT_SCOPE)
        endif()
    endif()
endfunction()


# Explicit program list, emitted by extract_autotools.py from a Makefile.am
# bin_PROGRAMS line (used for Spin and QuantumDots, whose Programs/
# directories hold sources upstream deliberately does not build):
#
#    diagham_add_programs(PREFIX Spin
#        PROGRAMS
#            Cobalt Cobalt.cc
#            PeriodicQuantumDot2D Periodic2DQuantumDot.cc
#    )
#
# PREFIX keeps target names unique across modules (src/Programs,
# Spin/src/Programs and QuantumDots/src/Programs all share the leaf name
# "Programs"); the output executable keeps the upstream binary name.
function(diagham_add_programs)
    if(NOT DIAGHAM_BUILD_PROGRAMS)
        return()
    endif()
    cmake_parse_arguments(ARG "" "PREFIX" "PROGRAMS" ${ARGN})
    list(LENGTH ARG_PROGRAMS n)
    math(EXPR odd "${n} % 2")
    if(NOT odd EQUAL 0)
        message(FATAL_ERROR "diagham_add_programs: PROGRAMS must be name/source pairs")
    endif()
    while(ARG_PROGRAMS)
        list(POP_FRONT ARG_PROGRAMS prog_name src)
        _diagham_program_excluded(${prog_name} excluded)
        if(excluded)
            continue()
        endif()
        set(target_name "${ARG_PREFIX}_${prog_name}")
        diagham_add_program(${target_name} ${src})
        set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${prog_name})
    endwhile()
endfunction()


# Discover and add every .cc in a Programs/ subdirectory as its own binary.
# This is more robust than enumerating them: adding a new program means
# just dropping a .cc file in the directory, no Makefile.am editing.
# That's one concrete win of the CMake migration.
#
# Some upstream program names collide across directories (e.g. QHEBosonsDelta
# exists in both FQHE/src/Programs/FQHEOnSphere and FQHE/src/Programs/FQHEOnDisk
# as genuinely different programs that happen to share a binary name). In
# autotools this means whichever installs second silently clobbers the first.
# CMake disallows duplicate target names, so we make the *target* name unique
# by prefixing with the leaf directory, while keeping the *output executable*
# name unchanged to match the autotools binary layout.
function(diagham_add_programs_in_directory)
    if(NOT DIAGHAM_BUILD_PROGRAMS)
        return()
    endif()

    get_filename_component(leaf ${CMAKE_CURRENT_SOURCE_DIR} NAME)

    file(GLOB program_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cc)
    foreach(src ${program_sources})
        get_filename_component(prog_name ${src} NAME_WE)

        # Skip files known to be broken upstream where patching is deferred
        # (see patches/PATCHES.md for the rationale on each excluded file).
        _diagham_program_excluded(${prog_name} excluded)
        if(excluded)
            continue()
        endif()

        set(target_name "${leaf}_${prog_name}")
        diagham_add_program(${target_name} ${src})
        set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${prog_name})
    endforeach()
endfunction()
