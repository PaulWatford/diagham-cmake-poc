# ============================================================================
#  DiagHam: CMake build system (modernisation proof-of-concept)
# ============================================================================
#
#  This CMakeLists.txt is a proof-of-concept replacement for DiagHam's
#  legacy autotools build.
#
#  Scope of this PoC:
#    - All static libraries under Base/src/ and src/ (the core libraries)
#    - All 34 generic executables under src/Programs/
#    - FQHE module: libraries and programs across FQHEOnSphere, FQHEOnDisk,
#      FQHEOnTorus, FQHEOnCylinder, FQHEOnLattice
#    - FTI module: libraries and Hubbard / FCI / FTI executables (FTI's
#      hamiltonians and operators are derived from FQHE's, so FTI needs
#      FQHE built and linked)
#
#    - Spin module (spin chains, SU(N) spins, PEPS/MPS tools)
#    - QuantumDots module (single-particle quantum dot / well programs)
#    - Every autotools configure flag except the defunct --enable-m64 (see
#      CONFIGURE_FLAGS.md for the one-to-one mapping)
#    - ctest suite (tests/), install + find_package(DiagHam) export,
#      CMakePresets.json
#
#  Out of scope (deferred to a real production CMake rollout):
#    - Anyons module: referenced by configure.ac but has no sources in the
#      tree (its AC_CONFIG_FILES entry is commented out upstream)
#    - Spack/EasyBuild packaging metadata
#
#  Author: Paul Watford
# ============================================================================

cmake_minimum_required(VERSION 3.16)

project(DiagHam
    VERSION 0.01
    DESCRIPTION "Exact Diagonalization toolkit for strongly correlated systems"
    LANGUAGES CXX C
)

# Out-of-source builds only (AGENTS.md lists an in-source build as Damage).
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed. Use: cmake -S . -B build  (or a preset)")
endif()

# cmake/ holds our own Find*.cmake modules, including FindMKL.cmake taken
# (with attribution) from Gunnar Möller & Lars Schonenberg's BDMC_UFL
# project, which already solves MKL-as-LAPACK/BLAS/FFTW3-provider for this
# group's own machines. See cmake/FindMKL.cmake header for provenance.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Optional Kent-cluster site defaults, modelled on Gunnar's TCMDefaults.cmake
# pattern (off by default; see cmake/KentDefaults.cmake -- it's a scaffold,
# not filled in with real Kent paths yet).
option(DIAGHAM_KENT_DEFAULTS "Use Kent-cluster site defaults (see cmake/KentDefaults.cmake)" OFF)
if(DIAGHAM_KENT_DEFAULTS)
    include(KentDefaults)
endif()

# ----------------------------------------------------------------------------
# Build configuration
# ----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)   # gnu++11, as upstream's g++ default (M_PI on non-glibc)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Build-type validation + per-compiler warning flags, adapted from Gunnar's
# BDMC_UFL DefaultHeader.cmake (attribution and rationale in the file
# itself). Included here, after the default above, because it validates
# whatever CMAKE_BUILD_TYPE is already set rather than choosing one.
include(CompilerDefaults)

# Flags stated explicitly rather than inheriting CMake's Release (-O3 -DNDEBUG):
# autotools' default is " -O2 -Wall " with no NDEBUG (DiagHam's asserts stay live);
# --enable-debug appends "-g -fPIC" to that, which RelWithDebInfo mirrors.
set(CMAKE_CXX_FLAGS_RELEASE        "-O2"    CACHE STRING "Release flags (autotools default -O2)" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING "Optimised with symbols (autotools --enable-debug)" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHCHECKS  "-O2 -g" CACHE STRING "RelWithChecks (was accepted but undefined)" FORCE)

# ----------------------------------------------------------------------------
# Options, mirroring the autotools AC_ARG_ENABLE / AC_ARG_WITH flags.
# CONFIGURE_FLAGS.md maps every configure.ac flag to its CMake equivalent.
# ----------------------------------------------------------------------------
option(DIAGHAM_BUILD_FQHE     "Build the FQHE module (sphere/disk/torus/cylinder/lattice)" ON)
option(DIAGHAM_BUILD_FTI      "Build the FTI module (Hubbard models, FCI, FTI; requires FQHE)" ON)
option(DIAGHAM_BUILD_SPIN     "Build the Spin module (spin chains, SU(N) spins, PEPS/MPS)" ON)
option(DIAGHAM_BUILD_QUANTUMDOTS "Build the QuantumDots module" ON)
option(DIAGHAM_BUILD_TESTS    "Register the ctest suite in tests/ (needs DIAGHAM_BUILD_PROGRAMS)" ON)

# FTI's hamiltonians and operators inherit from FQHE's, so FTI cannot build
# or link without FQHE. Enforce the dependency at configure time.
if(DIAGHAM_BUILD_FTI AND NOT DIAGHAM_BUILD_FQHE)
    message(FATAL_ERROR
        "DIAGHAM_BUILD_FTI requires DIAGHAM_BUILD_FQHE. "
        "FTI links against QHEHamiltonian, QHEHilbertSpace, QHEOperator, "
        "QHEArchitectureOperation, QHEMainTask, QHEFiles, QHEFunctionBasis "
        "and QHEQuantumNumber. Re-configure with -DDIAGHAM_BUILD_FQHE=ON.")
endif()
option(DIAGHAM_USE_SMP        "Enable pthread-based shared-memory parallelism" ON)
option(DIAGHAM_USE_DEBUG      "Define __DEBUG__ for runtime sanity checks" ON)
option(DIAGHAM_USE_LAPACK     "Enable LAPACK support (optional)" OFF)
option(DIAGHAM_USE_MKL        "When LAPACK is enabled, get BLAS/LAPACK (and FFTW3) from Intel MKL instead of plain find_package(LAPACK) -- matches Gunnar's BDMC_UFL precedent" OFF)
option(DIAGHAM_USE_MPI        "Enable MPI support (optional)" OFF)
option(DIAGHAM_BUILD_PROGRAMS "Build executables (not just libraries)" ON)
option(DIAGHAM_USE_LAPACK_ONLY "Use LAPACK routines instead of DiagHam's own (--enable-lapack-only; requires LAPACK)" OFF)
option(DIAGHAM_USE_BZ2        "Enable bzip2 support (--enable-bz2)" OFF)
option(DIAGHAM_USE_GMP        "Enable GMP arbitrary-precision support (--enable-gmp)" OFF)
option(DIAGHAM_USE_MPACK      "Enable MPACK arbitrary-precision LAPACK (--enable-mpack; requires LAPACK)" OFF)
option(DIAGHAM_USE_FFTW       "Enable FFTW3 support (--enable-fftw)" OFF)
option(DIAGHAM_USE_SCALAPACK  "Enable ScaLAPACK support (--enable-scalapack; requires MPI and LAPACK)" OFF)
option(DIAGHAM_USE_PROFILE    "Build with gprof instrumentation, -g -pg (--enable-profile)" OFF)
set(DIAGHAM_COMMAND_LOG "" CACHE STRING
    "Global logfile where DiagHam commands are stored (--with-cmd-log); empty disables it")

# Library overrides, the equivalents of --with-<pkg>-libs="...". When set,
# they replace what find_package would pick, exactly like the configure
# flags do. Each is a ;-list (or space-separated string) of link items.
set(DIAGHAM_LAPACK_LIBS    "" CACHE STRING "Override LAPACK+BLAS link line (--with-lapack-libs/--with-blas-libs)")
set(DIAGHAM_GSL_LIBS       "" CACHE STRING "Override GSL link line (--with-gsl-libs)")
set(DIAGHAM_GMP_LIBS       "" CACHE STRING "Override GMP link line (--with-gmp-libs)")
set(DIAGHAM_MPACK_LIBS     "mlapack_gmp;mblas_gmp" CACHE STRING "MPACK link line (--with-mpack-libs)")
set(DIAGHAM_SCALAPACK_LIBS "" CACHE STRING "Override ScaLAPACK link line (--with-scalapack-libs)")
set(DIAGHAM_FORTRAN2C_LIBS "" CACHE STRING "Extra Fortran runtime libraries (--with-fortran2c-libs)")

if(DIAGHAM_USE_LAPACK_ONLY AND NOT DIAGHAM_USE_LAPACK)
    message(FATAL_ERROR "DIAGHAM_USE_LAPACK_ONLY requires DIAGHAM_USE_LAPACK=ON")
endif()
if(DIAGHAM_USE_MPACK AND NOT DIAGHAM_USE_LAPACK)
    message(FATAL_ERROR "DIAGHAM_USE_MPACK requires DIAGHAM_USE_LAPACK=ON (configure.ac appends MPACK to the LAPACK link line)")
endif()
if(DIAGHAM_USE_SCALAPACK AND NOT (DIAGHAM_USE_MPI AND DIAGHAM_USE_LAPACK))
    message(FATAL_ERROR "DIAGHAM_USE_SCALAPACK requires DIAGHAM_USE_MPI=ON and DIAGHAM_USE_LAPACK=ON (same rule as configure.ac)")
endif()
if(DIAGHAM_BUILD_TESTS AND NOT DIAGHAM_BUILD_PROGRAMS)
    message(STATUS "DiagHam: DIAGHAM_BUILD_PROGRAMS=OFF, so the ctest suite is not registered")
    set(DIAGHAM_BUILD_TESTS OFF)
endif()

# Every external library any DiagHam program links against is attached to
# this one interface target; diagham_add_program() links it, and the
# exported DiagHam::DiagHam target forwards it to downstream projects.
add_library(diagham_external_deps INTERFACE)
target_link_libraries(diagham_external_deps INTERFACE m)

# ----------------------------------------------------------------------------
# Compile-time macros for src/config.h, corresponds to autoconf's config_ac.h
# ----------------------------------------------------------------------------
set(DIAGHAM_COMPILE_DEFS "")

# Endianness
#
# REMOVED 22/09 (audited by a second Claude session, same bug class as
# the HAVE_LAPACK/HAVE_MPI fix below): this used to inject __LITTLEENDIAN__
# / __BIGENDIAN__ directly via compile definitions. Confirmed by reading
# cmake/config_ac.h.in: HAVE_LITTLEENDIAN is unconditionally `#define`d
# there (not gated by any real endianness detection upstream -- this PoC
# is x86/x64-only regardless), so src/config.h's own
# `#ifdef HAVE_LITTLEENDIAN -> #define __LITTLEENDIAN__` chain already
# defines it on every supported platform. The direct injection here was
# fully redundant with that and produced a macro-redefinition warning on
# every translation unit (confirmed: real Vector-library build showed
# these warnings before this fix). No HAVE_BIGENDIAN path exists to
# preserve, so nothing is lost by removing the TestBigEndian block too.
include(TestBigEndian)
test_big_endian(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
    message(WARNING "DiagHam: big-endian target detected, but this PoC's "
        "config header (cmake/config_ac.h.in) unconditionally sets "
        "HAVE_LITTLEENDIAN with no HAVE_BIGENDIAN path -- big-endian is "
        "not actually supported here despite this check existing.")
endif()

# __64_BITS__ is NOT injected here: upstream src/config.h derives it from the
# compiler's own macros (__x86_64__, __LP64__, __aarch64__); injecting it from
# CMAKE_SIZEOF_VOID_P duplicated that on Linux and was wrong on LLP64 (Win64).

# __int128 support, gcc/clang have it, MSVC does not
#
# REMOVED 22/09 (audited by a second Claude session): the __128_BIT_LONGLONG__
# injection here was removed for the same reason as endianness above.
# Confirmed by reading upstream src/config.h: on x86_64 with a
# non-Intel compiler (i.e. GCC/Clang, what this PoC actually builds
# with) it unconditionally does
#   typedef int int128_t __attribute__((__mode__(TI)));
#   #define __128_BIT_LONGLONG__
#   #define LONGLONG int128_t
# with no HAVE_* gate at all, so the CMake-side injection was a pure
# duplicate, producing a macro-redefinition warning. The __int128
# detection below (CheckTypeSize) is kept only as a build-time sanity
# check that the compiler actually supports it -- it no longer defines
# anything, since config.h always will on the platforms this PoC targets.
include(CheckTypeSize)
check_type_size("__int128" SIZEOF_INT128 LANGUAGE CXX)
if(NOT HAVE_SIZEOF_INT128)
    message(WARNING "DiagHam: compiler lacks __int128 support, but "
        "upstream src/config.h unconditionally assumes it on x86_64 with "
        "a non-Intel compiler -- the build may fail or misbehave.")
endif()

# sstream style (modern C++ <sstream> rather than ancient strstream)
#
# REMOVED 22/09 (audited by a second Claude session): src/config.h
# unconditionally does `#define __SSTREAM_STYLE__` itself, with no
# HAVE_* gate -- this injection was a pure duplicate causing the same
# macro-redefinition warning as the two removals above.

# Machine precision constant, matches autotools' MACHINE_PRECISION
list(APPEND DIAGHAM_COMPILE_DEFS "MACHINE_PRECISION=1e-14")

# Cluster architecture and Hilbert space "feature" macros, these come
# from src/config.h via the HAVE_CONFIG_H -> config.h -> config_ac.h chain.
# We do NOT redefine them here, because some source files (e.g.
# DelocalizedRealVector.cc) use them in an #ifdef without including
# config.h first, so they're effectively dead code in the autotools
# build, and we preserve that behaviour to match output byte-for-byte.
#
# A production migration would either:
#   (a) fix the offending files to #include "config.h" properly, then
#       enable USE_CLUSTER_ARCHITECTURE here, or
#   (b) remove the dead-code paths entirely.
# This PoC chooses neither, it just matches what autotools produces.

# FLAGGED, NOT FIXED, 22/09 (audited by a second Claude session): __DEBUG__
# and __SMP__ below are the same bug class as the three removed above,
# but a full fix is more involved, so this is documented rather than
# changed this pass. Confirmed by reading upstream src/config.h: it
# unconditionally does `#define __DEBUG__` and `#define __SMP__`, with
# no HAVE_* gate on either -- so DIAGHAM_USE_DEBUG=OFF or
# DIAGHAM_USE_SMP=OFF does NOT actually turn off either macro; config.h
# defines both regardless, and the injections below just duplicate one
# of the two when the option is ON (another macro-redefinition warning),
# doing nothing when the option is OFF. DIAGHAM_USE_SMP still meaningfully
# controls whether Threads::Threads gets linked, so it isn't a no-op
# option overall, just a no-op for the macro it claims to toggle. A real
# fix needs a new patch to upstream src/config.h to gate these behind
# HAVE_DEBUG/HAVE_SMP the way HAVE_LAPACK/HAVE_MPI already are, which is
# a bigger change than this pass's scope (documentation + the three safe
# removals above only).
if(DIAGHAM_USE_DEBUG)
    # (__DEBUG__ itself is defined unconditionally by src/config.h; nothing to inject)
endif()

# Pthread (DiagHam's misleadingly-named "SMP" mode is actually pthreads,
# not OpenMP, see briefing document section 2)
if(DIAGHAM_USE_SMP)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    # (__SMP__ itself is defined unconditionally by src/config.h; nothing to inject)
    target_link_libraries(diagham_external_deps INTERFACE Threads::Threads)
endif()

# MPI
#
# BUG FIX (same class as the LAPACK fix below): was injecting __MPI__
# directly, bypassing src/config.h's "#ifdef HAVE_MPI -> #define __MPI__"
# block. Dormant today only because DIAGHAM_USE_MPI defaults OFF, but it
# would misbehave the moment someone turns MPI on, the same way LAPACK did.
if(DIAGHAM_USE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    set(HAVE_MPI 1)
    # Compile, not just link, against MPI: sources with __MPI__ blocks
    # include <mpi.h>. (--with-mpi-cxx maps to CMAKE_CXX_COMPILER=mpicxx,
    # which also works; this path works with a plain compiler too.)
    include_directories(${MPI_CXX_INCLUDE_DIRS})
    target_link_libraries(diagham_external_deps INTERFACE MPI::MPI_CXX)
endif()

# LAPACK
#
# BUG FIX (see patches/PATCHES.md): this used to inject the compile
# definition __LAPACK__ directly, bypassing src/config.h's own
# "#ifdef HAVE_LAPACK -> #define __LAPACK__ -> typedef doublecomplex..."
# block. That block is where upstream DiagHam defines `doublecomplex`
# (src/config.h, guarded on HAVE_LAPACK). Injecting __LAPACK__ directly
# turned on LAPACK-consuming code (e.g. HermitianMatrix.h's
# `#ifdef __LAPACK__` member block, which uses doublecomplex) WITHOUT
# ever defining HAVE_LAPACK, so the typedef never happened -> the
# ~5000 "doublecomplex undeclared" compile errors. This was a PoC
# translation bug, not an upstream DiagHam defect: LAPACK is disabled
# by default upstream too (configure.ac: HAVE_PKG_LAPACK=no unless
# --enable-lapack), and doublecomplex is only ever needed once LAPACK
# is actually on. The fix: set HAVE_LAPACK for configure_file(), the
# same way HAVE_FQHE/HAVE_FTI are set above, and let src/config.h do
# what it already does correctly.
#
# Provider choice (new): DIAGHAM_USE_MKL switches between two ways of
# satisfying LAPACK, both producing the same DIAGHAM_LAPACK_LIBRARIES
# variable that DiagHamHelpers.cmake links against:
#   OFF (default) - plain find_package(LAPACK), matching what upstream
#                    DiagHam's own configure.ac does (--enable-lapack
#                    links -llapack -lblas directly).
#   ON             - find_package(MKL) via cmake/FindMKL.cmake (taken from
#                    Gunnar's BDMC_UFL), which supplies BLAS+LAPACK+FFTW3
#                    together. This is his own group's established pattern
#                    (see BDMC_UFL/CMakeLists.txt: USE_MKL defaults ON
#                    there), offered here as an option, not a replacement,
#                    since DiagHam's own default path is plain LAPACK.
if(DIAGHAM_USE_LAPACK)
    set(HAVE_LAPACK 1)
    if(DIAGHAM_USE_MKL)
        set(MKL_MULTI_THREADED 0 CACHE BOOL "Use MKL library in multi-threaded mode")
        find_package(MKL REQUIRED COMPONENTS FFTW3)
        if(MKL_FOUND)
            set(HAVE_INTELMKL 1)
            include_directories("${MKL_INCLUDE_DIRS}" "${MKL_FFTW3_INCLUDE_DIRS}")
            set(DIAGHAM_LAPACK_LIBRARIES ${MKL_LIBRARIES})
            message(STATUS "LAPACK/BLAS provider: Intel MKL (${MKL_VERSION})")
        endif()
    elseif(DIAGHAM_LAPACK_LIBS)
        separate_arguments(DIAGHAM_LAPACK_LIBRARIES UNIX_COMMAND "${DIAGHAM_LAPACK_LIBS}")
        message(STATUS "LAPACK/BLAS provider: user-supplied DIAGHAM_LAPACK_LIBS")
    else()
        find_package(LAPACK REQUIRED)
        set(DIAGHAM_LAPACK_LIBRARIES ${LAPACK_LIBRARIES})
        message(STATUS "LAPACK/BLAS provider: plain find_package(LAPACK)")
    endif()
    target_link_libraries(diagham_external_deps INTERFACE ${DIAGHAM_LAPACK_LIBRARIES})
    if(DIAGHAM_USE_LAPACK_ONLY)
        set(HAVE_LAPACK_ONLY 1)
    endif()
endif()

# GSL
#
# Highest-priority unimplemented library per the parity audit (20 files
# use it, more than any other optional library -- see patches/PATCHES.md).
# Unlike Gunnar's BDMC_UFL, DiagHam's own GSL usage never touches GSL's
# cblas interface (checked: DiagHam includes gsl_errno/gsl_integration/
# gsl_randist/gsl_rng/gsl_roots/gsl_sf_airy/gsl_sf_gamma/gsl_sf_hyperg/
# gsl_spline -- root finding, integration, random numbers, special
# functions, splines; no gsl_blas.h, no cblas anywhere), so BDMC_UFL's
# EXTERNAL_CBLAS/mkl_cblas.h coexistence patch (cmake/gsl_external_cblas.patch,
# copied in alongside the other Find*.cmake modules for future reference)
# solves a problem DiagHam doesn't have. Wired simply instead: plain
# find_package(GSL), no cblas patching, whether or not DIAGHAM_USE_MKL
# is also on.
option(DIAGHAM_USE_GSL "Enable GSL support (optional)" OFF)
if(DIAGHAM_USE_GSL)
    set(HAVE_GSL 1)
    if(DIAGHAM_GSL_LIBS)
        separate_arguments(_gsl_libs UNIX_COMMAND "${DIAGHAM_GSL_LIBS}")
        target_link_libraries(diagham_external_deps INTERFACE ${_gsl_libs})
    else()
        find_package(GSL 2.0 REQUIRED)
        include_directories(${GSL_INCLUDE_DIRS})
        target_link_libraries(diagham_external_deps INTERFACE ${GSL_LIBRARIES})
        message(STATUS "GSL: ${GSL_VERSION} (${GSL_INCLUDE_DIRS})")
    endif()
endif()

# bzip2 (--enable-bz2)
if(DIAGHAM_USE_BZ2)
    find_package(BZip2 REQUIRED)
    set(HAVE_BZ2 1)
    target_link_libraries(diagham_external_deps INTERFACE BZip2::BZip2)
    include_directories(${BZIP2_INCLUDE_DIRS})
endif()

# GMP (--enable-gmp). configure.ac links "-lgmp -lgmpxx"; DiagHam uses the
# C++ interface (gmpxx.h), so gmpxx is required, not optional.
if(DIAGHAM_USE_GMP)
    set(HAVE_GMP 1)
    if(DIAGHAM_GMP_LIBS)
        separate_arguments(_gmp_libs UNIX_COMMAND "${DIAGHAM_GMP_LIBS}")
    else()
        find_package(GMP REQUIRED)
        find_library(GMPXX_LIBRARIES NAMES gmpxx HINTS ${GMP_LIBRARIES_DIR} REQUIRED)
        include_directories(${GMP_INCLUDE_DIR})
        set(_gmp_libs ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
    endif()
    target_link_libraries(diagham_external_deps INTERFACE ${_gmp_libs})
endif()

# MPACK (--enable-mpack): multiple-precision LAPACK, appended to the LAPACK
# link line exactly as configure.ac does. No upstream find module exists
# for it, so DIAGHAM_MPACK_LIBS is the link line (default matches
# configure.ac's "-lmlapack_gmp -lmblas_gmp").
if(DIAGHAM_USE_MPACK)
    set(HAVE_MPACK 1)
    separate_arguments(_mpack_libs UNIX_COMMAND "${DIAGHAM_MPACK_LIBS}")
    target_link_libraries(diagham_external_deps INTERFACE ${_mpack_libs})
endif()

# FFTW3 (--enable-fftw). MKL already provides the FFTW3 interface when
# DIAGHAM_USE_MKL is on, so only look for libfftw3 otherwise.
if(DIAGHAM_USE_FFTW)
    set(HAVE_FFTW 1)
    if(NOT (DIAGHAM_USE_LAPACK AND DIAGHAM_USE_MKL))
        find_package(FFTW3 REQUIRED)
        include_directories(${FFTW3_INCLUDE_DIR})
        target_link_libraries(diagham_external_deps INTERFACE ${FFTW3_LIBRARIES})
    endif()
endif()

# ScaLAPACK (--enable-scalapack). Prefers a CMake package config (as shipped
# by Netlib ScaLAPACK >= 2.1 and most HPC module stacks), then a plain
# library search, then the DIAGHAM_SCALAPACK_LIBS override.
if(DIAGHAM_USE_SCALAPACK)
    set(HAVE_SCALAPACK 1)
    if(DIAGHAM_SCALAPACK_LIBS)
        separate_arguments(_scalapack_libs UNIX_COMMAND "${DIAGHAM_SCALAPACK_LIBS}")
    else()
        find_package(scalapack CONFIG QUIET)
        if(TARGET scalapack)
            set(_scalapack_libs scalapack)
        else()
            find_library(SCALAPACK_LIBRARY NAMES scalapack scalapack-openmpi scalapack-mpich REQUIRED)
            set(_scalapack_libs ${SCALAPACK_LIBRARY})
        endif()
    endif()
    target_link_libraries(diagham_external_deps INTERFACE ${_scalapack_libs})
endif()

if(DIAGHAM_FORTRAN2C_LIBS)
    separate_arguments(_f2c_libs UNIX_COMMAND "${DIAGHAM_FORTRAN2C_LIBS}")
    target_link_libraries(diagham_external_deps INTERFACE ${_f2c_libs})
endif()

# Profiling (--enable-profile adds -g -pg to CXXFLAGS)
if(DIAGHAM_USE_PROFILE)
    add_compile_options(-g -pg)
    target_link_options(diagham_external_deps INTERFACE -pg)
endif()

# Global command log (--with-cmd-log=FILE)
if(DIAGHAM_COMMAND_LOG)
    set(HAVE_GLOBAL_COMMAND_LOG 1)
    set(GLOBAL_COMMAND_LOG "${DIAGHAM_COMMAND_LOG}")
endif()

# Apply globally to every target in this project
add_compile_definitions(${DIAGHAM_COMPILE_DEFS})

# DiagHam expects HAVE_CONFIG_H so it can include "config.h"
add_compile_definitions(HAVE_CONFIG_H)

# Suppress noise from the legacy codebase: many warnings are spurious in
# 24-year-old C++98-ish code. A production migration would clean these.
add_compile_options(-Wno-deprecated -Wno-write-strings -Wno-narrowing -Wno-format)

# ----------------------------------------------------------------------------
# Include directories, every source file uses paths like
#   #include "Architecture/AbstractArchitecture.h"
# so src/ and Base/src/ must both be on the search path.
# ----------------------------------------------------------------------------
include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/Base/src
    ${CMAKE_BINARY_DIR}/src              # for the generated config.h
)

if(DIAGHAM_BUILD_FQHE)
    # FQHE sources include their own headers via paths like
    #   #include "Architecture/ArchitectureOperation/FQHEDiskQuasiholePropagatorOperation.h"
    # so FQHE/src must be on the search path.
    include_directories(${CMAKE_SOURCE_DIR}/FQHE/src)
endif()

if(DIAGHAM_BUILD_FTI)
    # FTI sources include their own headers via paths like
    #   #include "Architecture/ArchitectureOperation/FTIComputeBandStructureOperation.h"
    # and reach into FQHE for the base classes their hamiltonians derive from.
    #
    # NOTE: this is a real fix relative to the autotools build, which sets
    # only -I src and -I Base/src on FTI compilations. That's a long-standing
    # bug: `make -C FTI/src` fails out of the box on a clean autotools setup.
    # The CMake migration is the right moment to fix this.
    include_directories(${CMAKE_SOURCE_DIR}/FTI/src)
endif()

if(DIAGHAM_BUILD_QUANTUMDOTS)
    include_directories(${CMAKE_SOURCE_DIR}/QuantumDots/src)
endif()

if(DIAGHAM_BUILD_SPIN)
    include_directories(${CMAKE_SOURCE_DIR}/Spin/src)
endif()

# ----------------------------------------------------------------------------
# Generate the config.h header that DiagHam includes everywhere
#
# We surface the module-presence macros HAVE_FQHE / HAVE_FTI here, mirroring
# what autotools does on `--enable-fqhe` / `--enable-fti`. These macros gate
# FTI-aware constructor declarations in shared FQHE headers, so they must be
# defined whenever the corresponding modules are built (otherwise downstream
# FCI programs see only the copy-constructor and fail to compile).
# ----------------------------------------------------------------------------
if(DIAGHAM_BUILD_FQHE)
    set(HAVE_FQHE 1)
endif()
if(DIAGHAM_BUILD_FTI)
    set(HAVE_FTI 1)
endif()

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/src)
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/config_ac.h.in
    ${CMAKE_BINARY_DIR}/src/config_ac.h
    @ONLY
)

# ----------------------------------------------------------------------------
# A helper that wraps add_library and add_executable so source files can be
# auto-discovered while preserving the autotools library granularity.
# ----------------------------------------------------------------------------
include(${CMAKE_SOURCE_DIR}/cmake/DiagHamHelpers.cmake)

# ----------------------------------------------------------------------------
# Apply upstream-bug patches and set the exclusion list. See patches/PATCHES.md
# for what each patch does and why each excluded file was excluded.
# ----------------------------------------------------------------------------
# Upstream fixes are ordinary commits on this branch (see patches/PATCHES.md for
# the audit trail); nothing is patched at configure time. One legacy program is
# excluded from the build (dead duplicate of FQHETorusFermionsWithSpin, DEFERRED.md).
set(DIAGHAM_UPSTREAM_EXCLUDED_PROGRAMS "QHEFermionsTorusWithSpin" CACHE INTERNAL "Programs excluded from the build")

# ----------------------------------------------------------------------------
# Build the libraries
# ----------------------------------------------------------------------------
add_subdirectory(Base/src)
add_subdirectory(src)

# Same order as the autotools top-level SUBDIRS:
#   Base src QuantumDots Spin FQHE Anyons FTI
if(DIAGHAM_BUILD_QUANTUMDOTS)
    add_subdirectory(QuantumDots/src)
endif()

if(DIAGHAM_BUILD_SPIN)
    add_subdirectory(Spin/src)
endif()

if(DIAGHAM_BUILD_FQHE)
    add_subdirectory(FQHE/src)
endif()

if(DIAGHAM_BUILD_FTI)
    add_subdirectory(FTI/src)
endif()

# ----------------------------------------------------------------------------
# Install rules, find_package(DiagHam) export, and the ctest suite
# ----------------------------------------------------------------------------
include(DiagHamInstall)

if(DIAGHAM_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ----------------------------------------------------------------------------
# Summary message
# ----------------------------------------------------------------------------
message(STATUS "")
message(STATUS "DiagHam CMake configuration:")
message(STATUS "  Kent defaults:  ${DIAGHAM_KENT_DEFAULTS}")
message(STATUS "  Build type:     ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ standard:   ${CMAKE_CXX_STANDARD}")
message(STATUS "  SMP (pthreads): ${DIAGHAM_USE_SMP}")
message(STATUS "  LAPACK:         ${DIAGHAM_USE_LAPACK}")
message(STATUS "  MPI:            ${DIAGHAM_USE_MPI}")
message(STATUS "  LAPACK only:    ${DIAGHAM_USE_LAPACK_ONLY}")
message(STATUS "  ScaLAPACK:      ${DIAGHAM_USE_SCALAPACK}")
message(STATUS "  GSL:            ${DIAGHAM_USE_GSL}")
message(STATUS "  GMP:            ${DIAGHAM_USE_GMP}")
message(STATUS "  MPACK:          ${DIAGHAM_USE_MPACK}")
message(STATUS "  FFTW:           ${DIAGHAM_USE_FFTW}")
message(STATUS "  bzip2:          ${DIAGHAM_USE_BZ2}")
message(STATUS "  Profile:        ${DIAGHAM_USE_PROFILE}")
message(STATUS "  Command log:    ${DIAGHAM_COMMAND_LOG}")
message(STATUS "  FQHE module:    ${DIAGHAM_BUILD_FQHE}")
message(STATUS "  FTI module:     ${DIAGHAM_BUILD_FTI}")
message(STATUS "  Spin module:    ${DIAGHAM_BUILD_SPIN}")
message(STATUS "  QuantumDots:    ${DIAGHAM_BUILD_QUANTUMDOTS}")
message(STATUS "  Build programs: ${DIAGHAM_BUILD_PROGRAMS}")
message(STATUS "  Tests (ctest):  ${DIAGHAM_BUILD_TESTS}")
message(STATUS "")
