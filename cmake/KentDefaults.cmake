# ============================================================================
#  cmake/KentDefaults.cmake
# ============================================================================
#
#  Site-defaults file for building DiagHam on the University of Kent's own
#  machines/cluster, modelled directly on Gunnar Moeller & Lars Schonenberg's
#  cmake/TCMDefaults.cmake from BDMC_UFL (their equivalent file for the TCM
#  group's cluster at Cambridge). Toggled the same way TCMDefaults.cmake is:
#  a cache flag (there TCM_DEFAULTS, here DIAGHAM_KENT_DEFAULTS), off by
#  default, included conditionally.
#
#  IMPORTANT -- this is a SCAFFOLD, not a finished site-defaults file.
#  Every path below is a placeholder. Unlike TCM_DEFAULTS (which has
#  Gunnar's own real, working cluster paths), I have no actual knowledge
#  of Kent's HPC layout, module system, or where MKL/NAG/HDF5 live on
#  Kent's machines -- inventing plausible-looking paths here would be
#  actively worse than leaving this unfinished, since a wrong hardcoded
#  path fails silently at link time, not at configure time. The FATAL_ERROR
#  guard below is deliberate: it means this file can be merged now (the
#  toggle mechanism, matching Gunnar's own pattern, is genuinely useful
#  to have in place) without anyone being able to silently build against
#  placeholder paths -- it has to be filled in first, loudly.
#
#  To finish this: ask Gunnar (he'll know Kent's setup, or who to ask --
#  likely Kent's e-Research/RCS team) for the real values, matching each
#  TODO below to what TCMDefaults.cmake sets for Cambridge's TCM cluster.
# ============================================================================

if(NOT DIAGHAM_KENT_MKL_ROOT)
    message(FATAL_ERROR
        "DIAGHAM_KENT_DEFAULTS is ON but this file hasn't been filled in yet. "
        "cmake/KentDefaults.cmake is a scaffold (see its header comment) -- "
        "ask Gunnar for Kent's actual MKL/compiler/MPI/HDF5 paths before "
        "enabling this option, or configure without DIAGHAM_KENT_DEFAULTS "
        "and pass paths manually via the usual find_package hints.")
endif()

# TODO (ask Gunnar / Kent e-Research): MKL installation root on Kent's
# cluster, equivalent to TCMDefaults.cmake's MKL_ROOT ("/misc/shared/mkl"
# for TCM).
set(DIAGHAM_KENT_MKL_ROOT "" CACHE PATH "MKL root directory on Kent's cluster")
if(DIAGHAM_KENT_MKL_ROOT)
    set(MKL_ROOT ${DIAGHAM_KENT_MKL_ROOT} CACHE PATH "MKL root directory" FORCE)
endif()

# TODO: compiler root, equivalent to TCMDefaults.cmake's INTEL_ROOT, if
# Kent's cluster uses Intel compilers/MKL via a similar module path scheme.
set(DIAGHAM_KENT_COMPILER_ROOT "" CACHE PATH "Compiler installation root on Kent's cluster")

# TODO: MPI compiler wrapper path, equivalent to TCMDefaults.cmake's
# MPI_CXX_COMPILER ("/usr/local/shared/MPI/OpenMPI-1.10.1/intel64/bin/mpiCC"
# for TCM). Only relevant once DIAGHAM_USE_MPI is exercised for real.
set(DIAGHAM_KENT_MPI_CXX_COMPILER "" CACHE PATH "Path to Kent's MPI C++ compiler wrapper")
if(DIAGHAM_KENT_MPI_CXX_COMPILER)
    set(MPI_CXX_COMPILER ${DIAGHAM_KENT_MPI_CXX_COMPILER} CACHE PATH "Path to the MPI c++ compiler." FORCE)
endif()

# HDF5 static linking preference -- TCMDefaults.cmake sets this to prefer
# static (1). Kept as a real default here since it's a policy choice, not
# a hardcoded path, and static linking is a reasonable default on a
# shared cluster filesystem either way.
set(HDF5_USE_STATIC_LIBRARIES 1 CACHE BOOL "Prefer static linking of hdf5 libraries")

message(STATUS "Kent site defaults loaded (DIAGHAM_KENT_DEFAULTS=ON)")
