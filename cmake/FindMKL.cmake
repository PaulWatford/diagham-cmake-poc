# BDMC - Bold Diagrammatic Monte Carlo Project
# Copyright (C) 2013-2014 Gunnar Möller & Lars Schonenberg
#
# Find the Intel MKL libraries
#
# Author: Lars Schonenberg
# Adapted from: https://github.com/hanjianwei/cmake-modules/blob/master/FindMKL.cmake
# Date added: 09/04/2014


# This script locates the intel MKL library.
#
# Currently, only the 64 bit version of the library using the 32 BIT INTERFACE (lp64) is supported.
# The script has been tested on Linux and Mac OS X, but not for Windows
#
# The Intel MKL libaries supports linking in various configurations,
# these can be controlled by setting the following variables prior to
# calling this script:
#
#   MKL_STATIC        :   use static linking
#   MKL_MULTI_THREADED:   use multi-threading
#   MKL_SDL           :   Single Dynamic Library interface
#                         (postpone configurations options till runtime, NOT RECOMMENDED)
#
# In additionan, the MKL libaries include additional support for clusters,
# the cluster libraries can be including by requiring the following components:
#
#   CLUSTER              : Search for cluster components of MKL library.
#   FFTW3                : Search for the FFTW3 interface to the MKL FFT library.
#
# This module defines the following variables for external use:
#
#   MKL_FOUND            : True if MKL include dirs and libraries are detected
#   MKL_INCLUDE_DIRS     : set when MKL_INCLUDE_DIR found
#   MKL_LIBRARIES        : the library to link against.
#   MKL_CLUSTER_FOUND    : True if the optional cluster components are also detected
#   MKL_CLUSTER_LIBRARIES: The cluster libraries to link against



include(FindPackageHandleStandardArgs)

# Set dummy vars for additional components to search for
foreach(component ${MKL_FIND_COMPONENTS})
    set("MKL_${component}_SEARCH" TRUE)
endforeach()


######################### MKL Root #######################

# Determine the location of the MKL root directory
if(IS_DIRECTORY "${MKL_ROOT}/include")
  # Path already specified, no further work to be done.
elseif(DEFINED ENV{MKLROOT})
  # Use path specified by environment variable
  set(MKL_ROOT $ENV{MKLROOT} CACHE PATH "MKL root directory")
elseif(IS_DIRECTORY "/opt/intel/mkl/include")
  # Use default installation directory
  set(MKL_ROOT "/opt/intel/mkl" CACHE PATH "MKL root directory")
else()
  # MKL root not found, display warning
  message("Failed to locate MKL root directory, please specify manually if necesarry")
  set(MKL_ROOT "" CACHE PATH "MKL root directory")
endif()


######################### Intel root #######################

# Determine the location of the MKL root directory
if(IS_DIRECTORY "${INTEL_ROOT}/lib")
  # Path already specified, no further work to be done.
elseif(IS_DIRECTORY "/opt/intel/lib")
  # Use default installation directory
  set(INTEL_ROOT "/opt/intel" CACHE PATH "Intel root directory")
else()
  # Intel root not found, display warning
    message("Failed to locate Intel root directory, please specify manually if necesarry")
    set(INTEL_ROOT "" CACHE PATH "Intel root directory")
endif()


######################### Include directories #######################

# Find MKL include dir
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    # Windows based systems need both Intel and MKL include dirs
    find_path(INTEL_INCLUDE_DIR omp.h
        PATHS ${INTEL_ROOT}/include/ ENV INCLUDE)
    find_path(MKL_INCLUDE_DIR mkl.h
        PATHS ${MKL_ROOT}/include/ ENV INCLUDE)
    set(MKL_INCLUDE_DIR ${MKL_INCLUDE_DIR} ${INTEL_INCLUDE_DIR})
else()
    # Unix based system only need MKL include dir
    find_path(MKL_INCLUDE_DIR mkl.h
        PATHS ${MKL_ROOT}/include/ ENV INCLUDE)
endif()


######################### MKL Libraries #######################

# Handle suffix
set(_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    if(MKL_STATIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES _dll.lib)
    endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    if(MKL_STATIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .dylib)
    endif()
else()
    if(MKL_STATIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
    endif()
endif()


# MKL is composed by four layers: Interface, Threading, Computational and RTL
if(MKL_SDL)
    find_library(MKL_LIBRARY mkl_rt
        PATHS ${MKL_ROOT}/lib/ ${MKL_ROOT}/lib/intel64/ ENV LIBRARY_PATH)
else()
    ######################### Interface layer #######################
    if(WIN32)
        set(MKL_INTERFACE_LIBNAME mkl_intel_lp64_c)
    else()
        set(MKL_INTERFACE_LIBNAME mkl_intel_lp64)
    endif()

    find_library(MKL_INTERFACE_LIBRARY ${MKL_INTERFACE_LIBNAME}
        PATHS ${MKL_ROOT}/lib/ ${MKL_ROOT}/lib/intel64/ ENV LIBRARY_PATH)

    ####################### Computational layer #####################
    find_library(MKL_CORE_LIBRARY mkl_core
        PATHS ${MKL_ROOT}/lib/ ${MKL_ROOT}/lib/intel64/ ENV LIBRARY_PATH)

    ###################### Threading & RTL layers ###################
    if(MKL_MULTI_THREADED)
        # Use multi-threaded MKL thread library
        set(MKL_THREADING_LIBNAME mkl_intel_thread)
        find_library(MKL_THREADING_LIBRARY ${MKL_THREADING_LIBNAME}
            PATHS ${MKL_ROOT}/lib/ ${MKL_ROOT}/lib/intel64/ ENV LIBRARY_PATH)

        # Also locate the RTL library (i.e. OpenMP interface)
        if(WIN32)
            set(MKL_RTL_LIBNAME iomp5md)
        else()
            set(MKL_RTL_LIBNAME iomp5)
        endif()
        find_library(MKL_RTL_LIBRARY ${MKL_RTL_LIBNAME}
            PATHS ${INTEL_ROOT}/lib/ ${INTEL_ROOT}/lib/intel64/ ENV LIBRARY_PATH)

        # Set MKL_LIBRARY variable to include all necesarry libraries
        set(MKL_LIBRARY ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_RTL_LIBRARY})
    else()
        # Use single threaded MKL thread library
        set(MKL_THREADING_LIBNAME mkl_sequential)
        find_library(MKL_THREADING_LIBRARY ${MKL_THREADING_LIBNAME}
            PATHS ${MKL_ROOT}/lib/ ${MKL_ROOT}/lib/intel64/ ENV LIBRARY_PATH)

        # No further work to be done

        # Set MKL_LIBRARY variable to include all necesarry libraries
        set(MKL_LIBRARY ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY})
    endif()
endif()


######################### Cluster components #######################

if(MKL_CLUSTER_SEARCH)
    if(MKL_SDL)
        message(FATAL_ERROR "Not possible to use cluster components with Intel MKL SDL interface")
    else()
        find_library(MKL_BLACS_LIBRARY mkl_blacs_intelmpi_lp64
            PATHS ${MKL_ROOT}/lib/ ${MKL_ROOT}/lib/intel64/ ENV LIBRARY_PATH)
        find_library(MKL_FFT_LIBRARY mkl_cdft_core
            PATHS ${MKL_ROOT}/lib/ ${MKL_ROOT}/lib/intel64/ ENV LIBRARY_PATH)
        find_library(MKL_SCALAPACK_LIBRARY mkl_scalapack_core
            PATHS ${MKL_ROOT}/lib/ ${MKL_ROOT}/lib/intel64/ ENV LIBRARY_PATH)
        set(MKL_CLUSTER_LIBRARY ${MKL_BLACS_LIBRARY} ${MKL_FFT_LIBRARY} ${MKL_SCALAPACK_LIBRARY})
        find_package_handle_standard_args(MKL_CLUSTER
            REQUIRED_VARS MKL_CLUSTER_LIBRARY)
    endif()
endif()


######################### FFTW3 Interface #######################

if(MKL_FFTW3_SEARCH)
    find_path(MKL_FFTW3_INCLUDE_DIR fftw3.h
        PATHS ${MKL_ROOT}/include ENV INCLUDE
        PATH_SUFFIXES fftw)
    # Manually set a variable to track if FFTW3 component was found
    if (MKL_FFTW3_INCLUDE_DIR)
        set(MKL_FFTW3_FOUND TRUE)
    else()
        set(MKL_FFTW3_FOUND FALSE)
    endif()
endif()


######################### Finalize #######################

#restore original suffixes
set(CMAKE_FIND_LIBRARY_SUFFIXES ${_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

#check if package is detected
find_package_handle_standard_args(MKL
    REQUIRED_VARS MKL_INCLUDE_DIR MKL_LIBRARY
    HANDLE_COMPONENTS)

#if package is detected, set return variables
if(MKL_FOUND)
    set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
    set(MKL_LIBRARIES ${MKL_LIBRARY})

    if(MKL_CLUSTER_FOUND)
        set(MKL_CLUSTER_LIBRARIES ${MKL_CLUSTER_LIBRARY})
    endif()

    MESSAGE(STATUS "Using MKL Library under root: ${MKL_ROOT}" )
    if(MKL_FFTW3_FOUND)
        set(MKL_FFTW3_INCLUDE_DIRS ${MKL_FFTW3_INCLUDE_DIR})
        MESSAGE(STATUS "FFTW3 component of MKL Library with header files at: ${MKL_FFTW3_INCLUDE_DIR}" )

    endif()

endif()
