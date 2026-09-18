# BDMC - Bold Diagrammatic Monte Carlo Project
# Copyright (C) 2013-2014 Gunnar Möller & Lars Schonenberg
#
# Find the Intel NAG libraries
#
# Author: Gunnar Möller
# Adapted from: https://github.com/hanjianwei/cmake-modules/blob/master/FindNAG.cmake
# Date added: 10/29/2014


# This script locates the NAG C library.
#
#
# The NAG libaries supports linking in various configurations,
# these can be controlled by setting the following variables prior to
# calling this script:
#
#   NAG_STATIC        :   use static linking
#   NAG_VECTORIZED    :   use veclib framework or MKL
#   NAG_INT64         :   NAG library uses 64 bit integers
#   NAG_NEEDS_IFCORE  :   also link with ifcore
#
# This module defines the following variables for external use:
#
#   NAG_FOUND            : True if NAG include dirs and libraries are detected
#   NAG_INCLUDE_DIRS     : set when NAG_INCLUDE_DIR found
#   NAG_LIBRARIES        : the library to link against.


include(FindPackageHandleStandardArgs)

MESSAGE(STATUS "NAG option NAG_STATIC: ${NAG_STATIC}" )
MESSAGE(STATUS "NAG option NAG_VECTORIZED: ${NAG_VECTORIZED}" )


######################### NAG Root #######################

# Determine the location of the NAG root directory
# Path specified by environment variable
if(DEFINED ENV{NAGROOT})
    set(NAG_ROOT $ENV{NAGROOT} CACHE PATH "NAG root directory")
# Default installation directory on Darwin
elseif(IS_DIRECTORY "/opt/NAG/clmi623dgl/include/")
    set(NAG_ROOT "/opt/NAG/clmi623dgl/" CACHE PATH "NAG root directory")
else()
    # Try default locations for Cambridge TCM group
    # with intel compiler
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        if(IS_DIRECTORY "/misc/shared/nag/C/Mark24_icc64/include/")
            set(NAG_ROOT "/misc/shared/nag/C/Mark24_icc64/" CACHE PATH "NAG root directory")
        # NAG root not found, display warning
        else()
            message("Failed to locate NAG root directory, please specify manually if necessary, using environment variable NAGROOT")
            set(NAG_ROOT "" CACHE PATH "NAG root directory")
        endif()
    elseif(CMAKE_COMPILER_IS_GNUCXX) # and alternatively, for gcc
        if(IS_DIRECTORY "/misc/shared/nag/C/Mark23_gcc64/include/")
            set(NAG_ROOT "/misc/shared/nag/C/Mark23_gcc64/" CACHE PATH "NAG root directory")
        # NAG root not found, display warning
        else()
            message("Failed to locate NAG root directory, please specify manually if necessary, using environment variable NAGROOT")
            set(NAG_ROOT "" CACHE PATH "NAG root directory")
        endif()
    endif()
endif()


######################### Include directories #######################

# Find NAG include dir
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    find_path(NAG_INCLUDE_DIR nag.h
        PATHS ${NAG_ROOT}/include/ ENV INCLUDE)
    set(NAG_INCLUDE_DIR ${NAG_INCLUDE_DIR})
else()
    find_path(NAG_INCLUDE_DIR nag.h
        PATHS ${NAG_ROOT}/include/ ENV INCLUDE)
    set(NAG_INCLUDE_DIR ${NAG_INCLUDE_DIR})
endif()


######################### NAG Libraries #######################

# Handle suffix
set(_NAG_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    SET(CMAKE_FIND_LIBRARY_PREFIXES "")
    if(NAG_STATIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES _dll.lib)
    endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    SET(CMAKE_FIND_LIBRARY_PREFIXES "lib")
    if(NAG_STATIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .dylib)
    endif()
else()
    SET(CMAKE_FIND_LIBRARY_PREFIXES "lib")
    if(NAG_STATIC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .so.*)
    endif()
endif()


# NAG can use native nag routines or rely on an external vector framework
# initialize optional veclib variabla
set(VECLIB_LIBRARY "")

# reset any previous NAG_LIBRARY value, overriding previous settings in the cache
unset(NAG_LIBRARY CACHE)

if(NAG_VECTORIZED)
    # check if we have the MKL library
    if(MKL_FOUND)
       # try to find nagc_mkl, then.
       find_library(NAG_LIBRARY nagc_mkl
            PATHS ${NAG_ROOT}/lib/ ${NAG_ROOT}/lib/shared/ ENV LIBRARY_PATH)
       if (NOT NAG_LIBRARY)
           # did not find nagc_mkl library
           # instead, try nagc_vl and also link the veclib framework
           find_library(NAG_LIBRARY nagc_vl
               PATHS ${NAG_ROOT}/lib/ ${NAG_ROOT}/lib/shared/ ENV LIBRARY_PATH)
           find_library(VECLIB_LIBRARY veclib ENV LIBRARY_PATH)
           if (NOT NAG_LIBRARY)
               # did not find MKL library, still
               MESSAGE(WARNING "Could not find any vectorized NAG libraries" )
           else()
               MESSAGE(STATUS "Using NAG Library with veclib support: ${NAG_LIBRARY}" )
           endif()
       else()
         MESSAGE(STATUS "Using NAG Library with MKL support: ${NAG_LIBRARY}" )
       endif()
    else()
       # otherwise, use nagc_vl and also link the veclib framework
       find_library(NAG_LIBRARY nagc_vl
           PATHS ${NAG_ROOT}/lib/ ${NAG_ROOT}/lib/shared/ ENV LIBRARY_PATH)
       find_library(VECLIB_LIBRARY veclib ENV LIBRARY_PATH)
       if (NOT NAG_LIBRARY)
           # did not find nagc_vl
           MESSAGE(WARNING "Could not find vectorized NAG library nagc_vl" )
       else()
           MESSAGE(STATUS "Using NAG Library with veclib support: ${NAG_LIBRARY}" )
       endif()
    endif()
else()
    # try to find nagc_nag for the native implementation
    find_library(NAG_LIBRARY nagc_nag
        PATHS ${NAG_ROOT}/lib/ ${NAG_ROOT}/lib/shared/ ENV LIBRARY_PATH)
    if (NOT NAG_LIBRARY)
        # did not find nagc_nag
        MESSAGE(WARNING "Could not find NAG libraries" )
    else()
        MESSAGE(STATUS "Using NAG Library with native vector routines: ${NAG_LIBRARY}" )
    endif()
endif()

#restore original suffixes, the search for standard libraries
set(CMAKE_FIND_LIBRARY_SUFFIXES ${_NAG_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

# search for pthread library
find_library(PTHREAD_LIBRARY pthread)

if (NOT PTHREAD_LIBRARY)
    # did not find pthread
    MESSAGE(WARNING "Could not find pthread libraries" )
endif()

if (NAG_NEEDS_IFCORE)
    find_library(IFCORE_LIBRARY ifcore
        ${INTEL_ROOT}/lib/ ${INTEL_ROOT}/lib/intel64/ ENV LIBRARY_PATH)
    if (NOT IFCORE_LIBRARY)
        # did not find ifcore
        MESSAGE(WARNING "Could not find ifcore libraries" )
    endif()
endif()


#check if package is detected
find_package_handle_standard_args(NAG DEFAULT_MSG
    NAG_LIBRARY NAG_INCLUDE_DIR)



# Set NAG_LIBRARY variable to include the required veclib and pthread library, as well
set(NAG_LIBRARY ${NAG_LIBRARY} ${VECLIB_LIBRARY} ${PTHREAD_LIBRARY} ${IFCORE_LIBRARY})


MESSAGE(STATUS "Final linker options for NAG library: ${NAG_LIBRARY}" )

######################### Finalize #######################

#if package is detected, set return variables
if(NAG_FOUND)
    set(NAG_INCLUDE_DIRS ${NAG_INCLUDE_DIR})
    set(NAG_LIBRARIES ${NAG_LIBRARY})
    set(HAVE_LIBNAG "yes")
    # NAG wants a compiler flag "NAG_INT" to allow for use of 32-bit integers (64 bit appears to be prone to bugs)
    # could automake choice of 32bit vs 64bit by probing NAG variable NAG_IFMT, which is either "d" (32 bit) or "ld" (64 bit)
    if (NAG_INT64)
       MESSAGE(STATUS "Assuming NAG Library uses 64 bit integers.")
    else()
       # 32 bit integers require a compiler flag
       set(NAG_INT "yes")
    endif()
    # BDMC code uses compiler flag "HAVE_LIBNAG"
    # add_definitions(-DHAVE_LIBNAG) # - now provided in nag_headers.h
    MESSAGE(STATUS "Using NAG Library under root: ${NAG_ROOT}" )
endif()
