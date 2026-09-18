# Try to find the FFTW3 libraries
# FFTW3_FOUND - system has FFTW3 lib
# FFTW3_INCLUDE_DIR - the FFTW3 include directory
# FFTW3_LIBRARIES_DIR - Directory where the FFTW3 libraries are located
# FFTW3_LIBRARIES - the FFTW3 libraries

include(FindPackageHandleStandardArgs)

if(FFTW3_INCLUDE_DIR)
  set(FFTW3_in_cache TRUE)
else()
  set(FFTW3_in_cache FALSE)
endif()
if(NOT FFTW3_LIBRARIES)
  set(FFTW3_in_cache FALSE)
endif()

# Is it already configured?
if (NOT FFTW3_in_cache)

  find_path(FFTW3_INCLUDE_DIR
            NAMES fftw3.h
            HINTS ENV FFTW3_INC_DIR
                  ENV FFTW3_DIR
            PATH_SUFFIXES include
                DOC "The directory containing the FFTW3 header files"
           )

  find_library(FFTW3_LIBRARIES NAMES fftw3
    HINTS ENV FFTW3_LIB_DIR
          ENV FFTW3_DIR
    PATH_SUFFIXES lib
    DOC "Path to the FFTW3 library"
    )

  if ( FFTW3_LIBRARIES )
    get_filename_component(FFTW3_LIBRARIES_DIR ${FFTW3_LIBRARIES} PATH CACHE )
  endif()

  # Attempt to load a user-defined configuration for FFTW3 if couldn't be found
  if ( NOT FFTW3_INCLUDE_DIR OR NOT FFTW3_LIBRARIES_DIR )
    include( FFTW3Config OPTIONAL )
  endif()

endif()

find_package_handle_standard_args(FFTW3 "DEFAULT_MSG" FFTW3_LIBRARIES FFTW3_INCLUDE_DIR)
