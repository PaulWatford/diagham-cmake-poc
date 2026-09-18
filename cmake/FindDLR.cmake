# Try to find the DLR libraries
# DLR_FOUND - system has DLR lib
# DLR_INCLUDE_DIR - the DLR include directory
# DLR_LIBRARIES_DIR - Directory where the DLR libraries are located
# DLR_LIBRARIES - the DLR libraries

include(FindPackageHandleStandardArgs)

if(DLR_INCLUDE_DIR)
  set(DLR_in_cache TRUE)
else()
  set(DLR_in_cache FALSE)
endif()
if(NOT DLR_LIBRARIES)
  set(DLR_in_cache FALSE)
endif()

# Is it already configured?
if (NOT DLR_in_cache)

  find_path(DLR_INCLUDE_DIR 
	   NAMES dlr_c.h
           HINTS /usr/local/include/dlr_c include/dlr_c
		 END DLR_DIR
                DOC "The directory containing the DLR header files"
           )

  find_library(DLR_LIBRARIES NAMES libdlr.dylib libdlr_c.dylib libdlr.so libdlr_c.so libdlr.a libdlr_c.a
              HINTS /usr/local/lib lib
    DOC "Path to the DLR library"
    )

  if ( DLR_LIBRARIES )
    get_filename_component(DLR_LIBRARIES_DIR ${DLR_LIBRARIES} PATH CACHE )
  endif()

  # Attempt to load a user-defined configuration for GMP if couldn't be found
  if ( NOT DLR_INCLUDE_DIR OR NOT DLR_LIBRARIES_DIR )
    include( DLRConfig OPTIONAL )
  endif()

endif()

find_package_handle_standard_args(DLR "DEFAULT_MSG" DLR_LIBRARIES DLR_INCLUDE_DIR)
