# Try to find the MPC libraries
# MPC_FOUND - system has MPC lib
# MPC_INCLUDE_DIR - the MPC include directory
# MPC_LIBRARIES_DIR - Directory where the MPC libraries are located
# MPC_LIBRARIES - the MPC libraries

include(FindPackageHandleStandardArgs)

if(MPC_INCLUDE_DIR)
  set(MPC_in_cache TRUE)
else()
  set(MPC_in_cache FALSE)
endif()
if(NOT MPC_LIBRARIES)
  set(MPC_in_cache FALSE)
endif()

# Is it already configured?
if (NOT MPC_in_cache)

  find_path(MPC_INCLUDE_DIR
            NAMES mpc.h
            HINTS ENV MPC_INC_DIR
                  ENV MPC_DIR
            PATH_SUFFIXES include
                DOC "The directory containing the MPC header files"
           )

  find_library(MPC_LIBRARIES NAMES mpc libmpc-3
    HINTS ENV MPC_LIB_DIR
          ENV MPC_DIR
    PATH_SUFFIXES lib
    DOC "Path to the MPC library"
    )

  if ( MPC_LIBRARIES )
    get_filename_component(MPC_LIBRARIES_DIR ${MPC_LIBRARIES} PATH CACHE )
  endif()

  # Attempt to load a user-defined configuration for MPC if couldn't be found
  if ( NOT MPC_INCLUDE_DIR OR NOT MPC_LIBRARIES_DIR )
    include( MPCConfig OPTIONAL )
  endif()

endif()

find_package_handle_standard_args(MPC "DEFAULT_MSG" MPC_LIBRARIES MPC_INCLUDE_DIR)
