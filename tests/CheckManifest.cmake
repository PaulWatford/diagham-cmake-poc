# CheckManifest.cmake: the list of registered tests must equal tests/manifest.txt.
#
# A green ctest run that has silently lost tests is the failure mode this
# guards against (AGENTS.md "Locked artefacts"). Invoked by ctest as
#   cmake -DBUILD_DIR=... -DMANIFEST=... -P CheckManifest.cmake
# Tests whose registration depends on the environment rather than on the
# source (the python.* cross-checks need numpy) are compared only when
# present on both sides.
foreach(var BUILD_DIR MANIFEST)
    if(NOT DEFINED ${var})
        message(FATAL_ERROR "CheckManifest.cmake: ${var} not set")
    endif()
endforeach()

execute_process(COMMAND ${CMAKE_CTEST_COMMAND} -N
    WORKING_DIRECTORY "${BUILD_DIR}"
    OUTPUT_VARIABLE listing RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "ctest -N failed in ${BUILD_DIR}")
endif()
string(REGEX MATCHALL "Test +#[0-9]+: [^\n]+" lines "${listing}")
set(actual "")
foreach(line ${lines})
    string(REGEX REPLACE "Test +#[0-9]+: " "" name "${line}")
    string(STRIP "${name}" name)
    if(NOT name MATCHES "^python\\.")
        list(APPEND actual "${name}")
    endif()
endforeach()
list(SORT actual)

file(STRINGS "${MANIFEST}" expected)
list(FILTER expected EXCLUDE REGEX "^python\\.")
list(FILTER expected EXCLUDE REGEX "^$")
list(SORT expected)
# tests registered only when an optional library is on (LAPACK, GSL): allowed to be present or absent
set(optional "")
get_filename_component(manifest_dir "${MANIFEST}" DIRECTORY)
if(EXISTS "${manifest_dir}/manifest-optional.txt")
    file(STRINGS "${manifest_dir}/manifest-optional.txt" optional)
    list(FILTER optional EXCLUDE REGEX "^#")
    list(FILTER optional EXCLUDE REGEX "^$")
endif()

set(missing ${expected})
list(REMOVE_ITEM missing ${actual})
set(extra ${actual})
list(REMOVE_ITEM extra ${expected})
if(optional)
    list(REMOVE_ITEM missing ${optional})
    list(REMOVE_ITEM extra ${optional})
endif()
list(LENGTH missing n_missing)
list(LENGTH extra n_extra)
list(LENGTH actual n_actual)
if(n_missing GREATER 0 OR n_extra GREATER 0)
    message(FATAL_ERROR "test manifest mismatch: ${n_missing} missing (${missing}), ${n_extra} not in manifest (${extra}). "
        "If the change is intended, regenerate tests/manifest.txt: python3 tests/coverage.py ${BUILD_DIR} --manifest tests/manifest.txt")
endif()
message(STATUS "manifest OK: ${n_actual} tests")
