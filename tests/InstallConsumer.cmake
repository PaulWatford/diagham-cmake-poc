# InstallConsumer.cmake: install DiagHam's Development component into a
# scratch prefix, then configure, build and run tests/consumer against it
# using nothing but find_package(DiagHam).
#
# Inputs: BUILD_DIR, CONSUMER_SOURCE_DIR, WORK_DIR, CXX_COMPILER, GENERATOR

file(REMOVE_RECURSE "${WORK_DIR}")
set(prefix "${WORK_DIR}/prefix")

function(run)
    execute_process(COMMAND ${ARGN} RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
        string(REPLACE ";" " " cmd "${ARGN}")
        message(FATAL_ERROR "command failed (${result}): ${cmd}")
    endif()
endfunction()

run(${CMAKE_COMMAND} --install "${BUILD_DIR}" --component Development --prefix "${prefix}")
run(${CMAKE_COMMAND} -S "${CONSUMER_SOURCE_DIR}" -B "${WORK_DIR}/build"
    -G "${GENERATOR}"
    "-DCMAKE_PREFIX_PATH=${prefix}"
    "-DCMAKE_CXX_COMPILER=${CXX_COMPILER}"
    -DCMAKE_BUILD_TYPE=Release)
run(${CMAKE_COMMAND} --build "${WORK_DIR}/build")
run("${WORK_DIR}/build/diagham_consumer")
