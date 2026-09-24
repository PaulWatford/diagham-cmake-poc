# RunAndCheck.cmake: run one DiagHam program, then check its spectrum file.
#
# Invoked by ctest as `cmake -D... -P RunAndCheck.cmake` (see
# tests/CMakeLists.txt). Every test gets a fresh, empty working directory so
# stale output from an earlier run can never satisfy a check.
#
# Inputs:
#   PROGRAM       absolute path of the DiagHam executable
#   PROGRAM_ARGS  ;-list of its arguments
#   WORK_DIR      scratch directory (wiped first)
#   OUTPUT_GLOB   glob, relative to WORK_DIR, matching exactly one output file
#   CHECKER       absolute path of check_spectrum
#   CHECK_ARGS    ;-list of check_spectrum arguments; the token @OUTPUT@ is
#                 replaced by the matched output file
#
# Optional, for comparing two programs against each other:
#   REF_PROGRAM, REF_PROGRAM_ARGS, REF_OUTPUT_GLOB
#                 a second program, run in WORK_DIR/reference; the token
#                 @REFERENCE@ in CHECK_ARGS is replaced by its output file

foreach(var PROGRAM WORK_DIR OUTPUT_GLOB CHECKER CHECK_ARGS)
    if(NOT DEFINED ${var})
        message(FATAL_ERROR "RunAndCheck.cmake: ${var} not set")
    endif()
endforeach()

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")

string(REPLACE "|" ";" CHECK_ARGS "${CHECK_ARGS}")

# Runs one program in its own directory; sets <out_var> to its one output file.
function(run_program program args_string dir glob out_var)
    string(REPLACE "|" ";" args "${args_string}")
    file(MAKE_DIRECTORY "${dir}")
    execute_process(
        COMMAND "${program}" ${args}
        WORKING_DIRECTORY "${dir}"
        RESULT_VARIABLE run_result
        OUTPUT_FILE "${dir}/program.log"
        ERROR_FILE "${dir}/program.err"
    )
    if(NOT run_result EQUAL 0)
        file(READ "${dir}/program.err" err)
        message(FATAL_ERROR "${program} exited with ${run_result}\n${err}\n(full log: ${dir}/program.log)")
    endif()
    file(GLOB outputs "${dir}/${glob}")
    list(LENGTH outputs n_outputs)
    if(NOT n_outputs EQUAL 1)
        file(GLOB all_files RELATIVE "${dir}" "${dir}/*")
        message(FATAL_ERROR "expected exactly one file matching '${glob}' in ${dir}, found ${n_outputs}; "
            "directory holds: ${all_files}")
    endif()
    set(${out_var} "${outputs}" PARENT_SCOPE)
endfunction()

run_program("${PROGRAM}" "${PROGRAM_ARGS}" "${WORK_DIR}" "${OUTPUT_GLOB}" output)
string(REPLACE "@OUTPUT@" "${output}" CHECK_ARGS "${CHECK_ARGS}")

if(DEFINED REF_PROGRAM)
    run_program("${REF_PROGRAM}" "${REF_PROGRAM_ARGS}" "${WORK_DIR}/reference" "${REF_OUTPUT_GLOB}" reference)
    string(REPLACE "@REFERENCE@" "${reference}" CHECK_ARGS "${CHECK_ARGS}")
endif()
execute_process(
    COMMAND "${CHECKER}" ${CHECK_ARGS}
    RESULT_VARIABLE check_result
)
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "spectrum check failed (${CHECKER} ${CHECK_ARGS})")
endif()
