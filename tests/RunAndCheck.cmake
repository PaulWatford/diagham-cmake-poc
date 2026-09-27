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
#
# Optional, for chains (eigenstate -> Jack polynomial -> overlap):
#   PRE_STEPS     programs run in WORK_DIR before PROGRAM, in order, each
#                 as "executable|arg|arg", steps separated by "^^"; every
#                 step must exit 0. Their stdout goes to WORK_DIR/step<i>.log.

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
    # OUTPUT_GLOB may name program.log itself: programs that print their result
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

if(DEFINED PRE_STEPS)
    string(REPLACE "^^" ";" steps "${PRE_STEPS}")
    set(i 0)
    foreach(step ${steps})
        math(EXPR i "${i} + 1")
        string(REPLACE "|" ";" step_cmd "${step}")
        execute_process(
            COMMAND ${step_cmd}
            WORKING_DIRECTORY "${WORK_DIR}"
            RESULT_VARIABLE step_result
            OUTPUT_FILE "${WORK_DIR}/step${i}.log"
            ERROR_FILE "${WORK_DIR}/step${i}.err")
        if(NOT step_result EQUAL 0)
            file(READ "${WORK_DIR}/step${i}.err" err)
            message(FATAL_ERROR "preparatory step ${i} (${step_cmd}) exited with ${step_result}\n${err}")
        endif()
    endforeach()
endif()

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
