# RunAndMatch.cmake: run one DiagHam program and check what it prints.
#
# Invoked by ctest as `cmake -D... -P RunAndMatch.cmake` (see
# diagham_dimension_test in tests/CMakeLists.txt). Fresh working directory
# every time.
#
# Inputs:
#   PROGRAM       absolute path of the executable
#   PROGRAM_ARGS  |-separated arguments
#   WORK_DIR      scratch directory (wiped first)
#   EXPECT_FILE   expected-output file (tests/data/dimensions/*.txt):
#                   lines starting with '#' are comments;
#                   a line "COUNT <n>" means: the number of stdout lines
#                   matching COUNT_REGEX must be n;
#                   every other line must occur verbatim as a whole line
#                   of stdout (trailing whitespace ignored)
#   COUNT_REGEX   regex a printed basis state matches (only with COUNT)
foreach(var PROGRAM WORK_DIR EXPECT_FILE)
    if(NOT DEFINED ${var})
        message(FATAL_ERROR "RunAndMatch.cmake: ${var} not set")
    endif()
endforeach()
file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
string(REPLACE "|" ";" args "${PROGRAM_ARGS}")
execute_process(
    COMMAND "${PROGRAM}" ${args}
    WORKING_DIRECTORY "${WORK_DIR}"
    RESULT_VARIABLE rc
    OUTPUT_FILE "${WORK_DIR}/program.log"
    ERROR_FILE "${WORK_DIR}/program.err")
if(NOT rc EQUAL 0)
    file(READ "${WORK_DIR}/program.err" err)
    message(FATAL_ERROR "${PROGRAM} exited with ${rc}\n${err}")
endif()
file(STRINGS "${WORK_DIR}/program.log" printed)
set(printed_trimmed "")
foreach(line ${printed})
    string(REGEX REPLACE "[ \t]+$" "" line "${line}")
    list(APPEND printed_trimmed "${line}")
endforeach()

file(STRINGS "${EXPECT_FILE}" expected)
set(failures "")
foreach(line ${expected})
    if(line MATCHES "^#" OR line STREQUAL "")
        continue()
    endif()
    if(line MATCHES "^COUNT ([0-9]+)$")
        set(want ${CMAKE_MATCH_1})
        if(NOT DEFINED COUNT_REGEX)
            message(FATAL_ERROR "COUNT line needs COUNT_REGEX")
        endif()
        set(got 0)
        foreach(p ${printed_trimmed})
            if(p MATCHES "${COUNT_REGEX}")
                math(EXPR got "${got} + 1")
            endif()
        endforeach()
        if(NOT got EQUAL want)
            list(APPEND failures "expected ${want} basis states (lines matching '${COUNT_REGEX}'), the program printed ${got}")
        endif()
    else()
        string(REGEX REPLACE "[ \t]+$" "" line "${line}")
        list(FIND printed_trimmed "${line}" idx)
        if(idx EQUAL -1)
            list(APPEND failures "line not printed: '${line}'")
        endif()
    endif()
endforeach()
if(failures)
    string(REPLACE ";" "\n  " failures "${failures}")
    message(FATAL_ERROR "output mismatch (${WORK_DIR}/program.log):\n  ${failures}")
endif()
