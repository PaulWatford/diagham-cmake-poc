# ============================================================================
# ApplyUpstreamPatches.cmake: apply upstream-bug patches at configure time
# ============================================================================
#
# CMake migration of DiagHam surfaced ten files in the FQHE and FTI modules
# that do not compile against the canonical DiagHam codebase. The autotools
# build hides these failures: most of the broken files are unreachable from
# the default build, and the few that are reachable produce errors that get
# lost in the noise.
#
# This module applies a numbered patch series under ${CMAKE_SOURCE_DIR}/patches
# at configure time. A sentinel file in the build tree records the patch
# fingerprint so re-runs skip work that has already been done.
#
# What each patch does is documented in patches/PATCHES.md, with explicit
# audit trails (sibling-file comparisons, header inspections, and physics
# notes) that justify every change.
#
# Files NOT patched (excluded from the build instead, see DiagHamHelpers.cmake)
# are listed in DIAGHAM_UPSTREAM_EXCLUDED_PROGRAMS and documented in PATCHES.md.
#
# ----------------------------------------------------------------------------
# Apply strategy (revised 22/09, see patches/PATCHES.md "audited by a second
# Claude session" entry for the full story -- the first version of this fix
# was wrong):
#
#   1. `patch --binary -p1 --forward` on the patch file as committed.
#   2. If that fails, byte-convert a copy of the SAME patch file to CRLF and
#      retry with `patch --binary -p1 --forward` again.
#
# Why: every patch in this series is stored in git as pure LF (confirmed by
# inspecting the git blobs directly, independent of whatever line endings a
# particular checkout happens to have on disk). Ten of the eleven patches
# target upstream files that are themselves LF, so step 1 matches them
# byte-for-byte. Patch 09 targets exactly one upstream file
# (FTI/src/Programs/FCI/FCIHofstadterCorrelation.cc) that is CRLF in the
# canonical DiagHam tree itself -- an LF patch can never match CRLF context,
# under any tool, no matter how it's invoked. Step 2 exists for exactly that
# file: convert the patch's own line endings to match its target, and use
# --binary so GNU patch doesn't strip the CRs back off before comparing
# (that stripping is documented, deliberate GNU patch behaviour, not a bug
# in any particular patch version -- an earlier version of this comment
# wrongly blamed a "GNU patch 2.7.6 CRLF-detection bug").
#
# A `git apply` fallback was tried first and looked like it worked, but only
# because it was verified against a Windows checkout's CRLF working-tree
# copies of the patches, not the actual LF content git has stored. Against
# the real git-tracked patches, `git apply` fails on patch 09 for the same
# reason `patch` alone does (LF patch, CRLF target) -- it was never a fix,
# it was untested against the actual repository content. This version is
# verified against fresh clones of both the upstream mirror AND this
# repository's real git-tracked patch files, not a local working copy.
# ----------------------------------------------------------------------------

set(DIAGHAM_PATCH_DIR ${CMAKE_SOURCE_DIR}/patches)
# The sentinel lives in the SOURCE tree, because that is what the patches
# modify: every build directory configured from one tree (one per preset,
# say) shares the patched sources. With the sentinel in the build directory,
# a second build directory re-applied the series to already-patched
# sources and failed at configure time.
set(DIAGHAM_PATCH_SENTINEL ${CMAKE_SOURCE_DIR}/.diagham_upstream_patches_applied)

if(NOT EXISTS ${DIAGHAM_PATCH_DIR})
    message(STATUS "DiagHam: no upstream patch directory found, skipping")
    set(DIAGHAM_UPSTREAM_EXCLUDED_PROGRAMS "" CACHE INTERNAL "")
    return()
endif()

file(GLOB patch_files ${DIAGHAM_PATCH_DIR}/[0-9]*.patch)
list(SORT patch_files)

if(NOT patch_files)
    message(STATUS "DiagHam: no upstream patches to apply")
else()
    # Compute a content fingerprint of the patches so we re-apply on any change
    set(patch_fingerprint "")
    foreach(p ${patch_files})
        file(MD5 ${p} hash)
        string(APPEND patch_fingerprint "${hash}\n")
    endforeach()

    set(need_apply TRUE)
    if(EXISTS ${DIAGHAM_PATCH_SENTINEL})
        file(READ ${DIAGHAM_PATCH_SENTINEL} existing_fingerprint)
        if("${existing_fingerprint}" STREQUAL "${patch_fingerprint}")
            set(need_apply FALSE)
            message(STATUS "DiagHam: upstream patches already applied (sentinel matches)")
        else()
            message(FATAL_ERROR
                "DiagHam: this source tree was patched with a different patch series "
                "(${DIAGHAM_PATCH_SENTINEL} does not match patches/). Patches cannot be "
                "re-applied on top of each other: start again from a pristine upstream "
                "tree (e.g. `git checkout -- . && git clean -fd` in an upstream clone, "
                "then scripts_cmake/overlay.py).")
        endif()
    endif()

    if(need_apply)
        find_program(PATCH_EXECUTABLE patch REQUIRED)
        find_package(Python3 COMPONENTS Interpreter QUIET)

        # Runs `patch --binary -p1 <extra args> -i <file>` in the source tree.
        function(_diagham_run_patch patch_file result_var error_var)
            execute_process(
                COMMAND ${PATCH_EXECUTABLE} --binary -p1 --silent ${ARGN} -i ${patch_file}
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                RESULT_VARIABLE result
                OUTPUT_VARIABLE output
                ERROR_VARIABLE error
            )
            set(${result_var} ${result} PARENT_SCOPE)
            set(${error_var} "${output}${error}" PARENT_SCOPE)
        endfunction()

        foreach(p ${patch_files})
            get_filename_component(pname ${p} NAME)

            # Candidate forms of this patch: as committed (LF) and, for the
            # one patch whose target is CRLF upstream, a CRLF-converted copy
            # (see the file header for why this is targeted, not global).
            set(candidates ${p})
            if(Python3_EXECUTABLE)
                set(crlf_copy ${CMAKE_BINARY_DIR}/_diagham_patch_crlf_retry_${pname})
                execute_process(
                    COMMAND ${Python3_EXECUTABLE} -c
                        "import sys; data=open(sys.argv[1],'rb').read(); data=data.replace(b'\\r\\n', b'\\n').replace(b'\\n', b'\\r\\n'); open(sys.argv[2],'wb').write(data)"
                        ${p} ${crlf_copy}
                    RESULT_VARIABLE convert_result
                )
                if(convert_result EQUAL 0)
                    list(APPEND candidates ${crlf_copy})
                endif()
            endif()

            set(status "")
            set(errors "")
            foreach(candidate ${candidates})
                # Dry run first, so a candidate that only half-matches never
                # leaves a partially patched file behind.
                _diagham_run_patch(${candidate} patch_result patch_error --forward --dry-run)
                if(patch_result EQUAL 0)
                    _diagham_run_patch(${candidate} patch_result patch_error --forward)
                endif()
                if(patch_result EQUAL 0)
                    if(candidate STREQUAL p)
                        set(status "applied")
                    else()
                        set(status "applied (CRLF-converted retry -- its target file is CRLF upstream)")
                    endif()
                    break()
                endif()
                string(APPEND errors "  ${candidate}: ${patch_error}\n")
            endforeach()

            if(status)
                message(STATUS "DiagHam: ${status}: ${pname}")
            else()
                message(FATAL_ERROR
                    "DiagHam: failed to apply ${pname}\n${errors}"
                    "  (If no CRLF retry is listed, no Python3 interpreter was found -- "
                    "if this patch needs it, install python3.)")
            endif()
        endforeach()

        file(WRITE ${DIAGHAM_PATCH_SENTINEL} "${patch_fingerprint}")
    endif()
endif()

# Files excluded from the build because they cannot be patched without
# physics-domain decisions that belong to the maintainer. See PATCHES.md.
set(DIAGHAM_UPSTREAM_EXCLUDED_PROGRAMS
    "QHEFermionsTorusWithSpin"         # legacy duplicate of FQHETorusFermionsWithSpin (see DEFERRED.md)
    CACHE INTERNAL "Programs excluded due to upstream issues"
)
