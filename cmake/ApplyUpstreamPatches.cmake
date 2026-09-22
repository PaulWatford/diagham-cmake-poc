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
#   1. `patch --binary -p1` on the patch file as committed.
#   2. If that fails, byte-convert a copy of the SAME patch file to CRLF and
#      retry with `patch --binary -p1` again.
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
set(DIAGHAM_PATCH_SENTINEL ${CMAKE_BINARY_DIR}/.diagham_upstream_patches_applied)

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
        endif()
    endif()

    if(need_apply)
        find_program(PATCH_EXECUTABLE patch REQUIRED)
        find_package(Python3 COMPONENTS Interpreter QUIET)

        foreach(p ${patch_files})
            get_filename_component(pname ${p} NAME)
            execute_process(
                COMMAND ${PATCH_EXECUTABLE} --binary -p1 --silent -i ${p}
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                RESULT_VARIABLE patch_result
                OUTPUT_VARIABLE patch_output
                ERROR_VARIABLE patch_error
            )

            if(NOT patch_result EQUAL 0 AND Python3_EXECUTABLE)
                # This specific patch's LF line endings didn't match its
                # target file's CRLF content. Byte-convert a scratch copy of
                # the patch to CRLF and retry with --binary (still needed:
                # --binary stops `patch` re-stripping the CRs we just added).
                # See the file header for why this is targeted at one patch,
                # not applied to all of them.
                set(crlf_copy ${CMAKE_BINARY_DIR}/_diagham_patch_crlf_retry.patch)
                execute_process(
                    COMMAND ${Python3_EXECUTABLE} -c
                        "import sys; data=open(sys.argv[1],'rb').read(); data=data.replace(b'\\r\\n', b'\\n').replace(b'\\n', b'\\r\\n'); open(sys.argv[2],'wb').write(data)"
                        ${p} ${crlf_copy}
                    RESULT_VARIABLE convert_result
                )
                if(convert_result EQUAL 0)
                    execute_process(
                        COMMAND ${PATCH_EXECUTABLE} --binary -p1 --silent -i ${crlf_copy}
                        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                        RESULT_VARIABLE retry_result
                        OUTPUT_VARIABLE retry_output
                        ERROR_VARIABLE retry_error
                    )
                    if(retry_result EQUAL 0)
                        set(patch_result 0)
                        set(pname "${pname} (CRLF-converted retry -- its target file is CRLF upstream)")
                    else()
                        set(patch_error
                            "as-committed (LF): ${patch_error}\nCRLF retry: ${retry_error}")
                    endif()
                endif()
            endif()

            if(patch_result EQUAL 0)
                message(STATUS "DiagHam: applied ${pname}")
            else()
                message(FATAL_ERROR
                    "DiagHam: failed to apply ${pname}\n"
                    "  result: ${patch_result}\n"
                    "  output: ${patch_output}\n"
                    "  error:  ${patch_error}\n"
                    "  (No Python3 interpreter found for the CRLF-retry path -- "
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
