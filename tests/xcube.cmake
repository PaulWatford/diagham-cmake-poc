# ----------------------------------------------------------------------------
# X-cube model entanglement entropy (label: physics, entanglement, xcube) against
# tests/oracles/xcube.py: the state is a uniform superposition over a GF(2) subspace, so
# S(A)/ln 2 = rank P_A V + rank P_B V - rank V (nothing shared with DiagHam; the oracle also checks
# one case against an explicit reduced density matrix in mpmath). The oracle agrees with the program
# on all eight cases below in the default and the LAPACK build. A region has three spins per site,
# so there is no 1x1 block and U42 is not reached. The 3x3x3 system needs 2^48 components in the
# program (the GF(2) rank is 46), the same defect as U43.
# ----------------------------------------------------------------------------
set(xcube_data ${DIAGHAM_TEST_DATA}/xcube)
set(xcube Spin_XCubeEntanglementEntropy)
function(diagham_xcube_test name nx ny nz ax ay az)
    cmake_parse_arguments(ARG "WILL_FAIL" "" "LABELS" ${ARGN})
    if(ARG_WILL_FAIL)
        set(wf WILL_FAIL)
    else()
        set(wf "")
    endif()
    diagham_physics_test(${name}
        PROGRAM ${xcube}
        ARGS -x ${nx} -y ${ny} -z ${nz} --nbra-sitex ${ax} --nbra-sitey ${ay} --nbra-sitez ${az} --export-entspectrum
        OUTPUT "xcube_entspectrum_*.dat"
        CHECK line @OUTPUT@ "# Entangement entropy" ${xcube_data}/xcube_${nx}x${ny}x${nz}_A${ax}x${ay}x${az}.txt 1e-9
        ${wf}
        LABELS entanglement spin xcube ${ARG_LABELS})
endfunction()

diagham_xcube_test(physics.entanglement.xcube.2x2x2_region_1x1x1 2 2 2 1 1 1)
diagham_xcube_test(physics.entanglement.xcube.2x2x2_region_1x1x2 2 2 2 1 1 2)
diagham_xcube_test(physics.entanglement.xcube.2x2x2_region_1x2x2 2 2 2 1 2 2)
diagham_xcube_test(physics.entanglement.xcube.2x2x2_region_2x1x1 2 2 2 2 1 1)
diagham_xcube_test(physics.entanglement.xcube.2x2x3_region_1x1x2 2 2 3 1 1 2)
diagham_xcube_test(physics.entanglement.xcube.2x3x2_region_1x2x1 2 3 2 1 2 1)
diagham_xcube_test(physics.entanglement.xcube.3x2x2_region_2x1x1 3 2 2 2 1 1)
diagham_xcube_test(physics.entanglement.xcube.3x2x2_region_1x1x1 3 2 2 1 1 1)
# U43: the array is sized 2^(number of Z terms) = 2^48 for 3x3x3 (rank 46); oracle entropy of a 2x2x2 region is 18
diagham_xcube_test(knownbug.xcube.ground_state_array_2_pow_nterms.3x3x3_region_2x2x2 3 3 3 2 2 2 WILL_FAIL LABELS known-bug)

if(Python3_Interpreter_FOUND)
    add_test(NAME selftest.xcube_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/xcube.py --check ${xcube_data})
    set_tests_properties(selftest.xcube_oracle PROPERTIES LABELS "selftest;python" TIMEOUT 300)
endif()
