# ----------------------------------------------------------------------------
# Haah cubic code entanglement entropy (label: physics, entanglement, haah_code) against
# tests/oracles/haah.py: the program's state is a uniform superposition over a GF(2) subspace,
# so S(A)/ln 2 = rank P_A V + rank P_B V - rank V, computed with nothing shared with DiagHam and
# cross-checked there against an explicit dense reduced density matrix diagonalized with mpmath.
# Cases with a one-site region give a 1x1 density-matrix block, which crashes the internal
# diagonaliser (U42), so they are physics tests only with LAPACK and known-bug tests without.
# The 4x4x4 system needs 2^63 components in the program (U43).
# ----------------------------------------------------------------------------
set(haah_data ${DIAGHAM_TEST_DATA}/haah)
set(haah Spin_HaahCodeEntropy)
function(diagham_haah_test name nx ny nz ax ay az)
    cmake_parse_arguments(ARG "WILL_FAIL" "DATA" "LABELS;EXTRA" ${ARGN})
    if(ARG_DATA)
        set(data ${ARG_DATA})
    else()
        set(data haah_${nx}x${ny}x${nz}_A${ax}x${ay}x${az})
    endif()
    if(ARG_WILL_FAIL)
        set(wf WILL_FAIL)
    else()
        set(wf "")
    endif()
    diagham_physics_test(${name}
        PROGRAM ${haah}
        ARGS -x ${nx} -y ${ny} -z ${nz} --nbra-sitex ${ax} --nbra-sitey ${ay} --nbra-sitez ${az} --export-entspectrum ${ARG_EXTRA}
        OUTPUT "haahcode_entspectrum_*.dat"
        CHECK line @OUTPUT@ "# Entangement entropy" ${haah_data}/${data}.txt 1e-9
        ${wf}
        LABELS entanglement spin haah_code ${ARG_LABELS})
endfunction()

diagham_haah_test(physics.entanglement.haah_code.2x2x2_region_1x1x2 2 2 2 1 1 2)
diagham_haah_test(physics.entanglement.haah_code.2x2x2_region_1x1x2_low_memory 2 2 2 1 1 2 EXTRA --low-memory DATA haah_2x2x2_A1x1x2)
diagham_haah_test(physics.entanglement.haah_code.2x2x2_region_2x1x1 2 2 2 2 1 1)
diagham_haah_test(physics.entanglement.haah_code.2x2x2_region_1x2x2 2 2 2 1 2 2)
diagham_haah_test(physics.entanglement.haah_code.2x2x3_region_1x1x2 2 2 3 1 1 2)
diagham_haah_test(physics.entanglement.haah_code.2x3x2_region_1x2x1 2 3 2 1 2 1)
diagham_haah_test(physics.entanglement.haah_code.3x2x2_region_2x1x1 3 2 2 2 1 1)
diagham_haah_test(physics.entanglement.haah_code.4x2x2_region_2x2x1 4 2 2 2 2 1)

if(DIAGHAM_USE_LAPACK)
    diagham_haah_test(physics.entanglement.haah_code.2x2x2_region_1x1x1 2 2 2 1 1 1)
else()
    # U42: RealSymmetricMatrix::Householder writes outside its arrays for a 1x1 matrix
    diagham_haah_test(knownbug.haah_code.householder_1x1_crash.2x2x2_region_1x1x1 2 2 2 1 1 1 WILL_FAIL LABELS known-bug DATA haah_2x2x2_A1x1x1)
endif()
# U43: the ground-state array is sized 2^(number of Z terms) = 2^63 for 4x4x4 (the GF(2) rank is 57);
# the oracle entropy of a 2x2x2 cube is 14 (6 a^2 - 6 a + 2)
diagham_physics_test(knownbug.haah_code.ground_state_array_2_pow_nterms.4x4x4_region_2x2x2
    PROGRAM ${haah}
    ARGS -x 4 -y 4 -z 4 --nbra-sitex 2 --nbra-sitey 2 --nbra-sitez 2 --export-entspectrum
    OUTPUT "haahcode_entspectrum_*.dat"
    CHECK line @OUTPUT@ "# Entangement entropy" ${haah_data}/haah_4x4x4_A2x2x2.txt 1e-9
    WILL_FAIL LABELS entanglement spin haah_code known-bug)

if(Python3_Interpreter_FOUND)
    add_test(NAME selftest.haah_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/haah.py --check ${haah_data})
    set_tests_properties(selftest.haah_oracle PROPERTIES LABELS "selftest;python" TIMEOUT 300)
endif()
