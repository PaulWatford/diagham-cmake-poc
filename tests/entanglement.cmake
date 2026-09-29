# ----------------------------------------------------------------------------
# Entanglement entropies and spectra (label: physics, entanglement) against
# tests/oracles/entanglement.py:
#  - toric code: entropy of a region from GF(2) ranks of the vertex-mask span (independent linear
#    algebra on the state the program builds), compared with "# Entangement entropy = S * log 2";
#  - torus particle entanglement spectrum of the Laughlin 1/3 state (N = 4, Nphi = 12, three ground
#    states): the number of non-zero levels per Ky sector equals the cyclic admissible counting;
#  - spin-1/2 chain Bell pair (ln 2, two levels 1/2) and product state (0), which need the LAPACK
#    build (U40: without it the program aborts on the ungated --use-svd option).
# ----------------------------------------------------------------------------
set(ent_data ${DIAGHAM_TEST_DATA}/entanglement)

# --- toric code -----------------------------------------------------------------------
set(toric Spin_2DToricCodeEntanglementEntropy)
function(diagham_toric_test name case)
    cmake_parse_arguments(ARG "" "" "ARGS" ${ARGN})
    diagham_physics_test(physics.entanglement.toric_code.${name}
        PROGRAM ${toric}
        ARGS ${ARG_ARGS} --export-entspectrum
        OUTPUT "2dtoriccode_entspectrum_*.dat"
        CHECK line @OUTPUT@ "# Entangement entropy" ${ent_data}/toric_${case}.txt 1e-9
        LABELS entanglement spin toric_code)
endfunction()
diagham_toric_test(3x3_region_2x2 3x3_A2x2 ARGS -x 3 -y 3 --nbra-sitex 2 --nbra-sitey 2)
diagham_toric_test(3x3_region_1x1 3x3_A1x1 ARGS -x 3 -y 3 --nbra-sitex 1 --nbra-sitey 1)
diagham_toric_test(4x3_region_2x2 4x3_A2x2 ARGS -x 4 -y 3 --nbra-sitex 2 --nbra-sitey 2)
diagham_toric_test(3x4_region_2x3 3x4_A2x3 ARGS -x 3 -y 4 --nbra-sitex 2 --nbra-sitey 3)
diagham_toric_test(3x3_wrapping_column 3x3_A1x3 ARGS -x 3 -y 3 --nbra-sitex 1 --nbra-sitey 3)
diagham_toric_test(3x3_region_2x2_other_parity_sector 3x3_A2x2 ARGS -x 3 -y 3 --nbra-sitex 2 --nbra-sitey 2 --gs-parity 1)
diagham_toric_test(3x3_region_2x2_low_memory 3x3_A2x2 ARGS -x 3 -y 3 --nbra-sitex 2 --nbra-sitey 2 --low-memory)

# --- torus particle entanglement spectrum, Laughlin 1/3 -------------------------------------
# counts: torus_pes_laughlin_n4_2s12_na<N_A>.txt lines "ky count"; the density-matrix file has columns N_A Ky lambda
set(tf2_pes FQHEOnTorus_FQHETorusFermionsTwoBodyGeneric)
set(torus_pes FQHEOnTorus_FQHETorusEntanglementEntropyParticlePartition)
function(diagham_torus_pes_test name)
    cmake_parse_arguments(ARG "" "NA;KY;COUNT" "" ${ARGN})
    set(filters 0 ${ARG_NA})
    if(NOT "${ARG_KY}" STREQUAL "")
        list(APPEND filters 1 ${ARG_KY})
    endif()
    diagham_physics_test(physics.entanglement.torus.laughlin_1_3.${name}
        PRE_STEPS
            "${tf2_pes}|-p|4|-l|12|--interaction-file|${DIAGHAM_TEST_DATA}/v1_pseudopotential.dat|--interaction-name|v1|--redundant-kymomenta|--full-diag|5000|--eigenstate|-n|1"
        PROGRAM ${torus_pes}
        ARGS --degenerated-groundstate ${ent_data}/torus_laughlin_n4_2s12_ground_states.txt -p 4 -l 12 --max-na 2 --density-matrix pes.dat
        OUTPUT "pes.dat"
        CHECK nonzero @OUTPUT@ 2 1e-10 ${ARG_COUNT} ${filters}
        TIMEOUT 300
        LABELS fqhe torus entanglement)
endfunction()
function(diagham_torus_pes_count na ky out)
    file(STRINGS ${ent_data}/torus_pes_laughlin_n4_2s12_na${na}.txt lines REGEX "^${ky} ")
    string(REGEX MATCH "^${ky} ([0-9]+)" _ "${lines}")
    set(${out} ${CMAKE_MATCH_1} PARENT_SCOPE)
endfunction()
diagham_torus_pes_count(1 0 c)
diagham_torus_pes_test(one_particle_one_level_per_sector NA 1 KY 0 COUNT ${c})
diagham_torus_pes_test(one_particle_twelve_levels NA 1 COUNT 12)
diagham_torus_pes_count(2 0 c)
diagham_torus_pes_test(two_particles_ky0_admissible NA 2 KY 0 COUNT ${c})
diagham_torus_pes_count(2 1 c)
diagham_torus_pes_test(two_particles_ky1_admissible NA 2 KY 1 COUNT ${c})

# --- spin-1/2 chain: Bell pair and product state ---------------------------------------------
# hand-written binary vectors (int32 dimension + doubles): (|01> - |10>)/sqrt 2 in n_2_sz_0.0, |1111> in n_4_sz_4.0
set(spin_ent Spin_SpinChainEntanglementEntropy)
function(diagham_spin_entropy_tests prefix)
    cmake_parse_arguments(SE "WILL_FAIL" "" "LABELS" ${ARGN})
    set(EXTRA ${SE_LABELS})
    if(SE_WILL_FAIL)
        set(WF WILL_FAIL)
    else()
        set(WF "")
    endif()
    diagham_physics_test(${prefix}.bell_pair_entropy_ln2
        PROGRAM ${spin_ent}
        ARGS -s 1 ${ent_data}/spin_1_2_n_2_sz_0.0.vec -o bell.ent
        OUTPUT "bell.ent"
        CHECK min-abs @OUTPUT@ 1 0.6931471805599453 1e-12
        ${WF} LABELS entanglement spin ${EXTRA})
    diagham_physics_test(${prefix}.product_state_entropy_zero
        PROGRAM ${spin_ent}
        ARGS -s 1 ${ent_data}/spin_1_2_n_4_sz_4.0.vec -o product.ent
        OUTPUT "product.ent"
        CHECK min-abs @OUTPUT@ 1 0 1e-12
        ${WF} LABELS entanglement spin ${EXTRA})
endfunction()
if(DIAGHAM_USE_LAPACK)
    diagham_spin_entropy_tests(physics.entanglement.spin_chain)
    diagham_physics_test(physics.entanglement.spin_chain.bell_pair_two_levels_one_half
        PROGRAM ${spin_ent}
        ARGS -s 1 ${ent_data}/spin_1_2_n_2_sz_0.0.vec -o bell.ent
        OUTPUT "bell.full.ent"
        CHECK nonzero @OUTPUT@ 2 1e-10 2 0 1
        LABELS entanglement spin )
    diagham_physics_test(physics.entanglement.spin_chain.product_state_one_level
        PROGRAM ${spin_ent}
        ARGS -s 1 ${ent_data}/spin_1_2_n_4_sz_4.0.vec -o product.ent
        OUTPUT "product.full.ent"
        CHECK nonzero @OUTPUT@ 2 1e-10 1 0 1
        LABELS entanglement spin )
else()
    # U40: use-svd is read unconditionally but only registered with LAPACK, so without it the program
    # stops with "Option 'use-svd' was requested, but is not implemented!". WILL_FAIL until gated.
    diagham_spin_entropy_tests(knownbug.spin_chain_entanglement.use_svd_ungated WILL_FAIL LABELS known-bug)
endif()

if(Python3_Interpreter_FOUND)
    add_test(NAME selftest.entanglement_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/entanglement.py --check ${ent_data})
    set_tests_properties(selftest.entanglement_oracle PROPERTIES LABELS "selftest;python" TIMEOUT 300)
endif()
