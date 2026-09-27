# ----------------------------------------------------------------------------
# FQHE on the sphere: goldens with answers from tests/oracles/fqhe_sphere.py
# (label: physics, fqhe, sphere). Expected values live in
# tests/data/fqhe_sphere/ and are re-derived by selftest.fqhe_sphere_oracle.
#
#  1. Zero-mode counts: the (k+1)-body hard-core zero modes are counted by
#     (k,r)-admissible configurations (Read-Rezayi / Bernevig-Haldane);
#     zero_modes.txt gives the count per case, read here so no number is
#     duplicated.
#  2. Coulomb pseudopotentials: CoulombPseudopotentials vs the Wigner-3j/6j
#     closed form, all 2S+1 values to 1e-13.
#  3. Model states from Jack polynomials vs exact diagonalisation: the
#     overlap must be 1 (Laughlin fermions and bosons, Moore-Read bosons).
#  4. Particle entanglement spectrum of the Laughlin state: level counts
#     per N_A (and per Lz_A) equal the quasihole counting.
# ----------------------------------------------------------------------------
set(sphere_data ${DIAGHAM_TEST_DATA}/fqhe_sphere)
set(v0 ${DIAGHAM_TEST_DATA}/v0_pseudopotential.dat)

# count for a case, from zero_modes.txt ("name 2Lz count # ...")
function(diagham_zero_mode_count case out_two_lz out_count)
    file(STRINGS ${sphere_data}/zero_modes.txt lines REGEX "^${case} ")
    list(LENGTH lines n)
    if(NOT n EQUAL 1)
        message(FATAL_ERROR "zero_modes.txt: case ${case} not found exactly once")
    endif()
    string(REGEX MATCH "^${case} (-?[0-9]+) ([0-9]+)" _ "${lines}")
    set(${out_two_lz} ${CMAKE_MATCH_1} PARENT_SCOPE)
    set(${out_count} ${CMAKE_MATCH_2} PARENT_SCOPE)
endfunction()

# diagham_zero_mode_test(<name> CASE <case> PROGRAM <target> ARGS <args before Lz options> OUTPUT <glob>)
function(diagham_zero_mode_test name)
    cmake_parse_arguments(ARG "" "CASE;PROGRAM;OUTPUT" "ARGS" ${ARGN})
    diagham_zero_mode_count(${ARG_CASE} two_lz count)
    diagham_physics_test(${name}
        PROGRAM ${ARG_PROGRAM}
        ARGS ${ARG_ARGS} --initial-lz ${two_lz} --nbr-lz 1 --full-diag 20000
        OUTPUT "${ARG_OUTPUT}"
        CHECK count @OUTPUT@ -1 0 1e-10 ${count}
        LABELS fqhe sphere)
endfunction()

# --- 1. zero modes -------------------------------------------------------------
set(bos2 FQHEOnSphere_FQHESphereBosonsTwoBodyGeneric)
set(bosN FQHEOnSphere_FQHESphereBosonsNBodyHardCore)
set(fer2 FQHEOnSphere_FQHESphereFermionsTwoBodyGeneric)
set(ferN FQHEOnSphere_QHEFermionsNBodyHardCore)
diagham_zero_mode_test(physics.fqhe.sphere.bosons_laughlin_1_2.unique_zero_mode CASE bosons_laughlin_1_2_n6_2s10
    PROGRAM ${bos2} ARGS -p 6 -l 10 --interaction-file ${v0} --interaction-name v0 OUTPUT "bosons_v0_n_6_2s_10_lz.dat")
diagham_zero_mode_test(physics.fqhe.sphere.bosons_laughlin_1_2.two_quasihole_count CASE bosons_laughlin_1_2_n6_2s12
    PROGRAM ${bos2} ARGS -p 6 -l 12 --interaction-file ${v0} --interaction-name v0 OUTPUT "bosons_v0_n_6_2s_12_lz.dat")
diagham_zero_mode_test(physics.fqhe.sphere.bosons_laughlin_1_2.one_quasihole_half_integer_lz CASE bosons_laughlin_1_2_n7_2s13
    PROGRAM ${bos2} ARGS -p 7 -l 13 --interaction-file ${v0} --interaction-name v0 OUTPUT "bosons_v0_n_7_2s_13_lz.dat")
diagham_zero_mode_test(physics.fqhe.sphere.bosons_mooreread.unique_zero_mode CASE bosons_mooreread_n6_2s4
    PROGRAM ${bosN} ARGS -p 6 -l 4 --nbr-nbody 3 OUTPUT "bosons_hardcore_nbody_3_n_6_2s_4_lz.dat")
diagham_zero_mode_test(physics.fqhe.sphere.bosons_mooreread.two_quasihole_count CASE bosons_mooreread_n6_2s5
    PROGRAM ${bosN} ARGS -p 6 -l 5 --nbr-nbody 3 OUTPUT "bosons_hardcore_nbody_3_n_6_2s_5_lz.dat")
diagham_zero_mode_test(physics.fqhe.sphere.bosons_mooreread.unique_zero_mode_n8 CASE bosons_mooreread_n8_2s6
    PROGRAM ${bosN} ARGS -p 8 -l 6 --nbr-nbody 3 OUTPUT "bosons_hardcore_nbody_3_n_8_2s_6_lz.dat")
diagham_zero_mode_test(physics.fqhe.sphere.bosons_readrezayi3.quasihole_count CASE bosons_readrezayi3_n6_2s3
    PROGRAM ${bosN} ARGS -p 6 -l 3 --nbr-nbody 4 OUTPUT "bosons_hardcore_nbody_4_n_6_2s_3_lz.dat")
diagham_zero_mode_test(physics.fqhe.sphere.bosons_readrezayi3.unique_zero_mode CASE bosons_readrezayi3_n9_2s4
    PROGRAM ${bosN} ARGS -p 9 -l 4 --nbr-nbody 4 OUTPUT "bosons_hardcore_nbody_4_n_9_2s_4_lz.dat")
diagham_zero_mode_test(physics.fqhe.sphere.laughlin_1_3.three_quasihole_count CASE fermions_laughlin_1_3_n6_2s18
    PROGRAM ${fer2} ARGS -p 6 -l 18 --interaction-file ${v1} --interaction-name v1 OUTPUT "fermions_v1_n_6_2s_18_lz.dat")
diagham_zero_mode_test(physics.fqhe.sphere.fermions_mooreread.unique_zero_mode CASE fermions_mooreread_n6_2s9
    PROGRAM ${ferN} ARGS -p 6 -l 9 --nbr-nbody 3 OUTPUT "fermions_hardcore_nbody_3_n_6_2s_9_lz.dat")
diagham_zero_mode_test(physics.fqhe.sphere.fermions_mooreread.two_quasihole_count CASE fermions_mooreread_n6_2s10
    PROGRAM ${ferN} ARGS -p 6 -l 10 --nbr-nbody 3 OUTPUT "fermions_hardcore_nbody_3_n_6_2s_10_lz.dat")

# --- 2. Coulomb pseudopotentials -----------------------------------------------------
foreach(two_s 6 10 15)
    diagham_physics_test(physics.fqhe.sphere.coulomb_pseudopotentials.2s${two_s}
        PROGRAM FQHEOnSphere_CoulombPseudopotentials
        ARGS -s ${two_s}
        OUTPUT "pseudopotential_coulomb_l_0_2s_${two_s}.dat"
        CHECK line @OUTPUT@ Pseudopotentials ${sphere_data}/coulomb_2s${two_s}.txt 1e-13
        LABELS fqhe sphere)
endforeach()

# --- 3. Jack polynomial vs exact diagonalisation ------------------------------------------
# Each chain: ED eigenstate -> Jack generator (squeezed basis) -> full basis -> overlap.
# GenericOverlap --quiet prints the squared overlap alone; program.log is checked.
set(jack FQHEOnSphere_FQHESphereJackGenerator)
set(conv FQHEOnSphere_FQHESphereConvertHaldaneBasis)
set(overlap Programs_GenericOverlap)
diagham_physics_test(physics.fqhe.sphere.laughlin_1_3.jack_equals_ed
    PRE_STEPS
        "${fer2}|-p|6|-l|15|--interaction-file|${v1}|--interaction-name|v1|--nbr-lz|1|--full-diag|2000|--eigenstate|-n|1"
        "${jack}|-a|-2|--fermion|-n|--reference-file|${sphere_data}/root_laughlin_1_3_n6_2s15.dat|-o|fermions_haldane_jack_n_6_2s_15_lz_0.0.vec"
        "${conv}|fermions_haldane_jack_n_6_2s_15_lz_0.0.vec|--reference-file|${sphere_data}/root_laughlin_1_3_n6_2s15.dat|-o|fermions_jack_n_6_2s_15_lz_0.0.vec"
    PROGRAM ${overlap}
    ARGS --quiet fermions_jack_n_6_2s_15_lz_0.0.vec fermions_v1_n_6_2s_15_lz_0.0.vec
    OUTPUT "program.log"
    CHECK min-abs @OUTPUT@ 0 1 1e-12
    LABELS fqhe sphere)
diagham_physics_test(physics.fqhe.sphere.bosons_laughlin_1_2.jack_equals_ed
    PRE_STEPS
        "${bos2}|-p|6|-l|10|--interaction-file|${v0}|--interaction-name|v0|--nbr-lz|1|--full-diag|5000|--eigenstate|-n|1"
        "${jack}|-a|-2|-n|--reference-file|${sphere_data}/root_laughlin_1_2_n6_2s10.dat|-o|bosons_haldane_jack_n_6_2s_10_lz_0.0.vec"
        "${conv}|bosons_haldane_jack_n_6_2s_10_lz_0.0.vec|--reference-file|${sphere_data}/root_laughlin_1_2_n6_2s10.dat|-o|bosons_jack_n_6_2s_10_lz_0.0.vec"
    PROGRAM ${overlap}
    ARGS --quiet bosons_jack_n_6_2s_10_lz_0.0.vec bosons_v0_n_6_2s_10_lz_0.0.vec
    OUTPUT "program.log"
    CHECK min-abs @OUTPUT@ 0 1 1e-12
    LABELS fqhe sphere)
diagham_physics_test(physics.fqhe.sphere.bosons_mooreread.jack_equals_ed
    PRE_STEPS
        "${bosN}|-p|6|-l|4|--nbr-nbody|3|--initial-lz|0|--nbr-lz|1|--full-diag|5000|--eigenstate|-n|1"
        "${jack}|-a|-3|-n|--reference-file|${sphere_data}/root_mooreread_n6_2s4.dat|-o|bosons_haldane_jack_n_6_2s_4_lz_0.0.vec"
        "${conv}|bosons_haldane_jack_n_6_2s_4_lz_0.0.vec|--reference-file|${sphere_data}/root_mooreread_n6_2s4.dat|-o|bosons_jack_n_6_2s_4_lz_0.0.vec"
    PROGRAM ${overlap}
    ARGS --quiet bosons_jack_n_6_2s_4_lz_0.0.vec bosons_hardcore_nbody_3_n_6_2s_4_lz_0.0.vec
    OUTPUT "program.log"
    CHECK min-abs @OUTPUT@ 0 1 1e-12
    LABELS fqhe sphere)

# --- 4. particle entanglement spectrum of the Laughlin state -------------------------------
# pes_laughlin_n6_2s15.txt: "name N_A 2Lz|all count"; the density-matrix file has columns N_A 2Lz lambda
set(pes FQHEOnSphere_FQHESphereFermionEntanglementEntropyParticlePartition)
file(STRINGS ${sphere_data}/pes_laughlin_n6_2s15.txt pes_lines REGEX "^pes_")
foreach(line ${pes_lines})
    string(REGEX MATCH "^([a-z0-9_-]+) ([0-9]+) ([a-z0-9-]+) ([0-9]+)" _ "${line}")
    set(case ${CMAKE_MATCH_1})
    set(na ${CMAKE_MATCH_2})
    set(lz ${CMAKE_MATCH_3})
    set(count ${CMAKE_MATCH_4})
    set(filters 0 ${na})
    if(NOT lz STREQUAL "all")
        list(APPEND filters 1 ${lz})
    endif()
    diagham_physics_test(physics.fqhe.sphere.laughlin_1_3.${case}
        PRE_STEPS
            "${fer2}|-p|6|-l|15|--interaction-file|${v1}|--interaction-name|v1|--nbr-lz|1|--full-diag|2000|--eigenstate|-n|1"
        PROGRAM ${pes}
        ARGS fermions_v1_n_6_2s_15_lz_0.0.vec --max-na ${na} --density-matrix pes.dat
        OUTPUT "pes.dat"
        CHECK nonzero @OUTPUT@ 2 1e-10 ${count} ${filters}
        LABELS fqhe sphere)
endforeach()

if(Python3_Interpreter_FOUND)
    add_test(NAME selftest.fqhe_sphere_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/fqhe_sphere.py --check ${sphere_data})
    set_tests_properties(selftest.fqhe_sphere_oracle PROPERTIES LABELS "selftest;python" TIMEOUT 300)
endif()
