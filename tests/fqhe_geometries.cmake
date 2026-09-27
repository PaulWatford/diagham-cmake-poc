# ----------------------------------------------------------------------------
# FQHE on the torus, the cylinder and the disk: zero-mode counts against
# tests/oracles/fqhe_geometries.py (label: physics, fqhe, <geometry>).
# The count for each case is read from tests/data/fqhe_geometries/zero_modes.txt
# ("name count # ..."), re-derived by selftest.fqhe_geometries_oracle.
#
# Torus programs run every momentum sector (--redundant-kymomenta, or the
# magnetic-translation programs) and the count is the total; cylinder
# programs list the sectors 2Ky >= 0; disk programs run one Lz.
# ----------------------------------------------------------------------------
set(geom_data ${DIAGHAM_TEST_DATA}/fqhe_geometries)

function(diagham_geometry_count case out_count)
    file(STRINGS ${geom_data}/zero_modes.txt lines REGEX "^${case} ")
    list(LENGTH lines n)
    if(NOT n EQUAL 1)
        message(FATAL_ERROR "fqhe_geometries/zero_modes.txt: case ${case} not found exactly once")
    endif()
    string(REGEX MATCH "^${case} ([0-9]+)" _ "${lines}")
    set(${out_count} ${CMAKE_MATCH_1} PARENT_SCOPE)
endfunction()

# diagham_geometry_test(<name> CASE <case> GEOMETRY <label> PROGRAM <target> ARGS <args...> OUTPUT <glob>)
function(diagham_geometry_test name)
    cmake_parse_arguments(ARG "" "CASE;GEOMETRY;PROGRAM;OUTPUT;TIMEOUT" "ARGS" ${ARGN})
    diagham_geometry_count(${ARG_CASE} count)
    if(NOT ARG_TIMEOUT)
        set(ARG_TIMEOUT 300)
    endif()
    diagham_physics_test(${name}
        PROGRAM ${ARG_PROGRAM}
        ARGS ${ARG_ARGS}
        OUTPUT "${ARG_OUTPUT}"
        CHECK count @OUTPUT@ -1 0 1e-10 ${count}
        TIMEOUT ${ARG_TIMEOUT}
        LABELS fqhe ${ARG_GEOMETRY})
endfunction()

# --- torus ------------------------------------------------------------------------
set(tb2 FQHEOnTorus_FQHETorusBosonsTwoBodyGeneric)
set(tf2 FQHEOnTorus_FQHETorusFermionsTwoBodyGeneric)
set(tbN FQHEOnTorus_FQHETorusBosonsWithTranslationsNBodyHardCore)
set(tfN FQHEOnTorus_FQHETorusFermionsWithTranslationsNBodyHollowCore)
diagham_geometry_test(physics.fqhe.torus.bosons_laughlin_1_2.twofold_degeneracy CASE torus_bosons_laughlin_1_2_n4_nphi8 GEOMETRY torus
    PROGRAM ${tb2} ARGS -p 4 -l 8 --interaction-file ${v0} --interaction-name v0 --redundant-kymomenta --full-diag 5000
    OUTPUT "bosons_torus_kysym_v0_n_4_2s_8_ratio_1.000000.dat")
diagham_geometry_test(physics.fqhe.torus.bosons_laughlin_1_2.twofold_degeneracy_n5 CASE torus_bosons_laughlin_1_2_n5_nphi10 GEOMETRY torus
    PROGRAM ${tb2} ARGS -p 5 -l 10 --interaction-file ${v0} --interaction-name v0 --redundant-kymomenta --full-diag 5000
    OUTPUT "bosons_torus_kysym_v0_n_5_2s_10_ratio_1.000000.dat")
diagham_geometry_test(physics.fqhe.torus.bosons_laughlin_1_2.one_quasihole_count CASE torus_bosons_laughlin_1_2_n4_nphi9 GEOMETRY torus
    PROGRAM ${tb2} ARGS -p 4 -l 9 --interaction-file ${v0} --interaction-name v0 --redundant-kymomenta --full-diag 5000
    OUTPUT "bosons_torus_kysym_v0_n_4_2s_9_ratio_1.000000.dat")
diagham_geometry_test(physics.fqhe.torus.laughlin_1_3.one_quasihole_count CASE torus_fermions_laughlin_1_3_n4_nphi13 GEOMETRY torus
    PROGRAM ${tf2} ARGS -p 4 -l 13 --interaction-file ${v1} --interaction-name v1 --redundant-kymomenta --full-diag 5000
    OUTPUT "fermions_torus_kysym_v1_n_4_2s_13_ratio_1.000000.dat")
diagham_geometry_test(physics.fqhe.torus.laughlin_1_3.two_quasihole_count CASE torus_fermions_laughlin_1_3_n4_nphi14 GEOMETRY torus
    PROGRAM ${tf2} ARGS -p 4 -l 14 --interaction-file ${v1} --interaction-name v1 --redundant-kymomenta --full-diag 5000
    OUTPUT "fermions_torus_kysym_v1_n_4_2s_14_ratio_1.000000.dat")
diagham_geometry_test(physics.fqhe.torus.bosons_mooreread.threefold_degeneracy CASE torus_bosons_mooreread_n4_nphi4 GEOMETRY torus
    PROGRAM ${tbN} ARGS -p 4 -l 4 --nbr-nbody 3 --all-points --full-diag 5000
    OUTPUT "bosons_torus_3body_hardcore_n_4_2s_4_ratio_1.000000.dat")
diagham_geometry_test(physics.fqhe.torus.bosons_mooreread.threefold_degeneracy_n6 CASE torus_bosons_mooreread_n6_nphi6 GEOMETRY torus
    PROGRAM ${tbN} ARGS -p 6 -l 6 --nbr-nbody 3 --all-points --full-diag 5000
    OUTPUT "bosons_torus_3body_hardcore_n_6_2s_6_ratio_1.000000.dat")
diagham_geometry_test(physics.fqhe.torus.fermions_mooreread.sixfold_degeneracy CASE torus_fermions_mooreread_n4_nphi8 GEOMETRY torus
    PROGRAM ${tfN} ARGS -p 4 -l 8 --nbr-nbody 3 --all-points --full-diag 5000
    OUTPUT "fermions_torus_3body_hollowcore_n_4_2s_8_ratio_1.000000.dat")

# --- cylinder ----------------------------------------------------------------------------
set(cbd FQHEOnCylinder_FQHECylinderBosonsDeltaInteraction)
set(cb3 FQHEOnCylinder_FQHECylinderBosonsThreeBodyDeltaInteraction)
diagham_geometry_test(physics.fqhe.cylinder.bosons_laughlin_1_2.unique_zero_mode CASE cylinder_bosons_laughlin_1_2_n4_l6 GEOMETRY cylinder
    PROGRAM ${cbd} ARGS -p 4 -l 6 --full-diag 5000 OUTPUT "bosons_cylinder_ky_delta_n_4_2s_6_ratio_1.000000.dat")
diagham_geometry_test(physics.fqhe.cylinder.bosons_laughlin_1_2.unique_zero_mode_n5 CASE cylinder_bosons_laughlin_1_2_n5_l8 GEOMETRY cylinder
    PROGRAM ${cbd} ARGS -p 5 -l 8 --full-diag 5000 OUTPUT "bosons_cylinder_ky_delta_n_5_2s_8_ratio_1.000000.dat")
diagham_geometry_test(physics.fqhe.cylinder.bosons_laughlin_1_2.one_quasihole_count CASE cylinder_bosons_laughlin_1_2_n4_l7 GEOMETRY cylinder
    PROGRAM ${cbd} ARGS -p 4 -l 7 --full-diag 5000 OUTPUT "bosons_cylinder_ky_delta_n_4_2s_7_ratio_1.000000.dat")
diagham_geometry_test(physics.fqhe.cylinder.bosons_laughlin_1_2.two_quasihole_count CASE cylinder_bosons_laughlin_1_2_n4_l8 GEOMETRY cylinder
    PROGRAM ${cbd} ARGS -p 4 -l 8 --full-diag 5000 OUTPUT "bosons_cylinder_ky_delta_n_4_2s_8_ratio_1.000000.dat")
diagham_geometry_test(physics.fqhe.cylinder.bosons_mooreread.unique_zero_mode CASE cylinder_bosons_mooreread_n4_l2 GEOMETRY cylinder
    PROGRAM ${cb3} ARGS -p 4 -l 2 --full-diag 5000 OUTPUT "bosons_cylinder_ky_3b_n_4_2s_2_ratio_1.000000.dat")
diagham_geometry_test(physics.fqhe.cylinder.bosons_mooreread.unique_zero_mode_n6 CASE cylinder_bosons_mooreread_n6_l4 GEOMETRY cylinder
    PROGRAM ${cb3} ARGS -p 6 -l 4 --full-diag 5000 OUTPUT "bosons_cylinder_ky_3b_n_6_2s_4_ratio_1.000000.dat")
diagham_geometry_test(physics.fqhe.cylinder.bosons_mooreread.two_quasihole_count CASE cylinder_bosons_mooreread_n6_l5 GEOMETRY cylinder
    PROGRAM ${cb3} ARGS -p 6 -l 5 --full-diag 5000 OUTPUT "bosons_cylinder_ky_3b_n_6_2s_5_ratio_1.000000.dat")

# --- disk (bosons only: FQHEDiskFermionsTwoBodyGeneric hangs, defect U29) ----------------
set(db2 FQHEOnDisk_FQHEDiskBosonsTwoBodyGeneric)
foreach(lz 12 13 14 15 16)
    diagham_geometry_test(physics.fqhe.disk.bosons_laughlin_1_2.n4_lz${lz} CASE disk_bosons_laughlin_1_2_n4_lz${lz} GEOMETRY disk
        PROGRAM ${db2} ARGS -p 4 --interaction-file ${v0} --interaction-name v0 --minimum-momentum ${lz} --maximum-momentum ${lz} --full-diag 5000
        OUTPUT "bosons_disk_v0_n_4_lz_${lz}.dat")
endforeach()
foreach(lz 20 22)
    diagham_geometry_test(physics.fqhe.disk.bosons_laughlin_1_2.n5_lz${lz} CASE disk_bosons_laughlin_1_2_n5_lz${lz} GEOMETRY disk
        PROGRAM ${db2} ARGS -p 5 --interaction-file ${v0} --interaction-name v0 --minimum-momentum ${lz} --maximum-momentum ${lz} --full-diag 5000
        OUTPUT "bosons_disk_v0_n_5_lz_${lz}.dat")
endforeach()

if(Python3_Interpreter_FOUND)
    add_test(NAME selftest.fqhe_geometries_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/fqhe_geometries.py --check ${geom_data})
    set_tests_properties(selftest.fqhe_geometries_oracle PROPERTIES LABELS "selftest;python" TIMEOUT 600)
endif()
