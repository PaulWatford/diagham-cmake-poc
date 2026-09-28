# ----------------------------------------------------------------------------
# Independent oracles, algorithm consensus, and self-tests of the harness.
#
#  1. tests/oracles/sphere_ed.py: a from-scratch exact diagonalisation of
#     pseudopotential Hamiltonians in the lowest Landau level (Clebsch-Gordan
#     pair projectors, explicit Fock basis, numpy). The program's whole
#     spectrum of an Lz sector must match it (label: physics, sphere, ed).
#  1b. tests/oracles/torus_ed.py: the same on the torus (Landau gauge, Coulomb
#     or pseudopotentials, fermions or bosons; label: physics, fqhe, torus, ed).
#  1c. tests/oracles/fci_bands.py: the checkerboard lattice bands from the
#     published Bloch Hamiltonian (label: physics, fci, checkerboard, bands).
#  2. tests/oracles/spin_ed.py: dense numpy diagonalisation of XXZ / J1-J2
#     chains with a field, open or periodic (label: physics, spin, ed).
#  3. crosscheck: two algorithms on the same Hamiltonian must agree --
#     Lanczos (-n 4, reorthogonalised) against full diagonalisation, and
#     LAPACK against DiagHam's own diagonaliser when the build has LAPACK.
#  4. selftest: the checker and the runner must reject wrong values, missing
#     files, NaN and missing output (WILL_FAIL tests).
# ----------------------------------------------------------------------------
set(ed_data ${DIAGHAM_TEST_DATA}/sphere_ed)
set(spin_ed_data ${DIAGHAM_TEST_DATA}/spin_ed)
set(selftest_data ${DIAGHAM_TEST_DATA}/selftest)

# --- 1. sphere pseudopotential ED -------------------------------------------------------
# name, program, N, 2S, 2Lz (the sector), file glob; the pseudopotential file is <name>_pp.dat
foreach(case
        "fermions_n4_2s9_coulomb|${fer2}|4|9|0"
        "fermions_n5_2s12_coulomb|${fer2}|5|12|0"
        "fermions_n4_2s9_v1_v3|${fer2}|4|9|0"
        "fermions_n4_2s11_coulomb_lz2|${fer2}|4|11|4"
        "bosons_n4_2s6_coulomb|${bos2}|4|6|0"
        "bosons_n5_2s6_v0_v2|${bos2}|5|6|0"
        "bosons_n6_2s8_coulomb|${bos2}|6|8|0")
    string(REPLACE "|" ";" c "${case}")
    list(GET c 0 name)
    list(GET c 1 prog)
    list(GET c 2 n)
    list(GET c 3 two_s)
    list(GET c 4 two_lz)
    if(prog MATCHES "Bosons")
        set(stat bosons)
    else()
        set(stat fermions)
    endif()
    diagham_physics_test(physics.fqhe.sphere.ed.${name}
        PROGRAM ${prog}
        ARGS -p ${n} -l ${two_s} --interaction-file ${ed_data}/${name}_pp.dat --interaction-name ed --initial-lz ${two_lz} --nbr-lz 1 --full-diag 100000
        OUTPUT "${stat}_ed_n_${n}_2s_${two_s}_lz.dat"
        CHECK spectrum @OUTPUT@ -1 ${ed_data}/${name}_spectrum.dat -1 1e-9
        LABELS fqhe sphere ed)
endforeach()

# --- 1b. torus ED -------------------------------------------------------------------------------
# tests/oracles/torus_ed.py: Landau-gauge lowest-Landau-level exact diagonalisation on the torus,
# Coulomb or Haldane pseudopotentials, fermions or bosons, one Ky sector; the program's whole
# spectrum of that sector must match. Tolerance 1e-8 rather than 1e-9: the oracle's momentum sums
# are converged to 1e-14 but DiagHam truncates its own at about 1e-10 (its degenerate pairs differ
# in the tenth digit).
# name | program | N | Nphi | ratio | Ky | interaction (coulomb, or the pseudopotential file stem)
set(torus_ed_data ${DIAGHAM_TEST_DATA}/torus_ed)
set(tfc FQHEOnTorus_FQHETorusFermionsCoulomb)
set(tbc FQHEOnTorus_FQHETorusBosonsCoulomb)
set(tfg FQHEOnTorus_FQHETorusFermionsTwoBodyGeneric)
set(tbg FQHEOnTorus_FQHETorusBosonsTwoBodyGeneric)
foreach(case
        "fermions_n3_nphi9_coulomb_ky0|${tfc}|3|9|1.0|0|coulomb"
        "fermions_n4_nphi8_coulomb_ky0|${tfc}|4|8|1.0|0|coulomb"
        "fermions_n3_nphi9_coulomb_r2_ky0|${tfc}|3|9|2.0|0|coulomb"
        "fermions_n4_nphi12_coulomb_ky0|${tfc}|4|12|1.0|0|coulomb"
        "bosons_n3_nphi6_coulomb_ky0|${tbc}|3|6|1.0|0|coulomb"
        "fermions_n3_nphi9_v1_ky0|${tfg}|3|9|1.0|0|v1"
        "fermions_n4_nphi12_v1v3_ky0|${tfg}|4|12|1.0|0|v1v3"
        "bosons_n4_nphi8_v0v2_ky0|${tbg}|4|8|1.0|0|v0v2")
    string(REPLACE "|" ";" c "${case}")
    list(GET c 0 name)
    list(GET c 1 prog)
    list(GET c 2 n)
    list(GET c 3 nphi)
    list(GET c 4 ratio)
    list(GET c 5 ky)
    list(GET c 6 inter)
    if(prog MATCHES "Bosons")
        set(stat bosons)
    else()
        set(stat fermions)
    endif()
    if(inter STREQUAL "coulomb")
        set(extra "")
    else()
        set(extra --interaction-file ${torus_ed_data}/${name}_pp.dat --interaction-name ${inter})
    endif()
    diagham_physics_test(physics.fqhe.torus.ed.${name}
        PROGRAM ${prog}
        ARGS -p ${n} -l ${nphi} --ratio ${ratio} --ky-momentum ${ky} --full-diag 100000 ${extra}
        OUTPUT "${stat}_torus_kysym_${inter}_n_${n}_2s_${nphi}_ratio_*.dat"
        CHECK spectrum @OUTPUT@ -1 ${torus_ed_data}/${name}_spectrum.dat -1 1e-8
        LABELS fqhe torus ed)
endforeach()

# --- 1c. FCI band structures ---------------------------------------------------------------------
# tests/oracles/fci_bands.py: the two bands of the checkerboard lattice model from the published
# Bloch Hamiltonian (Sun, Gu, Katsura, Das Sarma 2011) at every lattice momentum. The program's
# --export-onebodytext file (kx ky E_0 E_1 ...) must match column by column; one test per band.
# name | Nx | Ny | extra program options
set(fci_bands_data ${DIAGHAM_TEST_DATA}/fci_bands)
set(fci_cb FCI_FCICheckerboardLatticeModel)
foreach(case
        "checkerboard_3x3|3|3|"
        "checkerboard_4x4|4|4|"
        "checkerboard_4x3_t2_0.2_tpp_0.1|4|3|--t2;0.2;--tpp;0.1")
    string(REPLACE "|" ";" c "${case}")
    list(GET c 0 name)
    list(GET c 1 nx)
    list(GET c 2 ny)
    list(LENGTH c len)
    set(extra "")
    if(len GREATER 3)
        list(SUBLIST c 3 -1 extra)
    endif()
    foreach(band 0 1)
        math(EXPR col "${band} + 2")
        diagham_physics_test(physics.fci.checkerboard.bands.${name}.band${band}
            PROGRAM ${fci_cb}
            ARGS -p 2 -x ${nx} -y ${ny} --singleparticle-spectrum --export-onebodytext ${extra}
            OUTPUT "fermions_checkerboardlattice_n_2_ns_*_x_${nx}_y_${ny}_tightbinding.dat"
            CHECK spectrum @OUTPUT@ ${col} ${fci_bands_data}/${name}.dat ${col} 1e-12
            LABELS fci checkerboard bands)
    endforeach()
endforeach()

# --- 2. spin chain ED ------------------------------------------------------------------------
function(diagham_spin_ed_test name case)
    cmake_parse_arguments(ARG "" "PROGRAM;OUTPUT" "ARGS" ${ARGN})
    file(STRINGS ${spin_ed_data}/values.txt lines REGEX "^${case} ")
    list(LENGTH lines n)
    if(NOT n EQUAL 1)
        message(FATAL_ERROR "spin_ed/values.txt: case ${case} not found exactly once")
    endif()
    string(REGEX MATCH "^${case} ([-0-9.e+]+) ([-0-9.e+]+)" _ "${lines}")
    diagham_physics_test(${name}
        PROGRAM ${ARG_PROGRAM} ARGS ${ARG_ARGS} OUTPUT "${ARG_OUTPUT}"
        CHECK min-abs @OUTPUT@ -1 ${CMAKE_MATCH_1} ${CMAKE_MATCH_2}
        LABELS spin ed)
endfunction()
set(gos Spin_GenericOpenSpinChain)
diagham_spin_ed_test(physics.spin.ed.open_heisenberg_L8 open_heisenberg_L8
    PROGRAM ${gos} ARGS -p 8 --nbr-sz 1 --full-diag 5000 OUTPUT "spin_1_2_openchain_n_8_j_1.000000.dat")
diagham_spin_ed_test(physics.spin.ed.open_spin1_heisenberg_L6 open_spin1_heisenberg_L6
    PROGRAM ${gos} ARGS -s 2 -p 6 --nbr-sz 1 --full-diag 5000 OUTPUT "spin_1_openchain_n_6_j_1.000000.dat")
diagham_spin_ed_test(physics.spin.ed.open_xxz_L8_field open_xxz_L8_jz1.5_hz0.3
    PROGRAM ${gos} ARGS -p 8 --nbr-sz 1 -z 0.5 --hz-value 0.3 --full-diag 5000 OUTPUT "spin_1_2_openchain_n_8_*.dat")
diagham_spin_ed_test(physics.spin.ed.ring_xxz_L8 ring_xxz_L8_jz1.5
    PROGRAM ${gps} ARGS -p 8 --nbr-sz 1 --nn-coupling 0.5 --full-diag 5000 OUTPUT "spin_1_2_periodicchain_n_8.dat")
diagham_spin_ed_test(physics.spin.ed.ring_j1j2_L8 ring_j1j2_L8_j2_0.4
    PROGRAM ${gps} ARGS -p 8 --nbr-sz 1 --nnn-coupling 0.4 --full-diag 5000 OUTPUT "spin_1_2_periodicchain_n_8.dat")
diagham_spin_ed_test(physics.spin.ed.ring_xxz_L10 ring_xxz_L10_jz0.5
    PROGRAM ${gps} ARGS -p 10 --nbr-sz 1 --nn-coupling -0.5 --full-diag 5000 OUTPUT "spin_1_2_periodicchain_n_10.dat")

# --- 3. algorithm consensus (crosscheck) ------------------------------------------------------
# Lanczos with reorthogonalisation (-n 4) vs full diagonalisation: same program, same
# Hamiltonian (Coulomb, N=6, 2S=15, Lz=0, 338 states), the 4 lowest eigenvalues to 1e-10.
diagham_physics_test(crosscheck.fqhe.sphere.lanczos_vs_full_diagonalisation
    PROGRAM ${fer2}
    ARGS -p 6 -l 15 --interaction-file ${ed_data}/fermions_n5_2s12_coulomb_pp.dat --interaction-name c15 --nbr-lz 1 --full-diag 100 -n 4
    OUTPUT "fermions_c15_n_6_2s_15_lz.dat"
    REF_PROGRAM ${fer2}
    REF_ARGS -p 6 -l 15 --interaction-file ${ed_data}/fermions_n5_2s12_coulomb_pp.dat --interaction-name c15 --nbr-lz 1 --full-diag 100000
    REF_OUTPUT "fermions_c15_n_6_2s_15_lz.dat"
    CHECK lowest @OUTPUT@ -1 @REFERENCE@ -1 4 1e-10
    LABELS fqhe sphere)
# (the 2S=12 pseudopotential file has 13 entries, enough for the 2S=15 run's odd m up to 11;
#  this is a consensus check, the interaction only needs to be the same on both sides)
if(DIAGHAM_USE_LAPACK)
    diagham_physics_test(crosscheck.fqhe.sphere.lapack_vs_internal_diagonaliser
        PROGRAM ${fer2}
        ARGS -p 6 -l 15 --interaction-file ${ed_data}/fermions_n5_2s12_coulomb_pp.dat --interaction-name c15 --nbr-lz 1 --full-diag 100000 --use-lapack
        OUTPUT "fermions_c15_n_6_2s_15_lz.dat"
        REF_PROGRAM ${fer2}
        REF_ARGS -p 6 -l 15 --interaction-file ${ed_data}/fermions_n5_2s12_coulomb_pp.dat --interaction-name c15 --nbr-lz 1 --full-diag 100000
        REF_OUTPUT "fermions_c15_n_6_2s_15_lz.dat"
        CHECK spectrum @OUTPUT@ -1 @REFERENCE@ -1 1e-10
        LABELS fqhe sphere)
endif()

# --- 4. self-tests: the harness must fail when it should ------------------------------------------
add_test(NAME selftest.checker.accepts_correct_value
    COMMAND check_spectrum min-abs ${selftest_data}/spectrum.dat -1 -1.5 1e-12)
add_test(NAME selftest.checker.rejects_wrong_value
    COMMAND check_spectrum min-abs ${selftest_data}/spectrum.dat -1 -1.4 1e-6)
add_test(NAME selftest.checker.rejects_wrong_count
    COMMAND check_spectrum count ${selftest_data}/spectrum.dat -1 -1.5 1e-10 2)
add_test(NAME selftest.checker.rejects_missing_file
    COMMAND check_spectrum min-abs ${selftest_data}/does_not_exist.dat -1 -1.5 1e-12)
add_test(NAME selftest.checker.rejects_nan
    COMMAND check_spectrum min-abs ${selftest_data}/spectrum_nan.dat -1 -1.5 1e-12)
add_test(NAME selftest.checker.rejects_short_spectrum
    COMMAND check_spectrum spectrum ${selftest_data}/spectrum.dat -1 ${selftest_data}/spectrum_short.dat -1 1e-12)
set_tests_properties(selftest.checker.rejects_wrong_value selftest.checker.rejects_wrong_count
    selftest.checker.rejects_missing_file selftest.checker.rejects_nan selftest.checker.rejects_short_spectrum
    PROPERTIES WILL_FAIL TRUE)
set_tests_properties(selftest.checker.accepts_correct_value selftest.checker.rejects_wrong_value
    selftest.checker.rejects_wrong_count selftest.checker.rejects_missing_file selftest.checker.rejects_nan
    selftest.checker.rejects_short_spectrum PROPERTIES LABELS "selftest" TIMEOUT 30)
# the runner: a program that writes no spectrum file must fail the test (stale output can never pass)
diagham_physics_test(selftest.runner.rejects_missing_output
    PROGRAM ${fer2} ARGS --help OUTPUT "no_such_output_*.dat"
    CHECK min-abs @OUTPUT@ -1 0 1
    WILL_FAIL LABELS selftest)
if(TEST selftest.runner.rejects_missing_output)
    set_tests_properties(selftest.runner.rejects_missing_output PROPERTIES LABELS "selftest")
endif()

# --- oracle self-tests ------------------------------------------------------------------------------
if(Python3_Interpreter_FOUND AND numpy_missing EQUAL 0)
    add_test(NAME selftest.sphere_ed_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/sphere_ed.py --check ${ed_data})
    add_test(NAME selftest.spin_ed_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/spin_ed.py --check ${spin_ed_data})
    add_test(NAME selftest.torus_ed_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/torus_ed.py --check ${torus_ed_data})
    add_test(NAME selftest.fci_bands_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/fci_bands.py --check ${fci_bands_data})
    set_tests_properties(selftest.sphere_ed_oracle selftest.spin_ed_oracle selftest.torus_ed_oracle
        selftest.fci_bands_oracle PROPERTIES LABELS "selftest;python" TIMEOUT 900)
endif()
