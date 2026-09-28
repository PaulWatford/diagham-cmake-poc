# ----------------------------------------------------------------------------
# Spin chains and Hubbard: ground energies against tests/oracles/spin_hubbard.py
# (label: physics, spin | fti hubbard). Values and tolerances come from
# tests/data/spin_hubbard/values.txt ("name value tol # ..."), re-derived by
# selftest.spin_hubbard_oracle (needs numpy for three of them).
# ----------------------------------------------------------------------------
set(sh_data ${DIAGHAM_TEST_DATA}/spin_hubbard)

function(diagham_sh_value case out_value out_tol)
    file(STRINGS ${sh_data}/values.txt lines REGEX "^${case} ")
    list(LENGTH lines n)
    if(NOT n EQUAL 1)
        message(FATAL_ERROR "spin_hubbard/values.txt: case ${case} not found exactly once")
    endif()
    string(REGEX MATCH "^${case} ([-0-9.e+]+) ([-0-9.e+]+)" _ "${lines}")
    set(${out_value} ${CMAKE_MATCH_1} PARENT_SCOPE)
    set(${out_tol} ${CMAKE_MATCH_2} PARENT_SCOPE)
endfunction()

# diagham_ground_energy_test(<name> CASE <case> PROGRAM <target> ARGS <args...> OUTPUT <glob> LABELS <labels...>)
function(diagham_ground_energy_test name)
    cmake_parse_arguments(ARG "" "CASE;PROGRAM;OUTPUT" "ARGS;LABELS" ${ARGN})
    diagham_sh_value(${ARG_CASE} value tol)
    diagham_physics_test(${name}
        PROGRAM ${ARG_PROGRAM}
        ARGS ${ARG_ARGS}
        OUTPUT "${ARG_OUTPUT}"
        CHECK min-abs @OUTPUT@ -1 ${value} ${tol}
        LABELS ${ARG_LABELS})
endfunction()

# --- spin chains ---------------------------------------------------------------
set(gps Spin_GenericPeriodicSpinChain)
diagham_ground_energy_test(physics.spin.heisenberg_ring_L8 CASE heisenberg_ring_L8
    PROGRAM ${gps} ARGS -p 8 --nbr-sz 1 --full-diag 5000 OUTPUT "spin_1_2_periodicchain_n_8.dat" LABELS spin)
diagham_ground_energy_test(physics.spin.spin1_heisenberg_ring_L6 CASE spin1_heisenberg_ring_L6
    PROGRAM ${gps} ARGS -s 2 -p 6 --nbr-sz 1 --full-diag 5000 OUTPUT "spin_1_periodicchain_n_6.dat" LABELS spin)
# --nn-coupling -1 removes the Sz Sz term (the ZZ coupling is j + nn-coupling): the XX chain
foreach(L 6 8 10)
    diagham_ground_energy_test(physics.spin.xx_ring_L${L} CASE xx_ring_L${L}
        PROGRAM ${gps} ARGS -p ${L} --nbr-sz 1 --nn-coupling -1 --full-diag 5000 OUTPUT "spin_1_2_periodicchain_n_${L}.dat" LABELS spin)
endforeach()
foreach(L 6 8)
    diagham_ground_energy_test(physics.spin.aklt_ring_L${L} CASE aklt_ring_L${L}
        PROGRAM Spin_PeriodicSpinChainAKLT ARGS -p ${L} --nbr-sz 1 --full-diag 5000 OUTPUT "spin_1_periodicaklt_n_${L}.dat" LABELS spin)
endforeach()
# open AKLT chain as a sum of projectors: the two spin-1/2 edge spins give exactly four zero-energy
# states (singlet + triplet); the program lists the 2Sz >= 0 sectors, so 3 appear (2 at Sz=0, 1 at Sz=1)
diagham_physics_test(physics.spin.aklt_open_chain_L6.edge_state_count
    PROGRAM Spin_SpinChainAKLT
    ARGS -p 6 --projector-normalization --full-diag 5000
    OUTPUT "spin_1_open_projnormaklt_n_6.dat"
    CHECK count @OUTPUT@ -1 0 1e-10 3
    LABELS spin)
foreach(L 6 8 10)
    diagham_ground_energy_test(physics.spin.haldane_shastry_L${L} CASE haldane_shastry_L${L}
        PROGRAM Spin_HaldaneShastrySpinChain ARGS -p ${L} --nbr-sz 1 --full-diag 5000 OUTPUT "spin_1_2_haldaneshastrychain_n_${L}.dat" LABELS spin)
endforeach()
diagham_ground_energy_test(physics.spin.tfim_ring_L8 CASE tfim_ring_L8_h0.5
    PROGRAM Spin_SpinChainXYZ ARGS -p 8 -x 1 -y 0 -z 0 -f 0.5 -b 1 --full-diag 5000
    OUTPUT "spin_1_2_x_1.000000_y_0.000000_z_0.000000_h_0.500000_b_1_n_8.dat" LABELS spin)

# --- Hubbard at U = 0: the tight-binding sum, by full diagonalisation ---------------------
foreach(case "2 4 8" "4 2 8" "3 2 6")
    separate_arguments(c UNIX_COMMAND "${case}")
    list(GET c 0 lx)
    list(GET c 1 ly)
    list(GET c 2 n)
    diagham_ground_energy_test(physics.hubbard.${lx}x${ly}.U0.tight_binding CASE hubbard_${lx}x${ly}_N${n}_U0
        PROGRAM ${hubbard} ARGS -p ${n} -x ${lx} -y ${ly} --u-potential 0 --full-diag 5000
        OUTPUT "fermions_hubbard_square_x_${lx}_y_${ly}_n_${n}_ns_*_t_1.000000_tp_0.000000_sz_0.dat" LABELS fti hubbard)
endforeach()

# --- U30, fixed here: plain Lanczos (-n 1, no reorthogonalisation) on the U = 0 2x4 case used
# to return -18.41 and -17.87 in two momentum sectors, BELOW the exact ground energy -16. On this
# degenerate spectrum the Krylov space of the start vector closes after 7 to 9 steps; the
# recurrence went on with the normalised round-off and produced Ritz values outside the spectrum.
# BasicLanczosAlgorithm and ComplexBasicLanczosAlgorithm now detect the closure and stop
# (docs/upstream-reports/, "Lanczos Krylov space closure"). The ground energy of every one of
# the eight momentum sectors must equal the tight-binding minimum of that sector, from
# tests/oracles/spin_hubbard.py (was the known-bug test knownbug.hubbard.2x4.U0.plain_lanczos_below_ground_state).
diagham_physics_test(physics.hubbard.2x4.U0.plain_lanczos_krylov_closure
    PROGRAM ${hubbard}
    ARGS -p 8 -x 2 -y 4 --u-potential 0 -n 1 --full-diag 100
    OUTPUT "fermions_hubbard_square_x_2_y_4_n_8_ns_8_t_1.000000_tp_0.000000_sz_0.dat"
    CHECK spectrum @OUTPUT@ -1 ${sh_data}/hubbard_2x4_N8_U0_sectors.dat -1 1e-8
    LABELS fti hubbard lanczos)

if(Python3_Interpreter_FOUND AND numpy_missing EQUAL 0)
    add_test(NAME selftest.spin_hubbard_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/spin_hubbard.py --check ${sh_data})
    set_tests_properties(selftest.spin_hubbard_oracle PROPERTIES LABELS "selftest;python" TIMEOUT 600)
endif()
