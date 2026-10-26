# ----------------------------------------------------------------------------
# Core-program gaps (session S4, 25 October 2026): physics tests for programs of the provisional core list
# that had only a start-up check. Each answer is known without DiagHam:
#  - VectorAscii2Binary / VectorBinary2Ascii: the documented binary vector layout, written and read by
#    tests/oracles/vector_io.py with Python's struct module (bit-exact one way, 14 printed digits the other);
#  - LevelStatistics: synthetic spectra whose spacing ratio is exact (1/2 for spacings 1, 2, 1, 2, ...; 1 for
#    equal spacings);
#  - FQHESphereL2Diagonalize: the L^2 spectrum of an Lz sector is L(L+1) once per multiplet, the multiplet
#    counts d(L) - d(L+1) from exact state counting (tests/oracles/sphere_l2.py);
#  - FCIHofstadterModel: against an independent Hofstadter exact diagonalisation (tests/oracles/hofstadter.py):
#    the lowest band of one particle on three lattices at two fluxes, and the full spectra of two and three
#    hard-core bosons in the real-space mode (-Y 1, the only unit-cell height it supports there).
# Two open findings met on the way, both the same on a pristine autotools build of r4493: U44 (--real-space
# with spinless fermions prints "case not implented yet" and crashes) and U45 (--real-space accepts -Y other
# than 1, which its wiki page says does not work, and then gives a wrong spectrum).
# ----------------------------------------------------------------------------
set(core_data ${DIAGHAM_TEST_DATA})

if(Python3_Interpreter_FOUND AND TARGET Programs_VectorAscii2Binary AND TARGET Programs_VectorBinary2Ascii)
    add_test(NAME physics.core.vector_ascii2binary.bit_exact
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/vector_io.py ascii2binary
            $<TARGET_FILE:Programs_VectorAscii2Binary> ${DIAGHAM_TEST_WORK}/physics.core.vector_ascii2binary.bit_exact)
    add_test(NAME physics.core.vector_binary2ascii.printed_digits
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/vector_io.py binary2ascii
            $<TARGET_FILE:Programs_VectorBinary2Ascii> ${DIAGHAM_TEST_WORK}/physics.core.vector_binary2ascii.printed_digits)
    set_tests_properties(physics.core.vector_ascii2binary.bit_exact physics.core.vector_binary2ascii.printed_digits
        PROPERTIES LABELS "physics;python;core" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
endif()

foreach(case "alternating_spacings_r_half|alternating" "uniform_spacings_r_one|uniform")
    string(REPLACE "|" ";" c "${case}")
    list(GET c 0 name)
    list(GET c 1 file)
    diagham_physics_test(physics.core.level_statistics.${name}
        PROGRAM Programs_LevelStatistics
        ARGS -s ${core_data}/levelstat/${file}.dat -o levelstat.out
        OUTPUT "program.log"
        CHECK line @OUTPUT@ <r> ${core_data}/levelstat/${file}_r.txt ${DIAGHAM_TEST_TOLERANCE}
        LABELS core)
endforeach()

foreach(case "fermions_n4_2s9_lz0|fermions|4|9|0" "fermions_n3_2s6_2lz2|fermions|3|6|2"
             "bosons_n4_2s6_lz0|bosons|4|6|0" "bosons_n3_2s5_2lz3|bosons|3|5|3")
    string(REPLACE "|" ";" c "${case}")
    list(GET c 0 name)
    list(GET c 1 stat)
    list(GET c 2 n)
    list(GET c 3 two_s)
    list(GET c 4 two_lz)
    diagham_physics_test(physics.fqhe.sphere.l2.${name}
        PROGRAM FQHEOnSphere_FQHESphereL2Diagonalize
        ARGS -p ${n} -l ${two_s} -z ${two_lz} --statistics ${stat}
        OUTPUT "${stat}_l2_n_${n}_2s_${two_s}_lz.dat"
        CHECK spectrum @OUTPUT@ -1 ${core_data}/sphere_l2/${name}.dat -1 ${DIAGHAM_TEST_TOLERANCE}
        LABELS fqhe sphere)
endforeach()

# one particle, lowest band (the program's default band projection); --redundant-kx/-ky list every momentum
foreach(case "4|2|1|3|1" "3|2|1|4|1" "2|3|1|4|1")
    string(REPLACE "|" ";" c "${case}")
    list(GET c 0 x)
    list(GET c 1 y)
    list(GET c 2 cx)
    list(GET c 3 cy)
    list(GET c 4 q)
    math(EXPR nk "${x} * ${y}")
    diagham_physics_test(physics.fci.hofstadter.one_particle_lowest_band.x${x}_y${y}_X${cx}_Y${cy}_q${q}
        PROGRAM FCI_FCIHofstadterModel
        ARGS -p 1 -x ${x} -y ${y} -X ${cx} -Y ${cy} -q ${q} --redundant-kx --redundant-ky --full-diag 5000
        OUTPUT "fermions_hofstadter_X_${cx}_Y_${cy}_q_${q}_n_1_x_${x}_y_${y}_*.dat"
        CHECK lowest @OUTPUT@ -1 ${core_data}/hofstadter/fermions_n1_x${x}_y${y}_X${cx}_Y${cy}_q${q}.dat -1 ${nk} ${DIAGHAM_TEST_TOLERANCE}
        LABELS fci hofstadter)
endforeach()

# real space, hard-core bosons, -Y 1: the full spectrum (276 and 2024 levels) to the agreed tolerance
foreach(case "2|4|2|3|1|1" "3|3|2|4|1|1")
    string(REPLACE "|" ";" c "${case}")
    list(GET c 0 n)
    list(GET c 1 x)
    list(GET c 2 y)
    list(GET c 3 cx)
    list(GET c 4 cy)
    list(GET c 5 q)
    diagham_physics_test(physics.fci.hofstadter.realspace_hardcore_bosons.n${n}_x${x}_y${y}_X${cx}_Y${cy}_q${q}
        PROGRAM FCI_FCIHofstadterModel
        ARGS --real-space --redundant-kx --redundant-ky -p ${n} -x ${x} -y ${y} -X ${cx} -Y ${cy} -q ${q} --boson --hardcore --full-diag 5000
        OUTPUT "bosons_realspace_gutzwiller_hofstadter_*.dat"
        CHECK spectrum @OUTPUT@ -1 ${core_data}/hofstadter/hardcore_bosons_n${n}_x${x}_y${y}_X${cx}_Y${cy}_q${q}.dat -1 ${DIAGHAM_TEST_TOLERANCE}
        LABELS fci hofstadter)
endforeach()

# U44: --real-space with spinless fermions prints "case not implented yet" and then segfaults (exit 139), with
# or without interaction and whatever -Y; the same on a pristine autotools build of r4493. WILL_FAIL until fixed.
diagham_physics_test(openfinding.fci.hofstadter.realspace_fermions_segfault
    PROGRAM FCI_FCIHofstadterModel
    ARGS --real-space -p 2 -x 4 -y 2 -X 3 -Y 1 -q 1 --full-diag 5000
    OUTPUT "*.dat"
    CHECK count @OUTPUT@ -1 0 1e-10 0
    WILL_FAIL LABELS fci hofstadter)

# U45: the program's wiki page says the code works only with -Y 1; with -Y 3 or 4 --real-space does not refuse
# the input but gives a spectrum whose ground state lies below 2 e_0, the rigorous bound from the program's own
# one-particle spectrum (-6.129 against -5.464 for the first case). These tests carry what it should give;
# the program's output is byte-identical on a pristine autotools build of r4493.
foreach(case "2|4|2|1|3|1" "3|3|2|1|4|1")
    string(REPLACE "|" ";" c "${case}")
    list(GET c 0 n)
    list(GET c 1 x)
    list(GET c 2 y)
    list(GET c 3 cx)
    list(GET c 4 cy)
    list(GET c 5 q)
    diagham_physics_test(openfinding.fci.hofstadter.realspace_unit_cell_y_not_1_accepted.n${n}_x${x}_y${y}
        PROGRAM FCI_FCIHofstadterModel
        ARGS --real-space --redundant-kx --redundant-ky -p ${n} -x ${x} -y ${y} -X ${cx} -Y ${cy} -q ${q} --boson --hardcore --full-diag 5000
        OUTPUT "bosons_realspace_gutzwiller_hofstadter_*.dat"
        CHECK spectrum @OUTPUT@ -1 ${core_data}/hofstadter/hardcore_bosons_n${n}_x${x}_y${y}_X${cx}_Y${cy}_q${q}.dat -1 ${DIAGHAM_TEST_TOLERANCE}
        WILL_FAIL LABELS fci hofstadter)
endforeach()

if(Python3_Interpreter_FOUND)
    add_test(NAME selftest.hofstadter_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/hofstadter.py --check ${core_data}/hofstadter)
    add_test(NAME selftest.sphere_l2_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/sphere_l2.py --check ${core_data}/sphere_l2)
    set_tests_properties(selftest.hofstadter_oracle selftest.sphere_l2_oracle PROPERTIES LABELS "selftest;python" TIMEOUT 300)
endif()

# ----------------------------------------------------------------------------
# Core programs without an independent test (session S9, 26 October 2026): programs that came out core under every
# weighting of the core-program assessment, taken two to four per session. Answers from tests/oracles/core_tools.py:
#  - GenericHamiltonianDiagonalization: matrices given by closed formulas, written with 17 digits, against
#    numpy.linalg.eigvalsh (full diagonalisation real and complex; Lanczos on a 600 x 600 banded matrix whose lowest
#    eigenvalues are well separated);
#  - FTIShowBasis: the n-body basis of a momentum sector against an exact enumeration, state by state, no tolerance;
#  - FQHESphereBosonEntanglementEntropy: entropies and reduced density matrix eigenvalues of every orbital cut of a
#    unique sphere ground state, against a Schmidt decomposition of the ground state of the same Hamiltonian built
#    from scratch (tests/oracles/sphere_ed.py). The ground state is computed by FQHESphereBosonsTwoBodyGeneric.
# Open finding U46, met on the way: DiagHam's built-in full diagonalisation skips Householder steps whose remainder
# is below 1e-7 in norm, so a saved eigenvector can be wrong at about 1e-8 (the bosonic Laughlin state below);
# eigenvalues are affected only at second order. --use-lapack avoids it. The same on a pristine autotools build.
# ----------------------------------------------------------------------------
set(ct_data ${DIAGHAM_TEST_DATA}/core_tools)
set(ghd Programs_GenericHamiltonianDiagonalization)
diagham_physics_test(physics.core.generic_hamiltonian_diagonalization.real_dense_40
    PROGRAM ${ghd}
    ARGS --hamiltonian ${ct_data}/real_dense_40.txt -o real_dense_40.dat
    OUTPUT "real_dense_40.dat"
    CHECK spectrum @OUTPUT@ -1 ${ct_data}/real_dense_40_spectrum.dat -1 ${DIAGHAM_TEST_TOLERANCE}
    LABELS core)
diagham_physics_test(physics.core.generic_hamiltonian_diagonalization.complex_dense_30
    PROGRAM ${ghd}
    ARGS --hamiltonian ${ct_data}/complex_dense_30.txt -c -o complex_dense_30.dat
    OUTPUT "complex_dense_30.dat"
    CHECK spectrum @OUTPUT@ -1 ${ct_data}/complex_dense_30_spectrum.dat -1 ${DIAGHAM_TEST_TOLERANCE}
    LABELS core)
diagham_physics_test(physics.core.generic_hamiltonian_diagonalization.lanczos_sparse_600_lowest_10
    PROGRAM ${ghd}
    ARGS --hamiltonian ${ct_data}/real_sparse_600.txt -n 10 -o real_sparse_600.dat
    OUTPUT "real_sparse_600.dat"
    CHECK lowest @OUTPUT@ -1 ${ct_data}/real_sparse_600_spectrum.dat -1 10 ${DIAGHAM_TEST_TOLERANCE}
    LABELS core)

if(Python3_Interpreter_FOUND AND TARGET FTI_FTIShowBasis)
    foreach(case "fermions_n4_x3_y3_k1_2|-p;4;-x;3;-y;3;--kx;1;--ky;2"
                 "bosons_n3_x3_y2_k1_0|-p;3;-x;3;-y;2;--kx;1;--ky;0;--boson"
                 "bosons_n4_x4_y2_k2_1|-p;4;-x;4;-y;2;--kx;2;--ky;1;--boson"
                 "two_bands_fermions_n3_x2_y2_k0_1|-p;3;-x;2;-y;2;-s;2;--kx;0;--ky;1"
                 "two_bands_bosons_n3_x3_y2_k2_1|-p;3;-x;3;-y;2;-s;2;--kx;2;--ky;1;--boson")
        string(FIND "${case}" "|" bar)
        string(SUBSTRING "${case}" 0 ${bar} name)
        math(EXPR bar "${bar} + 1")
        string(SUBSTRING "${case}" ${bar} -1 args)
        add_test(NAME physics.fti.show_basis.${name}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/core_tools.py showbasis $<TARGET_FILE:FTI_FTIShowBasis> ${args})
        set_tests_properties(physics.fti.show_basis.${name} PROPERTIES LABELS "physics;python;fti" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endforeach()
endif()

if(Python3_Interpreter_FOUND AND TARGET FQHEOnSphere_FQHESphereBosonEntanglementEntropy AND TARGET FQHEOnSphere_FQHESphereBosonsTwoBodyGeneric)
    foreach(case bosons_v0_v2_n5_2s6 bosons_coulomb_n4_2s6)
        add_test(NAME physics.fqhe.sphere.boson_entanglement.${case}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/core_tools.py entanglement
                $<TARGET_FILE:FQHEOnSphere_FQHESphereBosonsTwoBodyGeneric> $<TARGET_FILE:FQHEOnSphere_FQHESphereBosonEntanglementEntropy>
                ${ct_data} ${case} ${DIAGHAM_TEST_WORK}/physics.fqhe.sphere.boson_entanglement.${case} ${DIAGHAM_TEST_TOLERANCE})
        set_tests_properties(physics.fqhe.sphere.boson_entanglement.${case} PROPERTIES LABELS "physics;python;fqhe;sphere" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endforeach()
    # U46: the bosonic Laughlin state (V0 only, N = 4, 2S = 6) saved by the built-in full diagonalisation is off by
    # 1.2e-8 (|Hv - Ev|), so its entropies are off by 3e-10 (the exact l_a = 1 eigenvalues are 3/7 and 4/7). WILL_FAIL until fixed.
    add_test(NAME openfinding.full_diagonalisation.householder_cutoff.laughlin_eigenvector
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/core_tools.py entanglement
            $<TARGET_FILE:FQHEOnSphere_FQHESphereBosonsTwoBodyGeneric> $<TARGET_FILE:FQHEOnSphere_FQHESphereBosonEntanglementEntropy>
            ${ct_data} bosons_laughlin_n4_2s6 ${DIAGHAM_TEST_WORK}/openfinding.full_diagonalisation.householder_cutoff.laughlin_eigenvector
            ${DIAGHAM_TEST_TOLERANCE})
    set_tests_properties(openfinding.full_diagonalisation.householder_cutoff.laughlin_eigenvector PROPERTIES
        WILL_FAIL TRUE LABELS "open-finding;python;fqhe;sphere" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
endif()

if(Python3_Interpreter_FOUND)
    add_test(NAME selftest.core_tools_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/core_tools.py --check ${ct_data})
    set_tests_properties(selftest.core_tools_oracle PROPERTIES LABELS "selftest;python" TIMEOUT 300)
endif()
