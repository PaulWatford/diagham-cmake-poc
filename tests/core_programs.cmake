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

# ----------------------------------------------------------------------------
# Session S9, second part (26 October 2026): three more programs core under every weighting.
#  - FQHESphereUnnormalizeState: the bosonic Laughlin state prod (z_i - z_j)^2, whose monomial coefficients are
#    exact integers: normalised by the program it must be an exact zero mode of V0 (tests/oracles/sphere_ed.py),
#    and the numpy zero mode, unnormalised by the program, must be proportional to the integers;
#  - FQHESphereWithSU2SpinRotation: a fully polarised state rotated about y against its exact image (the program
#    rotates by exp(+i t S_y); that sense is the only convention taken from it);
#  - HubbardKitaevHeisenbergModel on the wiki's single hexagon: the Hubbard ring (t, U) and the half-filled
#    Gutzwiller spin model (J1, J2, which multiply Pauli matrices, U47). The built-in full diagonalisation moves
#    degenerate eigenvalues of these spectra by up to 3.5e-9 (U46), so the physics tests run with --use-lapack in
#    builds with LAPACK and are open-finding tests without it.
# The first two run on vectors written by tests/oracles/core_tools.py in DiagHam's basis order.
# ----------------------------------------------------------------------------
set(fsb FQHEOnSphere_FQHESphereShowBasis)
if(Python3_Interpreter_FOUND AND TARGET ${fsb} AND TARGET FQHEOnSphere_FQHESphereUnnormalizeState)
    foreach(n 3 4 5)
        add_test(NAME physics.fqhe.sphere.unnormalize_state.laughlin_bosons_n${n}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/core_tools.py unnormalize
                $<TARGET_FILE:${fsb}> $<TARGET_FILE:FQHEOnSphere_FQHESphereUnnormalizeState> ${n}
                ${DIAGHAM_TEST_WORK}/physics.fqhe.sphere.unnormalize_state.laughlin_bosons_n${n} ${DIAGHAM_TEST_TOLERANCE})
        set_tests_properties(physics.fqhe.sphere.unnormalize_state.laughlin_bosons_n${n} PROPERTIES
            LABELS "physics;python;fqhe;sphere" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endforeach()
endif()
if(Python3_Interpreter_FOUND AND TARGET ${fsb} AND TARGET FQHEOnSphere_FQHESphereWithSU2SpinRotation)
    foreach(case "n3_2s4_angle0.5|3|4|0.5" "n3_2s4_angle0.3|3|4|0.3" "n4_2s6_angle0.25|4|6|0.25")
        string(REPLACE "|" ";" c "${case}")
        list(GET c 0 name)
        list(GET c 1 n)
        list(GET c 2 two_s)
        list(GET c 3 angle)
        add_test(NAME physics.fqhe.sphere.su2_spin_rotation.polarized_bosons_${name}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/core_tools.py su2rotation
                $<TARGET_FILE:${fsb}> $<TARGET_FILE:FQHEOnSphere_FQHESphereWithSU2SpinRotation> ${n} ${two_s} ${angle}
                ${DIAGHAM_TEST_WORK}/physics.fqhe.sphere.su2_spin_rotation.polarized_bosons_${name} ${DIAGHAM_TEST_TOLERANCE})
        set_tests_properties(physics.fqhe.sphere.su2_spin_rotation.polarized_bosons_${name} PROPERTIES
            LABELS "physics;python;fqhe;sphere" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endforeach()
endif()

set(kh HubbardModels_HubbardKitaevHeisenbergModel)
foreach(case "hubbard_ring_n4_u2|hex|-p;4;--isotropic-t;1;--anisotropic-t;0;--j1;0;--j2;0;--u-potential;2|fermions_kitaev_heisenberg_hex_n_4_*.dat"
             "spin_j1_0.3_j2_0.5|hexg|-p;6;--gutzwiller;--isotropic-t;0;--anisotropic-t;0;--j1;0.3;--j2;0.5|fermions_kitaev_heisenberg_gutzwiller_hexg_n_6_*.dat")
    string(REPLACE "|" ";" c "${case}")
    list(GET c 0 name)
    list(GET c 1 geo)
    list(LENGTH c len)
    math(EXPR last "${len} - 1")
    list(GET c ${last} output)
    list(SUBLIST c 2 -1 args)
    list(REMOVE_AT args -1)
    set(common ${args} --geometry-file ${ct_data}/kh_hexagon.dat --geometry-name ${geo} --full-diag 5000)
    if(DIAGHAM_USE_LAPACK)
        diagham_physics_test(physics.hubbard.kitaev_heisenberg.hexagon_${name}
            PROGRAM ${kh}
            ARGS ${common} --use-lapack
            OUTPUT "${output}"
            CHECK spectrum @OUTPUT@ -1 ${ct_data}/kh_${name}_spectrum.dat -1 ${DIAGHAM_TEST_TOLERANCE}
            LABELS fti hubbard)
    else()
        # U46: the built-in full diagonalisation moves degenerate eigenvalues of this spectrum by up to 3.5e-9
        diagham_physics_test(openfinding.full_diagonalisation.householder_cutoff.kitaev_heisenberg_hexagon_${name}
            PROGRAM ${kh}
            ARGS ${common}
            OUTPUT "${output}"
            CHECK spectrum @OUTPUT@ -1 ${ct_data}/kh_${name}_spectrum.dat -1 ${DIAGHAM_TEST_TOLERANCE}
            WILL_FAIL LABELS fti hubbard open-finding)
    endif()
endforeach()

# ----------------------------------------------------------------------------
# Overnight run (26 October 2026): FTIGenericInteractionFromFileTwoBands, core under every weighting. Generic one-body
# energies and a generic momentum-conserving interaction written to its input files, against an exact diagonalisation
# of the same Hamiltonian, sector by sector (tests/oracles/fti_generic.py). With --3band-convention the band labels
# mean the same in both files; in the default convention the interaction's band labels act on the opposite bands
# (U48), and the default test pins that behaviour.
# ----------------------------------------------------------------------------
set(fg_data ${DIAGHAM_TEST_DATA}/fti_generic)
if(Python3_Interpreter_FOUND AND TARGET FTI_FTIGenericInteractionFromFileTwoBands)
    foreach(case "fermions_x3_y2_n2|3band" "fermions_x3_y2_n3|3band" "fermions_x2_y2_n3|3band" "fermions_x3_y2_n2|default")
        string(REPLACE "|" ";" c "${case}")
        list(GET c 0 name)
        list(GET c 1 conv)
        set(tname physics.fti.generic_interaction_two_bands.${name}_${conv}_convention)
        add_test(NAME ${tname}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/fti_generic.py run
                $<TARGET_FILE:FTI_FTIGenericInteractionFromFileTwoBands> ${name} ${conv} ${fg_data} ${DIAGHAM_TEST_WORK}/${tname}
                ${DIAGHAM_TEST_TOLERANCE})
        set_tests_properties(${tname} PROPERTIES LABELS "physics;python;fti" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endforeach()
endif()
if(Python3_Interpreter_FOUND)
    add_test(NAME selftest.fti_generic_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/fti_generic.py --check ${fg_data})
    set_tests_properties(selftest.fti_generic_oracle PROPERTIES LABELS "selftest;python" TIMEOUT 300)
endif()

# Overnight run, program 2: FTIEntanglementEntropyParticlePartition (core under every weighting). Complex states with
# generic amplitudes, written in FTIShowBasis order, against the particle-partition reduced density matrix built
# from scratch (tests/oracles/core_tools.py); one case is a three-state --degenerated-groundstate mixture.
if(Python3_Interpreter_FOUND AND TARGET FTI_FTIShowBasis AND TARGET FTI_FTIEntanglementEntropyParticlePartition)
    foreach(case "n3_x3_y3_k1_2|3;3;3;1,2" "n4_x3_y3_k0_1|4;3;3;0,1" "n4_x4_y2_k1_1|4;4;2;1,1" "n3_x3_y3_three_sector_mixture|3;3;3;1,2;2,1;0,0")
        string(FIND "${case}" "|" bar)
        string(SUBSTRING "${case}" 0 ${bar} name)
        math(EXPR bar "${bar} + 1")
        string(SUBSTRING "${case}" ${bar} -1 args)
        set(tname physics.fti.particle_entanglement.fermions_${name})
        add_test(NAME ${tname}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/core_tools.py ftientanglement
                $<TARGET_FILE:FTI_FTIShowBasis> $<TARGET_FILE:FTI_FTIEntanglementEntropyParticlePartition>
                ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE} ${args})
        set_tests_properties(${tname} PROPERTIES LABELS "physics;python;fti;entanglement" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endforeach()
endif()

# Overnight run, program 3: FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing (core under every weighting),
# and FQHESphereQuasiholeMatrixElements, whose output it reads. tests/oracles/quasiholes_trs.py builds the Laughlin 1/3
# quasihole spaces as V1 zero modes and compares basis-free invariants of the matrix elements, then every Lz sector of
# the two-layer problem (potentials, pairing, charging energy). U49: with the default --memory the pairing program
# crashes after the first sector that uses fast multiplication, so the physics tests run with --memory 0.
set(qh_prog FQHEOnSphere_FQHESphereQuasiholeMatrixElements)
set(qh_pair FQHEOnSphere_FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing)
if(Python3_Interpreter_FOUND AND TARGET ${qh_prog})
    foreach(nphi 6 8)
        set(tname physics.fqhe.sphere.quasihole_matrix_elements.laughlin_1_3_nphi_${nphi})
        add_test(NAME ${tname}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/quasiholes_trs.py matrixelements
                $<TARGET_FILE:${qh_prog}> ${nphi} ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE})
        set_tests_properties(${tname} PROPERTIES LABELS "physics;python;fqhe;sphere" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endforeach()
    if(TARGET ${qh_pair})
        foreach(case pairing_2s6_sz0 pairing_2s5_sz1 fixed_n2_2s6_sz0)
            set(tname physics.fqhe.sphere.quasiholes_time_reversal_pairing.${case})
            add_test(NAME ${tname}
                COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/quasiholes_trs.py pairing
                    $<TARGET_FILE:${qh_prog}> $<TARGET_FILE:${qh_pair}> ${case} ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE})
            set_tests_properties(${tname} PROPERTIES LABELS "physics;python;fqhe;sphere" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
        endforeach()
        set(tname openfinding.fqhe.sphere.quasiholes_time_reversal_pairing.default_memory_crash)
        add_test(NAME ${tname}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/quasiholes_trs.py pairing
                $<TARGET_FILE:${qh_prog}> $<TARGET_FILE:${qh_pair}> pairing_2s6_sz0_default_memory ${DIAGHAM_TEST_WORK}/${tname}
                ${DIAGHAM_TEST_TOLERANCE})
        set_tests_properties(${tname} PROPERTIES WILL_FAIL TRUE LABELS "open-finding;python;fqhe;sphere" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endif()
endif()

# Overnight run, program 4: FQHESphereCP2GenericInteraction (core). tests/oracles/cp2.py writes the SU(3) Clebsch-Gordan
# files the program needs (from the two-particle Casimir; the wiki points to a download instead) and diagonalises
# H = sum_k V_k sum_pairs P_k in every sector the program writes. Fermions and bosons, N_phi = 2 and 3.
if(Python3_Interpreter_FOUND AND TARGET FQHEOnSphere_FQHESphereCP2GenericInteraction)
    foreach(case fermions_n3_nphi2 fermions_n3_nphi3 bosons_n3_nphi3 bosons_n4_nphi2)
        set(tname physics.fqhe.cp2.generic_interaction.${case})
        add_test(NAME ${tname}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/cp2.py run
                $<TARGET_FILE:FQHEOnSphere_FQHESphereCP2GenericInteraction> ${case} ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE})
        set_tests_properties(${tname} PROPERTIES LABELS "physics;python;fqhe;cp2" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endforeach()
endif()

# Overnight run, program 5: the real-space partition coefficients (review tier), against closed forms
# (tests/oracles/realspace.py): incomplete beta for a sphere cap, erfc for a cylinder cut, incomplete gamma for a disk.
# The sphere and disk programs need GSL; without it they write all-zero weights and exit 0 (U50), so they are physics
# tests in builds with GSL and open-finding tests without.
set(rs ${CMAKE_CURRENT_SOURCE_DIR}/oracles/realspace.py)
if(Python3_Interpreter_FOUND)
    set(rs_cases)
    if(TARGET FQHEOnCylinder_FQHECylinderRealSpacePartitionCoefficients)
        list(APPEND rs_cases "cylinder|FQHEOnCylinder_FQHECylinderRealSpacePartitionCoefficients|cylinder_2s17_perimeter10,cylinder_2s12_perimeter8")
    endif()
    if(TARGET FQHEOnSphere_FQHESphereRealSpacePartitionCoefficients)
        list(APPEND rs_cases "sphere|FQHEOnSphere_FQHESphereRealSpacePartitionCoefficients|sphere_2s15_hemisphere,sphere_2s10_theta0.3,sphere_2s7_theta0.8")
    endif()
    if(TARGET FQHEOnDisk_FQHEDiskRealSpacePartitionCoefficients)
        list(APPEND rs_cases "disk|FQHEOnDisk_FQHEDiskRealSpacePartitionCoefficients|disk_2s27_radius5,disk_2s10_radius2.5")
    endif()
    foreach(entry ${rs_cases})
        string(REPLACE "|" ";" e "${entry}")
        list(GET e 0 geo)
        list(GET e 1 target)
        list(GET e 2 cases)
        string(REPLACE "," ";" cases "${cases}")
        foreach(case ${cases})
            if(geo STREQUAL "cylinder" OR DIAGHAM_USE_GSL)
                set(tname physics.fqhe.${geo}.realspace_partition_coefficients.${case})
                set(wf FALSE)
                set(lab "physics;python;fqhe;${geo}")
            else()
                if(NOT case MATCHES "^(sphere_2s15_hemisphere|disk_2s10_radius2.5)$")
                    continue()
                endif()
                set(tname openfinding.fqhe.${geo}.realspace_partition_coefficients_zero_without_gsl.${case})
                set(wf TRUE)
                set(lab "open-finding;python;fqhe;${geo}")
            endif()
            add_test(NAME ${tname} COMMAND ${Python3_EXECUTABLE} ${rs} run $<TARGET_FILE:${target}> ${case} ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE})
            set_tests_properties(${tname} PROPERTIES WILL_FAIL ${wf} LABELS "${lab}" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
        endforeach()
    endforeach()
endif()

# Overnight run, program 6: MatrixExtractColumns, SpinChainComputeS2 and SpinChainMultipleEntanglementSpectra (review and
# optional tiers). Inputs written by tests/oracles/core_tools.py and spin_tools.py in DiagHam's layouts and basis order.
# U51: without LAPACK SpinChainMultipleEntanglementSpectra crashes (SVD returns a null pointer); U52: SpinChainComputeS2
# --degenerate-states uses only the first state when the file name has no momentum.
if(Python3_Interpreter_FOUND AND TARGET Programs_MatrixExtractColumns)
    foreach(case "real_7x5_columns_1_to_3|real;7;5;1;3" "complex_6x4_all_columns|complex;6;4;0;4")
        string(FIND "${case}" "|" bar)
        string(SUBSTRING "${case}" 0 ${bar} name)
        math(EXPR bar "${bar} + 1")
        string(SUBSTRING "${case}" ${bar} -1 args)
        set(tname physics.core.matrix_extract_columns.${name})
        add_test(NAME ${tname} COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/core_tools.py matrixcolumns
            $<TARGET_FILE:Programs_MatrixExtractColumns> ${DIAGHAM_TEST_WORK}/${tname} ${args})
        set_tests_properties(${tname} PROPERTIES LABELS "physics;python;core" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endforeach()
endif()
set(st ${CMAKE_CURRENT_SOURCE_DIR}/oracles/spin_tools.py)
if(Python3_Interpreter_FOUND AND TARGET Spin_GenericSpinChainShowBasis AND TARGET Spin_SpinChainComputeS2)
    foreach(case "spin_1_2_n8_sz0|8;1;0" "spin_1_2_n7_sz1|7;1;1" "spin_1_n6_sz2|6;2;2")
        string(FIND "${case}" "|" bar)
        string(SUBSTRING "${case}" 0 ${bar} name)
        math(EXPR bar "${bar} + 1")
        string(SUBSTRING "${case}" ${bar} -1 args)
        set(tname physics.spin.compute_s2.${name})
        add_test(NAME ${tname} COMMAND ${Python3_EXECUTABLE} ${st} s2 $<TARGET_FILE:Spin_GenericSpinChainShowBasis>
            $<TARGET_FILE:Spin_SpinChainComputeS2> ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE} ${args})
        set_tests_properties(${tname} PROPERTIES LABELS "physics;python;spin" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endforeach()
    set(tname openfinding.spin.compute_s2.degenerate_states_without_momentum_first_only)
    add_test(NAME ${tname} COMMAND ${Python3_EXECUTABLE} ${st} s2 $<TARGET_FILE:Spin_GenericSpinChainShowBasis>
        $<TARGET_FILE:Spin_SpinChainComputeS2> ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE} 6 1 0 3)
    set_tests_properties(${tname} PROPERTIES WILL_FAIL TRUE LABELS "open-finding;python;spin" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
endif()
if(Python3_Interpreter_FOUND AND TARGET Spin_GenericSpinChainShowBasis AND TARGET Spin_SpinChainMultipleEntanglementSpectra)
    if(DIAGHAM_USE_LAPACK)
        foreach(case "n8_sz0|8;0;4" "n10_sz2|10;2;3" "n9_sz1|9;1;2")
            string(FIND "${case}" "|" bar)
            string(SUBSTRING "${case}" 0 ${bar} name)
            math(EXPR bar "${bar} + 1")
            string(SUBSTRING "${case}" ${bar} -1 args)
            set(tname physics.spin.multiple_entanglement_spectra.spin_1_2_${name})
            add_test(NAME ${tname} COMMAND ${Python3_EXECUTABLE} ${st} multient $<TARGET_FILE:Spin_GenericSpinChainShowBasis>
                $<TARGET_FILE:Spin_SpinChainMultipleEntanglementSpectra> ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE} ${args})
            set_tests_properties(${tname} PROPERTIES LABELS "physics;python;spin;entanglement" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
        endforeach()
    else()
        set(tname openfinding.spin.multiple_entanglement_spectra.svd_null_without_lapack.spin_1_2_n8_sz0)
        add_test(NAME ${tname} COMMAND ${Python3_EXECUTABLE} ${st} multient $<TARGET_FILE:Spin_GenericSpinChainShowBasis>
            $<TARGET_FILE:Spin_SpinChainMultipleEntanglementSpectra> ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE} 8 0 4)
        set_tests_properties(${tname} PROPERTIES WILL_FAIL TRUE LABELS "open-finding;python;spin" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endif()
endif()

# Overnight run, program 7: FQHESphereMPSCreateState (review tier) with the Laughlin MPS, against the exact integer
# expansion of prod (z_i - z_j)^m (tests/oracles/mps_laughlin.py), in the monomial basis and, with -n, in the
# normalised sphere basis. Bosons need --boson-truncation N and a P truncation that reaches the whole squeezed space.
# U53: with fermions (no --boson) the Laughlin matrix constructor calls ComputeSiteDependentMatrices, which exits 0
# for the fermionic convention, so no state is written.
set(ml ${CMAKE_CURRENT_SOURCE_DIR}/oracles/mps_laughlin.py)
if(Python3_Interpreter_FOUND)
    add_test(NAME selftest.mps_laughlin_oracle COMMAND ${Python3_EXECUTABLE} ${ml} --check)
    set_tests_properties(selftest.mps_laughlin_oracle PROPERTIES LABELS "selftest;python" TIMEOUT 300)
endif()
if(Python3_Interpreter_FOUND AND TARGET FQHEOnSphere_FQHESphereMPSCreateState)
    foreach(case "boson_laughlin_m2_n5|5;2;6" "boson_laughlin_m2_n6|6;2;9" "boson_laughlin_m4_n4|4;4;8" "boson_laughlin_m6_n3|3;6;8"
                 "normalized_sphere.boson_laughlin_m2_n5|5;2;6;sphere" "normalized_sphere.boson_laughlin_m4_n4|4;4;8;sphere")
        string(FIND "${case}" "|" bar)
        string(SUBSTRING "${case}" 0 ${bar} name)
        math(EXPR bar "${bar} + 1")
        string(SUBSTRING "${case}" ${bar} -1 args)
        set(tname physics.fqhe.sphere.mps_create_state.${name})
        add_test(NAME ${tname} COMMAND ${Python3_EXECUTABLE} ${ml} run $<TARGET_FILE:FQHEOnSphere_FQHESphereMPSCreateState>
            ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE} ${args})
        set_tests_properties(${tname} PROPERTIES LABELS "physics;python;fqhe;sphere;mps" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endforeach()
    set(tname openfinding.fqhe.sphere.mps_create_state.fermion_laughlin_exits_without_state.m3_n4)
    add_test(NAME ${tname} COMMAND ${Python3_EXECUTABLE} ${ml} run $<TARGET_FILE:FQHEOnSphere_FQHESphereMPSCreateState>
        ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE} 4 3 8 fermion)
    set_tests_properties(${tname} PROPERTIES WILL_FAIL TRUE LABELS "open-finding;python;fqhe;sphere;mps" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
endif()

# Overnight run, program 8: FQHESphereMPSEntanglementSpectrumParticlePartition --orbital-es (review tier), the bosonic
# Laughlin MPS on a finite cylinder against the Schmidt decomposition of the exact state (tests/oracles/mps_oes.py).
# Without LAPACK the program diagonalises 1x1 overlap matrices with the built-in Householder and crashes (U42).
set(tname_base fqhe.cylinder.mps_orbital_entanglement_spectrum)
if(Python3_Interpreter_FOUND AND TARGET FQHEOnSphere_FQHESphereMPSEntanglementSpectrumParticlePartition)
    set(oes_target $<TARGET_FILE:FQHEOnSphere_FQHESphereMPSEntanglementSpectrumParticlePartition>)
    if(DIAGHAM_USE_LAPACK)
        foreach(case "boson_laughlin_m2_n5_perimeter8|5;2;6;8" "boson_laughlin_m2_n6_perimeter10|6;2;9;10"
                     "boson_laughlin_m4_n4_perimeter9|4;4;8;9" "boson_laughlin_m6_n3_perimeter8|3;6;8;8")
            string(FIND "${case}" "|" bar)
            string(SUBSTRING "${case}" 0 ${bar} name)
            math(EXPR bar "${bar} + 1")
            string(SUBSTRING "${case}" ${bar} -1 args)
            set(tname physics.${tname_base}.${name})
            add_test(NAME ${tname} COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/mps_oes.py run ${oes_target}
                ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE} ${args})
            set_tests_properties(${tname} PROPERTIES LABELS "physics;python;fqhe;cylinder;mps;entanglement" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
        endforeach()
    else()
        set(tname openfinding.${tname_base}.householder_1x1_crash.boson_laughlin_m2_n5_perimeter8)
        add_test(NAME ${tname} COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/mps_oes.py run ${oes_target}
            ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE} 5 2 6 8)
        set_tests_properties(${tname} PROPERTIES WILL_FAIL TRUE LABELS "open-finding;python;fqhe;cylinder;mps" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endif()
endif()

# Overnight run, program 9: FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity and ...ConvertStates (review tier).
# The pairing program's lowest eigenstates of one sector (its spectrum is tested above), then every <i|n_{m,sigma}|j>,
# the integrated charges and the charge imbalance against the same states built by tests/oracles/quasiholes_trs.py;
# for ConvertStates, fixed-N states converted to the space of all N and then measured. U42: without LAPACK the
# built-in diagonaliser that patch 03 falls back to crashes on the 1x1 matrices of a single input state. U54: with an
# odd number of orbitals the charge imbalance counts the central orbital in Q_R, unlike the wiki's definition.
set(qh_dens FQHEOnSphere_FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity)
set(qh_conv FQHEOnSphere_FQHESphereQuasiholesWithSpinTimeReversalSymmetryConvertStates)
if(Python3_Interpreter_FOUND AND TARGET ${qh_prog} AND TARGET ${qh_pair} AND TARGET ${qh_dens} AND TARGET ${qh_conv})
    set(qhd_cases
        "physics.fqhe.sphere.quasiholes_time_reversal_density.pairing_2s6_sz0_lz0_three_states|pairing_2s6_sz0_lz0_three_states|FALSE"
        "physics.fqhe.sphere.quasiholes_time_reversal_density.pairing_2s5_sz1_lz1_two_states|pairing_2s5_sz1_lz1_two_states|FALSE"
        "openfinding.fqhe.sphere.quasiholes_time_reversal_density.imbalance_central_orbital.pairing_2s6_sz0_lz0|pairing_2s6_sz0_lz0_imbalance_central_orbital|TRUE"
        "physics.fqhe.sphere.quasiholes_time_reversal_convert_states.fixed_n4_2s6_sz0_lz0_three_states|convert_fixed_n4_2s6_sz0_lz0_three_states|FALSE")
    if(DIAGHAM_USE_LAPACK)
        list(APPEND qhd_cases "physics.fqhe.sphere.quasiholes_time_reversal_density.pairing_2s5_sz1_lz1_single_state|pairing_2s5_sz1_lz1_single_state|FALSE")
    else()
        list(APPEND qhd_cases "openfinding.fqhe.sphere.quasiholes_time_reversal_density.householder_1x1_crash.pairing_2s5_sz1_lz1_single_state|pairing_2s5_sz1_lz1_single_state|TRUE")
    endif()
    foreach(entry ${qhd_cases})
        string(REPLACE "|" ";" e "${entry}")
        list(GET e 0 tname)
        list(GET e 1 case)
        list(GET e 2 wf)
        if(wf)
            set(lab "open-finding;python;fqhe;sphere")
        else()
            set(lab "physics;python;fqhe;sphere")
        endif()
        add_test(NAME ${tname}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/quasiholes_trs.py density
                $<TARGET_FILE:${qh_prog}> $<TARGET_FILE:${qh_pair}> $<TARGET_FILE:${qh_dens}> ${case} ${DIAGHAM_TEST_WORK}/${tname}
                ${DIAGHAM_TEST_TOLERANCE} $<TARGET_FILE:${qh_conv}>)
        set_tests_properties(${tname} PROPERTIES WILL_FAIL ${wf} LABELS "${lab}" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endforeach()
endif()

# Overnight run, program 10: FQHECylinderConfiningPotentialCoefficients (review tier), against the integral of the
# potential with each orbital's Gaussian density (tests/oracles/confining.py, closed form by recursion). The step
# potential (power 0) needs no GSL. Powers >= 1 need GSL: without it the program writes zeros and exits 0 (U50); with it
# the orbitals centred before the offset get a sign wrong (U55) and the left potential is put on the right (U56), so the
# physics tests keep every orbital centre beyond a right offset.
set(cf ${CMAKE_CURRENT_SOURCE_DIR}/oracles/confining.py)
if(Python3_Interpreter_FOUND)
    add_test(NAME selftest.confining_oracle COMMAND ${Python3_EXECUTABLE} ${cf} --check)
    set_tests_properties(selftest.confining_oracle PROPERTIES LABELS "selftest;python" TIMEOUT 300)
endif()
if(Python3_Interpreter_FOUND AND TARGET FQHEOnCylinder_FQHECylinderConfiningPotentialCoefficients)
    set(cf_cases
        "physics|step.right_2s24_perimeter10_x1.3|-s,24,--cylinder-perimeter,10.0,--confining-rightpower,0,--confining-rightoffset,1.3"
        "physics|step.right_2s15_r1.5_v3_x-0.7_spinful|-s,15,-r,1.5,--confining-rightpower,0,--confining-rightstrength,3,--confining-rightoffset,-0.7,--spinful")
    if(DIAGHAM_USE_GSL)
        list(APPEND cf_cases
            "physics|power2.right_2s24_perimeter10_x-9|-s,24,--cylinder-perimeter,10.0,--confining-rightpower,2,--confining-rightoffset,-9"
            "physics|power3.right_2s17_r0.8_v0.4_x-8|-s,17,-r,0.8,--confining-rightpower,3,--confining-rightstrength,0.4,--confining-rightoffset,-8"
            "physics|power1.right_2s20_perimeter12_v2.5_x-7_spinful|-s,20,--cylinder-perimeter,12,--confining-rightpower,1,--confining-rightstrength,2.5,--confining-rightoffset,-7,--spinful"
            "openfinding|sign_before_offset.wiki_example|-s,24,--cylinder-perimeter,10.0,--confining-rightpower,2,--spinful,--confining-rightoffset,4.5"
            "openfinding|left_potential_on_the_right.power2_2s24_perimeter10_x-4.5|-s,24,--cylinder-perimeter,10.0,--confining-rightstrength,0,--confining-leftstrength,1,--confining-leftpower,2,--confining-leftoffset,-4.5")
    else()
        list(APPEND cf_cases
            "openfinding|zero_without_gsl.power2_right_2s24_perimeter10_x-9|-s,24,--cylinder-perimeter,10.0,--confining-rightpower,2,--confining-rightoffset,-9")
    endif()
    foreach(case ${cf_cases})
        string(FIND "${case}" "|" bar)
        string(SUBSTRING "${case}" 0 ${bar} kind)
        math(EXPR bar "${bar} + 1")
        string(SUBSTRING "${case}" ${bar} -1 rest)
        string(FIND "${rest}" "|" bar)
        string(SUBSTRING "${rest}" 0 ${bar} name)
        math(EXPR bar "${bar} + 1")
        string(SUBSTRING "${rest}" ${bar} -1 args)
        string(REPLACE "," ";" args "${args}")
        set(tname ${kind}.fqhe.cylinder.confining_potential_coefficients.${name})
        add_test(NAME ${tname} COMMAND ${Python3_EXECUTABLE} ${cf} run $<TARGET_FILE:FQHEOnCylinder_FQHECylinderConfiningPotentialCoefficients>
            ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE} ${args})
        if(kind STREQUAL "openfinding")
            set_tests_properties(${tname} PROPERTIES WILL_FAIL TRUE LABELS "open-finding;python;fqhe;cylinder" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
        else()
            set_tests_properties(${tname} PROPERTIES LABELS "physics;python;fqhe;cylinder" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
        endif()
    endforeach()
endif()

# Overnight run, program 11: FQHECylinderDensity (review tier). Generic fermion states written in FQHESphereShowBasis
# order (n-body or Haldane basis); orbital occupations, the density along the axis and the integrated charge against
# Gaussian orbital densities (tests/oracles/cylinder_density.py).
if(Python3_Interpreter_FOUND AND TARGET FQHEOnSphere_FQHESphereShowBasis AND TARGET FQHEOnCylinder_FQHECylinderDensity)
    foreach(case "n4_2s9_lz2_perimeter7.5|4;9;2;7.5;1" "n5_2s12_lz0_perimeter9|5;12;0;9;1" "three_states.n4_2s9_lz0_perimeter8|4;9;0;8;3"
                 "haldane.laughlin_n4_2s9_perimeter9|4;9;0;9;1;1001001001" "haldane.two_states.laughlin_n5_2s12_perimeter10|5;12;0;10;2;1001001001001")
        string(FIND "${case}" "|" bar)
        string(SUBSTRING "${case}" 0 ${bar} name)
        math(EXPR bar "${bar} + 1")
        string(SUBSTRING "${case}" ${bar} -1 args)
        set(tname physics.fqhe.cylinder.density.${name})
        add_test(NAME ${tname} COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/cylinder_density.py run
            $<TARGET_FILE:FQHEOnSphere_FQHESphereShowBasis> $<TARGET_FILE:FQHEOnCylinder_FQHECylinderDensity>
            ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE} ${args})
        set_tests_properties(${tname} PROPERTIES LABELS "physics;python;fqhe;cylinder" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endforeach()
endif()

# Overnight run, program 12: FQHECylinderQuasiholeMatrixElements and
# FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing (review tier), against the V1 zero modes carried
# to cylinder orbitals (tests/oracles/cylinder_quasiholes.py): basis-free invariants of the matrix elements, then the
# spectrum of the two-layer problem with one-body terms (diagonal and off-diagonal, with phases), pairing and charging
# energy. U57: a confining-file line with m > n couples orbitals m and 2m - n. U58: a fixed particle number with
# off-diagonal terms dereferences a null Hilbert space.
set(cqh ${CMAKE_CURRENT_SOURCE_DIR}/oracles/cylinder_quasiholes.py)
set(cqh_prog FQHEOnCylinder_FQHECylinderQuasiholeMatrixElements)
set(cqh_jj FQHEOnCylinder_FQHECylinderQuasiholesWithSpinTimeReversalSymmetryAndJosephsonPairing)
if(Python3_Interpreter_FOUND AND TARGET ${cqh_prog})
    foreach(case "nphi5_perimeter8|5;8" "nphi6_perimeter6|6;6" "nphi7_perimeter5.5|7;5.5")
        string(FIND "${case}" "|" bar)
        string(SUBSTRING "${case}" 0 ${bar} name)
        math(EXPR bar "${bar} + 1")
        string(SUBSTRING "${case}" ${bar} -1 args)
        list(GET args 0 nphi)
        list(GET args 1 per)
        set(tname physics.fqhe.cylinder.quasihole_matrix_elements.${name})
        add_test(NAME ${tname} COMMAND ${Python3_EXECUTABLE} ${cqh} matrixelements $<TARGET_FILE:${cqh_prog}> ${nphi} ${per}
            ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE})
        set_tests_properties(${tname} PROPERTIES LABELS "physics;python;fqhe;cylinder" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
    endforeach()
    if(TARGET ${cqh_jj})
        foreach(entry "physics|diagonal_2s6_perimeter6_sz0|diagonal_2s6_perimeter6_sz0"
                      "physics|diagonal_2s5_perimeter5.5_sz1|diagonal_2s5_perimeter5.5_sz1"
                      "physics|offdiagonal_phases_2s6_perimeter6_sz0|offdiagonal_phases_2s6_perimeter6_sz0"
                      "openfinding|lower_triangle_line.both_triangles_2s6_perimeter6_sz0|both_triangles_2s6_perimeter6_sz0"
                      "openfinding|fixed_n_offdiagonal_crash.offdiagonal_fixed_n4_2s6_perimeter6_sz0|offdiagonal_fixed_n4_2s6_perimeter6_sz0")
            string(REPLACE "|" ";" e "${entry}")
            list(GET e 0 kind)
            list(GET e 1 name)
            list(GET e 2 case)
            set(tname ${kind}.fqhe.cylinder.quasiholes_time_reversal_josephson.${name})
            add_test(NAME ${tname} COMMAND ${Python3_EXECUTABLE} ${cqh} josephson $<TARGET_FILE:${cqh_prog}> $<TARGET_FILE:${cqh_jj}>
                ${case} ${DIAGHAM_TEST_WORK}/${tname} ${DIAGHAM_TEST_TOLERANCE})
            if(kind STREQUAL "openfinding")
                set_tests_properties(${tname} PROPERTIES WILL_FAIL TRUE LABELS "open-finding;python;fqhe;cylinder" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
            else()
                set_tests_properties(${tname} PROPERTIES LABELS "physics;python;fqhe;cylinder" TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
            endif()
        endforeach()
    endif()
endif()
