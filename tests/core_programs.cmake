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
