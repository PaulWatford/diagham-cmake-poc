# ============================================================================
# DiagHam ctest suite
# ============================================================================
#
# Labels (select with `ctest -L <label>`, or `ctest --preset physics`):
#
#   physics    the expected answer is known independently of DiagHam
#              (analytic, exact counting, an independent solver); a program
#              is run in a fresh directory and its spectrum file is checked
#              (tests/check_spectrum.cc). Every golden says in a comment why
#              the answer is known.
#   crosscheck two DiagHam programs, or two algorithms, must agree
#              (diagham_physics_test with REF_PROGRAM); no independent
#              answer is involved
#   regression a spectrum saved from an earlier build must be reproduced;
#              this detects change, it does not establish correctness
#   smoke      every built program is started with --help and must exit 0
#              (catches link-time and static-initialisation breakage; it is
#              not a physics check and docs/reference/test-coverage.md does
#              not count it as one)
#   selftest   the harness itself (the test manifest, the checker)
#   install    the Development install component is usable through
#              find_package(DiagHam) from a separate project
#   python     independent Python exact diagonalisation cross-check
#              (registered only when Python 3 with numpy is available)
#   known-bug  reproducers for documented upstream defects, registered with
#              WILL_FAIL so they pass while the defect is present and turn
#              red -- prompting their promotion to a physics test -- the
#              moment it is fixed
#
# Tests for a module are registered only when that module is built.
# ============================================================================

add_executable(check_spectrum check_spectrum.cc)
set(DIAGHAM_TEST_DATA ${CMAKE_CURRENT_SOURCE_DIR}/data)
set(DIAGHAM_TEST_WORK ${CMAKE_CURRENT_BINARY_DIR}/work)

# diagham_physics_test(<name>
#     PROGRAM <program target>
#     ARGS <program arguments...>
#     OUTPUT <glob matching the one spectrum file to check>
#     CHECK <check_spectrum arguments; @OUTPUT@ is the matched file>
#     [REF_PROGRAM <target> REF_ARGS <args...> REF_OUTPUT <glob>]
#         a second program whose output replaces @REFERENCE@ in CHECK, for
#         checking two independent implementations against each other
#     [LABELS <extra labels...>] [TIMEOUT <seconds>] [WILL_FAIL])
function(diagham_physics_test name)
    cmake_parse_arguments(ARG "WILL_FAIL" "PROGRAM;OUTPUT;TIMEOUT;REF_PROGRAM;REF_OUTPUT"
        "ARGS;CHECK;LABELS;REF_ARGS" ${ARGN})
    if(NOT TARGET ${ARG_PROGRAM})
        return()
    endif()
    set(reference_defs "")
    if(ARG_REF_PROGRAM)
        if(NOT TARGET ${ARG_REF_PROGRAM})
            return()
        endif()
        string(REPLACE ";" "|" ref_args "${ARG_REF_ARGS}")
        set(reference_defs
            -DREF_PROGRAM=$<TARGET_FILE:${ARG_REF_PROGRAM}>
            -DREF_PROGRAM_ARGS=${ref_args}
            -DREF_OUTPUT_GLOB=${ARG_REF_OUTPUT})
    endif()
    if(NOT ARG_TIMEOUT)
        set(ARG_TIMEOUT 300)
    endif()
    # "|" rather than ";" separates list items so they survive -D intact.
    string(REPLACE ";" "|" program_args "${ARG_ARGS}")
    string(REPLACE ";" "|" check_args "${ARG_CHECK}")
    add_test(NAME ${name}
        COMMAND ${CMAKE_COMMAND}
            -DPROGRAM=$<TARGET_FILE:${ARG_PROGRAM}>
            -DPROGRAM_ARGS=${program_args}
            -DWORK_DIR=${DIAGHAM_TEST_WORK}/${name}
            -DOUTPUT_GLOB=${ARG_OUTPUT}
            -DCHECKER=$<TARGET_FILE:check_spectrum>
            -DCHECK_ARGS=${check_args}
            ${reference_defs}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/RunAndCheck.cmake
    )
    set(labels physics ${ARG_LABELS})
    if(ARG_REF_PROGRAM)
        set(labels crosscheck ${ARG_LABELS})
    endif()
    if(ARG_WILL_FAIL)
        set(labels known-bug ${ARG_LABELS})
    endif()
    set_tests_properties(${name} PROPERTIES
        LABELS "${labels}"
        TIMEOUT ${ARG_TIMEOUT}
        WILL_FAIL ${ARG_WILL_FAIL}
    )
endfunction()

# ----------------------------------------------------------------------------
# FTI: Hubbard model on the square lattice (HubbardSquareLatticeModel)
#
# The locked golden (AGENTS.md): 2x2, half filling, U = 4t has the analytic
# ground state -4*sqrt(2). -5.6568542494923806 is its correctly rounded
# double (mpmath: -5.65685424949238019520...). GCC 13 -O3 reproduces it to
# 0 ulp; 4 ulp of headroom absorbs compiler/FMA differences without letting
# a real regression (which moves the 12th digit or earlier) through.
# ----------------------------------------------------------------------------
set(hubbard HubbardModels_HubbardSquareLatticeModel)

diagham_physics_test(physics.hubbard.2x2.U4.ground_state_is_minus_4sqrt2
    PROGRAM ${hubbard}
    ARGS -p 4 -x 2 -y 2 --u-potential 4 --nn-t 1.0 --full-diag 1000
    OUTPUT "*_u_4.000000_sz_0.dat"
    CHECK min @OUTPUT@ -1 -5.6568542494923806 4
    LABELS fti hubbard
)

# U = 0: free fermions, ground state is the filled band, exactly -8t.
diagham_physics_test(physics.hubbard.2x2.U0.ground_state_is_minus_8
    PROGRAM ${hubbard}
    ARGS -p 4 -x 2 -y 2 --u-potential 0 --nn-t 1.0 --full-diag 1000
    # (at U = 0 the program leaves "_u_..." out of the file name)
    OUTPUT "*_tp_0.000000_sz_0.dat"
    CHECK min-abs @OUTPUT@ -1 -8 1e-12
    LABELS fti hubbard
)

# Full 36-state 2x2 U = 4 spectrum against the saved reference in
# benchmarks/ (written with 14 significant digits, hence 1e-11).
diagham_physics_test(physics.hubbard.2x2.U4.full_spectrum
    PROGRAM ${hubbard}
    ARGS -p 4 -x 2 -y 2 --u-potential 4 --nn-t 1.0 --full-diag 1000
    OUTPUT "*_u_4.000000_sz_0.dat"
    CHECK spectrum @OUTPUT@ -1
          ${CMAKE_SOURCE_DIR}/benchmarks/fermions_hubbard_square_x_2_y_2_n_4_ns_4_t_1.000000_tp_0.000000_u_4.000000_sz_0.dat.txt
          -1 1e-11
    LABELS fti hubbard
)

# 2x4, U = 4, lowest state of each of the 8 momentum sectors (basis
# dimension 4900, Lanczos) against the saved reference.
diagham_physics_test(physics.hubbard.2x4.U4.sector_ground_states
    PROGRAM ${hubbard}
    ARGS -p 8 -x 2 -y 4 --u-potential 4 --nn-t 1.0 -n 1
    OUTPUT "*_u_4.000000_sz_0.dat"
    CHECK spectrum @OUTPUT@ -1
          ${CMAKE_SOURCE_DIR}/benchmarks/fermions_hubbard_square_x_2_y_4_n_8_ns_8_t_1.000000_tp_0.000000_u_4.000000_sz_0.dat.txt
          -1 1e-10
    LABELS fti hubbard
)

# ----------------------------------------------------------------------------
# FQHE: Laughlin states as exact zero modes of the V1 pseudopotential
# ----------------------------------------------------------------------------
set(v1 ${DIAGHAM_TEST_DATA}/v1_pseudopotential.dat)

# Sphere, N = 6, 2S = 15 = 3(N-1): the Laughlin state is the unique zero
# mode; in the Lz = 0 sector (dimension 338) exactly one eigenvalue is 0.
diagham_physics_test(physics.fqhe.sphere.laughlin_1_3.unique_zero_mode
    PROGRAM FQHEOnSphere_FQHESphereFermionsTwoBodyGeneric
    ARGS -p 6 -l 15 --interaction-file ${v1} --interaction-name v1 --nbr-lz 1 --full-diag 2000
    OUTPUT "fermions_v1_n_6_2s_15_lz.dat"
    CHECK count @OUTPUT@ -1 0 1e-10 1
    LABELS fqhe
)

# Sphere, N = 6, 2S = 17: two Laughlin quasiholes. Each carries angular
# momentum N/2 = 3 and they are bosonic, so the zero modes form
# L = 6, 4, 2, 0 multiplets (28 states in all), i.e. exactly 4 zero modes
# in the Lz = 0 sector.
diagham_physics_test(physics.fqhe.sphere.laughlin_1_3.two_quasihole_count
    PROGRAM FQHEOnSphere_FQHESphereFermionsTwoBodyGeneric
    ARGS -p 6 -l 17 --interaction-file ${v1} --interaction-name v1 --nbr-lz 1 --full-diag 5000
    OUTPUT "fermions_v1_n_6_2s_17_lz.dat"
    CHECK count @OUTPUT@ -1 0 1e-10 4
    LABELS fqhe
)

# Torus, N = 4, Nphi = 12, all Ky sectors (495 states): the Laughlin state
# is 3-fold degenerate (topological degeneracy q = 3), all at zero energy.
diagham_physics_test(physics.fqhe.torus.laughlin_1_3.threefold_degeneracy
    PROGRAM FQHEOnTorus_FQHETorusFermionsTwoBodyGeneric
    ARGS -p 4 -l 12 --interaction-file ${v1} --interaction-name v1 --redundant-kymomenta --full-diag 5000
    OUTPUT "fermions_torus_kysym_v1_n_4_2s_12_ratio_1.000000.dat"
    CHECK count @OUTPUT@ -1 0 1e-10 3
    LABELS fqhe
)

# ----------------------------------------------------------------------------
# Spin
# ----------------------------------------------------------------------------
# S = 1/2 Heisenberg ring H = J sum_i S_i . S_{i+1}, J = 1, Sz = 0 sector,
# all momenta (GenericPeriodicSpinChain). L = 4: E0 = -2 exactly. L = 6:
# E0 = -(2 + sqrt(13))/2 = -2.8027756377319946 in closed form; both are
# also reproduced by an independent numpy ED to 1e-15. L = 10 checks the Lanczos
# path against the numpy ED value -4.5154463544920365.
diagham_physics_test(physics.spin.heisenberg_ring_L4
    PROGRAM Spin_GenericPeriodicSpinChain
    ARGS -p 4 --nbr-sz 1
    OUTPUT "spin_1_2_periodicchain_n_4.dat"
    CHECK min-abs @OUTPUT@ -1 -2 1e-12
    LABELS spin
)
diagham_physics_test(physics.spin.heisenberg_ring_L6
    PROGRAM Spin_GenericPeriodicSpinChain
    ARGS -p 6 --nbr-sz 1
    OUTPUT "spin_1_2_periodicchain_n_6.dat"
    CHECK min-abs @OUTPUT@ -1 -2.8027756377319946 1e-12
    LABELS spin
)
diagham_physics_test(physics.spin.heisenberg_ring_L10
    PROGRAM Spin_GenericPeriodicSpinChain
    ARGS -p 10 --nbr-sz 1
    OUTPUT "spin_1_2_periodicchain_n_10.dat"
    CHECK min-abs @OUTPUT@ -1 -4.5154463544920365 1e-10
    LABELS spin
)

# ----------------------------------------------------------------------------
# FQHE: spinful (SU(2)) torus Coulomb, FQHETorusFermionsWithSpin
#
# Cross-checks guarding the spinful torus Coulomb fix
# (docs/explanation/torus-su2-coulomb-defect.md). Each compares the program
# with another DiagHam implementation, full spectrum:
#
#  - fully polarised (total spin = N) must reproduce the spinless Coulomb
#    spectrum (FQHETorusFermionsCoulomb). Before patch 14, N = 3 was off by
#    up to 0.25 because the basis ordering broke FindStateIndex.
#  - unpolarised (Sz = 0) must reproduce FQHETorusFermionsWithSpinAndTranslations,
#    which uses different Hilbert-space and Hamiltonian classes. Before
#    patch 14 the up-down interaction dropped every same-orbital term and
#    the N = 4 ground state came out at -1.2734 instead of -0.7624.
# ----------------------------------------------------------------------------
diagham_physics_test(crosscheck.fqhe.torus.su2_coulomb.polarized_N3_matches_spinless
    PROGRAM FQHEOnTorus_FQHETorusFermionsWithSpin
    ARGS --nbr-particles 3 --max-momentum 9 --total-spin 3 --ratio 1.0 --full-diag 9000 --redundantYMomenta
    OUTPUT "fermions_torus_su2_coulomb_n_3_2s_9_Sz_3_ratio_1.000000.dat"
    REF_PROGRAM FQHEOnTorus_FQHETorusFermionsCoulomb
    REF_ARGS --nbr-particles 3 --max-momentum 9 --ratio 1.0 --full-diag 9000
    REF_OUTPUT "fermions_torus_kysym_coulomb_n_3_2s_9_ratio_1.000000.dat"
    CHECK spectrum @OUTPUT@ -1 @REFERENCE@ -1 1e-10
    LABELS fqhe
)
diagham_physics_test(crosscheck.fqhe.torus.su2_coulomb.unpolarized_N3_matches_translations
    PROGRAM FQHEOnTorus_FQHETorusFermionsWithSpin
    ARGS --nbr-particles 3 --max-momentum 9 --total-spin 1 --ratio 1.0 --full-diag 9000 --redundantYMomenta
    OUTPUT "fermions_torus_su2_coulomb_n_3_2s_9_Sz_1_ratio_1.000000.dat"
    REF_PROGRAM FQHEOnTorus_FQHETorusFermionsWithSpinAndTranslations
    REF_ARGS --nbr-particles 3 --max-momentum 9 --total-spin 1 --ratio 1.0 --full-diag 9000 --all-points
    REF_OUTPUT "*.dat"
    CHECK spectrum @OUTPUT@ -1 @REFERENCE@ -1 1e-8
    LABELS fqhe
)

# ----------------------------------------------------------------------------
# Independent Python ED cross-check: DiagHam vs benchmarks/hubbard_ed.py
# ----------------------------------------------------------------------------
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND AND TARGET ${hubbard})
    execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy"
        RESULT_VARIABLE numpy_missing OUTPUT_QUIET ERROR_QUIET)
    if(numpy_missing EQUAL 0)
        add_test(NAME python.hubbard_ed_cross_check
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/hubbard_cross_check.py
                $<TARGET_FILE:${hubbard}> ${DIAGHAM_TEST_WORK}/python.hubbard_ed_cross_check)
        set_tests_properties(python.hubbard_ed_cross_check PROPERTIES
            LABELS "python;physics" TIMEOUT 600)
    else()
        message(STATUS "DiagHam tests: numpy not found, python.hubbard_ed_cross_check not registered")
    endif()
endif()

# ----------------------------------------------------------------------------
# Smoke: every program answers --help
#
# Three upstream programs have no option parser at all and are left out:
#   QHEBosons           reads argv[1] as a particle count, so "--help"
#                       becomes 0 bosons and it aborts
#                       (std::bad_array_new_length)
#   MultipleSpinChain   same pattern (argv[1] = spins per chain), segfaults
#   TestDiagHamVectors  a vector-routine benchmark whose main() takes no
#                       arguments and runs for minutes
# ----------------------------------------------------------------------------
#   EvaluateBroadening  (QuantumDots/src/Tools/Analysis) tests for its required
#                       --input before it looks at --help, so --help returns -1
#                       (upstream quirk; register entry U26)
set(DIAGHAM_SMOKE_NO_OPTION_PARSER QHEBosons MultipleSpinChain TestDiagHamVectors EvaluateBroadening)

get_property(programs GLOBAL PROPERTY DIAGHAM_ALL_PROGRAMS)
foreach(program ${programs})
    get_target_property(output_name ${program} OUTPUT_NAME)
    if(NOT output_name)
        set(output_name ${program})
    endif()
    if(output_name IN_LIST DIAGHAM_SMOKE_NO_OPTION_PARSER)
        continue()
    endif()
    add_test(NAME smoke.${program} COMMAND ${program} --help)
    set_tests_properties(smoke.${program} PROPERTIES LABELS "smoke" TIMEOUT 30)
endforeach()

# ----------------------------------------------------------------------------
# Install: find_package(DiagHam) from a separate project
# ----------------------------------------------------------------------------
add_test(NAME install.find_package_consumer
    COMMAND ${CMAKE_COMMAND}
        -DBUILD_DIR=${CMAKE_BINARY_DIR}
        -DCONSUMER_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/consumer
        -DWORK_DIR=${DIAGHAM_TEST_WORK}/install.find_package_consumer
        -DCXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DGENERATOR=${CMAKE_GENERATOR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/InstallConsumer.cmake
)
set_tests_properties(install.find_package_consumer PROPERTIES
    LABELS "install" TIMEOUT 600)

# ----------------------------------------------------------------------------
# Self-test: the registered tests equal tests/manifest.txt (a green run that
# silently lost tests must not pass). The manifest describes the full build;
# it is checked only when every module is built.
# ----------------------------------------------------------------------------
if(DIAGHAM_BUILD_FQHE AND DIAGHAM_BUILD_FTI AND DIAGHAM_BUILD_SPIN AND DIAGHAM_BUILD_QUANTUMDOTS)
    add_test(NAME selftest.manifest
        COMMAND ${CMAKE_COMMAND}
            -DBUILD_DIR=${CMAKE_BINARY_DIR}
            -DMANIFEST=${CMAKE_CURRENT_SOURCE_DIR}/manifest.txt
            -P ${CMAKE_CURRENT_SOURCE_DIR}/CheckManifest.cmake)
    set_tests_properties(selftest.manifest PROPERTIES LABELS "selftest" TIMEOUT 120)
endif()
