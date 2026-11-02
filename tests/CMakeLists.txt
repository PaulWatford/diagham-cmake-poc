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
#   open-finding  reproducers for documented upstream findings, registered with
#              WILL_FAIL so they pass while the finding is present and turn
#              red -- prompting their promotion to a physics test -- the
#              moment it is fixed
#
# Tests for a module are registered only when that module is built.
# ============================================================================

add_executable(check_spectrum check_spectrum.cc)
set(DIAGHAM_TEST_DATA ${CMAKE_CURRENT_SOURCE_DIR}/data)
set(DIAGHAM_TEST_WORK ${CMAKE_CURRENT_BINARY_DIR}/work)

# Numerical tolerance of the value comparisons (team decision of 20 October 2026: relative, about 1e-12 to
# 1e-13, configurable). "rel:X" means |a - b| <= X * max(|a|, |b|, 1) (tests/check_spectrum.cc). Tests that
# cannot reach it keep their own tolerance with the reason beside it; see docs/reference/tests.md.
set(DIAGHAM_TEST_TOLERANCE "rel:1e-12" CACHE STRING
    "Tolerance of the physics value comparisons: rel:X (relative) or X (absolute); default rel:1e-12")
# Time budget: every DiagHam test is designed to finish in about 30 s run on its own (measured: see
# docs/reference/tests.md). The ctest TIMEOUT below is larger, for parallel runs on a loaded machine.
set(DIAGHAM_TEST_TIMEOUT 120 CACHE STRING "ctest TIMEOUT (seconds) of each DiagHam program test")
# BLAS threads of the Python oracles (numpy). Without a limit each parallel ctest job starts a thread pool sized for
# the whole machine, which oversubscribes it; the matrices are small, so one thread is as fast. A value already in
# the environment at configure time is kept. The two self-tests that do gain from threads get
# DIAGHAM_TEST_BLAS_THREADS_HEAVY threads and reserve as many ctest slots (PROCESSORS). Applied at the end of this file.
if(DEFINED ENV{OMP_NUM_THREADS})
    set(_blas_threads_default "$ENV{OMP_NUM_THREADS}")
else()
    set(_blas_threads_default 1)
endif()
set(DIAGHAM_TEST_BLAS_THREADS ${_blas_threads_default} CACHE STRING "BLAS threads of each Python oracle test")
set(DIAGHAM_TEST_BLAS_THREADS_HEAVY 4 CACHE STRING
    "BLAS threads (and ctest PROCESSORS) of selftest.hofstadter_oracle and selftest.spin_ed_oracle")

# diagham_physics_test(<name>
#     PROGRAM <program target>
#     ARGS <program arguments...>
#     OUTPUT <glob matching the one spectrum file to check>
#     CHECK <check_spectrum arguments; @OUTPUT@ is the matched file>
#     [REF_PROGRAM <target> REF_ARGS <args...> REF_OUTPUT <glob>]
#         a second program whose output replaces @REFERENCE@ in CHECK, for
#         checking two independent implementations against each other
#     [PRE_STEPS <step>...]   programs run first in the work directory, each
#         step "target|arg|arg..." (the target name is replaced by its file)
#     [LABELS <extra labels...>] [TIMEOUT <seconds>] [WILL_FAIL]
#     [KNOWN_FROM <where the expected answer is known from>]
#     [REASON <why the tolerance is looser than the policy>])
# KNOWN_FROM (test property KNOWN_FROM) starts with its kind: "oracle tests/oracles/<file>.py", "analytic: ...",
# "exact count: ...", "independent program: <target>", "paper: ..." or "saved reference: <file>, cross-checked by
# <test>". Left out, it is "independent program: <REF_PROGRAM>" for a crosscheck and "oracle tests/oracles/<dir>.py"
# for a CHECK that reads tests/data/<dir>/ when that oracle exists. REASON becomes TOLERANCE_REASON. Both are
# indexed by tests/claims.py (docs/reference/claims.md) and enforced by selftest.claims (tests/check_claims.py).
function(diagham_physics_test name)
    cmake_parse_arguments(ARG "WILL_FAIL" "PROGRAM;OUTPUT;TIMEOUT;REF_PROGRAM;REF_OUTPUT;KNOWN_FROM;REASON"
        "ARGS;CHECK;LABELS;REF_ARGS;PRE_STEPS" ${ARGN})
    set(pre_defs "")
    if(ARG_PRE_STEPS)
        set(steps "")
        foreach(step ${ARG_PRE_STEPS})
            string(REPLACE "|" ";" parts "${step}")
            list(GET parts 0 target)
            if(NOT TARGET ${target})
                return()
            endif()
            list(REMOVE_AT parts 0)
            string(REPLACE ";" "|" rest "${parts}")
            list(APPEND steps "$<TARGET_FILE:${target}>|${rest}")
        endforeach()
        string(REPLACE ";" "^^" steps "${steps}")
        set(pre_defs -DPRE_STEPS=${steps})
    endif()
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
        set(ARG_TIMEOUT ${DIAGHAM_TEST_TIMEOUT})
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
            ${pre_defs}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/RunAndCheck.cmake
    )
    set(labels physics ${ARG_LABELS})
    if(ARG_REF_PROGRAM)
        set(labels crosscheck ${ARG_LABELS})
    endif()
    if(ARG_WILL_FAIL)
        set(labels open-finding ${ARG_LABELS})
    endif()
    set_tests_properties(${name} PROPERTIES
        LABELS "${labels}"
        TIMEOUT ${ARG_TIMEOUT}
        WILL_FAIL ${ARG_WILL_FAIL}
    )
    if(NOT ARG_KNOWN_FROM)
        if(ARG_REF_PROGRAM)
            set(ARG_KNOWN_FROM "independent program: ${ARG_REF_PROGRAM}")
        elseif(check_args MATCHES "/tests/data/([A-Za-z0-9_-]+)/")
            # ${CMAKE_MATCH_1} is expanded before if() evaluates, so the match and the EXISTS are two steps.
            set(data_dir ${CMAKE_MATCH_1})
            if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/oracles/${data_dir}.py)
                set(ARG_KNOWN_FROM "oracle tests/oracles/${data_dir}.py")
            endif()
        endif()
    endif()
    if(ARG_KNOWN_FROM)
        set_tests_properties(${name} PROPERTIES KNOWN_FROM "${ARG_KNOWN_FROM}")
    endif()
    if(ARG_REASON)
        set_tests_properties(${name} PROPERTIES TOLERANCE_REASON "${ARG_REASON}")
    endif()
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

# Results are written with upstream's 14 significant digits, so the default build is checked to 1e-12 (the team's
# agreed precision); with DIAGHAM_FULL_PRECISION_OUTPUT (17 digits) the same value is checked to 4 ulp.
if(DIAGHAM_FULL_PRECISION_OUTPUT)
    set(hubbard_u4_check min @OUTPUT@ -1 -5.6568542494923806 4)
else()
    set(hubbard_u4_check min-abs @OUTPUT@ -1 -5.6568542494923806 1e-12)
endif()
diagham_physics_test(physics.hubbard.2x2.U4.ground_state_is_minus_4sqrt2
    KNOWN_FROM "analytic: 2x2 Hubbard, half filling, U = 4t, Sz = 0: E0 = -4*sqrt(2), correctly rounded with mpmath"
    PROGRAM ${hubbard}
    ARGS -p 4 -x 2 -y 2 --u-potential 4 --nn-t 1.0 --full-diag 1000
    OUTPUT "*_u_4.000000_sz_0.dat"
    CHECK ${hubbard_u4_check}
    LABELS fti hubbard
)

# U = 0: free fermions, ground state is the filled band, exactly -8t.
diagham_physics_test(physics.hubbard.2x2.U0.ground_state_is_minus_8
    KNOWN_FROM "analytic: at U = 0 the 2x2 lattice is free fermions and the filled band gives exactly -8t"
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
    KNOWN_FROM "saved reference: benchmarks/fermions_hubbard_square_x_2_y_2_n_4_ns_4_t_1.000000_tp_0.000000_u_4.000000_sz_0.dat.txt, cross-checked by benchmarks/hubbard_ed.py (benchmarks/BENCHMARK.md test 5, outside ctest)"
    PROGRAM ${hubbard}
    ARGS -p 4 -x 2 -y 2 --u-potential 4 --nn-t 1.0 --full-diag 1000
    OUTPUT "*_u_4.000000_sz_0.dat"
    CHECK spectrum @OUTPUT@ -1
          ${CMAKE_SOURCE_DIR}/benchmarks/fermions_hubbard_square_x_2_y_2_n_4_ns_4_t_1.000000_tp_0.000000_u_4.000000_sz_0.dat.txt
          -1 ${DIAGHAM_TEST_TOLERANCE}
    LABELS fti hubbard
)

# 2x4, U = 4, lowest state of each of the 8 momentum sectors (basis
# dimension 4900, Lanczos) against the saved reference.
diagham_physics_test(physics.hubbard.2x4.U4.sector_ground_states
    KNOWN_FROM "saved reference: benchmarks/fermions_hubbard_square_x_2_y_4_n_8_ns_8_t_1.000000_tp_0.000000_u_4.000000_sz_0.dat.txt, cross-checked by benchmarks/hubbard_ed.py (benchmarks/BENCHMARK.md test 5, outside ctest)"
    PROGRAM ${hubbard}
    ARGS -p 8 -x 2 -y 4 --u-potential 4 --nn-t 1.0 -n 1
    OUTPUT "*_u_4.000000_sz_0.dat"
    CHECK spectrum @OUTPUT@ -1
          ${CMAKE_SOURCE_DIR}/benchmarks/fermions_hubbard_square_x_2_y_4_n_8_ns_8_t_1.000000_tp_0.000000_u_4.000000_sz_0.dat.txt
          -1 ${DIAGHAM_TEST_TOLERANCE}
    LABELS fti hubbard
)

# ----------------------------------------------------------------------------
# FQHE: Laughlin states as exact zero modes of the V1 pseudopotential
# ----------------------------------------------------------------------------
set(v1 ${DIAGHAM_TEST_DATA}/v1_pseudopotential.dat)

# Sphere, N = 6, 2S = 15 = 3(N-1): the Laughlin state is the unique zero
# mode; in the Lz = 0 sector (dimension 338) exactly one eigenvalue is 0.
diagham_physics_test(physics.fqhe.sphere.laughlin_1_3.unique_zero_mode
    KNOWN_FROM "exact count: at 2S = 3(N-1) the Laughlin state is the unique V1 zero mode, one state at Lz = 0"
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
    KNOWN_FROM "exact count: two bosonic quasiholes of angular momentum N/2 = 3 form L = 6, 4, 2, 0, four zero modes at Lz = 0"
    PROGRAM FQHEOnSphere_FQHESphereFermionsTwoBodyGeneric
    ARGS -p 6 -l 17 --interaction-file ${v1} --interaction-name v1 --nbr-lz 1 --full-diag 5000
    OUTPUT "fermions_v1_n_6_2s_17_lz.dat"
    CHECK count @OUTPUT@ -1 0 1e-10 4
    LABELS fqhe
)

# Torus, N = 4, Nphi = 12, all Ky sectors (495 states): the Laughlin state
# is 3-fold degenerate (topological degeneracy q = 3), all at zero energy.
diagham_physics_test(physics.fqhe.torus.laughlin_1_3.threefold_degeneracy
    KNOWN_FROM "exact count: the Laughlin 1/3 state on the torus is threefold degenerate (q = 3)"
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
# also reproduced by an independent numpy ED to 1e-15. L = 10 (the Lanczos path) is in
# tests/spin_hubbard.cmake, against tests/oracles/spin_hubbard.py.
diagham_physics_test(physics.spin.heisenberg_ring_L4
    KNOWN_FROM "analytic: the 4-site spin-1/2 Heisenberg ring has E0 = -2J"
    PROGRAM Spin_GenericPeriodicSpinChain
    ARGS -p 4 --nbr-sz 1
    OUTPUT "spin_1_2_periodicchain_n_4.dat"
    CHECK min-abs @OUTPUT@ -1 -2 1e-12
    LABELS spin
)
diagham_physics_test(physics.spin.heisenberg_ring_L6
    KNOWN_FROM "analytic: the 6-site spin-1/2 Heisenberg ring has E0 = -(2 + sqrt(13))/2"
    PROGRAM Spin_GenericPeriodicSpinChain
    ARGS -p 6 --nbr-sz 1
    OUTPUT "spin_1_2_periodicchain_n_6.dat"
    CHECK min-abs @OUTPUT@ -1 -2.8027756377319946 1e-12
    LABELS spin
)

# ----------------------------------------------------------------------------
# FQHE: spinful (SU(2)) torus Coulomb, FQHETorusFermionsWithSpin
#
# Cross-checks guarding the spinful torus Coulomb fix
# (docs/explanation/torus-su2-coulomb-finding.md). Each compares the program
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
    CHECK spectrum @OUTPUT@ -1 @REFERENCE@ -1 ${DIAGHAM_TEST_TOLERANCE}
    LABELS fqhe
)
diagham_physics_test(crosscheck.fqhe.torus.su2_coulomb.unpolarized_N3_matches_translations
    REASON "DiagHam truncates its momentum sums at about 1e-10, so its degenerate pairs differ in the tenth digit (docs/reference/tests.md, tolerance table)"
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

include(${CMAKE_CURRENT_SOURCE_DIR}/dimensions.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/fqhe_sphere.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/fqhe_geometries.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/spin_hubbard.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/independent.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/entanglement.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/haah.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/xcube.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/core_programs.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/regression.cmake)

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
# Which locked list applies: every program and every module -> manifest.txt, exactly; the core set alone
# (FQHE, FTI and Spin on; QuantumDots has no core program) -> manifest-core-programs.txt, exactly; any other
# selection (custom, or core plus extra programs) -> SUBSET: every registered test must be in manifest.txt, so a
# renamed or stray test still fails, but tests of programs not built may be absent.
set(_manifest_mode "")
if(NOT DIAGHAM_PROGRAM_SELECTION AND DIAGHAM_BUILD_FQHE AND DIAGHAM_BUILD_FTI AND DIAGHAM_BUILD_SPIN AND DIAGHAM_BUILD_QUANTUMDOTS)
    set(_manifest_mode -DMANIFEST=${CMAKE_CURRENT_SOURCE_DIR}/manifest.txt)
elseif(DIAGHAM_PROGRAM_SET STREQUAL "core" AND NOT DIAGHAM_PROGRAMS AND DIAGHAM_BUILD_FQHE AND DIAGHAM_BUILD_FTI AND DIAGHAM_BUILD_SPIN)
    set(_manifest_mode -DMANIFEST=${CMAKE_CURRENT_SOURCE_DIR}/manifest-core-programs.txt)
elseif(DIAGHAM_PROGRAM_SELECTION)
    set(_manifest_mode -DMANIFEST=${CMAKE_CURRENT_SOURCE_DIR}/manifest.txt -DSUBSET=ON)
endif()
if(_manifest_mode)
    add_test(NAME selftest.manifest
        COMMAND ${CMAKE_COMMAND} -DBUILD_DIR=${CMAKE_BINARY_DIR} ${_manifest_mode}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/CheckManifest.cmake)
    set_tests_properties(selftest.manifest PROPERTIES LABELS "selftest" TIMEOUT 120)
endif()

# Self-test: the claim rules (tests/check_claims.py): every physics and crosscheck test says where its answer is known
# from (KNOWN_FROM, required since 2 November 2026), loose tolerances carry a TOLERANCE_REASON, labels are in
# tests/labels.txt, WILL_FAIL goes with the label open-finding, and every open-finding test is in the register.
if(Python3_FOUND)
    add_test(NAME selftest.claims
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_claims.py ${CMAKE_BINARY_DIR}
            --known-from-required)
    set_tests_properties(selftest.claims PROPERTIES LABELS "selftest" TIMEOUT 120
        ENVIRONMENT "CTEST_COMMAND=${CMAKE_CTEST_COMMAND}")
endif()

# Where the answer of each Python oracle test is known from (KNOWN_FROM, see diagham_physics_test): the oracle the
# test runs, by test name. selftest.claims checks that the oracle named here is the script in the test's command.
set(_python_known_from
    "^python\\.hubbard_ed_cross_check$=tests/hubbard_cross_check.py"
    "^physics\\.fqhe\\.cylinder\\.confining_potential_coefficients\\.=tests/oracles/confining.py"
    "^physics\\.core\\.matrix_extract_columns\\.=tests/oracles/core_tools.py"
    "^physics\\.fqhe\\.sphere\\.(boson_entanglement|su2_spin_rotation|unnormalize_state)\\.=tests/oracles/core_tools.py"
    "^physics\\.fti\\.(particle_entanglement|show_basis)\\.=tests/oracles/core_tools.py"
    "^physics\\.fqhe\\.cp2\\.=tests/oracles/cp2.py"
    "^physics\\.fqhe\\.cylinder\\.density\\.=tests/oracles/cylinder_density.py"
    "^physics\\.fqhe\\.cylinder\\.(quasihole_matrix_elements|quasiholes_time_reversal_josephson)\\.=tests/oracles/cylinder_quasiholes.py"
    "^physics\\.fti\\.generic_interaction_two_bands\\.=tests/oracles/fti_generic.py"
    "^physics\\.hubbard\\.superconductor_order_parameter\\.=tests/oracles/hubbard_order_parameter.py"
    "^physics\\.(fti\\.min_nbr_singlet_projection|hubbard\\.square_lattice)\\.=tests/oracles/hubbard_square.py"
    "^physics\\.fqhe\\.sphere\\.mps_create_state\\.=tests/oracles/mps_laughlin.py"
    "^physics\\.fqhe\\.sphere\\.(quasihole_matrix_elements|quasiholes_time_reversal_[a-z_]+)\\.=tests/oracles/quasiholes_trs.py"
    "^physics\\.fqhe\\.(cylinder|sphere|disk)\\.realspace_partition_coefficients\\.=tests/oracles/realspace.py"
    "^physics\\.spin\\.(compute_s2|multiple_entanglement_spectra)\\.=tests/oracles/spin_tools.py"
    "^physics\\.fqhe\\.cylinder\\.mps_orbital_entanglement_spectrum\\.=tests/oracles/mps_oes.py"
    "^physics\\.fqhe\\.torus\\.mps_create_state\\.=tests/oracles/torus_mps.py"
    "^physics\\.fqhe\\.cylinder\\.quasiholes_time_reversal_generate_basis\\.=tests/oracles/trs_generate_basis.py"
    "^physics\\.core\\.vector_(ascii2binary|binary2ascii)\\.=tests/oracles/vector_io.py")
get_property(_all_tests DIRECTORY PROPERTY TESTS)
foreach(_t IN LISTS _all_tests)
    get_test_property(${_t} LABELS _labels)
    get_test_property(${_t} KNOWN_FROM _known)
    if(_labels AND "python" IN_LIST _labels AND NOT _known)
        foreach(_entry IN LISTS _python_known_from)
            string(FIND "${_entry}" "=" _eq REVERSE)
            string(SUBSTRING "${_entry}" 0 ${_eq} _re)
            math(EXPR _eq "${_eq} + 1")
            string(SUBSTRING "${_entry}" ${_eq} -1 _oracle)
            if(_t MATCHES "${_re}")
                set_tests_properties(${_t} PROPERTIES KNOWN_FROM "oracle ${_oracle}")
                break()
            endif()
        endforeach()
    endif()
endforeach()

# Thread limit of every Python oracle test (label "python"; see DIAGHAM_TEST_BLAS_THREADS above).
foreach(_t IN LISTS _all_tests)
    get_test_property(${_t} LABELS _labels)
    if(_labels AND "python" IN_LIST _labels)
        set(_n ${DIAGHAM_TEST_BLAS_THREADS})
        if(_t MATCHES "^selftest\\.(hofstadter|spin_ed)_oracle$")
            set(_n ${DIAGHAM_TEST_BLAS_THREADS_HEAVY})
            set_tests_properties(${_t} PROPERTIES PROCESSORS ${_n})
        endif()
        set_tests_properties(${_t} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=${_n};OPENBLAS_NUM_THREADS=${_n};MKL_NUM_THREADS=${_n}")
    endif()
endforeach()
