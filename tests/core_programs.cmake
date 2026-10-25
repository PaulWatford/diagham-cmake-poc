# ----------------------------------------------------------------------------
# Core-program gaps (session S4, 25 October 2026): physics tests for programs of the provisional core list
# that had only a start-up check. Each answer is known without DiagHam:
#  - VectorAscii2Binary / VectorBinary2Ascii: the documented binary vector layout, written and read by
#    tests/oracles/vector_io.py with Python's struct module (bit-exact one way, 14 printed digits the other);
#  - LevelStatistics: synthetic spectra whose spacing ratio is exact (1/2 for spacings 1, 2, 1, 2, ...; 1 for
#    equal spacings);
#  - FQHESphereL2Diagonalize: the L^2 spectrum of an Lz sector is L(L+1) once per multiplet, the multiplet
#    counts d(L) - d(L+1) from exact state counting (tests/oracles/sphere_l2.py);
#  - FCIHofstadterModel: the lowest band of one particle against an independent Hofstadter exact
#    diagonalisation (tests/oracles/hofstadter.py), on three lattices and two fluxes.
# Two open findings met on the way: U44 (--real-space with fermions prints "case not implented yet" and
# crashes) and U45 (--real-space hard-core bosons: ground state below the bound set by the program's own
# one-particle spectrum). Both are the same on a pristine autotools build of r4493.
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

# U44: --real-space with fermions prints "case not implented yet" and then segfaults (exit 139), with or
# without interaction; the same on a pristine autotools build of r4493. WILL_FAIL until it is fixed.
diagham_physics_test(openfinding.fci.hofstadter.realspace_fermions_segfault
    PROGRAM FCI_FCIHofstadterModel
    ARGS --real-space -p 2 -x 4 -y 2 -X 1 -Y 3 -q 1 --full-diag 5000
    OUTPUT "*.dat"
    CHECK count @OUTPUT@ -1 0 1e-10 0
    WILL_FAIL LABELS fci hofstadter)

# U45: --real-space hard-core bosons give a spectrum whose ground state lies below 2 e_0, the rigorous lower
# bound from the program's own one-particle spectrum (-6.129 against -5.464 here), and which matches the
# hard-core restriction of that one-particle Hamiltonian in neither Landau gauge nor any uniform flux.
# These tests carry the spectrum of the consistent model; byte-identical output on a pristine build of r4493.
foreach(case "2|4|2|1|3|1" "3|3|2|1|4|1")
    string(REPLACE "|" ";" c "${case}")
    list(GET c 0 n)
    list(GET c 1 x)
    list(GET c 2 y)
    list(GET c 3 cx)
    list(GET c 4 cy)
    list(GET c 5 q)
    diagham_physics_test(openfinding.fci.hofstadter.realspace_hardcore_bosons_below_bound.n${n}_x${x}_y${y}
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
