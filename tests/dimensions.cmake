# ----------------------------------------------------------------------------
# Hilbert-space dimension goldens (label: physics, dimension)
#
# The number of many-body basis states in a sector is exact combinatorics
# and is computed without DiagHam by tests/oracles/dimensions.py; the
# expected files under tests/data/dimensions/ are its output (each file
# says which formula). The program's printed dimensions must match line by
# line, or the number of basis states it prints must equal the count.
# selftest.dimension_oracle re-derives every expected file (needs Python 3).
#
# What this proves: the Hilbert-space classes enumerate the right states
# in the right sectors -- for the sphere (Lz and L), the torus with
# magnetic translations (k_x, k_y), the disk, single-band lattice models
# (kx, ky), the Hubbard model with and without Sz, and spin chains. It says
# nothing about matrix elements or eigenvalues; those are the other goldens.
# ----------------------------------------------------------------------------

# diagham_dimension_test(<name> PROGRAM <target> ARGS <args...> EXPECT <case>
#                        [COUNT_REGEX <regex>] [LABELS <extra>...])
function(diagham_dimension_test name)
    cmake_parse_arguments(ARG "" "PROGRAM;EXPECT;COUNT_REGEX;TIMEOUT" "ARGS;LABELS" ${ARGN})
    if(NOT TARGET ${ARG_PROGRAM})
        return()
    endif()
    if(NOT ARG_TIMEOUT)
        set(ARG_TIMEOUT 120)
    endif()
    string(REPLACE ";" "|" program_args "${ARG_ARGS}")
    set(extra "")
    if(ARG_COUNT_REGEX)
        set(extra -DCOUNT_REGEX=${ARG_COUNT_REGEX})
    endif()
    add_test(NAME ${name}
        COMMAND ${CMAKE_COMMAND}
            -DPROGRAM=$<TARGET_FILE:${ARG_PROGRAM}>
            -DPROGRAM_ARGS=${program_args}
            -DWORK_DIR=${DIAGHAM_TEST_WORK}/${name}
            -DEXPECT_FILE=${DIAGHAM_TEST_DATA}/dimensions/${ARG_EXPECT}.txt
            ${extra}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/RunAndMatch.cmake)
    set_tests_properties(${name} PROPERTIES LABELS "physics;dimension;${ARG_LABELS}" TIMEOUT ${ARG_TIMEOUT})
endfunction()

# --- sphere: FQHESphereGetDimension prints "Lz = ..." (states per 2Lz = 0,2,...)
#     and "L = ..." (per L); both lines must match exactly.
set(sph FQHEOnSphere_FQHESphereGetDimension)
diagham_dimension_test(physics.dimension.sphere.fermions_n4_2s9   PROGRAM ${sph} ARGS -n 4 -s 9 --fermion  EXPECT sphere_getdim_fermions_n4_2s9  LABELS fqhe)
diagham_dimension_test(physics.dimension.sphere.fermions_n6_2s15  PROGRAM ${sph} ARGS -n 6 -s 15 --fermion EXPECT sphere_getdim_fermions_n6_2s15 LABELS fqhe)
diagham_dimension_test(physics.dimension.sphere.bosons_n4_2s6     PROGRAM ${sph} ARGS -n 4 -s 6 --boson    EXPECT sphere_getdim_bosons_n4_2s6    LABELS fqhe)
diagham_dimension_test(physics.dimension.sphere.bosons_n5_2s8     PROGRAM ${sph} ARGS -n 5 -s 8 --boson    EXPECT sphere_getdim_bosons_n5_2s8    LABELS fqhe)
# --- sphere: FQHESphereShowBasis prints one occupation line per basis state
set(sphb FQHEOnSphere_FQHESphereShowBasis)
diagham_dimension_test(physics.dimension.sphere.showbasis_fermions_n3_2s6_lz0 PROGRAM ${sphb} ARGS -p 3 -l 6 -z 0 --fermion EXPECT sphere_showbasis_fermions_n3_2s6_lz0 COUNT_REGEX "^[01]( [01])*$" LABELS fqhe)
diagham_dimension_test(physics.dimension.sphere.showbasis_fermions_n4_2s9_lz0 PROGRAM ${sphb} ARGS -p 4 -l 9 -z 0 --fermion EXPECT sphere_showbasis_fermions_n4_2s9_lz0 COUNT_REGEX "^[01]( [01])*$" LABELS fqhe)
diagham_dimension_test(physics.dimension.sphere.showbasis_bosons_n3_2s4_lz0   PROGRAM ${sphb} ARGS -p 3 -l 4 -z 0 --boson   EXPECT sphere_showbasis_bosons_n3_2s4_lz0   COUNT_REGEX "^[0-9]+( [0-9]+)*$" LABELS fqhe)

# --- torus with magnetic translations: FQHETorusGetDimension prints one
#     " (k_x = a, k_y = b) : n" line per sector; every line must match.
set(tor FQHEOnTorus_FQHETorusGetDimension)
diagham_dimension_test(physics.dimension.torus.bosons_n3_nphi6    PROGRAM ${tor} ARGS -p 3 -q 6 --boson    EXPECT torus_getdim_bosons_n3_nphi6    LABELS fqhe)
diagham_dimension_test(physics.dimension.torus.bosons_n4_nphi8    PROGRAM ${tor} ARGS -p 4 -q 8 --boson    EXPECT torus_getdim_bosons_n4_nphi8    LABELS fqhe)
diagham_dimension_test(physics.dimension.torus.fermions_n3_nphi6  PROGRAM ${tor} ARGS -p 3 -q 6 --fermion  EXPECT torus_getdim_fermions_n3_nphi6  LABELS fqhe)
diagham_dimension_test(physics.dimension.torus.fermions_n4_nphi12 PROGRAM ${tor} ARGS -p 4 -q 12 --fermion EXPECT torus_getdim_fermions_n4_nphi12 LABELS fqhe)
set(torb FQHEOnTorus_FQHETorusShowBasis)
diagham_dimension_test(physics.dimension.torus.showbasis_fermions_n4_nphi12_k00 PROGRAM ${torb} ARGS -p 4 -l 12 -x 0 -y 0 --fermion EXPECT torus_showbasis_fermions_n4_nphi12_kx0_ky0 COUNT_REGEX "^[01]( [01])*$" LABELS fqhe)
diagham_dimension_test(physics.dimension.torus.showbasis_bosons_n3_nphi6_k00    PROGRAM ${torb} ARGS -p 3 -l 6 -x 0 -y 0 --boson    EXPECT torus_showbasis_bosons_n3_nphi6_kx0_ky0    COUNT_REGEX "^[0-9]+( [0-9]+)*$" LABELS fqhe)

# --- disk: FQHEDiskShowBasis, one line per state (bosons append key/lzmax text)
set(disk FQHEOnDisk_FQHEDiskShowBasis)
diagham_dimension_test(physics.dimension.disk.fermions_n3_lz6  PROGRAM ${disk} ARGS -p 3 -z 6 --fermion  EXPECT disk_showbasis_fermions_n3_lz6  COUNT_REGEX "^[01]( [01])*$"  LABELS fqhe)
diagham_dimension_test(physics.dimension.disk.fermions_n4_lz10 PROGRAM ${disk} ARGS -p 4 -z 10 --fermion EXPECT disk_showbasis_fermions_n4_lz10 COUNT_REGEX "^[01]( [01])*$"  LABELS fqhe)
diagham_dimension_test(physics.dimension.disk.bosons_n3_lz6    PROGRAM ${disk} ARGS -p 3 -z 6 --boson    EXPECT disk_showbasis_bosons_n3_lz6    COUNT_REGEX "^[0-9]+( [0-9]+)*  key =" LABELS fqhe)

# --- single-band lattice models: FTIGetDimension, all (kx,ky) sectors
#     (--no-inversion: by default the program prints only kx <= Nx/2)
set(fti FTI_FTIGetDimension)
diagham_dimension_test(physics.dimension.fti.fermions_n4_3x3 PROGRAM ${fti} ARGS -p 4 -x 3 -y 3 --no-inversion          EXPECT fti_getdim_fermions_n4_3x3 LABELS fti)
diagham_dimension_test(physics.dimension.fti.bosons_n3_2x3   PROGRAM ${fti} ARGS -p 3 -x 2 -y 3 --bosons --no-inversion EXPECT fti_getdim_bosons_n3_2x3   LABELS fti)
diagham_dimension_test(physics.dimension.fti.fermions_n3_4x2 PROGRAM ${fti} ARGS -p 3 -x 4 -y 2 --no-inversion          EXPECT fti_getdim_fermions_n3_4x2 LABELS fti)

# --- Hubbard: HubbardGetDimension prints "Hilbert space dimension = n"
set(hub HubbardModels_HubbardGetDimension)
diagham_dimension_test(physics.dimension.hubbard.n4_sites4     PROGRAM ${hub} ARGS -p 4 -x 4                              EXPECT hubbard_getdim_n4_sites4     LABELS fti hubbard)
diagham_dimension_test(physics.dimension.hubbard.n4_sites4_sz0 PROGRAM ${hub} ARGS -p 4 -x 4 --conserve-sz --total-spin 0 EXPECT hubbard_getdim_n4_sites4_sz0 LABELS fti hubbard)
diagham_dimension_test(physics.dimension.hubbard.n3_sites4_sz1 PROGRAM ${hub} ARGS -p 3 -x 4 --conserve-sz --total-spin 1 EXPECT hubbard_getdim_n3_sites4_sz1 LABELS fti hubbard)
diagham_dimension_test(physics.dimension.hubbard.n6_sites6_sz0 PROGRAM ${hub} ARGS -p 6 -x 6 --conserve-sz --total-spin 0 EXPECT hubbard_getdim_n6_sites6_sz0 LABELS fti hubbard)

# --- spin chains: GenericSpinChainShowBasis, one line per state
set(spn Spin_GenericSpinChainShowBasis)
diagham_dimension_test(physics.dimension.spin.s1_2_n6_sz0  PROGRAM ${spn} ARGS -s 1 -p 6 -z 0 EXPECT spin_showbasis_s1_2_n6_sz0  COUNT_REGEX "^[-+]( [-+])*$"       LABELS spin)
diagham_dimension_test(physics.dimension.spin.s1_2_n8_2sz2 PROGRAM ${spn} ARGS -s 1 -p 8 -z 2 EXPECT spin_showbasis_s1_2_n8_2sz2 COUNT_REGEX "^[-+]( [-+])*$"       LABELS spin)
diagham_dimension_test(physics.dimension.spin.s1_n4_sz0    PROGRAM ${spn} ARGS -s 2 -p 4 -z 0 EXPECT spin_showbasis_s1_n4_sz0    COUNT_REGEX "^-?[0-9]+( -?[0-9]+)*$"   LABELS spin)
diagham_dimension_test(physics.dimension.spin.s1_n6_sz0    PROGRAM ${spn} ARGS -s 2 -p 6 -z 0 EXPECT spin_showbasis_s1_n6_sz0    COUNT_REGEX "^-?[0-9]+( -?[0-9]+)*$"   LABELS spin)

# --- the oracle re-derives every committed expected file
if(Python3_Interpreter_FOUND)
    add_test(NAME selftest.dimension_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/dimensions.py --check ${DIAGHAM_TEST_DATA}/dimensions)
    set_tests_properties(selftest.dimension_oracle PROPERTIES LABELS "selftest;python" TIMEOUT 300)
endif()
