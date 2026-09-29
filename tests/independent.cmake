# ----------------------------------------------------------------------------
# Independent oracles, algorithm consensus, and self-tests of the harness.
#
#  1. tests/oracles/sphere_ed.py: a from-scratch exact diagonalisation of
#     pseudopotential Hamiltonians in the lowest Landau level (Clebsch-Gordan
#     pair projectors, explicit Fock basis, numpy). The program's whole
#     spectrum of an Lz sector must match it (label: physics, sphere, ed).
#  1a. the same oracle with SU(2), SU(3) and SU(4) spin (spectrum_su2, spectrum_sun),
#     including the Sz and Lz parity sectors of the symmetrised bases (label: physics,
#     fqhe, sphere, ed, spin).
#  1a'. tests/oracles/sphere_nbody.py: three- and four-body pseudopotentials on the sphere
#     from L-multiplet projectors (label: physics, fqhe, sphere, ed, nbody).
#  1b. tests/oracles/torus_ed.py: the same on the torus (Landau gauge, Coulomb
#     or pseudopotentials, fermions or bosons; label: physics, fqhe, torus, ed).
#  1b'. the same on the torus with SU(2)/SU(3)/SU(4) spin, bilayer Coulomb, Landau
#     level 1 and the magnetic-translation (Kx, Ky) sectors (spectrum_species).
#  1b''. tests/oracles/geometry_ed.py: two-body pseudopotentials on the cylinder and the
#     disk (label: physics, fqhe, cylinder / disk, ed).
#  1c. tests/oracles/fci_bands.py: the checkerboard lattice bands from the
#     published Bloch Hamiltonian (label: physics, fci, checkerboard, bands).
#  2. tests/oracles/spin_ed.py: dense numpy diagonalisation of XXZ / J1-J2
#     chains with a field, open or periodic (label: physics, spin, ed).
#  2b. tests/oracles/spin_models.py: dense ED of models given as operator lists
#     (XYZ, 2D Heisenberg and Ising, J1-J2, double triangle, generalised AKLT,
#     O'Brien-Fendley as coded, disorder, Potts; label: physics, spin, ed).
#  2c. tests/oracles/lattice_fermions.py: Fock-basis ED of Hubbard-family models (square,
#     Haldane honeycomb, SSH; label: physics, fti, hubbard, ed).
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
        "bosons_n6_2s8_coulomb|${bos2}|6|8|0"
        "fermions_n4_2s8_ll1_coulomb|${fer2}|4|10|0")
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

# --- 1a. sphere ED with SU(2), SU(3) and SU(4) spin -----------------------------------------------
# The same oracle with several species (tests/oracles/sphere_ed.py spectrum_su2 / spectrum_sun); the
# pseudopotential files carry the species-pair keys the programs read. Same-species and
# different-species pairs of relative angular momentum m both cost 2 V_m.
set(su2f FQHEOnSphere_FQHESphereFermionsWithSpin)
set(su2b FQHEOnSphere_FQHESphereBosonsWithSpin)
set(su3f FQHEOnSphere_FQHESphereFermionsWithSU3Spin)
set(su3b FQHEOnSphere_FQHESphereBosonsWithSU3Spin)
set(su4b FQHEOnSphere_FQHESphereBosonsWithSU4Spin)
function(diagham_spin_sphere_ed_test name case)
    cmake_parse_arguments(ARG "" "PROGRAM;OUTPUT" "ARGS" ${ARGN})
    diagham_physics_test(${name}
        PROGRAM ${ARG_PROGRAM}
        ARGS ${ARG_ARGS} --interaction-file ${ed_data}/${case}_pp.dat --interaction-name ed --full-diag 100000
        OUTPUT "${ARG_OUTPUT}"
        CHECK spectrum @OUTPUT@ -1 ${ed_data}/${case}_spectrum.dat -1 1e-9
        LABELS fqhe sphere ed spin)
endfunction()
diagham_spin_sphere_ed_test(physics.fqhe.sphere.ed.su2.fermions_n4_2s6_coulomb_sz0 fermions_su2_n4_2s6_coulomb_sz0
    PROGRAM ${su2f} ARGS -p 4 -l 6 -s 0 --initial-lz 0 --nbr-lz 1 OUTPUT "fermions_sphere_su2_ed_n_4_2s_6_sz_0_lz*.dat")
diagham_spin_sphere_ed_test(physics.fqhe.sphere.ed.su2.fermions_n4_2s6_aniso_sz0 fermions_su2_n4_2s6_aniso_sz0
    PROGRAM ${su2f} ARGS -p 4 -l 6 -s 0 --initial-lz 0 --nbr-lz 1 OUTPUT "fermions_sphere_su2_ed_n_4_2s_6_sz_0_lz*.dat")
diagham_spin_sphere_ed_test(physics.fqhe.sphere.ed.su2.fermions_n4_2s7_coulomb_sz2 fermions_su2_n4_2s7_coulomb_sz2
    PROGRAM ${su2f} ARGS -p 4 -l 7 -s 2 --initial-lz 0 --nbr-lz 1 OUTPUT "fermions_sphere_su2_ed_n_4_2s_7_sz_2_lz*.dat")
diagham_spin_sphere_ed_test(physics.fqhe.sphere.ed.su2.bosons_n4_2s4_v0_v2_sz0 bosons_su2_n4_2s4_v0_v2_sz0
    PROGRAM ${su2b} ARGS -p 4 -l 4 -s 0 --initial-lz 0 --nbr-lz 1 OUTPUT "bosons_sphere_su2_ed_n_4_2s_4_sz_0_lz*.dat")
diagham_spin_sphere_ed_test(physics.fqhe.sphere.ed.su2.bosons_n4_2s5_v0_v2_sz2 bosons_su2_n4_2s5_v0_v2_sz2
    PROGRAM ${su2b} ARGS -p 4 -l 5 -s 2 --initial-lz 0 --nbr-lz 1 OUTPUT "bosons_sphere_su2_ed_n_4_2s_5_sz_2_lz*.dat")
# the symmetrised bases: each parity sector of the Coulomb Sz=0 case against the oracle's projection
diagham_spin_sphere_ed_test(physics.fqhe.sphere.ed.su2.symmetrized_basis.sz_plus fermions_su2_n4_2s6_coulomb_sz0_szplus
    PROGRAM ${su2f} ARGS -p 4 -l 6 -s 0 --initial-lz 0 --nbr-lz 1 --szsymmetrized-basis OUTPUT "fermions_sphere_su2_ed_n_4_2s_6_sz_0_lz*.dat")
diagham_spin_sphere_ed_test(physics.fqhe.sphere.ed.su2.symmetrized_basis.sz_minus fermions_su2_n4_2s6_coulomb_sz0_szminus
    PROGRAM ${su2f} ARGS -p 4 -l 6 -s 0 --initial-lz 0 --nbr-lz 1 --szsymmetrized-basis --minus-szparity OUTPUT "fermions_sphere_su2_ed_n_4_2s_6_sz_0_lz*.dat")
diagham_spin_sphere_ed_test(physics.fqhe.sphere.ed.su2.symmetrized_basis.lz_plus fermions_su2_n4_2s6_coulomb_sz0_lzplus
    PROGRAM ${su2f} ARGS -p 4 -l 6 -s 0 --initial-lz 0 --nbr-lz 1 --lzsymmetrized-basis OUTPUT "fermions_sphere_su2_ed_n_4_2s_6_sz_0_lz*.dat")
diagham_spin_sphere_ed_test(physics.fqhe.sphere.ed.su2.symmetrized_basis.lz_minus fermions_su2_n4_2s6_coulomb_sz0_lzminus
    PROGRAM ${su2f} ARGS -p 4 -l 6 -s 0 --initial-lz 0 --nbr-lz 1 --lzsymmetrized-basis --minus-lzparity OUTPUT "fermions_sphere_su2_ed_n_4_2s_6_sz_0_lz*.dat")
# U32: with both --lzsymmetrized-basis and --szsymmetrized-basis the fermionic basis comes out empty
# (dimension 0 for every size tried; the bosonic program is fine). The oracle's even-even sector has
# 20 states. WILL_FAIL: passes while the defect is present, turns red when the basis is fixed.
diagham_physics_test(knownbug.fqhe.sphere.fermions_with_spin.lz_sz_symmetrized_basis_empty
    PROGRAM ${su2f}
    ARGS -p 4 -l 6 -s 0 --initial-lz 0 --nbr-lz 1 --lzsymmetrized-basis --szsymmetrized-basis --interaction-file ${ed_data}/fermions_su2_n4_2s6_coulomb_sz0_szplus_lzplus_pp.dat --interaction-name ed --full-diag 100000
    OUTPUT "fermions_sphere_su2_ed_n_4_2s_6_sz_0_lz*.dat"
    CHECK spectrum @OUTPUT@ -1 ${ed_data}/fermions_su2_n4_2s6_coulomb_sz0_szplus_lzplus_spectrum.dat -1 1e-9
    WILL_FAIL
    LABELS fqhe sphere spin)
# the legacy program QHEFermionsSphereWithSpin (-v, -w): its own normalisation, every Lz >= 0 sector
diagham_physics_test(physics.fqhe.sphere.ed.su2.legacy_fermions_n4_2s6_v1_w0.5
    PROGRAM FQHEOnSphere_QHEFermionsSphereWithSpin
    ARGS -p 4 -l 6 -s 0 -v 1 -w 0.5 --initial-lz 0 --nbr-lz 1 --full-diag 100000
    OUTPUT "fermions_sphere_spin_n_4_2S_6_Sz_0_V_*_lz*.dat"
    CHECK spectrum @OUTPUT@ -1 ${ed_data}/fermions_su2_n4_2s6_legacy_v1_w0.5_spectrum.dat -1 1e-9
    LABELS fqhe sphere ed spin)
# SU(3) and SU(4)
diagham_spin_sphere_ed_test(physics.fqhe.sphere.ed.su3.fermions_n3_2s4_coulomb fermions_su3_n3_2s4_coulomb
    PROGRAM ${su3f} ARGS -p 3 -l 4 -t 0 -y 0 --initial-lz 0 --nbr-lz 1 OUTPUT "fermions_sphere_su3_ed_n_3_2s_4_tz_0_y_0*.dat")
diagham_spin_sphere_ed_test(physics.fqhe.sphere.ed.su3.fermions_n4_2s5_generic_tz1_y1 fermions_su3_n4_2s5_generic_tz1_y1
    PROGRAM ${su3f} ARGS -p 4 -l 5 -t 1 -y 1 --initial-lz 0 --nbr-lz 1 OUTPUT "fermions_sphere_su3_ed_n_4_2s_5_tz_1_y_1*.dat")
diagham_spin_sphere_ed_test(physics.fqhe.sphere.ed.su3.bosons_n4_2s3_generic bosons_su3_n4_2s3_generic
    PROGRAM ${su3b} ARGS -p 4 -l 3 --nbr-n1 2 --nbr-n2 1 --nbr-n3 1 --initial-lz 0 --nbr-lz 1 OUTPUT "bosons_sphere_su3_ed_n_4_2s_3_tz_1_y_1*.dat")
diagham_spin_sphere_ed_test(physics.fqhe.sphere.ed.su3.bosons_n3_2s5_coulomb bosons_su3_n3_2s5_coulomb
    PROGRAM ${su3b} ARGS -p 3 -l 5 --nbr-n1 1 --nbr-n2 1 --nbr-n3 1 --initial-lz 1 --nbr-lz 1 OUTPUT "bosons_sphere_su3_ed_n_3_2s_5_tz_0_y_0*.dat")
diagham_spin_sphere_ed_test(physics.fqhe.sphere.ed.su4.bosons_n4_2s3_generic bosons_su4_n4_2s3_generic
    PROGRAM ${su4b} ARGS -p 4 -l 3 --nbr-n1 1 --nbr-n2 1 --nbr-n3 1 --nbr-n4 1 --initial-lz 0 --nbr-lz 1 OUTPUT "bosons_sphere_su4_ed_n_4_2s_3_sz_0_iz_0_pz_0*.dat")

# --- 1a'. k-body sphere ED --------------------------------------------------------------------------
# tests/oracles/sphere_nbody.py: three- and four-body pseudopotential Hamiltonians from the L-multiplet
# projectors of the k-particle space (no Clebsch-Gordan algebra). DiagHam's k-body pseudopotentials
# multiply projectors normalised per (k, m, 2S); those constants are read off the k-particle runs
# (the only thing taken from the programs), every N > k spectrum is a prediction.
set(nbody_data ${DIAGHAM_TEST_DATA}/sphere_nbody)
set(bos3 FQHEOnSphere_FQHESphereBosonsThreeBodyGeneric)
set(fer3 FQHEOnSphere_FQHESphereFermionsThreeBodyGeneric)
function(diagham_nbody_test name case)
    cmake_parse_arguments(ARG "" "PROGRAM;OUTPUT" "ARGS" ${ARGN})
    diagham_physics_test(${name}
        PROGRAM ${ARG_PROGRAM}
        ARGS ${ARG_ARGS} --initial-lz 0 --nbr-lz 1 --full-diag 100000
        OUTPUT "${ARG_OUTPUT}"
        CHECK spectrum @OUTPUT@ -1 ${nbody_data}/${case}_spectrum.dat -1 1e-9
        LABELS fqhe sphere ed nbody)
endfunction()
diagham_nbody_test(physics.fqhe.sphere.ed.nbody.bosons_n5_2s6_hardcore3 bosons_n5_2s6_hardcore3
    PROGRAM ${bosN} ARGS -p 5 -l 6 --nbr-nbody 3 OUTPUT "bosons_hardcore_nbody_3_n_5_2s_6_lz*.dat")
diagham_nbody_test(physics.fqhe.sphere.ed.nbody.bosons_n6_2s4_hardcore4 bosons_n6_2s4_hardcore4
    PROGRAM ${bosN} ARGS -p 6 -l 4 --nbr-nbody 4 OUTPUT "bosons_hardcore_nbody_4_n_6_2s_4_lz*.dat")
diagham_nbody_test(physics.fqhe.sphere.ed.nbody.fermions_n5_2s10_hardcore3 fermions_n5_2s10_hardcore3
    PROGRAM ${ferN} ARGS -p 5 -l 10 --nbr-nbody 3 OUTPUT "fermions_hardcore_nbody_3_n_5_2s_10_lz*.dat")
diagham_nbody_test(physics.fqhe.sphere.ed.nbody.bosons_n5_2s6_v0_v2_v3 bosons_n5_2s6_v0_v2_v3
    PROGRAM ${bos3} ARGS -p 5 -l 6 --interaction-file ${nbody_data}/bosons_n5_2s6_v0_v2_v3_pp.dat --interaction-name ed OUTPUT "bosons_ed_n_5_2s_6_lz*.dat")
diagham_nbody_test(physics.fqhe.sphere.ed.nbody.bosons_n5_2s6_v0_v2_v3_plus_twobody bosons_n5_2s6_v0_v2_v3_plus_twobody
    PROGRAM ${bos3} ARGS -p 5 -l 6 --interaction-file ${nbody_data}/bosons_n5_2s6_v0_v2_v3_plus_twobody_pp.dat --interaction-name ed OUTPUT "bosons_ed_n_5_2s_6_lz*.dat")
diagham_nbody_test(physics.fqhe.sphere.ed.nbody.bosons_n5_2s8_v0_v2 bosons_n5_2s8_v0_v2
    PROGRAM ${bos3} ARGS -p 5 -l 8 --interaction-file ${nbody_data}/bosons_n5_2s8_v0_v2_pp.dat --interaction-name ed OUTPUT "bosons_ed_n_5_2s_8_lz*.dat")
diagham_nbody_test(physics.fqhe.sphere.ed.nbody.fermions_n5_2s10_v3_v5 fermions_n5_2s10_v3_v5
    PROGRAM ${fer3} ARGS -p 5 -l 10 --interaction-file ${nbody_data}/fermions_n5_2s10_v3_v5_pp.dat --interaction-name ed OUTPUT "fermions_ed_n_5_2s_10_lz*.dat")
diagham_nbody_test(physics.fqhe.sphere.ed.nbody.fermions_n5_2s10_v3_v5_plus_twobody fermions_n5_2s10_v3_v5_plus_twobody
    PROGRAM ${fer3} ARGS -p 5 -l 10 --interaction-file ${nbody_data}/fermions_n5_2s10_v3_v5_plus_twobody_pp.dat --interaction-name ed OUTPUT "fermions_ed_n_5_2s_10_lz*.dat")

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

# --- 1b'. torus ED with spin, Landau level and magnetic-translation sectors ------------------------
# spectrum_species in tests/oracles/torus_ed.py: several species (SU(2), SU(3), SU(4); bilayer Coulomb
# 2 pi exp(-q d)/q between layers), Landau level n (form factor [L_n(q^2/2)]^2), and the (Kx, Ky)
# sectors of the translation programs. Pseudopotential files carry the species-pair keys and
# "Name = ed"; Coulomb programs take no file. Reference rows are "Ky E" or "Kx Ky E".
function(diagham_torus_species_test name case)
    cmake_parse_arguments(ARG "PP;KX" "PROGRAM;OUTPUT" "ARGS" ${ARGN})
    set(extra "")
    if(ARG_PP)
        set(extra --interaction-file ${torus_ed_data}/${case}_pp.dat)
    endif()
    if(ARG_KX)
        set(col 2)
    else()
        set(col -1)
    endif()
    diagham_physics_test(${name}
        PROGRAM ${ARG_PROGRAM}
        ARGS ${ARG_ARGS} ${extra} --full-diag 100000
        OUTPUT "${ARG_OUTPUT}"
        CHECK spectrum @OUTPUT@ ${col} ${torus_ed_data}/${case}_spectrum.dat ${col} 1e-8
        LABELS fqhe torus ed)
endfunction()
set(tfwt FQHEOnTorus_FQHETorusFermionsWithTranslations)
set(tbwt FQHEOnTorus_FQHETorusBosonsWithTranslations)
set(tfs FQHEOnTorus_FQHETorusFermionsWithSpin)
set(tfst FQHEOnTorus_FQHETorusFermionsWithSpinAndTranslations)
set(tfsg FQHEOnTorus_FQHETorusFermionsWithSpinTwoBodyGeneric)
set(tbsg FQHEOnTorus_FQHETorusBosonsWithSpinTwoBodyGeneric)
set(tbst FQHEOnTorus_FQHETorusBosonsWithSpinAndTranslations)
set(tb3g FQHEOnTorus_FQHETorusBosonsWithSU3SpinTwoBodyGeneric)
set(tb3t FQHEOnTorus_FQHETorusBosonsWithSU3SpinAndTranslations)
set(tb4g FQHEOnTorus_FQHETorusBosonsWithSU4SpinTwoBodyGeneric)
diagham_torus_species_test(physics.fqhe.torus.ed.fermions_n3_nphi9_coulomb_ll1_ky0 fermions_n3_nphi9_coulomb_ll1_ky0
    PROGRAM ${tfc} ARGS -p 3 -l 9 --ratio 1.0 --landau-level 1 --ky-momentum 0 OUTPUT "fermions_torus_kysym_coulomb_n_3_2s_9_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.fermions_n4_nphi8_v1v3_allky fermions_n4_nphi8_v1v3_allky PP
    PROGRAM FQHEOnTorus_FQHETorusFermionsTwoBodyGenericAllMomenta ARGS -p 4 -l 8 -r 1 --interaction-name ed OUTPUT "fermions_torus_noky_ed_n_4_2s_8_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.translations.fermions_n4_nphi8_v1v3_kx2_ky0 fermions_n4_nphi8_v1v3_kx2_ky0 PP KX
    PROGRAM ${tfwt} ARGS -p 4 -l 8 -R 1 -x 2 -y 0 OUTPUT "fermions_torus_ed_n_4_2s_8_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.translations.fermions_n4_nphi8_v1v3_kx1_ky2 fermions_n4_nphi8_v1v3_kx1_ky2 PP KX
    PROGRAM ${tfwt} ARGS -p 4 -l 8 -R 1 -x 1 -y 2 OUTPUT "fermions_torus_ed_n_4_2s_8_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.translations.bosons_n4_nphi6_v0v2_kx0_ky0 bosons_n4_nphi6_v0v2_kx0_ky0 PP KX
    PROGRAM ${tbwt} ARGS -p 4 -l 6 -R 1 -x 0 -y 0 OUTPUT "bosons_torus_ed_n_4_2s_6_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.translations.bosons_n4_nphi6_v0v2_kx1_ky2 bosons_n4_nphi6_v0v2_kx1_ky2 PP KX
    PROGRAM ${tbwt} ARGS -p 4 -l 6 -R 1 -x 1 -y 2 OUTPUT "bosons_torus_ed_n_4_2s_6_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.su2.fermions_n4_nphi8_generic_sz0_ky0 fermions_su2_n4_nphi8_generic_sz0_ky0 PP
    PROGRAM ${tfsg} ARGS -p 4 -l 8 -s 0 -r 1 -y 0 --interaction-name ed OUTPUT "fermions_torus_su2_kysym_ed_n_4_2s_8_sz_0_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.su2.fermions_n4_nphi8_generic_sz2_ky0 fermions_su2_n4_nphi8_generic_sz2_ky0 PP
    PROGRAM ${tfsg} ARGS -p 4 -l 8 -s 2 -r 1 -y 0 --interaction-name ed OUTPUT "fermions_torus_su2_kysym_ed_n_4_2s_8_sz_2_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.su2.bosons_n4_nphi6_generic_sz0_ky0 bosons_su2_n4_nphi6_generic_sz0_ky0 PP
    PROGRAM ${tbsg} ARGS -p 4 -l 6 -s 0 -r 1 -y 0 --interaction-name ed OUTPUT "bosons_torus_su2_kysym_ed_n_4_2s_6_sz_0_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.su2.fermions_n4_nphi8_coulomb_sz0_ky0 fermions_su2_n4_nphi8_coulomb_sz0_ky0
    PROGRAM ${tfs} ARGS -p 4 -l 8 -s 0 -R 1 -y 0 OUTPUT "fermions_torus_su2_coulomb_n_4_2s_8_Sz_0_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.su2.fermions_n4_nphi8_bilayer_d1_sz0_ky0 fermions_su2_n4_nphi8_bilayer_d1_sz0_ky0
    PROGRAM ${tfs} ARGS -p 4 -l 8 -s 0 -R 1 -y 0 -d 1.0 OUTPUT "fermions_torus_d_1.000000_coulomb_n_4_2s_8_Sz_0_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.su2.translations.fermions_n4_nphi8_bilayer_d0.5_kx0_ky0 fermions_su2_n4_nphi8_bilayer_d0.5_sz0_kx0_ky0 KX
    PROGRAM ${tfst} ARGS -p 4 -l 8 -s 0 -r 1 -d 0.5 -x 0 -y 0 OUTPUT "fermions_torus_su2_coulomb_n_4_2s_8_d_0.500000_sz_0_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.su2.translations.fermions_n4_nphi8_bilayer_d0.5_kx1_ky3 fermions_su2_n4_nphi8_bilayer_d0.5_sz0_kx1_ky3 KX
    PROGRAM ${tfst} ARGS -p 4 -l 8 -s 0 -r 1 -d 0.5 -x 1 -y 3 OUTPUT "fermions_torus_su2_coulomb_n_4_2s_8_d_0.500000_sz_0_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.su2.translations.bosons_n4_nphi6_generic_kx1_ky1 bosons_su2_n4_nphi6_generic_sz0_kx1_ky1 PP KX
    PROGRAM ${tbst} ARGS -p 4 -l 6 -s 0 -r 1 -x 1 -y 1 --interaction-name ed OUTPUT "bosons_torus_su2_ed_n_4_2s_6_sz_0_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.su3.bosons_n4_nphi6_generic_ky0 bosons_su3_n4_nphi6_generic_ky0 PP
    PROGRAM ${tb3g} ARGS -p 4 -l 6 --nbr-n1 2 --nbr-n2 1 --nbr-n3 1 -r 1 -y 0 --interaction-name ed OUTPUT "bosons_torus_su3_kysym_ed_n_4_2s_6_tz_1_y_1_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.su3.translations.bosons_n3_nphi6_generic_kx1_ky2 bosons_su3_n3_nphi6_generic_kx1_ky2 PP KX
    PROGRAM ${tb3t} ARGS -p 3 -l 6 --nbr-n1 1 --nbr-n2 1 --nbr-n3 1 -r 1 -x 1 -y 2 --interaction-name ed OUTPUT "bosons_torus_su3_ed_n_3_2s_6_tz_0_y_0_ratio_*.dat")
diagham_torus_species_test(physics.fqhe.torus.ed.su4.bosons_n4_nphi5_generic_ky0 bosons_su4_n4_nphi5_generic_ky0 PP
    PROGRAM ${tb4g} ARGS -p 4 -l 5 --nbr-n1 1 --nbr-n2 1 --nbr-n3 1 --nbr-n4 1 -r 1 -y 0 --interaction-name ed OUTPUT "bosons_torus_su4_kysym_ed_n_4_2s_5_sz_0_iz_0_pz_0_ratio_*.dat")

# --- 1b''. cylinder and disk ED ---------------------------------------------------------------------
# tests/oracles/geometry_ed.py: the Landau-gauge matrix element on a cylinder (no periodic images,
# integral over q_x; L = sqrt(2 pi r N_orb) as the programs define it; the programs label momentum
# sectors by twice the momentum, -y 2k --nbr-ky 1), and the symmetric-gauge pair projectors on the
# disk. The pseudopotential cylinder programs need GSL (U28), so those tests exist only with it and
# are listed in tests/manifest-optional.txt; the delta programs run in every build.
set(geo_data ${DIAGHAM_TEST_DATA}/geometry_ed)
function(diagham_geometry_test name case)
    cmake_parse_arguments(ARG "" "PROGRAM;OUTPUT;TOL" "ARGS;LABELS" ${ARGN})
    if(NOT ARG_TOL)
        set(ARG_TOL 1e-9)
    endif()
    diagham_physics_test(${name}
        PROGRAM ${ARG_PROGRAM}
        ARGS ${ARG_ARGS} --full-diag 100000
        OUTPUT "${ARG_OUTPUT}"
        CHECK spectrum @OUTPUT@ -1 ${geo_data}/${case}_spectrum.dat -1 ${ARG_TOL}
        LABELS fqhe ed ${ARG_LABELS})
endfunction()
set(cyf2 FQHEOnCylinder_FQHECylinderFermionsTwoBodyGeneric)
if(DIAGHAM_USE_GSL)
    foreach(case "fermions_n4_norb9_v1v3_r1_ky0|1|0" "fermions_n4_norb9_v1v3_r2_ky0|2|0" "fermions_n4_norb9_v1v3_r0.5_ky0|0.5|0" "fermions_n4_norb9_v1v3_r1_ky2|1|4")
        string(REPLACE "|" ";" c "${case}")
        list(GET c 0 name)
        list(GET c 1 ratio)
        list(GET c 2 y)
        diagham_geometry_test(physics.fqhe.cylinder.ed.${name} ${name}
            PROGRAM ${cyf2} ARGS -p 4 -l 8 -r ${ratio} -y ${y} --nbr-ky 1 --interaction-file ${geo_data}/${name}_pp.dat --interaction-name ed
            OUTPUT "fermions_cylinder_ky_ed_n_4_2s_8_ratio_*.dat" LABELS cylinder)
    endforeach()
endif()
diagham_geometry_test(physics.fqhe.cylinder.ed.bosons_n4_norb7_delta_r1_ky0 bosons_n4_norb7_delta_r1_ky0
    PROGRAM ${cbd} ARGS -p 4 -l 6 -r 1 -y 0 --nbr-ky 1 OUTPUT "bosons_cylinder_ky_delta_n_4_2s_6_ratio_*.dat" LABELS cylinder)
diagham_geometry_test(physics.fqhe.cylinder.ed.bosons_n4_norb7_delta_r2_ky1 bosons_n4_norb7_delta_r2_ky1
    PROGRAM ${cbd} ARGS -p 4 -l 6 -r 2 -y 2 --nbr-ky 1 OUTPUT "bosons_cylinder_ky_delta_n_4_2s_6_ratio_*.dat" LABELS cylinder)
diagham_geometry_test(physics.fqhe.cylinder.ed.fermions_n4_norb9_laplacian_delta_r1_ky0 fermions_n4_norb9_laplacian_delta_r1_ky0
    PROGRAM FQHEOnCylinder_FQHECylinderFermionsLaplacianDelta ARGS -p 4 -l 8 -r 1 -y 0 --nbr-ky 1 OUTPUT "fermions_cylinder_ky_delta_n_4_2s_8_ratio_*.dat" LABELS cylinder)
# (FQHECylinderFermionsCoulomb, U33: needs GSL too, now fails clearly without it; its Coulomb
#  matrix elements on a cylinder are not yet in the oracle)
# disk: the Coulomb pseudopotentials have the closed form Gamma(m+1/2)/(2 m!)
diagham_physics_test(physics.fqhe.disk.coulomb_pseudopotentials.2s8
    PROGRAM FQHEOnDisk_FQHEDiskCoulombPseudopotentials
    ARGS -s 8
    OUTPUT "pseudopotential_disk_coulomb_l_0_2s_8.dat"
    CHECK line @OUTPUT@ Pseudopotentials ${geo_data}/disk_coulomb_2s8.txt 1e-13
    LABELS fqhe disk)
# U34: the disk pseudopotential Hamiltonian is wrong for N >= 3 (zero modes right, nonzero energies off);
# the exact spectra are in the reference files (two independent routes agree). WILL_FAIL until fixed.
foreach(case "bosons_n3_lz6_v0|3|6" "bosons_n4_lz12_v0_v2|4|12" "bosons_n3_lz9_v2|3|9")
    string(REPLACE "|" ";" c "${case}")
    list(GET c 0 name)
    list(GET c 1 n)
    list(GET c 2 lz)
    diagham_physics_test(knownbug.fqhe.disk.bosons_two_body_generic.${name}
        PROGRAM ${db2}
        ARGS -p ${n} --minimum-momentum ${lz} --maximum-momentum ${lz} --interaction-file ${geo_data}/${name}_pp.dat --interaction-name ed --full-diag 100000
        OUTPUT "bosons_disk_ed_n_${n}_lz_${lz}.dat"
        CHECK spectrum @OUTPUT@ -1 ${geo_data}/${name}_spectrum.dat -1 1e-9
        WILL_FAIL
        LABELS fqhe disk)
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

# Further band-structure models: one test per listed band, comparing the sorted column of the
# program's export with the numpy eigenvalues of the Bloch matrix (tests/oracles/fci_bands.py).
# name | program | dims | output glob | first data column | bands | extra options
# (2D programs write kx ky E_0 ... in *_tightbinding.dat; the 3D ones write kx ky kz E_0 ... in the
# plain .dat, because their *_tightbinding.dat only holds the kz = 0 slice, see U36)
function(diagham_band_case name program nx ny nz glob col0 bands)
    foreach(band ${bands})
        math(EXPR col "${band} + ${col0}")
        diagham_physics_test(physics.${program}.bands.${name}.band${band}
            PROGRAM ${program}
            ARGS -p 2 -x ${nx} -y ${ny} ${ARGN}
            OUTPUT "${glob}"
            CHECK spectrum @OUTPUT@ ${col} ${fci_bands_data}/${name}.dat ${col} ${fci_band_tol}
            LABELS fci bands)
    endforeach()
endfunction()

set(fci_export --singleparticle-spectrum --export-onebodytext)
set(fci_band_tol 1e-10)
diagham_band_case(haldane_3x3 FCI_FCIHaldaneModel 3 3 0
    "fermions_singleband_haldane_n_2_x_3_y_3_*_tightbinding.dat" 2 "0;1" --single-band ${fci_export})
diagham_band_case(haldane_4x3_t2_0.4_phi_0.9_mus_0.3 FCI_FCIHaldaneModel 4 3 0
    "fermions_singleband_haldane_n_2_x_4_y_3_*_tightbinding.dat" 2 "0;1"
    --single-band --t2 0.4 --phi 0.9 --mu-s 0.3 ${fci_export})
diagham_band_case(kagome_3x3 FCI_FCIKagomeLatticeModel 3 3 0
    "fermions_threeband_kagomelattice_n_2_x_3_y_3_tightbinding.dat" 2 "0;1;2" --three-bands ${fci_export})
diagham_band_case(kagome_4x3_t2_0.1_l1_0.5 FCI_FCIKagomeLatticeModel 4 3 0
    "fermions_threeband_kagomelattice_n_2_x_4_y_3_*tightbinding.dat" 2 "0;1;2"
    --three-bands --t2 0.1 --l1 0.5 --l2 0 ${fci_export})
diagham_band_case(zhangqi_3x3 FCI_FCIZhangQiLatticeModel 3 3 0
    "fermions_singleband_zhangqi_n_2_x_3_y_3_tightbinding.dat" 2 "0;1" ${fci_export})
diagham_band_case(zhangqi_4x3_theta_0.15_mus_0.2 FCI_FCIZhangQiLatticeModel 4 3 0
    "fermions_singleband_zhangqi_n_2_x_4_y_3_*tightbinding.dat" 2 "0;1" --theta 0.15 --mu-s 0.2 ${fci_export})
diagham_band_case(bhz_3x3 FTI_FQSH2DBHZModel 3 3 0
    "fermions_quantumspihall2d_BHZlattice_n_2_x_3_y_3_*_tightbinding.dat" 2 "0;1;2;3" --nbr-bands 2 ${fci_export})
# the 3D plain .dat carries six significant digits only
set(fci_band_tol 1e-5)
diagham_band_case(simple_ti3d_3x3x2_m1.5 FTI_FTI3DSimpleTI 3 3 2
    "fermions_quantumspinhall3d_simpleti_fourbands_n_2_x_3_y_3_z_2_u_*_gz_0.000000.dat" 3 "0;3"
    -z 2 --four-bands --mass 1.5 ${fci_export})
diagham_band_case(simple_ti3d_3x2x3_m2.5 FTI_FTI3DSimpleTI 3 2 3
    "fermions_quantumspinhall3d_simpleti_fourbands_n_2_x_3_y_2_z_3_u_*_gz_0.000000.dat" 3 "0;3"
    -z 3 --four-bands --mass 2.5 ${fci_export})

# U36: FTI3DHopf writes its band file and then segfaults (exit 139, every size and --lambda tried;
# the same on a pristine autotools build of r4493). WILL_FAIL until it is fixed.
diagham_physics_test(knownbug.fti.hopf3d.segfault_after_band_file
    PROGRAM FTI_FTI3DHopf
    ARGS -p 2 -x 3 -y 3 -z 3 --singleparticle-spectrum
    OUTPUT "fermions_singleband_hopf_p_2_x_3_y_3_z_3_l_1.dat"
    CHECK count @OUTPUT@ -1 0 1e-10 0
    WILL_FAIL LABELS fci)

# U37: for the 3D tight-binding programs --export-onebodytext writes the two-dimensional layout
# (header "kx ky E_0 ...") and only the kz = 0 slice: 9 rows for a 3x3x2 grid instead of 18.
diagham_physics_test(knownbug.fti.simple_ti3d.export_onebodytext_kz0_slice_only
    PROGRAM FTI_FTI3DSimpleTI
    ARGS -p 2 -x 3 -y 3 -z 2 --four-bands --mass 1.5 ${fci_export}
    OUTPUT "fermions_quantumspinhall3d_simpleti_fourbands_n_2_x_3_y_3_z_2_*_tightbinding.dat"
    CHECK spectrum @OUTPUT@ 2 ${fci_bands_data}/simple_ti3d_3x3x2_m1.5.dat 3 1e-5
    WILL_FAIL LABELS fci)

# Chern number of the lowest Haldane band (--singleparticle-chernnumber prints it on a line of its
# own in the log): the Fukui-Hatsugai integer with DiagHam's orientation, at N = 64 where the
# program's small-angle formula is within 0.008 of it.
foreach(case "haldane_phi0.5_mus0.3|0.5|0.3" "haldane_phi-0.5_mus0.3|-0.5|0.3" "haldane_phi0.5_mus6|0.5|6")
    string(REPLACE "|" ";" c "${case}")
    list(GET c 0 name)
    list(GET c 1 phi)
    list(GET c 2 mus)
    diagham_physics_test(physics.fci.chern.${name}
        PROGRAM FCI_FCIHaldaneModel
        ARGS -p 2 -x 64 -y 64 --single-band --phi ${phi} --mu-s ${mus} --singleparticle-spectrum --singleparticle-chernnumber
        OUTPUT "program.log"
        CHECK numbers @OUTPUT@ ${fci_bands_data}/chern_${name}.txt 0.02
        LABELS fci bands)
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

# --- 1d. FCI many-body spectra ---------------------------------------------------------------------
# tests/oracles/fci_manybody.py: dense ED of the band-projected interacting problem (numpy), sector by sector
# in the total momentum, and the exact atomic limit. The programs' full spectra (--full-momentum) must match as
# sorted lists. Physics of the checkerboard checks: A sites at R, B sites at R + (1/2, 1/2), U on the four
# nearest-neighbour A-B bonds, projected onto the lowest band.
set(fci_mb_data ${DIAGHAM_TEST_DATA}/fci_manybody)
set(fci_al FCI_FCIAtomicLimitLatticeModel)

# flat band (eps = 0, U fixed at 1 by the program): only the projected interaction remains
diagham_physics_test(physics.fci.checkerboard.manybody.3x3_n2_flat
    PROGRAM ${fci_cb}
    ARGS -p 2 -x 3 -y 3 --single-band --flat-band --full-diag 5000 --full-momentum
    OUTPUT "fermions_singleband_checkerboardlattice_n_2_ns_18_x_3_y_3_*_gx_*.dat"
    CHECK spectrum @OUTPUT@ 2 ${fci_mb_data}/checkerboard_3x3_n2_flat.dat 2 1e-9
    LABELS fci checkerboard manybody)
diagham_physics_test(physics.fci.checkerboard.manybody.3x3_n3_flat
    PROGRAM ${fci_cb}
    ARGS -p 3 -x 3 -y 3 --single-band --flat-band --full-diag 5000 --full-momentum
    OUTPUT "fermions_singleband_checkerboardlattice_n_3_ns_18_x_3_y_3_*_gx_*.dat"
    CHECK spectrum @OUTPUT@ 2 ${fci_mb_data}/checkerboard_3x3_n3_flat.dat 2 1e-9
    LABELS fci checkerboard manybody)
diagham_physics_test(physics.fci.checkerboard.manybody.4x3_n2_flat_t2_0.2_tpp_0.1
    PROGRAM ${fci_cb}
    ARGS -p 2 -x 4 -y 3 --single-band --flat-band --t2 0.2 --tpp 0.1 --full-diag 5000 --full-momentum
    OUTPUT "fermions_singleband_checkerboardlattice_n_2_ns_24_x_4_y_3_*_gx_*.dat"
    CHECK spectrum @OUTPUT@ 2 ${fci_mb_data}/checkerboard_4x3_n2_flat_t2_0.2_tpp_0.1.dat 2 1e-9
    LABELS fci checkerboard manybody)

# atomic limit: bosons pay 4 U per same-cell pair (the convention of the program), fermions are inert
diagham_physics_test(physics.fci.atomiclimit.bosons_3x3_n2_u1.5
    PROGRAM ${fci_al}
    ARGS -p 2 -x 3 -y 3 --boson --single-band --u-potential 1.5 --full-diag 5000
    OUTPUT "bosons_singleband_atomiclimit_band_0_n_2_x_3_y_3_*_gx_*.dat"
    CHECK spectrum @OUTPUT@ 2 ${fci_mb_data}/atomic_boson_3x3_n2_u1.5.dat 2 1e-9
    LABELS fci atomiclimit manybody)
diagham_physics_test(physics.fci.atomiclimit.bosons_3x3_n3_u1
    PROGRAM ${fci_al}
    ARGS -p 3 -x 3 -y 3 --boson --single-band --u-potential 1 --full-diag 5000
    OUTPUT "bosons_singleband_atomiclimit_band_0_n_3_x_3_y_3_*_gx_*.dat"
    CHECK spectrum @OUTPUT@ 2 ${fci_mb_data}/atomic_boson_3x3_n3_u1.dat 2 1e-9
    LABELS fci atomiclimit manybody)
diagham_physics_test(physics.fci.atomiclimit.bosons_4x2_n2_u1
    PROGRAM ${fci_al}
    ARGS -p 2 -x 4 -y 2 --boson --single-band --u-potential 1 --full-diag 5000
    OUTPUT "bosons_singleband_atomiclimit_band_0_n_2_x_4_y_2_*_gx_*.dat"
    CHECK spectrum @OUTPUT@ 2 ${fci_mb_data}/atomic_boson_4x2_n2_u1.dat 2 1e-9
    LABELS fci atomiclimit manybody)
diagham_physics_test(physics.fci.atomiclimit.fermions_3x3_n2_inert
    PROGRAM ${fci_al}
    ARGS -p 2 -x 3 -y 3 --single-band --u-potential 1 --v-potential 1 --full-diag 5000
    OUTPUT "fermions_singleband_atomiclimit_band_0_n_2_x_3_y_3_*_gx_*.dat"
    CHECK spectrum @OUTPUT@ 2 ${fci_mb_data}/atomic_fermion_3x3_n2.dat 2 1e-9
    LABELS fci atomiclimit manybody)

# U39: the one-body energy of the single-band many-body Hamiltonian is wrong. Non-flat checkerboard: the
# band energy is overwritten by an interaction-derived term (EvaluateInteractionFactors, second assignment
# of OneBodyInteractionFactors), so U = 0 does not give the band-energy sum; Haldane and Kagome: with no
# interaction, N = 1 gives exactly half of the band energy. WILL_FAIL until fixed.
diagham_physics_test(knownbug.fci.checkerboard.manybody.3x3_n2_u0_band_energy_overwritten
    PROGRAM ${fci_cb}
    ARGS -p 2 -x 3 -y 3 --single-band --u-potential 0 --full-diag 5000 --full-momentum
    OUTPUT "fermions_singleband_checkerboardlattice_n_2_ns_18_x_3_y_3_*_gx_*.dat"
    CHECK spectrum @OUTPUT@ 2 ${fci_mb_data}/checkerboard_3x3_n2_u0.dat 2 1e-9
    WILL_FAIL LABELS fci checkerboard manybody known-bug)
diagham_physics_test(knownbug.fci.checkerboard.manybody.4x3_n2_u1_band_energy_overwritten
    PROGRAM ${fci_cb}
    ARGS -p 2 -x 4 -y 3 --single-band --t2 0.2 --tpp 0.1 --u-potential 1 --full-diag 5000 --full-momentum
    OUTPUT "fermions_singleband_checkerboardlattice_n_2_ns_24_x_4_y_3_*_gx_*.dat"
    CHECK spectrum @OUTPUT@ 2 ${fci_mb_data}/checkerboard_4x3_n2_u1_t2_0.2_tpp_0.1.dat 2 1e-9
    WILL_FAIL LABELS fci checkerboard manybody known-bug)
diagham_physics_test(knownbug.fci.haldane.manybody.3x3_n1_band_energy_halved
    PROGRAM FCI_FCIHaldaneModel
    ARGS -p 1 -x 3 -y 3 --single-band --full-diag 5000 --full-momentum
    OUTPUT "fermions_singleband_haldane_n_1_x_3_y_3_*_gy_0.000000.dat"
    CHECK min-abs @OUTPUT@ 2 -4.5350243553983 1e-9
    WILL_FAIL LABELS fci haldane manybody known-bug)
diagham_physics_test(knownbug.fci.kagome.manybody.3x3_n1_band_energy_halved
    PROGRAM FCI_FCIKagomeLatticeModel
    ARGS -p 1 -x 3 -y 3 --full-diag 5000 --full-momentum
    OUTPUT "fermions_singleband_kagomelattice_n_1_x_3_y_3_*_gx_*.dat"
    CHECK min-abs @OUTPUT@ 2 -2.8 1e-9
    WILL_FAIL LABELS fci kagome manybody known-bug)

# --- 2b. spin models from explicit operator lists ------------------------------------------------
# tests/oracles/spin_models.py: dense ED of a model given as operator terms, with Sz and momentum
# sectors; each program's convention (S or Pauli, signs, boundary terms) is recorded in the oracle.
set(sm_data ${DIAGHAM_TEST_DATA}/spin_models)
function(diagham_spin_model_test name case)
    cmake_parse_arguments(ARG "WILL_FAIL" "PROGRAM;OUTPUT;TOL" "ARGS;LABELS" ${ARGN})
    if(NOT ARG_TOL)
        set(ARG_TOL 1e-9)
    endif()
    set(wf "")
    if(ARG_WILL_FAIL)
        set(wf WILL_FAIL)
    endif()
    diagham_physics_test(${name}
        PROGRAM ${ARG_PROGRAM}
        ARGS ${ARG_ARGS} --full-diag 100000
        OUTPUT "${ARG_OUTPUT}"
        CHECK spectrum @OUTPUT@ -1 ${sm_data}/${case}_spectrum.dat -1 ${ARG_TOL}
        ${wf}
        LABELS spin ed ${ARG_LABELS})
endfunction()
diagham_spin_model_test(physics.spin.ed.xyz_open_L8_h0.2 xyz_open_L8_h0.2
    PROGRAM Spin_SpinChainXYZ ARGS -p 8 -x 1 -y 0.6 -z 0.3 -f 0.2 -b 0 --no-parity OUTPUT "spin_1_2_x_1.000000_y_0.600000_z_0.300000_h_0.200000_b_0_n_8.dat")
diagham_spin_model_test(physics.spin.ed.xyz_periodic_L8_parity_boundary xyz_periodic_L8
    PROGRAM Spin_SpinChainXYZ ARGS -p 8 -x 1 -y 0.6 -z 0.3 -b 1 --no-parity OUTPUT "spin_1_2_x_1.000000_y_0.600000_z_0.300000_b_1_n_8.dat")
diagham_spin_model_test(physics.spin.ed.fullgeneric_open_L6_three_fields fullgeneric_open_L6_fields
    PROGRAM Spin_FullGenericOpenSpinChain ARGS -s 1 -p 6 -x 1 -y 0.6 -z 0.3 --hx-value 0.2 --hy-value 0.1 --hz-value 0.3 OUTPUT "spin_1_2_openchain_n_6_jx_*.dat")
diagham_spin_model_test(physics.spin.ed.heisenberg2d_3x3_jz0.7 heisenberg2d_3x3_jz0.7_sz1
    PROGRAM Spin_TwoDimensionalHeisenbergModel ARGS -x 3 -y 3 --j-value 1 --jz-value 0.7 --nbr-sz 1 --disable-momentum --disable-inversion --disable-szsymmetry OUTPUT "spin_1_2_2dheisenberg_n_9_x_3_y_3_*.dat" LABELS lattice)
foreach(k 0 1 3)
    diagham_spin_model_test(physics.spin.ed.j1j2_L8_sz0_k${k} j1j2_L8_sz0_k${k}
        PROGRAM Spin_PeriodicSpinChainJ1J2 ARGS -s 1 -p 8 -1 1 -2 0.4 --nbr-sz 1 --momentum ${k} --disable-szsymmetry --disable-inversionsymmetry OUTPUT "spin_1_2_periodicj1j2_j1_1.000000_j2_0.400000_n_8.dat")
endforeach()
diagham_spin_model_test(physics.spin.ed.doubletriangle_L8_sz0 doubletriangle_L8_sz0
    PROGRAM Spin_DoubleTriangleSpinChain ARGS -s 1 -p 8 -j 1 -g 0.5 --djz1-value 0.2 --djz2-value 0.1 --nbr-sz 1 --no-translations OUTPUT "spin_1_2_doubletrianglechain_*_n_8.dat")
foreach(k 0 2)
    diagham_spin_model_test(physics.spin.ed.aklt_p3p4_spin2_L4_sz0_k${k} aklt_p3p4_spin2_L4_sz0_k${k}
        PROGRAM Spin_PeriodicSpinChainGeneralizedAKLT ARGS -s 4 -p 4 --nbr-sz 1 --momentum ${k} --disable-szsymmetry --disable-inversionsymmetry OUTPUT "spin_2_periodicaklt_*_n_4.dat" TOL 1e-7)
endforeach()
diagham_spin_model_test(physics.spin.ed.tfim2d_2x3_periodic tfim2d_2x3_periodic
    PROGRAM Spin_TwoDimensionalTransverseFieldIsingModel ARGS -x 2 -y 3 --jz-value 1 --hx-value 0.7 --hz-value 0.2 --use-periodic --disable-momentum --disable-inversion OUTPUT "spin_1_2_ising_transversefield_closed_n_6_*.dat" LABELS lattice)
diagham_spin_model_test(physics.spin.ed.tfim2d_3x3_periodic tfim2d_3x3_periodic
    PROGRAM Spin_TwoDimensionalTransverseFieldIsingModel ARGS -x 3 -y 3 --jz-value 1 --hx-value 0.7 --hz-value 0.2 --use-periodic --disable-momentum --disable-inversion OUTPUT "spin_1_2_ising_transversefield_closed_n_9_*.dat" LABELS lattice)
# U35: without --use-periodic the program still adds the wrap-around Ising bonds (its "open" spectrum equals the periodic one)
diagham_spin_model_test(knownbug.spin.tfim2d.open_boundaries_ignored tfim2d_2x3_open
    PROGRAM Spin_TwoDimensionalTransverseFieldIsingModel ARGS -x 2 -y 3 --jz-value 1 --hx-value 0.7 --hz-value 0.2 OUTPUT "spin_1_2_ising_transversefield_open_n_6_*.dat" WILL_FAIL LABELS lattice)
diagham_spin_model_test(physics.spin.ed.obrienfendley_spin1_L4_sz0 obrienfendley_spin1_L4_sz0
    PROGRAM Spin_PeriodicSpinChainOBrienFendley ARGS -s 2 -p 4 --nbr-sz 1 --disable-szsymmetry --disable-inversionsymmetry --disable-momentum OUTPUT "spin_1_periodicobrienfendley_nomomentum_n_4.dat")
diagham_spin_model_test(physics.spin.ed.potts3_open_L6 potts3_open_L6
    PROGRAM Spin_Potts3ChainModel ARGS -p 6 -j 1 -f 0.5 --nbr-q 0 OUTPUT "potts3_openchain_*_n_6.dat")
if(DIAGHAM_USE_GSL)
    diagham_spin_model_test(physics.spin.ed.disorder_open_L6_sz0 disorder_open_L6_sz0
        PROGRAM Spin_GenericOpenSpinChainWithDisorder ARGS -s 1 -p 6 -j 1 -z 0.2 --disorder-file ${sm_data}/disorder_open_L6.dat --nbr-sz 1 OUTPUT "spin_1_2_openchain_n_6.dat")
endif()

# --- 2c. lattice fermions (Hubbard-family programs) --------------------------------------------------
# tests/oracles/lattice_fermions.py: Fock-basis ED of spinful fermions on a Hermitian hopping matrix with an
# on-site U; the Haldane cluster is the inverse Fourier transform of DiagHam's own Bloch Hamiltonian
# (its --phi is in radians, its --mu-s sits on the A sublattice). Tolerance 1e-7: these programs
# reproduce the oracle to about 1e-8.
set(lf_data ${DIAGHAM_TEST_DATA}/lattice_fermions)
function(diagham_lattice_test name case)
    cmake_parse_arguments(ARG "" "PROGRAM;OUTPUT" "ARGS" ${ARGN})
    diagham_physics_test(${name}
        PROGRAM ${ARG_PROGRAM}
        ARGS ${ARG_ARGS} --full-diag 100000
        OUTPUT "${ARG_OUTPUT}"
        CHECK spectrum @OUTPUT@ -1 ${lf_data}/${case}_spectrum.dat -1 1e-7
        LABELS fti hubbard ed)
endfunction()
set(hald HubbardModels_HubbardHaldaneLatticeModel)
diagham_lattice_test(physics.hubbard.ed.square_2x3_n6_u2 square_2x3_n6_u2_sz0
    PROGRAM ${hubbard} ARGS -p 6 -x 2 -y 3 --u-potential 2 OUTPUT "fermions_hubbard_square_x_2_y_3_n_6_ns_6_*_u_2.000000_sz_0.dat")
diagham_lattice_test(physics.hubbard.ed.haldane_3x3_n2_u0 haldane_3x3_n2_u0_sz0
    PROGRAM ${hald} ARGS -p 2 -x 3 -y 3 --u-potential 0 --only-sz 0 OUTPUT "fermions_hubbard_haldane_x_3_y_3_n_2_ns_18_*_phi_0.333333_sz_0.dat")
diagham_lattice_test(physics.hubbard.ed.haldane_3x3_n2_u2 haldane_3x3_n2_u2_sz0
    PROGRAM ${hald} ARGS -p 2 -x 3 -y 3 --u-potential 2 --only-sz 0 OUTPUT "fermions_hubbard_haldane_x_3_y_3_n_2_ns_18_*_u_2.000000_sz_0.dat")
diagham_lattice_test(physics.hubbard.ed.haldane_3x3_n2_u2_mus0.1 haldane_3x3_n2_u2_mus0.1_sz0
    PROGRAM ${hald} ARGS -p 2 -x 3 -y 3 --u-potential 2 --mu-s 0.1 --only-sz 0 OUTPUT "fermions_hubbard_haldane_x_3_y_3_n_2_ns_18_*_u_2.000000_sz_0.dat")
diagham_lattice_test(physics.hubbard.ed.ssh_4cells_n4_delta0.3 ssh_4cells_n4_delta0.3
    PROGRAM HubbardModels_HubbardSSHModel ARGS -p 4 -x 4 --delta 0.3 --no-translation OUTPUT "fermions_ssh_x_4_n_4_ns_8_d_0.3.dat")
diagham_spin_model_test(physics.spin.ed.kitaev_heisenberg_honeycomb_2x2 kitaev_heisenberg_honeycomb_2x2_j1_jk0.5
    PROGRAM HubbardModels_HubbardExtendedKitaevHeisenbergHoneycombModel ARGS -x 2 -y 2 --spin --j1 1 --jK 0.5 OUTPUT "spin_kitaev_heisenberg_honeycomb_x_2_y_2_ns_8_j1_*_hz_0.000000.dat" TOL 1e-7 LABELS lattice)

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
    add_test(NAME selftest.sphere_nbody_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/sphere_nbody.py --check ${nbody_data})
    add_test(NAME selftest.geometry_ed_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/geometry_ed.py --check ${geo_data})
    add_test(NAME selftest.spin_models_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/spin_models.py --check ${sm_data})
    add_test(NAME selftest.lattice_fermions_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/lattice_fermions.py --check ${lf_data})
    add_test(NAME selftest.torus_ed_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/torus_ed.py --check ${torus_ed_data})
    add_test(NAME selftest.fci_bands_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/fci_bands.py --check ${fci_bands_data})
    add_test(NAME selftest.fci_manybody_oracle
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/oracles/fci_manybody.py --check ${fci_mb_data})
    set_tests_properties(selftest.sphere_ed_oracle selftest.spin_ed_oracle selftest.torus_ed_oracle
        selftest.fci_bands_oracle selftest.fci_manybody_oracle selftest.sphere_nbody_oracle selftest.geometry_ed_oracle
        selftest.spin_models_oracle selftest.lattice_fermions_oracle PROPERTIES LABELS "selftest;python" TIMEOUT 900)
endif()
