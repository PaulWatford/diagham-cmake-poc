#!/bin/bash
# Reproduce the registered defects on a pristine Subversion export of DiagHam trunk, built with
# DiagHam's own autotools (bootstrap.sh, configure, make): no CMake, no converted tree, nothing
# of this repository in the way except the reference spectra and the checker.
#
#   scripts_cmake/reproduce_on_pristine_trunk.sh [BUILD_DIR] [REV] [WORK]
#
#   BUILD_DIR  a CMake build of this repository (for tests/check_spectrum), default build/default
#   REV        trunk revision to export, default scripts_cmake/upstream-revision.txt
#   WORK       where to export and build, default $HOME/diagham/pristine-r<REV>
#
# Steps: (1) svn info on trunk (is REV still the head?), (2) svn export -r REV, (3) diff the export
# against this repository's upstream branch (the conversion is expected to differ only in empty
# directories and keyword expansions), (4) bootstrap + configure --enable-fqhe --enable-fti (no GSL,
# the default of both build systems) + make -k (upstream's own build failures U01-U03 stop a plain
# make), (5) run each reproducer and a positive control, print what came out next to what the
# register says. Needs svn, autotools, a C++ compiler; takes about half an hour for the build.
set -u
REPO=$(cd "$(dirname "$0")/.." && pwd)
BUILD=$(cd "${1:-$REPO/build/default}" && pwd)
REV=${2:-$(cat "$REPO/scripts_cmake/upstream-revision.txt")}
WORK=${3:-$HOME/diagham/pristine-r$REV}
URL=https://www.nick-ux.org/diagham/svn/DiagHam/trunk
CHECK=$BUILD/tests/check_spectrum
[ -x "$CHECK" ] || { echo "no check_spectrum in $BUILD; build this repository first" >&2; exit 2; }

echo "== 1. trunk head"
svn info "$URL" | grep -E "^(Revision|Last Changed Rev|Last Changed Date)"
echo "   this repository's upstream branch: r$REV"

FRESH=0
if [ ! -d "$WORK" ]; then
    echo "== 2. svn export -r $REV -> $WORK"
    svn export -q -r "$REV" "$URL" "$WORK" || exit 1
    FRESH=1
else
    echo "== 2. reusing $WORK"
fi

if [ "$FRESH" -eq 1 ]; then
    echo "== 3. export vs upstream branch (empty directories and keyword expansions are expected)"
    TREE=$(mktemp -d)
    (cd "$REPO" && git archive upstream | tar -x -C "$TREE")
    diff -rq "$WORK" "$TREE" | grep -vE "^Only in $TREE" | head -20
    rm -rf "$TREE"
else
    echo "== 3. tree diff skipped: $WORK was built in place (remove it for a fresh export and diff)"
fi

if [ ! -x "$WORK/FQHE/src/Programs/FQHEOnDisk/FQHEDiskBosonsTwoBodyGeneric" ]; then
    echo "== 4. autotools build (bootstrap, configure --enable-fqhe --enable-fti, make -k -j)"
    (cd "$WORK" && sh bootstrap.sh > bootstrap.log 2>&1 && ./configure --enable-fqhe --enable-fti > configure.log 2>&1 && make -k -j"$(nproc)" > make.log 2>&1)
    echo "   configure: $(grep -c "^checking" "$WORK/configure.log") checks; GSL: $(grep "using gsl" "$WORK/configure.log" | sed 's/.*\.\.\. //'); make errors (expected: the ungated LAPACK calls, U01-U03): $(grep -c "error:" "$WORK/make.log")"
else
    echo "== 4. reusing the build in $WORK"
fi

D=$WORK/FQHE/src/Programs/FQHEOnDisk; S=$WORK/FQHE/src/Programs/FQHEOnSphere
C=$WORK/FQHE/src/Programs/FQHEOnCylinder; H=$WORK/FTI/src/Programs/HubbardModels
DATA=$REPO/tests/data
R=$(mktemp -d); cd "$R"
for b in $D/FQHEDiskBosonsTwoBodyGeneric $D/FQHEDiskFermionsTwoBodyGeneric $S/FQHESphereFermionsWithSpin $S/FQHESphereFermionsTwoBodyGeneric $C/FQHECylinderFermionsCoulomb $H/HubbardSquareLatticeModel; do
    [ -x "$b" ] || { echo "missing $b (build failed there)"; exit 1; }
done

echo "== 5. reproducers on the pristine build"
echo "-- positive control: sphere Coulomb ED N=4 2S=9 Lz=0 against tests/data/sphere_ed (expect PASS)"
$S/FQHESphereFermionsTwoBodyGeneric -p 4 -l 9 --interaction-file "$DATA/sphere_ed/fermions_n4_2s9_coulomb_pp.dat" --interaction-name c --initial-lz 0 --nbr-lz 1 --full-diag 100000 > /dev/null 2>&1
$CHECK spectrum fermions_c_n_4_2s_9_lz*.dat -1 "$DATA/sphere_ed/fermions_n4_2s9_coulomb_spectrum.dat" -1 1e-9 | tail -1

echo "-- U34: disk bosons N=3 Lz=6 V0; exact 0 1.5 1.875 2.0625 2.25 3 6; the defect gives 0 1.5951 1.9825 2.0000 2.1174 2.8833 5.9842"
$D/FQHEDiskBosonsTwoBodyGeneric -p 3 --minimum-momentum 6 --maximum-momentum 6 --interaction-file "$DATA/geometry_ed/bosons_n3_lz6_v0_pp.dat" --interaction-name ed --full-diag 100000 > /dev/null 2>&1
grep -v '^#' bosons_disk_ed_n_3_lz_6.dat | awk '{printf "   %.6f", $2} END {print ""}'
echo "   N=2 Lz=5 V0 (expect a single level at 2: the pair states are right):"
$D/FQHEDiskBosonsTwoBodyGeneric -p 2 --minimum-momentum 5 --maximum-momentum 5 --interaction-file "$DATA/geometry_ed/bosons_n3_lz6_v0_pp.dat" --interaction-name ed --full-diag 100000 > /dev/null 2>&1
grep -v '^#' bosons_disk_ed_n_2_lz_5.dat | awk '{printf "   %.6f", $2} END {print ""}'

echo "-- U29: disk fermions N=3 Lz=6 (expect a hang: rc 124 after 20 s) and N=2 Lz=5 (expect a segfault: rc 139)"
printf 'Pseudopotentials = 0 1 0 0 0 0 0\n' > v1.dat
timeout 20 $D/FQHEDiskFermionsTwoBodyGeneric -p 3 --minimum-momentum 6 --maximum-momentum 6 --interaction-file v1.dat --interaction-name v1 --full-diag 100000 > f3.log 2>&1; echo "   N=3 rc=$?, last line: $(tail -1 f3.log)"
timeout 20 $D/FQHEDiskFermionsTwoBodyGeneric -p 2 --minimum-momentum 5 --maximum-momentum 5 --interaction-file v1.dat --interaction-name v1 --full-diag 100000 > f2.log 2>&1; echo "   N=2 Lz=5 rc=$?"

echo "-- U32: sphere fermions with spin, both symmetrised bases (expect dimension 0; one flag alone 25 and 28)"
SP=$DATA/sphere_ed/fermions_su2_n4_2s6_coulomb_sz0_pp.dat
for flags in "--lzsymmetrized-basis --szsymmetrized-basis" "--szsymmetrized-basis" "--lzsymmetrized-basis"; do
    echo "   $flags: $($S/FQHESphereFermionsWithSpin -p 4 -l 6 -s 0 --interaction-file "$SP" --interaction-name sp --initial-lz 0 --nbr-lz 1 --full-diag 100000 $flags 2>&1 | grep -i "dimension" | head -1)"
done

echo "-- U30: Hubbard 2x4 U=0, plain Lanczos -n 1 (exact ground energy -16 in sectors (0,0) and (1,2); the defect gives -17.87 and -18.41)"
$H/HubbardSquareLatticeModel -p 8 -x 2 -y 4 --u-potential 0 -n 1 --full-diag 100 > h.log 2>&1
grep -v '^#' fermions_hubbard_square_x_2_y_4_n_8_ns_8_t_1.000000_tp_0.000000_sz_0.dat | awk '{printf "   (%s,%s) %.4f", $1, $2, $4} END {print ""}'

echo "-- U33: cylinder fermions Coulomb without GSL (expect a segfault, rc 139)"
$C/FQHECylinderFermionsCoulomb -p 4 -l 8 -r 1 -y 0 --nbr-ky 1 --full-diag 100000 > cc.log 2>&1; echo "   rc=$?, last line: $(tail -1 cc.log)"
echo "== done; work directory $R"
