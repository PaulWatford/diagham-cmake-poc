#!/usr/bin/env bash
# Self-test of tests/check_claims.py (selftest.claims): one synthetic test per rule, in a throw-away CMake project
# whose tests are listed but never run. A clean project must pass; each broken test must be reported by its rule,
# and its corrected twin must not be. Run by scripts_cmake/ci.sh; needs cmake, ctest and python3.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
fails=0
ok()  { echo "ok   $1"; }
bad() { echo "FAIL $1" >&2; fails=$((fails + 1)); }

# A source tree with the three files the rules read.
make_tree() {
    mkdir -p "$1/tests/oracles" "$1/docs/reference"
    printf '# vocabulary\nphysics\ncrosscheck\npython\nselftest\nopen-finding\nspin\n' > "$1/tests/labels.txt"
    printf '| U1 | `openfinding.cited.exactly` and `openfinding.cited_by_prefix.*` |\n' \
        > "$1/docs/reference/open-findings.md"
    echo 'print(0)' > "$1/tests/oracles/right.py"
    echo 'print(0)' > "$1/tests/oracles/other.py"
}

# project <dir> <cmake body>: configure a project whose CMakeLists.txt is the body.
project() {
    mkdir -p "$1/src"
    printf 'cmake_minimum_required(VERSION 3.16)\nproject(claims_selftest NONE)\nenable_testing()\n%s\n' "$2" \
        > "$1/src/CMakeLists.txt"
    cmake -S "$1/src" -B "$1/build" > /dev/null
}

# t <name> <check args> <properties...>: a check_spectrum-style test (listed, never run)
t() {
    local name="$1" check="$2"; shift 2
    printf 'add_test(NAME %s COMMAND ${CMAKE_COMMAND} "-DCHECK_ARGS=%s" -P run.cmake)\n' "$name" "$check"
    printf 'set_tests_properties(%s PROPERTIES %s)\n' "$name" "$*"
}

tree="$work/tree"
make_tree "$tree"
O='KNOWN_FROM "oracle tests/oracles/right.py"'

# --- clean project: every rule satisfied -----------------------------------------------------------------------
clean="$(t physics.clean 'spectrum|@OUTPUT@|-1|ref.dat|-1|rel:1e-12' "LABELS physics $O")
$(t physics.min_four_ulp 'min|@OUTPUT@|-1|-1.5|4' "LABELS physics $O")
$(t physics.loose_with_reason 'spectrum|@OUTPUT@|-1|ref.dat|-1|1e-6' "LABELS physics $O TOLERANCE_REASON \"six digits\"")
$(t openfinding.cited.exactly 'min-abs|@OUTPUT@|-1|0|1e-12' 'LABELS open-finding WILL_FAIL ON')
$(t openfinding.cited_by_prefix.case 'min-abs|@OUTPUT@|-1|0|1e-12' 'LABELS open-finding WILL_FAIL ON')
$(t selftest.fails_on_purpose 'min-abs|@OUTPUT@|-1|0|1e-12' 'LABELS selftest WILL_FAIL ON')
add_test(NAME physics.python_right COMMAND python3 $tree/tests/oracles/right.py)
set_tests_properties(physics.python_right PROPERTIES LABELS \"physics;python\" $O)"
project "$work/clean" "$clean"
if out="$(python3 "$here/check_claims.py" "$work/clean/build" --source-dir "$tree")"; then
    ok "clean project passes"
else
    bad "clean project fails:"; echo "$out" >&2
fi

# --- one broken test per rule, each with a corrected twin above ---------------------------------------------------
broken="$clean
$(t physics.no_source 'spectrum|@OUTPUT@|-1|ref.dat|-1|rel:1e-12' 'LABELS physics')
$(t physics.unknown_kind 'spectrum|@OUTPUT@|-1|ref.dat|-1|rel:1e-12' 'LABELS physics KNOWN_FROM "it looked right"')
$(t physics.missing_oracle 'spectrum|@OUTPUT@|-1|ref.dat|-1|rel:1e-12' 'LABELS physics KNOWN_FROM "oracle tests/oracles/absent.py"')
$(t physics.loose_no_reason 'spectrum|@OUTPUT@|-1|ref.dat|-1|1e-6' "LABELS physics $O")
$(t physics.loose_relative 'lowest|@OUTPUT@|-1|ref.dat|-1|4|rel:1e-10' "LABELS physics $O")
$(t physics.min_five_ulp 'min|@OUTPUT@|-1|-1.5|5' "LABELS physics $O")
$(t physics.misspelt_label 'spectrum|@OUTPUT@|-1|ref.dat|-1|rel:1e-12' "LABELS \"physics;phyiscs\" $O")
$(t physics.will_fail_unlabelled 'spectrum|@OUTPUT@|-1|ref.dat|-1|rel:1e-12' "LABELS physics $O WILL_FAIL ON")
$(t openfinding.cited.exactly_without_will_fail 'min-abs|@OUTPUT@|-1|0|1e-12' 'LABELS open-finding')
$(t openfinding.not_in_register 'min-abs|@OUTPUT@|-1|0|1e-12' 'LABELS open-finding WILL_FAIL ON')
add_test(NAME physics.python_wrong COMMAND python3 $tree/tests/oracles/other.py)
set_tests_properties(physics.python_wrong PROPERTIES LABELS \"physics;python\" $O)"
project "$work/broken" "$broken"
set +e
out="$(python3 "$here/check_claims.py" "$work/broken/build" --source-dir "$tree")"; status=$?
out_req="$(KNOWN_FROM_REQUIRED=1 python3 "$here/check_claims.py" "$work/broken/build" --source-dir "$tree")"
set -e
[ "$status" -eq 1 ] && ok "broken project exits 1" || bad "broken project exits $status, not 1"

expect() {  # expect <rule> <test> <text>: a FAIL line for this test under this rule
    if grep -q "^FAIL $1: $2: " <<< "$out"; then ok "$1 rule reports $2"; else bad "$1 rule misses $2"; fi
}
expect known-from physics.unknown_kind
expect known-from physics.missing_oracle
expect known-from physics.python_wrong
expect tolerance physics.loose_no_reason
expect tolerance physics.loose_relative
expect tolerance physics.min_five_ulp
expect labels physics.misspelt_label
expect open-finding physics.will_fail_unlabelled
expect open-finding openfinding.cited.exactly_without_will_fail
expect register openfinding.not_in_register
# report-only by default, a failure when required
grep -q "^NOTE known-from: 1 .*physics.no_source" <<< "$out" && ok "missing KNOWN_FROM is a note by default" \
    || bad "missing KNOWN_FROM is not reported as a note"
grep -q "^FAIL known-from: 1 .*physics.no_source" <<< "$out_req" && ok "missing KNOWN_FROM fails when required" \
    || bad "missing KNOWN_FROM does not fail with KNOWN_FROM_REQUIRED=1"
# the corrected twins stay quiet
for twin in physics.clean physics.min_four_ulp physics.loose_with_reason physics.python_right \
            openfinding.cited.exactly openfinding.cited_by_prefix.case selftest.fails_on_purpose; do
    if grep -q "^FAIL [a-z-]*: ${twin}: " <<< "$out"; then bad "$twin reported"; else ok "$twin not reported"; fi
done

if [ "$fails" -ne 0 ]; then
    echo "== $fails claims self-test check(s) failed" >&2
    exit 1
fi
echo "== all claims self-test checks passed"
