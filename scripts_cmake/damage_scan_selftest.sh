#!/usr/bin/env bash
# damage_scan_selftest.sh: exercise scripts_cmake/damage_scan.py and the Guard scripts without GitHub.
#
#   scripts_cmake/damage_scan_selftest.sh            # from the repository root; needs git and python3
#
# Makes a throw-away clone of this repository in a temporary directory, commits synthetic changes on
# branches there, and checks that the scan reports each damage class and nothing for a clean change.
# Then runs the author gate and the report script on the fixtures in .github/scripts/fixtures with
# --dry-run (nothing is sent anywhere). The working tree you run it from is not touched.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
fail=0
say() { printf '\n== %s\n' "$*"; }
expect_rc() { # expect_rc WANT NAME CMD...
    local want="$1" name="$2"; shift 2
    set +e; "$@" > "$TMP/$name.out" 2>&1; local rc=$?; set -e
    if [ "$rc" -eq "$want" ]; then echo "ok   $name (exit $rc)"; else echo "FAIL $name: exit $rc, wanted $want"; sed -n 1,40p "$TMP/$name.out"; fail=1; fi
}
expect_in() { # expect_in NAME PATTERN
    if grep -q -- "$2" "$TMP/$1.out"; then echo "ok   $1 mentions '$2'"; else echo "FAIL $1 does not mention '$2'"; fail=1; fi
}

say "clone"
git clone -q --no-hardlinks --single-branch --branch "$(git -C "$ROOT" rev-parse --abbrev-ref HEAD)" "$ROOT" "$TMP/repo"
cd "$TMP/repo"
git config user.email selftest@example.invalid; git config user.name selftest
SCAN="python3 $ROOT/scripts_cmake/damage_scan.py"
BASE=$(git rev-parse HEAD)

say "clean change: a documentation line"
git checkout -q -b clean; echo "" >> README.md; git commit -qam "docs: blank line"
expect_rc 0 clean $SCAN "$BASE" clean
expect_in clean "Damage:\*\* none"

say "tests removed, count kept (names replaced)"
git checkout -q -b replaced "$BASE"
python3 - <<'EOF'
from pathlib import Path
p = Path("tests/manifest.txt"); lines = p.read_text().splitlines()
victim = next(l for l in lines if l.startswith("physics.fqhe.sphere."))
lines = [l for l in lines if l != victim] + ["physics.made_up.new_name"]
p.write_text("\n".join(lines) + "\n")
EOF
git commit -qam "tests: swap a name"
expect_rc 2 replaced $SCAN "$BASE" replaced
expect_in replaced "tests removed"; expect_in replaced "count is unchanged"

say "demo unregistered"
git checkout -q -b demo "$BASE"; sed -i '/^physics\.hubbard\.2x2\.U4\.ground_state_is_minus_4sqrt2$/d' tests/manifest.txt; git commit -qam "drop the locked golden"
expect_rc 2 demo $SCAN "$BASE" demo
expect_in demo "demo unregistered"

say "checks skipped: enable_testing removed, DIAGHAM_BUILD_TESTS OFF, a locked name made optional"
git checkout -q -b skipped "$BASE"
sed -i 's/^\(\s*\)enable_testing()/\1# enable_testing()/' CMakeLists.txt
sed -i 's/^\(option(DIAGHAM_BUILD_TESTS .*\) ON)/\1 OFF)/' CMakeLists.txt
echo "physics.hubbard.2x2.U0.ground_state_is_minus_8" >> tests/manifest-optional.txt
git commit -qam "skip checks"
expect_rc 2 skipped $SCAN "$BASE" skipped
expect_in skipped "enable_testing() removed"; expect_in skipped "DIAGHAM_BUILD_TESTS"; expect_in skipped "manifest-optional"

say "golden edited: a reference spectrum and an expected value"
git checkout -q -b golden "$BASE"
f=$(ls benchmarks/*.dat.txt | head -1); echo "0.0" >> "$f"
sed -i 's/-5\.6568542494923806/-5.6568542494923800/' tests/CMakeLists.txt
git commit -qam "edit goldens"
expect_rc 2 golden $SCAN "$BASE" golden
expect_in golden "golden edited"; expect_in golden "reference file"; expect_in golden "value or tolerance"

say "option default flipped and a module option removed"
git checkout -q -b option "$BASE"
sed -i 's/^\(option(DIAGHAM_BUILD_SPIN .*\) ON)/\1 OFF)/' CMakeLists.txt
sed -i '/^option(DIAGHAM_BUILD_QUANTUMDOTS /d' CMakeLists.txt
git commit -qam "flip options"
expect_rc 2 option $SCAN "$BASE" option
expect_in option "DIAGHAM_BUILD_SPIN"; expect_in option "DIAGHAM_BUILD_QUANTUMDOTS"

say "in-source build guard removed"
git checkout -q -b insource "$BASE"
python3 - <<'EOF'
from pathlib import Path
import re
p = Path("CMakeLists.txt"); s = p.read_text()
s = re.sub(r"if\(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR\).*?endif\(\)\n", "", s, count=1, flags=re.S)
p.write_text(s)
EOF
git commit -qam "allow in-source"
expect_rc 2 insource $SCAN "$BASE" insource
expect_in insource "in-source build"

say "patches changed"
git checkout -q -b patches "$BASE"; p=$(ls patches/*.patch | head -1); echo "# touched" >> "$p"; git commit -qam "touch a patch"
expect_rc 2 patches $SCAN "$BASE" patches
expect_in patches "patches changed"

say "test added (information only)"
git checkout -q -b added "$BASE"; echo "physics.selftest.new_name" >> tests/manifest.txt; git commit -qam "add a name"
expect_rc 0 added $SCAN "$BASE" added
expect_in added "1 added, 0 removed"

say "unknown commit"
expect_rc 1 unknown $SCAN "$BASE" no-such-ref

say "author gate on the fixtures (dry run; nothing sent)"
cd "$ROOT"
G="python3 $ROOT/.github/scripts/author_gate.py --repo-root $ROOT --dry-run"
expect_rc 1 gate_external_ai $G .github/scripts/fixtures/external_ai_physics.json
expect_in gate_external_ai "HELD_FOR_APPROVAL"; expect_in gate_external_ai "for maintainer review"; expect_in gate_external_ai "ai-assisted"
expect_rc 0 gate_approved $G .github/scripts/fixtures/approved_no_ai.json
expect_in gate_approved "| Hold | none |"
expect_rc 0 gate_released $G .github/scripts/fixtures/external_released_by_label.json
expect_in gate_released "released: @PaulWatford"
expect_rc 1 gate_bad_label $G .github/scripts/fixtures/bad_release_label.json
expect_in gate_bad_label "not declared"; expect_in gate_bad_label "remove labels \['approved-by-maintainer'\]"

say "report script on a damaged scan (dry run)"
$SCAN "$BASE" "$(git -C "$TMP/repo" rev-parse golden)" --repo "$TMP/repo" --markdown "$TMP/d.md" --json "$TMP/d.json" --quiet || true
R="python3 $ROOT/.github/scripts/guard_report.py $TMP/d.md $TMP/d.json --repo-root $ROOT --dry-run"
expect_rc 1 report_held $R .github/scripts/fixtures/external_ai_physics.json
expect_in report_held "add labels \['damage', 'held-for-approval'\]"
expect_rc 0 report_released $R .github/scripts/fixtures/external_released_by_label.json

say "workflow file parses as YAML"
expect_rc 0 yaml python3 -c "import yaml,sys; d=yaml.safe_load(open('$ROOT/.github/workflows/guard.yml')); assert set(d['jobs'])=={'author-gate','damage-scan'}; print('jobs:', list(d['jobs']))"

if [ "$fail" -eq 0 ]; then say "all checks passed"; else say "SOME CHECKS FAILED"; exit 1; fi
