#!/usr/bin/env python3
"""Write the claim index: what each physics and crosscheck test establishes and where its answer is known from.

Usage: tests/claims.py <build-dir> [--write docs/reference/claims.md] [--check docs/reference/claims.md] [--missing]

Reads `ctest --show-only=json-v1` from the build directory (through tests/coverage.py). Every test labelled
`physics` or `crosscheck` carries the test property KNOWN_FROM (set by diagham_physics_test and the other
registration helpers; see tests/CMakeLists.txt), whose first word is its kind:

  oracle tests/oracles/<file>.py     an independent Python calculation in the tree
  analytic: <one line>               a closed form
  exact count: <one line>            a dimension or degeneracy by counting
  independent program: <target>      another DiagHam program (the crosscheck case)
  paper: <reference>                 a published value
  saved reference: <file>, cross-checked by <test>   a stored spectrum whose own check is named

--write produces the page; --missing lists every physics or crosscheck test without KNOWN_FROM, with the line of
tests/ where its name (or the literal part of a loop-built name) is registered. tests/check_claims.py
(selftest.claims) uses the same parsing to enforce the rules.
"""
import argparse, os, re, sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import coverage  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
PAGE = "docs/reference/claims.md"
KINDS = ("oracle", "analytic", "exact count", "independent program", "paper", "saved reference")
GENERATED = {"test-coverage.md", "claims.md", "numbers.md"}


def kind_of(known_from):
    """The kind word of a KNOWN_FROM text, or None when it does not start with one."""
    k = (known_from or "").strip()
    for w in KINDS:
        if k.startswith(w + ":") or (w == "oracle" and k.startswith("oracle ")):
            return w
    return None


def records(build_dir):
    """One dict per test of the build, with what the claim index and the rules need."""
    build_abs = os.path.abspath(build_dir)
    out = []
    for t in coverage.ctest_json(build_dir)["tests"]:
        props = {p["name"]: p["value"] for p in t.get("properties", [])}
        cmd = t.get("command") or []
        r = {"name": t["name"], "labels": list(props.get("LABELS", [])), "command": cmd,
             "known_from": props.get("KNOWN_FROM", ""), "reason": props.get("TOLERANCE_REASON", ""),
             "will_fail": bool(props.get("WILL_FAIL", False)),
             "programs": sorted(coverage.programs_of(t, build_abs, {"check_spectrum"}))}
        check = next((a[len("-DCHECK_ARGS="):] for a in cmd if a.startswith("-DCHECK_ARGS=")), None)
        r["runner"] = ("python" if cmd and os.path.basename(cmd[0]).startswith("python")
                       else "check_spectrum" if check is not None
                       else "dimension" if any(a.endswith("RunAndMatch.cmake") for a in cmd) else "other")
        r["check"] = check.split("|") if check else []
        r["script"] = next((a for a in cmd[1:] if a.endswith(".py")), "") if r["runner"] == "python" else ""
        out.append(r)
    return out


def claim_tests(recs):
    return [r for r in recs if {"physics", "crosscheck"} & set(r["labels"])]


def checks_text(r):
    if r["runner"] == "check_spectrum" and r["check"]:
        mode = r["check"][0]
        col = r["check"][2] if len(r["check"]) > 2 and mode not in ("line", "numbers") else ""
        return f"`{mode}`" + (f" column {col}" if col not in ("", "@OUTPUT@") else "")
    if r["runner"] == "python":
        rel = os.path.relpath(r["script"], ROOT) if r["script"].startswith("/") else r["script"]
        mode = next((a for a in r["command"][2:3] if not a.startswith("/")), "")
        return f"`{rel}`" + (f" `{mode}`" if mode else "")
    if r["runner"] == "dimension":
        return "basis dimensions, line by line (`RunAndMatch.cmake`)"
    return "not read"


def registration_line(name):
    """tests/<file>:<line> where the name, or the literal prefix of a loop-built name, is registered."""
    files = [ROOT / "tests/CMakeLists.txt"] + sorted((ROOT / "tests").glob("*.cmake"))
    parts = name.split(".")
    for k in range(len(parts), 1, -1):
        pre = ".".join(parts[:k])
        pat = re.compile(re.escape(pre) + (r"\b" if k == len(parts) else r"\.?\$\{"))
        for f in files:
            for i, line in enumerate(f.read_text().splitlines(), 1):
                if pat.search(line):
                    return f"{f.relative_to(ROOT)}:{i}"
    return "not found"


def citations(names):
    """{test name: [markdown files that cite it in backticks]}, generated pages excepted."""
    files = [p for p in list((ROOT / "docs").rglob("*.md")) + list(ROOT.glob("*.md")) if p.name not in GENERATED]
    texts = {p: p.read_text(encoding="utf-8", errors="replace") for p in files}
    out = {}
    for n in names:
        key = f"`{n}`"
        out[n] = sorted(str(p.relative_to(ROOT)) for p, t in texts.items() if key in t)
    return out


def render(build_dir, recs=None):
    recs = recs or records(build_dir)
    claims = sorted(claim_tests(recs), key=lambda r: r["name"])
    cited = citations([r["name"] for r in claims])
    by_kind = Counter(kind_of(r["known_from"]) or ("no source" if not r["known_from"] else "unrecognised kind")
                      for r in claims)
    preset = Path(os.path.abspath(build_dir)).name
    lines = [
        "<!-- generated by tests/claims.py; do not edit -->",
        "# Claim index: what each physics test establishes and where its answer is known from",
        "",
        f"Source: `ctest --show-only=json-v1` of the `{preset}` build, the test property `KNOWN_FROM` of every "
        f"`physics` and `crosscheck` test; {len(claims)} tests. Regenerate with "
        f"`python3 tests/claims.py <build-dir> --write {PAGE}`; `scripts_cmake/ci.sh default` fails when it is stale.",
        "Purpose: one row per test that claims an independently known answer, saying what it checks and why the "
        "expected answer is known; `selftest.claims` (`tests/check_claims.py`) reports a row with no source, and fails "
        "the suite on a source of no recognised kind or an oracle that is not the one the test runs.",
        "",
        "| Known from (kind) | Tests |",
        "|---|---:|",
    ]
    for k in KINDS + ("unrecognised kind", "no source"):
        if by_kind.get(k):
            lines.append(f"| {k} | {by_kind[k]} |")
    lines += ["", "| Test | Program | Checks | Known from | Cited by |", "|---|---|---|---|---|"]
    for r in claims:
        # "<r>" would otherwise be read as an HTML tag by the renderer
        kf = (r["known_from"] or "**no source**").replace("|", "\\|").replace("<", "&lt;").replace(">", "&gt;")
        cite = ", ".join(cited[r["name"]]) or ""
        lines.append(f"| `{r['name']}` | {', '.join(r['programs']) or 'not read'} | {checks_text(r)} | {kf} | {cite} |")
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("build_dir")
    ap.add_argument("--write", help="write the page here")
    ap.add_argument("--check", metavar="PAGE", help="exit 1 when PAGE differs from what --write would produce")
    ap.add_argument("--missing", action="store_true", help="list the physics and crosscheck tests without KNOWN_FROM")
    a = ap.parse_args()
    recs = records(a.build_dir)
    if a.write:
        Path(a.write).write_text(render(a.build_dir, recs), encoding="utf-8")
    if a.check:
        if Path(a.check).read_text(encoding="utf-8") != render(a.build_dir, recs):
            print(f"{a.check} is stale: regenerate with python3 tests/claims.py build/default --write {PAGE}",
                  file=sys.stderr)
            return 1
        print(f"{a.check} is current")
    if a.missing:
        claims = claim_tests(recs)
        missing = sorted((r for r in claims if not r["known_from"]), key=lambda r: r["name"])
        for r in missing:
            print(f"{r['name']}  {registration_line(r['name'])}")
        print(f"{len(claims) - len(missing)} of {len(claims)} physics and crosscheck tests have KNOWN_FROM; "
              f"{len(missing)} missing")
    return 0


if __name__ == "__main__":
    sys.exit(main())
