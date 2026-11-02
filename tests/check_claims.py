#!/usr/bin/env python3
"""selftest.claims: the rules every registered test must keep (CLAUDE.md G4).

Usage: tests/check_claims.py <build-dir> [--known-from-required] [--source-dir DIR]

Reads `ctest --show-only=json-v1` (through tests/claims.py) and checks:

  1. known-from   every physics or crosscheck test has KNOWN_FROM, of a recognised kind; an "oracle <file>" source
                  exists in the tree and, for a Python test, is the script the test runs. A missing KNOWN_FROM is
                  reported but does not fail unless --known-from-required (or KNOWN_FROM_REQUIRED=1 in the
                  environment) is given; a wrong one always fails.
  2. tolerance    every check_spectrum comparison of a physics or crosscheck test is within the policy (relative
                  rel:1e-12 or tighter, at most 4 ulp for `min`, absolute at most 1e-9) or has TOLERANCE_REASON.
  3. labels       every label of every test is in tests/labels.txt.
  4. open-finding a test has WILL_FAIL exactly when it has the label open-finding; the one exception is a
                  `selftest` that fails on purpose (the checker rejecting a wrong value).
  5. register     every open-finding test is cited in docs/reference/open-findings.md, by name or by a cited
                  `prefix.*`.

Exit 0 when nothing fails, 1 otherwise. --source-dir points the file lookups (labels, register, oracles) at another
tree, for tests/claims_selftest.sh.
"""
import argparse, os, re, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import claims  # noqa: E402

# Position of the tolerance (or zero threshold) in the check_spectrum arguments, mode first; see check_spectrum.cc.
TOLERANCE_INDEX = {"min": 4, "min-abs": 4, "spectrum": 5, "count": 4, "lowest": 6, "numbers": 3, "line": 4,
                   "nonzero": 3}
MAX_REL, MAX_ABS, MAX_ULP = 1e-12, 1e-9, 4


def tolerance_problem(check):
    """None when the CHECK is within the policy, else a short description of why not."""
    if not check:
        return None
    mode = check[0]
    i = TOLERANCE_INDEX.get(mode)
    if i is None:
        return f"unknown check_spectrum mode {mode!r}"
    if len(check) <= i:
        return f"{mode}: no tolerance at argument {i}"
    tol = check[i]
    try:
        if mode == "min":
            return None if int(tol) <= MAX_ULP else f"min: {tol} ulp > {MAX_ULP}"
        if tol.startswith("rel:"):
            return None if float(tol[4:]) <= MAX_REL else f"{mode}: {tol} looser than rel:{MAX_REL:g}"
        return None if float(tol) <= MAX_ABS else f"{mode}: absolute {tol} looser than {MAX_ABS:g}"
    except ValueError:
        return f"{mode}: tolerance {tol!r} is not a number"


def register_names(text):
    cited = set(re.findall(r"`(openfinding\.[^`]+)`", text))
    exact = {c for c in cited if not c.endswith(".*")}
    prefixes = tuple(c[:-1] for c in cited if c.endswith(".*"))
    return exact, prefixes


def check(recs, src, known_from_required):
    fails, notes = [], []
    vocab = {l.strip() for l in (src / "tests/labels.txt").read_text().splitlines()
             if l.strip() and not l.startswith("#")}
    exact, prefixes = register_names((src / "docs/reference/open-findings.md").read_text(encoding="utf-8"))
    missing_known = []
    for r in recs:
        n, labels = r["name"], set(r["labels"])
        if {"physics", "crosscheck"} & labels:
            kf = r["known_from"]
            if not kf:
                missing_known.append(n)
            else:
                kind = claims.kind_of(kf)
                if kind is None:
                    fails.append(f"known-from: {n}: KNOWN_FROM {kf!r} does not start with a recognised kind")
                elif kind == "oracle":
                    oracle = kf[len("oracle "):].strip()
                    if not (src / oracle).is_file():
                        fails.append(f"known-from: {n}: oracle {oracle} is not in the tree")
                    elif r["runner"] == "python" and not r["script"].endswith("/" + oracle):
                        fails.append(f"known-from: {n}: KNOWN_FROM names {oracle} but the test runs {r['script']}")
            if r["runner"] == "check_spectrum" and not r["reason"]:
                p = tolerance_problem(r["check"])
                if p:
                    fails.append(f"tolerance: {n}: {p} and no TOLERANCE_REASON")
        for l in sorted(labels - vocab):
            fails.append(f"labels: {n}: label {l!r} is not in tests/labels.txt")
        # the checker's own self-tests (label selftest) fail on purpose to show that a wrong answer is rejected
        if r["will_fail"] != ("open-finding" in labels) and not (r["will_fail"] and "selftest" in labels):
            fails.append(f"open-finding: {n}: WILL_FAIL is {r['will_fail']} but the label open-finding is "
                         f"{'present' if 'open-finding' in labels else 'absent'}")
        if "open-finding" in labels and n not in exact and not n.startswith(prefixes):
            fails.append(f"register: {n}: not cited in docs/reference/open-findings.md")
    if missing_known:
        line = (f"known-from: {len(missing_known)} physics or crosscheck test(s) without KNOWN_FROM: "
                + ", ".join(missing_known))
        (fails if known_from_required else notes).append(line)
    return fails, notes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("build_dir")
    ap.add_argument("--known-from-required", action="store_true")
    ap.add_argument("--source-dir", default=str(claims.ROOT))
    a = ap.parse_args()
    required = a.known_from_required or os.environ.get("KNOWN_FROM_REQUIRED") == "1"
    recs = claims.records(a.build_dir)
    fails, notes = check(recs, Path(a.source_dir), required)
    n_claims = len(claims.claim_tests(recs))
    for f in fails:
        print("FAIL", f)
    for x in notes:
        print("NOTE", x, "(report-only until every test has one)")
    print(f"selftest.claims: {len(recs)} tests, {n_claims} physics or crosscheck; {len(fails)} failure(s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
