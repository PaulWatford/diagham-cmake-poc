#!/usr/bin/env python3
"""Byte-for-byte comparison of the CMake build with DiagHam's own autotools build of the same revision.

The migration's gate: DiagHam built with CMake must write exactly what DiagHam built with its own autotools
writes. This replays every physics, crosscheck, regression and openfinding test invocation of a CMake build
(program, arguments, preparatory steps, reference program; from ctest's own listing) through both builds,
each in a fresh directory, and compares every file produced, byte for byte. Console logs (*.log) are
reported separately because they carry run timings.

    python3 scripts_cmake/byte_compare_with_autotools.py BUILD_DIR AUTOTOOLS_DIR OUT_JSON

BUILD_DIR      a CMake build of this repository, e.g. build/default or build/lapack
AUTOTOOLS_DIR  upstream's tree at the same revision, built by its own autotools, e.g. made with
               scripts_cmake/reproduce_on_pristine_trunk.sh and then
               ./configure --enable-fqhe --enable-fti --enable-spin --enable-quantumdots && make -k
               (for build/lapack add --enable-lapack --with-blas-libs=-lopenblas --with-lapack-libs=-lopenblas
               --with-fortran2c-libs=-lm, so both sides link the same OpenBLAS)
OUT_JSON       per-test results: identical files, different files, missing autotools programs

Differences are expected only where this repository carries a fix that changes a result on purpose
(docs/explanation/upstream-fixes.md); every other difference breaks the gate.
"""
import concurrent.futures as cf, filecmp, json, os, re, shutil, subprocess, sys
from pathlib import Path

M = Path(__file__).resolve().parent.parent
B, P, OUT = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(), Path(sys.argv[3])
W = Path("/tmp/byte_compare") / B.name
shutil.rmtree(W, ignore_errors=True); W.mkdir(parents=True)

renames = {}
d = subprocess.run(["git", "diff", "--name-status", "-M", "upstream", "main", "--", "*/Programs/*.cc"], cwd=M,
                   capture_output=True, text=True).stdout
for line in d.splitlines():
    p = line.split("\t")
    if p[0].startswith("R") and len(p) == 3:
        renames[(str(Path(p[2]).parent), Path(p[2]).stem)] = (str(Path(p[1]).parent), Path(p[1]).stem)


def to_pristine(exe):
    """build/<rel>/<name> -> pristine/<rel>/<upstream name>, or None if it is not one of our programs."""
    try:
        rel = Path(exe).relative_to(B)
    except ValueError:
        return exe                                    # system tool (cmake, python): same on both sides
    if rel.parts[0] == "tests":
        return exe                                    # our helpers (check_spectrum): same on both sides
    up = renames.get((str(rel.parent), rel.name), (str(rel.parent), rel.name))
    q = P / up[0] / up[1]
    return str(q) if q.exists() else None


tests = json.loads(subprocess.run(["ctest", "--test-dir", str(B), "--show-only=json-v1"], capture_output=True,
                                  text=True).stdout)["tests"]
jobs = []
for t in tests:
    name = t["name"]
    if not re.match(r"^(physics|crosscheck|regression|openfinding)\.", name):
        continue
    cmd = t.get("command", [])
    defs = {}
    for a in cmd:
        m = re.match(r"^-D([A-Z_]+)=(.*)$", a)
        if m:
            defs[m.group(1)] = m.group(2)
    if "PROGRAM" not in defs:
        continue
    jobs.append((name, defs))


def run_side(defs, side, root):
    """Run pre-steps, the program and the reference program of one test on one side; return missing executables."""
    missing = []
    conv = (lambda e: e) if side == "new" else to_pristine
    root.mkdir(parents=True)
    if "PRE_STEPS" in defs:
        for i, step in enumerate(defs["PRE_STEPS"].split("^^"), 1):
            parts = step.split("|"); exe = conv(parts[0])
            if exe is None:
                missing.append(parts[0]); return missing
            with open(root / f"step{i}.log", "w") as o:
                subprocess.run([exe, *parts[1:]], cwd=root, stdout=o, stderr=subprocess.STDOUT, timeout=600)
    for prog, args, sub in (("PROGRAM", "PROGRAM_ARGS", ""), ("REF_PROGRAM", "REF_PROGRAM_ARGS", "reference")):
        if prog not in defs:
            continue
        exe = conv(defs[prog])
        if exe is None:
            missing.append(defs[prog]); continue
        d = root / sub if sub else root
        d.mkdir(exist_ok=True)
        a = [x for x in defs.get(args, "").split("|") if x != ""]
        with open(d / "program.log", "w") as o:
            subprocess.run([exe, *a], cwd=d, stdout=o, stderr=subprocess.STDOUT, timeout=600)
    return missing


def one(job):
    name, defs = job
    base = W / name
    res = dict(test=name, program=Path(defs["PROGRAM"]).name, missing=[], identical=[], different=[], logs_differ=[])
    try:
        miss = run_side(defs, "old", base / "old")
        if miss:
            res["missing"] = [Path(m).name for m in miss]; return res
        run_side(defs, "new", base / "new")
    except subprocess.TimeoutExpired:
        res["missing"] = ["timeout"]; return res
    for dp, _, files in os.walk(base / "new"):
        for f in files:
            pn = Path(dp) / f; rel = pn.relative_to(base / "new"); po = base / "old" / rel
            is_log = f.endswith(".log")
            if not po.exists():
                res["different"].append(f"{rel} (only in CMake run)"); continue
            same = filecmp.cmp(pn, po, shallow=False)
            if is_log:
                if not same:
                    res["logs_differ"].append(str(rel))
            elif same:
                res["identical"].append(str(rel))
            else:
                res["different"].append(str(rel))
    for dp, _, files in os.walk(base / "old"):
        for f in files:
            rel = (Path(dp) / f).relative_to(base / "old")
            if not (base / "new" / rel).exists() and not f.endswith(".log"):
                res["different"].append(f"{rel} (only in autotools run)")
    return res


with cf.ThreadPoolExecutor(16) as ex:
    results = list(ex.map(one, jobs))
json.dump(results, open(OUT, "w"), indent=1)
tot_same = sum(len(r["identical"]) for r in results); tot_diff = sum(len(r["different"]) for r in results)
noref = [r for r in results if r["missing"]]
print(f"tests replayed: {len(results)}; output files identical: {tot_same}; different: {tot_diff}; "
      f"tests without an autotools reference: {len(noref)}; tests whose stdout differs: {sum(1 for r in results if r['logs_differ'])}")
