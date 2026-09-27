#!/usr/bin/env python3
"""Regression references: run each case and save what the program printed/wrote.

    python3 tests/regression_reference.py <build-dir> --write            # (re)generate every reference
    python3 tests/regression_reference.py <build-dir> --write --only NAME

Cases are listed in tests/data/regression/cases.txt, one per line:

    name | target | output-glob | args...

`target` is the ctest program target (<leafdir>_<Program>); the executable
is located in the build directory. `output-glob` names the spectrum file
the program writes ("program.log" = its standard output; then only lines
made of numbers are compared). The program runs in an empty temporary
directory; its output is stored as tests/data/regression/<name>.dat with a
provenance header (program, arguments, DiagHam revision, commit, preset,
compiler, date).

A regression reference is NOT a golden: it records what this build
produced, so the tests built on it (regression.<name>) detect change, not
correctness. Regenerate only when a change is intended and understood, and
say so in the commit.
"""
import argparse, datetime, glob, os, shutil, subprocess, sys, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CASES = ROOT / "tests" / "data" / "regression" / "cases.txt"


def read_cases():
    out = []
    for line in CASES.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = [p.strip() for p in line.split("|")]
        name, target, glob_ = parts[0], parts[1], parts[2]
        args = parts[3].split() if len(parts) > 3 else []
        out.append((name, target, glob_, args))
    return out


def find_program(build_dir, target):
    leaf, _, prog = target.partition("_")
    hits = [p for p in Path(build_dir).rglob(prog) if p.is_file() and os.access(p, os.X_OK) and "CMakeFiles" not in p.parts]
    hits = sorted(hits, key=lambda p: (leaf not in p.parts, len(p.parts)))
    if not hits:
        raise SystemExit(f"program {target} not found under {build_dir}")
    return hits[0]


def provenance(build_dir, target, args):
    def git(*a):
        try:
            return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()
        except Exception:
            return "?"
    cache = Path(build_dir) / "CMakeCache.txt"
    preset = compiler = "?"
    if cache.exists():
        for l in cache.read_text().splitlines():
            if l.startswith("CMAKE_CXX_COMPILER:"):
                compiler = l.split("=", 1)[1]
    preset = Path(build_dir).name
    rev = (ROOT / "scripts_cmake" / "upstream-revision.txt").read_text().strip() if (ROOT / "scripts_cmake" / "upstream-revision.txt").exists() else "?"
    return [f"# regression reference: {target}",
            f"# args: {' '.join(args)}",
            f"# DiagHam r{rev}, commit {git('rev-parse', '--short', 'HEAD')}, preset {preset}, {compiler}, {datetime.date.today().isoformat()}",
            "# what this build produced -- detects change, not correctness (see tests/regression_reference.py)"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("build_dir")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--only")
    a = ap.parse_args()
    outdir = CASES.parent
    for name, target, glob_, args in read_cases():
        if a.only and name != a.only:
            continue
        prog = find_program(Path(a.build_dir).resolve(), target)
        args = [x.replace("@DATA@", str(outdir.resolve())) for x in args]
        work = Path(tempfile.mkdtemp(prefix="regref_"))
        with open(work / "program.log", "w") as log:
            rc = subprocess.run([str(prog), *args], cwd=work, stdout=log, stderr=subprocess.STDOUT).returncode
        if rc != 0:
            print(f"{name}: exit {rc}")
            continue
        files = glob.glob(str(work / glob_))
        if len(files) != 1:
            print(f"{name}: expected one file matching {glob_}, found {len(files)}")
            continue
        body = Path(files[0]).read_text()
        if glob_ == "program.log":   # keep only the numeric lines of stdout
            keep = []
            for l in body.splitlines():
                toks = l.split()
                try:
                    [float(t) for t in toks]
                    if toks:
                        keep.append(l)
                except ValueError:
                    pass
            body = "\n".join(keep) + "\n"
        text = "\n".join(provenance(a.build_dir, target, args)) + "\n" + body
        n = sum(1 for l in body.splitlines() if l.strip() and not l.startswith("#"))
        if a.write:
            (outdir / f"{name}.dat").write_text(text)
            print(f"{name}: {n} rows -> {outdir / (name + '.dat')}")
        else:
            print(f"{name}: {n} rows (not written)")
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
