#!/usr/bin/env python3
"""The upstream collaboration pack: one report and one Subversion-applicable diff per carried fix.

    python3 scripts_cmake/upstream_pack.py            # regenerate patches/upstream/ and docs/upstream-reports/
    python3 scripts_cmake/upstream_pack.py --check    # apply every diff in series to a pristine upstream tree; exit 1 on failure

Every fix to DiagHam sources on `main` is one commit with the trailers
`Upstream-Patch:` and `Upstream-Base:`. For each such commit (oldest first) this
writes
  patches/upstream/NN-<name>.diff       `git diff --no-prefix` of the DiagHam source files
                                        only (no CMake, tests, docs or stray *.orig), CRLF
                                        preserved, paths as they are upstream (two FQHEOnDisk
                                        programs were renamed here; see RENAMED); `svn patch`
                                        applies the series in order to a trunk working copy
                                        at the base revision
  docs/upstream-reports/NN-<name>.md    the report to send with it: what, why, how it was
                                        verified (the commit message, or the fix index row
                                        and the audit trail when the message is only a
                                        subject), files, how to apply
and an index. Hand-written reports for defects that are registered but not fixed
(docs/upstream-reports/open-*.md) are listed in the index as well.
"""
import re, shutil, subprocess, sys, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PATCHES = ROOT / "patches" / "upstream"
REPORTS = ROOT / "docs" / "upstream-reports"
SOURCE_DIRS = ("Base/", "src/", "FQHE/", "FTI/", "Spin/", "QuantumDots/")
# files renamed in this repository (name collisions at install); upstream keeps the old names
RENAMED = {
    "FQHE/src/Programs/FQHEOnDisk/FQHEDiskBosonsDelta.cc": "FQHE/src/Programs/FQHEOnDisk/QHEBosonsDelta.cc",
    "FQHE/src/Programs/FQHEOnDisk/FQHEDiskLaughlinMonteCarloOverlap.cc": "FQHE/src/Programs/FQHEOnDisk/QHEFermionsOverlap.cc",
}


def git(*args, binary=False):
    r = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, check=True)
    return r.stdout if binary else r.stdout.decode("utf-8", "replace")


def fixes():
    out = []
    log = git("log", "--reverse", "--format=%H%x00%s%x00%b%x00%(trailers:key=Upstream-Patch,valueonly)%x00%(trailers:key=Upstream-Base,valueonly)%x01",
              "^upstream", "main")
    for rec in log.split("\x01"):
        if not rec.strip():
            continue
        sha, subject, body, patch, base = rec.strip("\n").split("\x00")
        if not patch.strip():
            continue
        files = [f for f in git("diff-tree", "--no-commit-id", "--name-only", "-r", sha).split()
                 if f.startswith(SOURCE_DIRS) and not f.endswith("CMakeLists.txt") and not f.endswith((".orig", ".rej"))]
        if not files:
            continue
        body = re.sub(r"\n*Upstream-Patch:.*$", "", body, flags=re.S).strip()
        out.append(dict(sha=sha, subject=subject, body=body, patch=patch.strip(), base=base.strip(), files=files))
    return out


def slug(fix):
    name = re.sub(r"\.patch.*$", "", fix["patch"])
    name = re.sub(r"^\d+-", "", name)
    return re.sub(r"[^A-Za-z0-9]+", "-", name).strip("-")


def diff_of(fix):
    """the fix's diff for the files that exist under the same name upstream; renamed files are
    reported as hunks to apply by hand (their context differs by the rename commit)"""
    files = [f for f in fix["files"] if f not in RENAMED]
    if not files:
        return b""
    return git("diff", "--no-prefix", f"{fix['sha']}^", fix["sha"], "--", *files, binary=True)


def renamed_hunks(fix):
    files = [f for f in fix["files"] if f in RENAMED]
    if not files:
        return ""
    text = git("diff", "--no-prefix", "-U1", f"{fix['sha']}^", fix["sha"], "--", *files)
    for new, old in RENAMED.items():
        text = text.replace(new, old)
    return text


def index_rows():
    """rows of docs/explanation/upstream-fixes.md keyed by the (NN) patch number they cite"""
    rows = {}
    p = ROOT / "docs" / "explanation" / "upstream-fixes.md"
    if p.exists():
        for line in p.read_text().splitlines():
            m = re.match(r"^\| ([A-Z]) \| (.*?) \| ([^|]*) \| ([^|]*) \| ([^|]*) \|$", line)
            if m:
                for num in re.findall(r"\((\d\d)\)", m.group(2)):
                    rows[num] = dict(cls=m.group(1), what=m.group(2), kind=m.group(4).strip(), status=m.group(5).strip())
    return rows


def report(n, fix, rows):
    num = re.match(r"^(\d\d)-", fix["patch"])
    row = rows.get(num.group(1)) if num else None
    if fix["body"]:
        what = fix["body"]
    elif row:
        what = (f"The commit message carries only the subject; the fix index says: **kind** — {row['kind']}; "
                f"**upstream status** — {row['status']} (class {row['cls']}). The audit trail with the root cause and "
                f"the verification is class {row['cls']} of `patches/PATCHES.md`.")
    else:
        what = "The commit message carries only the subject; see `docs/explanation/upstream-fixes.md` and `patches/PATCHES.md`."
    renamed = [f for f in fix["files"] if f in RENAMED]
    note = ""
    if renamed:
        note = ("\n\n**Note**: " + "; ".join(f"`{Path(f).name}` is `{Path(RENAMED[f]).name}` upstream (renamed here because of a name collision at install)" for f in renamed)
                + ". Those files are not in the diff; the same change is given below as hunks with the upstream names, to apply by hand (their context lines differ by the rename).\n\n"
                "```diff\n" + renamed_hunks(fix).rstrip("\n") + "\n```")
    return f"""# {n:02d}. {fix['subject']}

Purpose: the report to send to the DiagHam authors with the diff `patches/upstream/{n:02d}-{slug(fix)}.diff`.
Source: commit `{fix['sha'][:9]}` on `main` (generated by `scripts_cmake/upstream_pack.py`, do not edit; edit the commit message or the fix index instead).

**Base**: {fix['base'] or 'SVN r4493'}. **Files**: {', '.join(f'`{f}`' for f in fix['files'])}.{note}

## What, why, how it was verified

{what}

## Applying to a Subversion working copy

The diffs form a series: apply `01` to `{n:02d}` in order (earlier ones may touch the same files).

```
cd DiagHam            # trunk working copy at the base revision
svn patch /path/to/{n:02d}-{slug(fix)}.diff
```

No `a/`/`b/` prefixes; `git apply -p0` accepts it too. Whether the fix has been taken upstream is tracked in `docs/explanation/upstream-fixes.md`.
"""


def check():
    fx = fixes()
    work = Path(tempfile.mkdtemp(prefix="upstream-check-"))
    subprocess.run(["git", "worktree", "add", "-q", str(work), "upstream"], cwd=ROOT, check=True)
    bad = 0
    try:
        for n, fix in enumerate(fx, 1):
            p = PATCHES / f"{n:02d}-{slug(fix)}.diff"
            if not p.exists():
                print(f"missing {p}"); bad += 1; continue
            r = subprocess.run(["git", "apply", "-p0", "--whitespace=nowarn", str(p)], cwd=work, capture_output=True, text=True)
            if r.returncode != 0:
                print(f"FAILS in series on upstream: {p.name}: {r.stderr.strip().splitlines()[0] if r.stderr.strip() else '?'}")
                bad += 1
            else:
                print(f"ok   {p.name}")
    finally:
        subprocess.run(["git", "worktree", "remove", "--force", str(work)], cwd=ROOT)
    print(f"{len(fx)} diffs, {bad} failures")
    sys.exit(1 if bad else 0)


def main():
    if subprocess.run(["git", "rev-parse", "--verify", "-q", "upstream"], cwd=ROOT, capture_output=True).returncode != 0:
        sys.exit("no `upstream` branch in this clone: git fetch origin upstream:upstream   (the pack is the difference between main and upstream)")
    if "--check" in sys.argv:
        check()
        return
    fx = fixes()
    rows = index_rows()
    PATCHES.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    for old in list(PATCHES.glob("*.diff")) + [p for p in REPORTS.glob("*.md") if re.match(r"\d\d-", p.name)]:
        old.unlink()
    index = ["# Reports for the DiagHam authors",
             "",
             "Purpose: everything needed to offer this repository's fixes and findings upstream — one report and one `svn patch`-applicable diff per fix, and a report per open defect.",
             "Source: generated from the commits on `main` that carry an `Upstream-Patch:` trailer (`scripts_cmake/upstream_pack.py`; `--check` applies the diffs in series to a pristine `upstream` tree); the open-defect reports are hand-written from `docs/reference/known-defects.md`. Nothing here is sent automatically: `docs/how-to/develop/send-a-fix-upstream.md` says how.",
             "",
             "## Fixes carried on `main` (diffs in `patches/upstream/`, a series to apply in order)",
             "",
             "| # | Fix | Files | Diff |",
             "|---|---|---|---|"]
    for n, fix in enumerate(fx, 1):
        s = slug(fix)
        (PATCHES / f"{n:02d}-{s}.diff").write_bytes(diff_of(fix))
        (REPORTS / f"{n:02d}-{s}.md").write_text(report(n, fix, rows))
        index.append(f"| {n:02d} | [{fix['subject'][:90]}]({n:02d}-{s}.md) | {len(fix['files'])} | `patches/upstream/{n:02d}-{s}.diff` |")
    opens = sorted(REPORTS.glob("open-*.md"))
    if opens:
        index += ["", "## Registered defects without a fix (reproducers)", ""]
        for p in opens:
            title = next((l[2:] for l in p.read_text().splitlines() if l.startswith("# ")), p.stem)
            index.append(f"- [{title}]({p.name})")
    index += ["", f"Generated from {len(fx)} fix commits; regenerate with `python3 scripts_cmake/upstream_pack.py`."]
    (REPORTS / "README.md").write_text("\n".join(index) + "\n")
    (PATCHES / "README.md").write_text(
        "# Subversion-applicable diffs of the fixes carried on `main`\n\n"
        "Generated by `scripts_cmake/upstream_pack.py` from the commits with an `Upstream-Patch:` trailer; do not edit. "
        "Each diff covers DiagHam source files only (no CMake, tests or docs), keeps the files' line endings, uses the "
        "upstream file names, and has no path prefixes; they form a series:\n\n"
        "    svn patch 01-....diff; svn patch 02-....diff; ...   # in a trunk working copy at the base revision\n\n"
        "`python3 scripts_cmake/upstream_pack.py --check` applies the series to a pristine `upstream` tree. "
        "The report to send with each diff is `docs/upstream-reports/NN-name.md`; the index with the upstream status of "
        "each fix is `docs/explanation/upstream-fixes.md`. `patches/PATCHES.md` is the frozen audit trail of the "
        "proof-of-concept era and is not regenerated.\n")
    print(f"{len(fx)} fixes -> {PATCHES.relative_to(ROOT)}/ and {REPORTS.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
