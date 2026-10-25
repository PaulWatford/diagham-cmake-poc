#!/usr/bin/env python3
"""The one version number must agree everywhere it is written.

    python3 scripts_cmake/check_version.py          # exit 1 on any mismatch

Source of truth: `project(DiagHam VERSION x.y.z)` in CMakeLists.txt (calendar
versioning: YEAR.MONTH.PATCH, e.g. 2026.9.0). It must equal the `version:` of
CITATION.cff and codemeta.json, and the newest release heading
of CHANGELOG.md (`## [x.y.z] — YYYY-MM-DD`); the release date must agree
between CITATION.cff, codemeta.json and that heading. The authors and the
licence of CITATION.cff must also be those of its preferred-citation and of
.zenodo.json (the Zenodo record's creators are the DiagHam authors). CI runs
this.
"""
import json, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
problems = []

cm = re.search(r"project\(DiagHam\s+VERSION\s+([0-9.]+)", (ROOT / "CMakeLists.txt").read_text())
version = cm.group(1) if cm else None
if not version:
    problems.append("CMakeLists.txt: no project(DiagHam VERSION ...)")

cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
cff_version = re.search(r"^version:\s*\"?([^\"\n]+)\"?", cff, re.M)
cff_date = re.search(r"^date-released:\s*\"?([0-9-]+)\"?", cff, re.M)
if not cff_version or cff_version.group(1).strip() != version:
    problems.append(f"CITATION.cff version {cff_version.group(1) if cff_version else None!r} != {version!r}")

meta = json.loads((ROOT / "codemeta.json").read_text())
if meta.get("version") != version:
    problems.append(f"codemeta.json version {meta.get('version')!r} != {version!r}")
if cff_date and meta.get("dateModified") != cff_date.group(1):
    problems.append(f"codemeta.json dateModified {meta.get('dateModified')!r} != CITATION.cff date-released {cff_date.group(1)!r}")

def cff_people(block):
    return re.findall(r"- family-names: (.+)\n\s+given-names: (.+)", block)

top, _, pref = cff.partition("\npreferred-citation:")
authors = cff_people(top.split("\nauthors:", 1)[1].split("\nlicense:", 1)[0]) if "\nauthors:" in top else []
cff_licence = re.search(r"^license:\s*(\S+)", top, re.M)
if not authors:
    problems.append("CITATION.cff: no authors")
if not pref:
    problems.append("CITATION.cff: no preferred-citation")
else:
    if cff_people(pref) != authors:
        problems.append("CITATION.cff: preferred-citation authors differ from authors")
    pl = re.search(r"^\s+license:\s*(\S+)", pref, re.M)
    if not pl or not cff_licence or pl.group(1) != cff_licence.group(1):
        problems.append("CITATION.cff: preferred-citation licence differs from license")
zen = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
if [c.get("name") for c in zen.get("creators", [])] != [f"{f}, {g}" for f, g in authors]:
    problems.append(".zenodo.json creators differ from CITATION.cff authors")
if cff_licence and zen.get("license") != cff_licence.group(1):
    problems.append(f".zenodo.json license {zen.get('license')!r} != CITATION.cff {cff_licence.group(1)!r}")

log = (ROOT / "CHANGELOG.md").read_text()
heads = re.findall(r"^## \[([^\]]+)\](?: — ([0-9-]+))?", log, re.M)
released = [(v, d) for v, d in heads if v != "Unreleased"]
if not released or released[0][0] != version:
    problems.append(f"CHANGELOG.md newest release heading {released[0][0] if released else None!r} != {version!r}")
elif cff_date and released[0][1] != cff_date.group(1):
    problems.append(f"CHANGELOG.md release date {released[0][1]!r} != CITATION.cff {cff_date.group(1)!r}")

if problems:
    print("version check FAILED:\n  " + "\n  ".join(problems))
    sys.exit(1)
print(f"version {version} ({cff_date.group(1) if cff_date else '?'}) consistent in CMakeLists.txt, CITATION.cff, codemeta.json, CHANGELOG.md; "
      f"{len(authors)} authors and licence consistent in CITATION.cff, its preferred-citation and .zenodo.json")
