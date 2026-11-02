#!/usr/bin/env python3
"""damage_scan.py: the Damage scan of the review pipeline (AGENTS.md; docs/explanation/guard-design.md).

Compares two commits of this repository and reports changes to the locked artefacts that a green
ctest run would not reveal. It reports; it never changes anything, and it never proposes a fix.

    python3 scripts_cmake/damage_scan.py BASE HEAD [--repo DIR] [--markdown FILE] [--json FILE] [--quiet]

Exit codes: 0 no damage found, 2 damage found, 1 the scan could not run (unknown commit, not a repository).

Damage classes (AGENTS.md, "Damage scan"):
  tests-removed      a test name listed in tests/manifest.txt or tests/manifest-core-programs.txt at BASE is not
                     listed at HEAD; reported even when the count is unchanged (names dropped, new names filling
                     the count back up: AGENTS.md "What an AI agent must never do here")
  checks-skipped     enable_testing() or add_test() calls removed; a registration file (tests/*.cmake) deleted;
                     a locked test name added to tests/manifest-optional.txt (which lets it be absent);
                     DIAGHAM_BUILD_TESTS no longer defaults ON
  demo-unregistered  a locked Hubbard golden test (physics.hubbard.2x2.*: the analytic -4*sqrt(2) ground state
                     and its companions) is gone
  golden-edited      a file under tests/data/ or a reference spectrum under benchmarks/ modified, renamed or
                     deleted; a line carrying an expected value or tolerance changed in tests/CMakeLists.txt or
                     tests/*.cmake
  option-default     an option() in CMakeLists.txt or cmake/*.cmake defaulted ON at BASE and OFF at HEAD, or a
                     DIAGHAM_BUILD_* option was removed (a DiagHam tree silently unbuilt)
  in-source-build    the in-source build guard in CMakeLists.txt is gone
  patches-changed    anything under patches/ changed (the frozen audit trail of the carried fixes)

Information, not damage: test names added; options added or now ON; WILL_FAIL registrations added.
The scan reads the two commits with git only; it does not run cmake, ctest or any program, and it does not
execute anything from either commit. Run it from a trusted checkout (CI runs the copy on the base branch).
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

LOCKED_TEST_PREFIXES = ("physics.hubbard.2x2.",)
MANIFESTS = ("tests/manifest.txt", "tests/manifest-core-programs.txt")
OPTIONAL_MANIFEST = "tests/manifest-optional.txt"
GOLDEN_DIRS = ("tests/data/", "benchmarks/")
GOLDEN_EXCLUDE = re.compile(r"(^|/)(__pycache__/|BENCHMARK\.md$|.*\.py$)")
PATCH_DIR = "patches/"
OPTION_RE = re.compile(r'^\s*option\(\s*([A-Za-z0-9_]+)\s+"(?:[^"\\]|\\.)*"\s+(ON|OFF)\s*\)', re.M)
INSOURCE_RE = re.compile(r"CMAKE_SOURCE_DIR\s+STREQUAL\s+CMAKE_BINARY_DIR")
VALUE_LINE_RE = re.compile(r"-?\d+\.\d+(?:[eE][-+]?\d+)?|\b(?:CHECK|EXPECT\w*|TOL\w*|rel:)")
WILL_FAIL_RE = re.compile(r"\bWILL_FAIL\b")
CLASS_ORDER = ["tests-removed", "demo-unregistered", "checks-skipped", "golden-edited",
               "option-default", "in-source-build", "patches-changed"]
CLASS_WORDS = {"tests-removed": "tests removed", "demo-unregistered": "demo unregistered",
               "checks-skipped": "checks skipped", "golden-edited": "golden edited",
               "option-default": "option default flipped", "in-source-build": "in-source build",
               "patches-changed": "patches changed"}


class Scan:
    def __init__(self, repo, base, head):
        self.repo = Path(repo)
        self.base = base
        self.head = head
        self.findings = []   # (class, detail)
        self.info = []
        self.tests_added = []
        self.tests_removed = []

    # -- git -----------------------------------------------------------------
    def git(self, *args, check=True):
        p = subprocess.run(["git", *args], cwd=self.repo, capture_output=True, text=True)
        if check and p.returncode != 0:
            raise RuntimeError(f"git {' '.join(args)}: {p.stderr.strip()}")
        return p.stdout

    def show(self, rev, path):
        p = subprocess.run(["git", "show", f"{rev}:{path}"], cwd=self.repo, capture_output=True, text=True)
        return p.stdout if p.returncode == 0 else None

    def tree(self, rev, prefix):
        out = self.git("ls-tree", "-r", "--name-only", rev, "--", prefix, check=False)
        return [l for l in out.splitlines() if l]

    def damage(self, cls, detail):
        self.findings.append((cls, detail))

    # -- helpers -------------------------------------------------------------
    @staticmethod
    def names(text):
        if text is None:
            return set()
        return {l.strip() for l in text.splitlines() if l.strip() and not l.lstrip().startswith("#")}

    @staticmethod
    def fmt_list(items, limit=50):
        items = list(items)
        shown = ", ".join(i if i.startswith("`") else f"`{i}`" for i in items[:limit])
        if len(items) > limit:
            shown += f", and {len(items) - limit} more"
        return shown

    def options(self, rev):
        files = ["CMakeLists.txt"] + [f for f in self.tree(rev, "cmake") if f.endswith(".cmake")]
        found = {}
        for f in files:
            text = self.show(rev, f)
            if text:
                for name, default in OPTION_RE.findall(text):
                    found.setdefault(name, (default, f))
        return found

    def registration_files(self):
        files = {"CMakeLists.txt", "tests/CMakeLists.txt"}
        for rev in (self.base, self.head):
            files.update(f for f in self.tree(rev, "tests") if f.endswith(".cmake"))
        return sorted(files)

    def count_calls(self, rev, path):
        text = self.show(rev, path)
        if text is None:
            return None
        text = re.sub(r"#[^\n]*", "", text)   # a commented-out call is a removed call
        return (len(re.findall(r"\benable_testing\s*\(", text)),
                len(re.findall(r"\badd_test\s*\(", text)),
                len(WILL_FAIL_RE.findall(text)))

    # -- checks --------------------------------------------------------------
    def check_manifests(self):
        for m in MANIFESTS:
            base_text, head_text = self.show(self.base, m), self.show(self.head, m)
            b, h = self.names(base_text), self.names(head_text)
            if base_text is not None and head_text is None:
                self.damage("tests-removed", f"`{m}` deleted ({len(b)} locked names)")
                continue
            removed, added = sorted(b - h), sorted(h - b)
            if m == MANIFESTS[0]:
                self.tests_removed, self.tests_added = removed, added
            if removed:
                note = " (the count is unchanged: names were replaced)" if len(b) == len(h) else ""
                self.damage("tests-removed", f"`{m}`: {len(removed)} name(s) gone{note}: " + self.fmt_list(removed))
            locked_gone = [n for n in removed if n.startswith(LOCKED_TEST_PREFIXES)]
            if locked_gone:
                self.damage("demo-unregistered",
                            f"locked Hubbard golden test(s) gone from `{m}`: " + self.fmt_list(locked_gone))
            if added and m == MANIFESTS[0]:
                self.info.append(f"{len(added)} test name(s) added to `{m}`")
        # names moved to the optional list become allowed to be absent
        ob = self.names(self.show(self.base, OPTIONAL_MANIFEST))
        oh = self.names(self.show(self.head, OPTIONAL_MANIFEST))
        locked_base = self.names(self.show(self.base, MANIFESTS[0]))
        moved = sorted((oh - ob) & locked_base)
        if moved:
            self.damage("checks-skipped", f"{len(moved)} locked test name(s) added to `{OPTIONAL_MANIFEST}` "
                        "(may now be absent): " + self.fmt_list(moved))

    def check_registrations(self):
        for f in self.registration_files():
            cb, ch = self.count_calls(self.base, f), self.count_calls(self.head, f)
            if cb is None:
                continue
            if ch is None:
                self.damage("checks-skipped", f"`{f}` deleted ({cb[1]} add_test call(s))")
                continue
            if ch[0] < cb[0]:
                self.damage("checks-skipped", f"`{f}`: enable_testing() removed ({cb[0]} -> {ch[0]})")
            if ch[1] < cb[1]:
                self.damage("checks-skipped", f"`{f}`: add_test() calls {cb[1]} -> {ch[1]}")
            if ch[2] > cb[2]:
                self.info.append(f"`{f}`: WILL_FAIL registrations {cb[2]} -> {ch[2]} (open-finding tests added; "
                                 "a physics test turned into one would look the same, so read the diff)")

    def check_goldens(self):
        status = self.git("diff", "--name-status", "-M", self.base, self.head).splitlines()
        golden_hits, patch_hits = set(), set()
        for line in status:
            parts = line.split("\t")
            st, path = parts[0], parts[-1]
            old = parts[1] if st.startswith("R") and len(parts) > 2 else path
            for p in {old, path}:
                if p.startswith(PATCH_DIR):
                    patch_hits.add(f"{st[0]} {p}")
            if old.startswith(GOLDEN_DIRS) and not GOLDEN_EXCLUDE.search(old) and st[0] in "MDR":
                golden_hits.add(f"{st[0]} {old}")
        if golden_hits:
            self.damage("golden-edited", f"{len(golden_hits)} reference file(s) modified, renamed or deleted: "
                        + self.fmt_list(sorted(golden_hits)))
        if patch_hits:
            self.damage("patches-changed", f"{len(patch_hits)} file(s) under `patches/`: "
                        + self.fmt_list(sorted(patch_hits)))
        # expected-value lines changed or removed in the registration files
        reg = [f for f in self.registration_files() if f.startswith("tests/")]
        diff = self.git("diff", "-U0", "--no-color", self.base, self.head, "--", *reg, check=False)
        cur, hits = None, []
        for line in diff.splitlines():
            if line.startswith("--- "):
                cur = None
            elif line.startswith("+++ b/"):
                cur = line[6:]
            elif line.startswith("-") and not line.startswith("---") and cur and VALUE_LINE_RE.search(line):
                s = line[1:].strip()
                if s and not s.startswith("#"):
                    hits.append(f"`{cur}`: `{s[:110]}`")
        if hits:
            self.damage("golden-edited", f"{len(hits)} line(s) carrying a value or tolerance changed or removed "
                        "in the test registrations: " + self.fmt_list(hits, 12))

    def check_options(self):
        ob, oh = self.options(self.base), self.options(self.head)
        for name, (default, f) in sorted(ob.items()):
            if name not in oh:
                if name.startswith("DIAGHAM_BUILD_"):
                    self.damage("option-default", f"`{name}` (was {default} in `{f}`) removed")
                else:
                    self.info.append(f"option `{name}` removed (was {default})")
            elif default == "ON" and oh[name][0] == "OFF":
                cls = "checks-skipped" if name == "DIAGHAM_BUILD_TESTS" else "option-default"
                self.damage(cls, f"`{name}` now defaults OFF (was ON) in `{oh[name][1]}`")
            elif default == "OFF" and oh[name][0] == "ON":
                self.info.append(f"option `{name}` now defaults ON (was OFF)")
        for name in sorted(set(oh) - set(ob)):
            self.info.append(f"option `{name}` added, default {oh[name][0]}")

    def check_insource(self):
        head = self.show(self.head, "CMakeLists.txt")
        if head is None:
            self.damage("in-source-build", "`CMakeLists.txt` missing at HEAD")
            return
        m = INSOURCE_RE.search(head)
        guarded = bool(m) and "FATAL_ERROR" in head[m.end():m.end() + 400]
        if not guarded:
            self.damage("in-source-build", "the in-source build guard (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR "
                        "-> FATAL_ERROR) is not in `CMakeLists.txt` at HEAD")

    # -- run -----------------------------------------------------------------
    def run(self):
        # resolve both ends to commit ids first: a ref that is also a path name (a branch called
        # "patches", say) would otherwise make git's revision/path parsing ambiguous
        self.base = self.git("rev-parse", "--verify", "--quiet", f"{self.base}^{{commit}}").strip()
        self.head = self.git("rev-parse", "--verify", "--quiet", f"{self.head}^{{commit}}").strip()
        self.check_manifests()
        self.check_registrations()
        self.check_goldens()
        self.check_options()
        self.check_insource()
        return self

    def classes(self):
        present = {c for c, _ in self.findings}
        return [c for c in CLASS_ORDER if c in present]

    def diffstat(self):
        short = self.git("diff", "--shortstat", self.base, self.head, check=False).strip()
        stat = self.git("diff", "--stat=100,80", self.base, self.head, check=False).splitlines()
        return short, stat[:40], max(0, len(stat) - 40)

    def markdown(self):
        b = self.git("rev-parse", "--short", self.base).strip()
        h = self.git("rev-parse", "--short", self.head).strip()
        out = [f"## Damage scan: `{b}` to `{h}`", ""]
        cls = self.classes()
        out.append("**Damage:** " + ("none" if not cls else "; ".join(CLASS_WORDS[c] for c in cls)))
        out.append("")
        if self.findings:
            out += ["| Class | Finding |", "|---|---|"]
            out += [f"| {c} | {d} |" for c, d in self.findings]
            out += ["", "A flagged change is held for a named maintainer. Nothing here is corrected automatically "
                    "(AGENTS.md: a golden is never rewritten to make a run green)."]
        out += ["", f"**Test names:** {len(self.tests_added)} added, {len(self.tests_removed)} removed "
                "(`tests/manifest.txt`)."]
        if self.tests_added:
            out += ["", "Added: " + self.fmt_list(self.tests_added, 30)]
        if self.tests_removed:
            out += ["", "Removed: " + self.fmt_list(self.tests_removed, 30)]
        if self.info:
            out += ["", "**Information (not damage):**"] + [f"- {i}" for i in self.info]
        short, stat, more = self.diffstat()
        out += ["", f"**Diffstat:** {short or 'no changes'}"]
        if stat:
            out += ["", "```"] + stat + ([f"... and {more} more lines"] if more else []) + ["```"]
        return "\n".join(out) + "\n"

    def as_json(self):
        short, _, _ = self.diffstat()
        return {"base": self.base, "head": self.head, "damage": bool(self.findings), "classes": self.classes(),
                "findings": [{"class": c, "detail": d} for c, d in self.findings],
                "tests_added": self.tests_added, "tests_removed": self.tests_removed,
                "info": self.info, "diffstat": short}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("base")
    ap.add_argument("head")
    ap.add_argument("--repo", default=".")
    ap.add_argument("--markdown", help="write the report here (it is also printed unless --quiet)")
    ap.add_argument("--json", help="write the machine-readable result here")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    try:
        scan = Scan(a.repo, a.base, a.head).run()
    except RuntimeError as e:
        print(f"damage_scan: {e}", file=sys.stderr)
        return 1
    md = scan.markdown()
    if a.markdown:
        Path(a.markdown).write_text(md, encoding="utf-8")
    if a.json:
        Path(a.json).write_text(json.dumps(scan.as_json(), indent=2), encoding="utf-8")
    if not a.quiet:
        sys.stdout.write(md)
    return 2 if scan.findings else 0


if __name__ == "__main__":
    sys.exit(main())
