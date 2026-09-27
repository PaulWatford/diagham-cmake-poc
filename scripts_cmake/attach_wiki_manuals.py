#!/usr/bin/env python3
"""Convert the archived DiagHam wiki program pages to Markdown manuals and attach them.

Usage: attach_wiki_manuals.py <wiki-raw-dir> <docs/reference/programs> <docs/drafts>

- A wiki page whose name is a built program becomes manuals/<Program>.md.
- Renamed programs (RENAMED map) are attached under the new name with a note;
  when the wiki also has a page under the new name, that page is the manual
  and the old one is parked as superseded.
- Stub pages (< 200 bytes of wikitext) and pages whose program cannot be
  identified are parked in drafts/manual-stubs/ with an INDEX.
- Each manual's documented --options are checked against the generated
  reference page (LAPACK build); options not present get a banner.
- Wiki links resolve to the manual when one exists, otherwise to the
  program's generated --help page; links in parked pages are rewritten to
  resolve from drafts/manual-stubs/.
- Every manual carries a provenance line.
"""
import re, sys
from pathlib import Path

WIKI_DATE = "2026-09-24"
RENAMED = {
    "FQHETopInsulatorEntanglementEntropyParticlePartition": "FTIEntanglementEntropyParticlePartition",
    "FQHETopInsulatorEntanglementSpectrum": "FTIEntanglementSpectrum",
    "FQHETopInsulatorGetDimension": "FTIGetDimension",
    "FQHETopInsulatorShowBasis": "FTIShowBasis",
    "FQHECheckerboardLatticeModel": "FCICheckerboardLatticeModel",
    "FQHEKagomeLatticeModel": "FCIKagomeLatticeModel",
    "FQHEQuantumSpinHallCheckerboardModelTwoBands": "FQSHCheckerboardModelTwoBands",
    "QHEBosonsTwoBodyGeneric": "FQHESphereBosonsTwoBodyGeneric",
    "QHEFermionsTwoBodyGeneric": "FQHESphereFermionsTwoBodyGeneric",
    "FQHEFermionsTorusWithTranslation": "FQHETorusFermionsWithTranslations",
    "FQHETorusFermionsWithTranslation": "FQHETorusFermionsWithTranslations",
}
UNCERTAIN = {"FQHETorusBosonsDelta": "program renamed upstream; the current equivalent was not identified with certainty"}
NOT_PROGRAMS = {"Articles", "Binary_matrices", "Binary_vectors", "Create_a_new_spin_Hilbert_space_with_a_discrete_symmetry",
                "Create_new_FQHE_code", "FCI_programs", "FQHE_programs", "Install", "Lanczos", "Laughlin_state", "Logo",
                "Long_double", "MPI", "Main_Page", "People", "Scalapack", "Spin_chain_programs"}
STUB_BYTES = 200

LINK_MANUAL = re.compile(r"\]\(([A-Za-z0-9_]+)\.md\)")
LINK_GEN = re.compile(r"\]\(\.\./([A-Za-z0-9_]+)/([A-Za-z0-9_]+)\.md\)")


def wikitext_to_md(t, manual_pages, gen_dirs):
    out = []
    in_code = False
    in_table = False
    rows = []
    # <math> blocks (possibly multi-line) are stashed first so literal dollars
    # elsewhere can be escaped without touching them
    maths = []

    def stash(m):
        maths.append(m.group(1).strip())
        return f"MATHPLACEHOLDER{len(maths) - 1}X"
    t = re.sub(r"<math>(.*?)</math>", stash, t.replace("\r\n", "\n"), flags=re.S)

    def wikilink(m):
        target, text = (m.group(1).split("|", 1) + [None])[:2]
        text = text or target
        tgt = target.replace(" ", "_")
        tgt = RENAMED.get(tgt, tgt)
        if tgt in manual_pages:                       # a manual exists for it
            return f"[{text}]({tgt}.md)"
        if tgt in gen_dirs:                           # no manual: the generated --help page
            return f"[{text}](../{gen_dirs[tgt]}/{tgt}.md)"
        return f"{text} *(wiki page {target})*"

    for line in t.split("\n"):
        if re.match(r"^=+$", line.strip()):
            out.append("---")
            continue
        # indented lines are preformatted in MediaWiki
        if line.startswith("    ") or line.startswith(" ") and line.strip() and not line.lstrip().startswith(("*", "#", "==", "{|", "|", "!")):
            if not in_code:
                out.append("```text")
                in_code = True
            out.append(line[4:] if line.startswith("    ") else line.lstrip())
            continue
        if in_code:
            out.append("```")
            in_code = False
        s = line
        if s.startswith("{|"):
            in_table = True
            rows = []
            continue
        if in_table:
            if s.startswith("|}"):
                in_table = False
                if rows:
                    out.append("| " + " | ".join(rows[0]) + " |")
                    out.append("|" + "---|" * len(rows[0]))
                    for r in rows[1:]:
                        out.append("| " + " | ".join(r) + " |")
                continue
            if s.startswith("|-"):
                rows.append([])
                continue
            if s.startswith(("|", "!")):
                cells = re.split(r"\s*\|\|\s*|\s*!!\s*", s[1:].strip())
                if not rows:
                    rows.append([])
                rows[-1].extend(cells)
            continue
        s = re.sub(r"<pre>|</pre>", "```", s)
        # literal dollars (shell variables in prose) must not become math delimiters
        s = s.replace("$", r"\$")
        s = re.sub(r"<br\s*/?>", "  ", s)
        m = re.match(r"^(=+)\s*(.*?)\s*\1\s*$", s)
        if m:
            out.append("#" * len(m.group(1)) + " " + m.group(2))
            continue
        s = re.sub(r"'''(.+?)'''", r"**\1**", s)
        s = re.sub(r"''(.+?)''", r"*\1*", s)
        s = re.sub(r"\[\[([^\]]+)\]\]", wikilink, s)
        s = re.sub(r"\[(https?://[^\s\]]+)\s+([^\]]+)\]", r"[\2](\1)", s)
        s = re.sub(r"\[(https?://[^\s\]]+)\]", r"<\1>", s)
        s = re.sub(r"^\*\*\*\s", "      - ", s)
        s = re.sub(r"^\*\*\s", "   - ", s)
        s = re.sub(r"^\*\s", "- ", s)
        s = re.sub(r"^#\s", "1. ", s)
        out.append(s)
    if in_code:
        out.append("```")
    res = "\n".join(out)
    for i, mth in enumerate(maths):
        res = res.replace(f"MATHPLACEHOLDER{i}X", ("\n$$\n" + mth + "\n$$\n") if "\n" in mth else "$" + mth + "$")
    return res.strip() + "\n"


def relink_for_drafts(text):
    """Links written for manuals/ must resolve from drafts/manual-stubs/."""
    text = LINK_GEN.sub(r"](../../reference/programs/\1/\2.md)", text)
    return LINK_MANUAL.sub(r"](../../reference/programs/manuals/\1.md)", text)


def options_in(text):
    return set(re.findall(r"--[a-z][a-z0-9-]+", text))


def main(wiki, progs, drafts):
    wiki, progs, drafts = Path(wiki), Path(progs), Path(drafts)
    manuals = progs / "manuals"
    manuals.mkdir(exist_ok=True)
    stubs = drafts / "manual-stubs"
    stubs.mkdir(parents=True, exist_ok=True)
    gen = {p.stem: p for p in progs.glob("*/*.md") if p.parent.name != "manuals"}
    gen_dirs = {name: path.parent.name for name, path in gen.items()}
    pages = sorted(wiki.glob("*.wiki"))
    wiki_names = {p.stem for p in pages}
    raw_of = {p.stem: p.read_text(encoding="utf-8", errors="replace") for p in pages}

    # classify first: which pages become manuals (so links can target them), which are parked
    def classify(name):
        raw = raw_of[name]
        target = RENAMED.get(name, name)
        if name in NOT_PROGRAMS:
            return "skip", target
        if name != target and target in wiki_names:
            return "superseded", target
        if name in UNCERTAIN or target not in gen:
            return "unmatched", target
        if len(raw.strip()) < STUB_BYTES:
            return "stub", target
        return "manual", target
    decisions = {n: classify(n) for n in raw_of}
    manual_pages = {t for n, (k, t) in decisions.items() if k == "manual"}

    for old in manuals.glob("*.md"):
        old.unlink()
    for old in stubs.glob("*.md"):
        old.unlink()
    attached = parked = renamed_n = banners = 0
    rows = []
    for name in sorted(raw_of):
        kind, target = decisions[name]
        raw = raw_of[name]
        if kind == "skip":
            continue
        md = wikitext_to_md(raw, manual_pages, gen_dirs)
        if kind == "superseded":
            (stubs / f"{name}.md").write_text(
                f"# {name} (superseded wiki page)\n\nSource: DiagHam wiki page `{name}`, as of {WIKI_DATE}. "
                f"Reason: the program is now `{target}` and the wiki has a page under that name, which is the attached manual; this older page is kept for comparison.\n\n"
                + relink_for_drafts(md))
            rows.append((name, f"superseded by the {target} page", len(raw)))
            parked += 1
            continue
        if kind == "unmatched":
            (stubs / f"{name}.md").write_text(
                f"# {name} (wiki page, not attached)\n\nSource: DiagHam wiki page `{name}`, as of {WIKI_DATE}. "
                f"Reason: {UNCERTAIN.get(name, 'no program of this name in the r4493 build')}.\n\n" + relink_for_drafts(md))
            rows.append((name, "no matching program / uncertain rename", len(raw)))
            parked += 1
            continue
        if kind == "stub":
            (stubs / f"{name}.md").write_text(
                f"# {name} (stub wiki page)\n\nSource: DiagHam wiki page `{name}`, as of {WIKI_DATE} ({len(raw)} bytes). "
                f"Reason: stub — the generated `--help` page is the only reference. To complete: describe purpose, inputs, outputs and a worked example.\n\n"
                + relink_for_drafts(md))
            rows.append((name, "stub", len(raw)))
            parked += 1
            continue
        help_text = gen[target].read_text()
        missing = sorted(options_in(raw) - options_in(help_text))
        head = [f"# {target} — manual", "",
                f"Source: DiagHam wiki page `{name}`, as of {WIKI_DATE} (archived copy). Changed: wikitext converted to Markdown; options checked against the program's current `--help` (r4493, LAPACK build)."
                + (f" The wiki page was named `{name}`; the program is now `{target}`." if name != target else ""), "",
                f"Generated option reference: [{target}](../{gen[target].parent.name}/{target}.md)", ""]
        if missing:
            head += ["> **Options in this manual that the current program does not list**: " + ", ".join(f"`{o}`" for o in missing)
                     + ". They may have been renamed, removed, or depend on a build option not enabled here (GMP, MPI). Trust `--help`.", ""]
            banners += 1
        (manuals / f"{target}.md").write_text("\n".join(head) + md)
        attached += 1
        if name != target:
            renamed_n += 1
        rows.append((name, "attached" + (f" as {target}" if name != target else "") + (f"; {len(missing)} option(s) flagged" if missing else ""), len(raw)))
    with open(stubs / "INDEX.md", "w") as f:
        f.write(f"# Parked wiki program pages\n\nSource: DiagHam wiki, archived {WIKI_DATE}. Each row: page, why it is parked, size. To complete a page: write purpose, inputs, outputs and an example, then move it to `docs/reference/programs/manuals/`.\n\n| Page | Reason | Bytes |\n|---|---|---:|\n")
        for n, why, sz in rows:
            if not why.startswith("attached"):
                f.write(f"| [{n}]({n}.md) | {why} | {sz} |\n")
    print(f"attached {attached} manuals ({renamed_n} under a new program name, {banners} with option banners); parked {parked} in {stubs}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
