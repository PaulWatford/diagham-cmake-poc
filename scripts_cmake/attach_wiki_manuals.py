#!/usr/bin/env python3
"""Convert the archived DiagHam wiki program pages to Markdown manuals and attach them.

Usage: attach_wiki_manuals.py <wiki-raw-dir> <docs/reference/programs> <docs/drafts>

- A wiki page whose name is a built program becomes manuals/<Program>.md.
- Renamed programs (RENAMED map) are attached under the new name with a note.
- Stub pages (< 200 bytes of wikitext) are parked in drafts/manual-stubs/.
- Each manual's documented --options are checked against the generated
  reference page (LAPACK build); options not present get a banner.
- Every manual carries a provenance line.
"""
import re, sys, os
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


def wikitext_to_md(t, program_pages):
    out = []
    in_code = False
    in_table = False
    # <math> blocks (possibly multi-line) are stashed first so literal dollars
    # elsewhere can be escaped without touching them
    maths = []
    def stash(m):
        maths.append(m.group(1).strip())
        return f"MATHPLACEHOLDER{len(maths) - 1}X"
    t = re.sub(r"<math>(.*?)</math>", stash, t.replace("\r\n", "\n"), flags=re.S)
    for line in t.split("\n"):
        if re.match(r"^=+$", line.strip()):
            out.append("---")
            continue
        # indented lines are preformatted in MediaWiki
        if line.startswith("    ") or line.startswith(" ") and line.strip() and not line.lstrip().startswith(("*", "#", "==", "{|", "|", "!")):
            if not in_code:
                out.append("```text"); in_code = True
            out.append(line[4:] if line.startswith("    ") else line.lstrip())
            continue
        if in_code:
            out.append("```"); in_code = False
        s = line
        if s.startswith("{|"):
            in_table = True; rows = []; continue
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
                rows.append([]); continue
            if s.startswith(("|", "!")):
                cells = re.split(r"\s*\|\|\s*|\s*!!\s*", s[1:].strip())
                if not rows: rows.append([])
                rows[-1].extend(cells)
            continue
        s = re.sub(r"<pre>|</pre>", "```", s)
        # literal dollars (shell variables in prose) must not become math delimiters
        s = s.replace("$", r"\$")
        s = re.sub(r"<br\s*/?>", "  ", s)
        m = re.match(r"^(=+)\s*(.*?)\s*\1\s*$", s)
        if m:
            out.append("#" * (len(m.group(1)) + 0) + " " + m.group(2)); continue
        s = re.sub(r"'''(.+?)'''", r"**\1**", s)
        s = re.sub(r"''(.+?)''", r"*\1*", s)
        def wikilink(m):
            target, text = (m.group(1).split("|", 1) + [None])[:2]
            text = text or target
            tgt = target.replace(" ", "_")
            if tgt in program_pages:
                return f"[{text}]({tgt}.md)"
            if tgt in RENAMED:
                return f"[{text}]({RENAMED[tgt]}.md)"
            return f"{text} *(wiki page {target})*"
        s = re.sub(r"\[\[([^\]]+)\]\]", wikilink, s)
        s = re.sub(r"\[(https?://[^\s\]]+)\s+([^\]]+)\]", r"[\2](\1)", s)
        s = re.sub(r"\[(https?://[^\s\]]+)\]", r"<\1>", s)
        s = re.sub(r"^\*\*\*\s", "      - ", s); s = re.sub(r"^\*\*\s", "   - ", s); s = re.sub(r"^\*\s", "- ", s)
        s = re.sub(r"^#\s", "1. ", s)
        out.append(s)
    if in_code:
        out.append("```")
    res = "\n".join(out)
    for i, mth in enumerate(maths):
        res = res.replace(f"MATHPLACEHOLDER{i}X", ("\n$$\n" + mth + "\n$$\n") if "\n" in mth else "$" + mth + "$")
    return res.strip() + "\n"


def options_in(text):
    return set(re.findall(r"--[a-z][a-z0-9-]+", text))


def main(wiki, progs, drafts):
    wiki, progs, drafts = Path(wiki), Path(progs), Path(drafts)
    manuals = progs / "manuals"; manuals.mkdir(exist_ok=True)
    stubs = drafts / "manual-stubs"; stubs.mkdir(parents=True, exist_ok=True)
    gen = {p.stem: p for p in progs.glob("*/*.md") if p.parent.name != "manuals"}
    program_pages = set()
    pages = sorted(wiki.glob("*.wiki"))
    for p in pages:
        n = p.stem
        if n in NOT_PROGRAMS: continue
        if n in gen or RENAMED.get(n) in gen: program_pages.add(RENAMED.get(n, n))
    def relink(text):
        # links written for manuals/ must resolve from drafts/manual-stubs/
        return re.sub(r"\]\(([A-Za-z0-9_]+)\.md\)", r"](../../reference/programs/manuals/\1.md)", text)
    for old in manuals.glob("*.md"):
        old.unlink()
    for old in stubs.glob("*.md"):
        old.unlink()
    attached = parked = renamed_n = banners = 0
    rows = []
    wiki_names = {p.stem for p in pages}
    for p in pages:
        name = p.stem
        if name in NOT_PROGRAMS: continue
        raw = p.read_text(encoding="utf-8", errors="replace")
        target = RENAMED.get(name, name)
        if name != target and target in wiki_names:
            # the wiki also has a page under the new name: that one is the manual; park the old one
            (stubs / f"{name}.md").write_text(f"# {name} (superseded wiki page)\n\nSource: DiagHam wiki page `{name}`, as of {WIKI_DATE}. "
                                             f"Reason: the program is now `{target}` and the wiki has a page under that name, which is the attached manual; this older page is kept for comparison.\n\n"
                                             + relink(wikitext_to_md(raw, program_pages)))
            rows.append((name, f"superseded by the {target} page", len(raw))); parked += 1; continue
        if name in UNCERTAIN or target not in gen:
            (stubs / f"{name}.md").write_text(f"# {name} (wiki page, not attached)\n\nSource: DiagHam wiki page `{name}`, as of {WIKI_DATE}. "
                                             f"Reason: {UNCERTAIN.get(name, 'no program of this name in the r4493 build')}.\n\n"
                                             + relink(wikitext_to_md(raw, program_pages)))
            rows.append((name, "no matching program / uncertain rename", len(raw))); parked += 1; continue
        if len(raw.strip()) < STUB_BYTES:
            (stubs / f"{name}.md").write_text(f"# {name} (stub wiki page)\n\nSource: DiagHam wiki page `{name}`, as of {WIKI_DATE} ({len(raw)} bytes). "
                                             f"Reason: stub — the generated `--help` page is the only reference. To complete: describe purpose, inputs, outputs and a worked example.\n\n"
                                             + relink(wikitext_to_md(raw, program_pages)))
            rows.append((name, "stub", len(raw))); parked += 1; continue
        md = wikitext_to_md(raw, program_pages)
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
        if name != target: renamed_n += 1
        rows.append((name, "attached" + (f" as {target}" if name != target else "") + (f"; {len(missing)} option(s) flagged" if missing else ""), len(raw)))
    with open(stubs / "INDEX.md", "w") as f:
        f.write(f"# Parked wiki program pages\n\nSource: DiagHam wiki, archived {WIKI_DATE}. Each row: page, why it is parked, size. To complete a page: write purpose, inputs, outputs and an example, then move it to `docs/reference/programs/manuals/`.\n\n| Page | Reason | Bytes |\n|---|---|---:|\n")
        for n, why, sz in rows:
            if not why.startswith("attached"):
                f.write(f"| {n} | {why} | {sz} |\n")
    print(f"attached {attached} manuals ({renamed_n} under a new program name, {banners} with option banners); parked {parked} in {stubs}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
