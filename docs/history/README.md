# docs/history — upstream files kept as found

Purpose: say what each preserved upstream file is and which living page replaces it.
Source: the DiagHam SVN tree at r4493 (2026-09-18); this index is new (2026-09-25).

These files are preserved unchanged from the DiagHam Subversion tree (r4493)
for the record. They are not maintained; the living equivalents are named.

| File | What it is | Living equivalent |
|---|---|---|
| `INSTALL` | the generic GNU autotools installation text (2003) | `docs/how-to/build/` and `README.md` |
| `ChangeLog-upstream` | DiagHam's ChangeLog, 2003-05-03 to 2005-04-26, generated from CVS by `scripts/docs/CreateChangelog.pl`; not updated since | `CHANGELOG.md`, and the Subversion commit log: `git log upstream` (4,477 trunk commits to r4493) |
| `TODO-2007` | the maintainers' to-do list dated 2007-11-29 (GM and NR) | none; kept for context ("Torus with spin" is still on it) |
| `NEWS` | one line ("Zlatko's first commit", 2009) | none |

The upstream `README` was empty and has been removed. `AUTHORS` and `COPYING`
stay at the repository root (AUTHORS updated from the commit log).

Source of the SVN history: https://www.nick-ux.org/diagham/svn/DiagHam
(trunk, branches, tags). The canonical DiagHam documentation is the wiki at
https://www.nick-ux.org/diagham/; its pages are being brought into `docs/`
with a provenance line each.

## Legacy documentation tooling in `scripts/docs/`

`ProgXMLDoc.pl`, `OldProgXMLDoc.pl` and `XMLDoc2{Html,Pdf,Php,Tex}.pl`
extracted program documentation from an XML description into several
formats, and `CreateChangelog.pl` generated the ChangeLog from CVS; none
has been used since the mid-2000s. They stay in `scripts/` untouched. The
living equivalent is the generated program reference
(`scripts_cmake/gen_program_reference.py` → `docs/reference/programs/`).
