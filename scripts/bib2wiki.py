#!/usr/bin/env python3
import argparse
import re
import logging
from typing import List, Optional, Dict
import bibtexparser

logger = logging.getLogger("bib2wiki")

# -------------------------------
# Utilities
# -------------------------------
def normalize_title(title: str) -> str:
    """Normalize titles for duplicate detection (strip LaTeX and punctuation)."""
    if not title:
        return ""
    t = title.lower()
    # remove LaTeX braces and math delimiters
    t = re.sub(r"[{}]", "", t)
    t = re.sub(r"\$.*?\$", "", t)
    # collapse spaces/punctuation
    t = re.sub(r"[^a-z0-9]+", " ", t)
    return t.strip()


def get_arxiv_id(entry: Dict) -> Optional[str]:
    """Extract canonical arXiv identifier from a bib entry."""
    arxiv_id = None

    if "eprint" in entry and entry.get("archiveprefix", "").lower() == "arxiv":
        arxiv_id = entry["eprint"]

    if not arxiv_id and "preprint" in entry and "arxiv" in entry["preprint"].lower():
        token = entry["preprint"].split()[-1]
        arxiv_id = token.split(":")[-1]

    if not arxiv_id and "url" in entry:
        m = re.search(r"arxiv\.org/(?:abs|pdf)/([^?#/]+)", entry["url"])
        if m:
            arxiv_id = m.group(1)

    # NEW: handle journal = arXiv.org
    if not arxiv_id and entry.get("journal", "").lower() in ["arxiv", "arxiv.org"]:
        if "eprint" in entry:
            arxiv_id = entry["eprint"]

    if not arxiv_id:
        return None

    # canonicalize: strip version and "arXiv:" prefix
    arxiv_id = re.sub(r"(?i)^arxiv:", "", arxiv_id)
    arxiv_id = re.sub(r"v\d+$", "", arxiv_id)
    return arxiv_id.strip()


def format_authors(author_field: str) -> str:
    if not author_field:
        return ""
    authors = [a.strip() for a in re.split(" and ", author_field)]
    return ", ".join(authors)


def format_entry(e: Dict, number: int, show_notes: bool = False) -> str:
    """Format one entry for MediaWiki output (APS style)."""
    authors = format_authors(e.get("author", ""))
    title = e.get("title", "")
    journal = e.get("journal", "")
    volume = e.get("volume", "")
    pages = e.get("pages", e.get("number", ""))
    year = e.get("year", "")
    doi = e.get("doi", "")
    url = e.get("url", "")

    ref = f"# {authors}, ''{title}''"

    if journal and journal.lower() not in ["arxiv", "arxiv.org"]:
        ref += f", {journal}"
        if volume:
            ref += f" '''{volume}'''"
        if pages:
            ref += f", {pages}"
        if year:
            ref += f" ({year})"
    else:
        if year:
            ref += f", arXiv.org ({year})"

    links = []
    if doi:
        links.append(f"[https://doi.org/{doi} doi]")
    elif url and "doi.org" in url:
        links.append(f"[{url} doi]")
    if url and "doi.org" not in url:
        links.append(f"[{url} pub]")
    if "_arxiv_id" in e:
        arxiv_id = e["_arxiv_id"]
        links.append(f"[https://arxiv.org/abs/{arxiv_id} arXiv:{arxiv_id}]")

    if links:
        ref += " " + " ".join(links)

    if show_notes and "note" in e:
        ref += f"\n: {e['note']}"

    return ref


# -------------------------------
# Processing
# -------------------------------
def process_bibtex(
    filename: str,
    merge_arxiv_with_published: bool = True,
    show_notes: bool = False,
    reverse: bool = True,
    year_headings: bool = True
) -> List[str]:
    with open(filename) as bibfile:
        db = bibtexparser.load(bibfile)

    published_entries: List[Dict] = []
    arxiv_entries: List[Dict] = []
    seen_ids = set()

    # Collect entries
    for e in db.entries:
        doi = e.get("doi", "").lower().strip()
        title_norm = normalize_title(e.get("title", ""))

        identifier = doi or title_norm
        if identifier in seen_ids:
            logger.info("Duplicate skipped: %s", e.get("title", "[untitled]"))
            continue
        seen_ids.add(identifier)

        arxiv_id = get_arxiv_id(e)
        if arxiv_id:
            e["_arxiv_id"] = arxiv_id
            arxiv_entries.append(e)
        else:
            # published paper (may lack DOI, but has journal/volume)
            if not doi and (e.get("journal") or e.get("volume")):
                logger.warning(
                    "No DOI for published entry: %s", e.get("title", "[untitled]")
                )
            published_entries.append(e)

    # Merge arXiv with published
    if merge_arxiv_with_published:
        still_unmatched = []
        for arx in arxiv_entries:
            matched = False
            arxiv_title = normalize_title(arx.get("title", ""))
            arxiv_id = arx.get("_arxiv_id")

            for pub in published_entries:
                if (
                    pub.get("doi")
                    and normalize_title(pub.get("title", "")) == arxiv_title
                ):
                    pub["_arxiv_id"] = arxiv_id
                    logger.info(
                        "Merged arXiv:%s into published '%s'",
                        arxiv_id,
                        pub.get("title", "[untitled]"),
                    )
                    matched = True
                    break

            if not matched:
                still_unmatched.append(arx)

        arxiv_entries = still_unmatched

    # Retained standalone arXiv warnings
    for arx in arxiv_entries:
        logger.warning(
            "Retaining standalone arXiv:%s (%s)",
            arx.get("_arxiv_id"),
            arx.get("title", "[untitled]"),
        )

    all_entries = published_entries + arxiv_entries

    # Sort by year ascending (oldest first, numbering starts at 1)
    def year_key(e):
        try:
            return int(e.get("year", 0))
        except ValueError:
            return 0

    # Sort by year
    all_entries.sort(key=year_key, reverse=reverse)

    output_lines = []
    last_year = None
    for i, e in enumerate(all_entries, start=1):
        if year_headings:
            year = e.get("year", "")
            if year != last_year:
                if year:
                    output_lines.append(f"== {year} ==")
                    last_year = year
                    line = format_entry(e, i, show_notes=show_notes)
                    output_lines.append(line)
                    output_lines.append("")  # small space between papers

    return output_lines


# -------------------------------
# CLI
# -------------------------------
def main():
    parser = argparse.ArgumentParser(description="Convert BibTeX to MediaWiki references.")
    parser.add_argument("bibfile", help="Input .bib file")
    parser.add_argument("--no-merge-arxiv", action="store_true", help="Do not merge arXiv with published")
    parser.add_argument("--show-notes", action="store_true", help="Show notes if present")
    parser.add_argument( "--reverse", action="store_false", help="Sort bibliography in reverse chronological order (newest first)")
    parser.add_argument( "--year-headings", action="store_false", help="Insert year headings in the output")

    parser.add_argument("--output", "-o", help="Output file (default: same name as input, with .wiki extension; use '-' for stdout)", default=None)

    parser.add_argument("--verbose", type=int, default=1, choices=[0, 1, 2], help="Verbosity level (0=warnings only, 1=info, 2=debug)")
    args = parser.parse_args()

    level = logging.WARNING
    if args.verbose == 1:
        level = logging.INFO
    elif args.verbose >= 2:
        level = logging.DEBUG
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")

    refs = process_bibtex(
        args.bibfile,
        merge_arxiv_with_published=not args.no_merge_arxiv,
        show_notes=args.show_notes,
        reverse=args.reverse,
        year_headings=args.year_headings
        )

    output_text = "\n\n".join(refs)

    if args.output == "-" or (args.output is None and args.bibfile == "-"):
        # force stdout
        print(output_text)
    else:
        if args.output:
            outfile = args.output
        else:
            outfile = re.sub(r"\.bib$", ".wiki", args.bibfile)
        with open(outfile, "w", encoding="utf-8") as f:
            f.write(output_text + "\n")
        logger.info("Wrote output to %s", outfile)


if __name__ == "__main__":
    main()
