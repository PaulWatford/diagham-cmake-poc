#!/usr/bin/env python3
import argparse
import re
import logging
from typing import List, Optional, Dict
import bibtexparser
import codecs

logger = logging.getLogger("bib2wiki")

# -------------------------------
# Utilities
# -------------------------------
# ---- robust LaTeX -> Unicode converter (tries latexcodec, falls back) ----
def latex_to_unicode(s: str) -> str:
    if not s:
        return s
    # try latexcodec first (recommended)
    try:
        import codecs
        return codecs.decode(s, "ulatex+utf8")
    except Exception:
        # fallback: some common accent mappings + compose using Unicode
        import unicodedata
        t = s

        # common named command replacements
        simple_map = {
            r"\\ss": "ß", r"\\ae": "æ", r"\\AE": "Æ", r"\\o": "ø", r"\\O": "Ø",
            r"\\aa": "å", r"\\AA": "Å", r"---": "—", r"--": "–", r"\\&": "&"
        }
        for k, v in simple_map.items():
            t = t.replace(k, v)

        # combining diacritic map: LaTeX accent -> combining mark
        accents = {
            '"': "\u0308",  # diaeresis
            "'": "\u0301",  # acute
            "`": "\u0300",  # grave
            "^": "\u0302",  # circumflex
            "~": "\u0303",  # tilde
            "=": "\u0304",  # macron
            "c": "\u0327",  # cedilla
            "v": "\u030C",  # caron
            "H": "\u030B",  # double acute
            "u": "\u0306",  # breve
            "k": "\u0328",  # ogonek
            "r": "\u030A",  # ring
            "d": "\u0323",  # dot under
        }

        # patterns like \"o  \"{o}  {\\"o}
        for acc, comb in accents.items():
            # forms: \<acc>{x}   or \<acc>x   or {\<acc>x}
            p1 = re.compile(r"\\" + re.escape(acc) + r"\{?([A-Za-z])\}?")
            t = p1.sub(lambda m: unicodedata.normalize("NFC", m.group(1) + comb), t)
            p2 = re.compile(r"\{\\" + re.escape(acc) + r"\s*\{?([A-Za-z])\}?\}")
            t = p2.sub(lambda m: unicodedata.normalize("NFC", m.group(1) + comb), t)

        # remove remaining LaTeX braces and stray backslashes
        t = re.sub(r"[{}]", "", t)
        t = re.sub(r"\\(?=[A-Za-z])", "", t)

        # final normalization
        return unicodedata.normalize("NFC", t)

# ---- helper to strip accents for canonical comparisons ----
def strip_accents(s: str) -> str:
    if not s:
        return s
    import unicodedata
    nfkd = unicodedata.normalize("NFKD", s)
    return "".join(ch for ch in nfkd if not unicodedata.combining(ch))

# ---- normalize title: convert accents first, then remove punctuation and lowercase ----
def normalize_title(title: str) -> str:
    if not title:
        return ""
    t = latex_to_unicode(title)
    t = t.lower()
    # remove math inline and display ($...$ and \[...\])
    t = re.sub(r"\$.*?\$", " ", t)
    t = re.sub(r"\\\[.*?\\\]", " ", t)
    # remove any remaining non-alphanumeric (keep ascii after stripping accents)
    t = strip_accents(t)
    t = re.sub(r"[^a-z0-9]+", " ", t)
    return t.strip()

# ---- normalize authors into a canonical list for comparison ----
def normalize_author_list(author_field: str) -> List[str]:
    if not author_field:
        return []
    # convert LaTeX first, then split on ' and '
    s = latex_to_unicode(author_field)
    parts = [p.strip() for p in re.split(r"\s+and\s+", s) if p.strip()]
    normalized = []
    for p in parts:
        p2 = re.sub(r"[{}]", "", p)  # remove braces
        p2 = p2.replace(",", " ")    # unify "Last, First" vs "First Last"
        p2 = p2.lower()
        # strip accents and non-letters/numbers
        p2 = strip_accents(p2)
        p2 = re.sub(r"[^a-z0-9 ]+", "", p2)
        p2 = " ".join(p2.split())
        if p2:
            normalized.append(p2)
    return normalized

# ---- prefer cleaner value: fewer LaTeX tokens/braces is cleaner ----
def prefer_cleaner(a: Optional[str], b: Optional[str]) -> Optional[str]:
    if not a:
        return b
    if not b:
        return a
    score = lambda s: (s.count("\\") + s.count("{") + s.count("}"))
    return a if score(a) <= score(b) else b

# ---- merge two bib entries with same identifier (DOI or title) ----
def merge_records(existing: Dict, new: Dict) -> Dict:
    core_fields = {
        "title", "journal", "volume", "number", "issue",
        "pages", "year", "month", "doi", "url", "publisher"
    }
    merged = dict(existing)  # shallow copy

    # canonicalize some fields on both sides first
    if "title" in merged:
        merged["title"] = latex_to_unicode(merged["title"])
    if "title" in new:
        new["title"] = latex_to_unicode(new["title"])
    if "author" in merged:
        merged["author"] = latex_to_unicode(merged["author"])
    if "author" in new:
        new["author"] = latex_to_unicode(new["author"])

    # Prefer cleaner values for core fields
    for k in core_fields:
        val_existing = merged.get(k, "")
        val_new = new.get(k, "")
        # if existing empty, take new; otherwise prefer cleaner
        if not val_existing and val_new:
            merged[k] = val_new
        elif val_existing and val_new and val_existing != val_new:
            merged[k] = prefer_cleaner(val_existing, val_new)

    # Authors: if normalized lists are equal, pick cleaner formatting
    a_exist_norm = normalize_author_list(merged.get("author", ""))
    a_new_norm = normalize_author_list(new.get("author", ""))
    if a_exist_norm and a_new_norm and a_exist_norm == a_new_norm:
     merged["author"] = prefer_cleaner(merged.get("author", ""), new.get("author", ""))
    else:
        # if they differ but DOIs same, prefer the entry that has a DOI/publisher info,
        # otherwise keep existing but log (caller should log)
        if new.get("doi") and not merged.get("doi"):
            # prefer new if it has DOI
            merged.update(new)
        # else keep merged's authors (do not overwrite silently)

    # Merge arXiv id if present
    if "_arxiv_id" not in merged and "_arxiv_id" in new:
        merged["_arxiv_id"] = new["_arxiv_id"]
    elif "_arxiv_id" in merged and "_arxiv_id" in new and merged["_arxiv_id"] != new["_arxiv_id"]:
        # keep the shorter canonical id
        merged["_arxiv_id"] = min(merged["_arxiv_id"], new["_arxiv_id"], key=len)

    # Merge non-core fields: note, abstract, keywords, file, affiliation
    extra_fields_to_append = {"note", "abstract", "keywords", "file", "affiliation"}
    for k, v in new.items():
        if not v:
            continue
        if k in core_fields or k in {"author", "title", "_arxiv_id"}:
            continue
        if k not in merged or not merged[k]:
            merged[k] = v
        else:
            # if both exist and differ for certain extra fields, append unique part
            if k in extra_fields_to_append and merged[k] != v:
                merged[k] = merged[k] + "; " + v

    return merged

def get_arxiv_id(entry: Dict) -> Optional[str]:
    """Extract canonical arXiv identifier from a bib entry."""
    arxiv_id = None

    # 1) eprint field
    if "eprint" in entry and entry.get("archiveprefix", "").lower() == "arxiv":
        arxiv_id = entry["eprint"]

    # 2) preprint field
    if not arxiv_id and "preprint" in entry and "arxiv" in entry["preprint"].lower():
        token = entry["preprint"].split()[-1]
        arxiv_id = token.split(":")[-1]

    # 3) url field
    if not arxiv_id and "url" in entry:
        m = re.search(r"arxiv\.org/(?:abs|pdf)/([^?#/]+)", entry["url"], re.I)
        if m:
            arxiv_id = m.group(1)

    # 4) journal field
    if not arxiv_id and "journal" in entry:
        m = re.search(r"arxiv\.org/(?:abs|pdf)/([^?#/]+)", entry["journal"], re.I)
        if not m:
            # handle journal="arXiv:2504.20139"
            m = re.match(r"arxiv[:/]\s*(\d{4}\.\d+)", entry["journal"], re.I)
        if m:
            arxiv_id = m.group(1)
            entry["journal"] = "arXiv.org"
            if "eprint" not in entry:
                entry["eprint"] = arxiv_id

    if not arxiv_id:
        return None

    # canonicalize: strip "arxiv:" prefix and version numbers
    arxiv_id = re.sub(r"(?i)^arxiv:", "", arxiv_id)
    arxiv_id = re.sub(r"v\d+$", "", arxiv_id).strip()

    # update entry's eprint if it still has prefix
    if "eprint" not in entry or entry["eprint"].lower().startswith("arxiv:"):
        entry["eprint"] = arxiv_id

    return arxiv_id


def year_key(e: Dict) -> int:
    """Return publication year, or infer from arXiv ID if missing."""
    y = e.get("year")
    if y:
        try:
            return int(y)
        except (TypeError, ValueError):
            pass

    arx_id = e.get("_arxiv_id")
    if arx_id:
        # strip leftover prefix
        arx_id = re.sub(r"(?i)^arxiv:", "", arx_id).strip()

        # short style: YYMM.NNNNN
        m = re.match(r"^(\d{2})\d{2}\.\d+", arx_id)
        if m:
            yy = int(m.group(1))
            return 2000 + yy if yy < 70 else 1900 + yy
        # legacy cond-mat style
        m2 = re.match(r".*/(\d{2})\d{4,}", arx_id)
        if m2:
            yy = int(m2.group(1))
            return 2000 + yy if yy < 70 else 1900 + yy

    return 0


def format_authors(authors_raw: str) -> str:
    """Format author field in APS style (first initials + last name)."""
    if not authors_raw:
        return ""

    # decode LaTeX accents
    authors_raw = latex_to_unicode(authors_raw)

    # split on 'and'
    parts = [p.strip() for p in re.split(r"\s+and\s+", authors_raw) if p.strip()]
    formatted = []

    for p in parts:
        # remove braces
        p = re.sub(r"[{}]", "", p).strip()

        if "," in p:
            # format "Last, First Middle"
            last, first = [s.strip() for s in p.split(",", 1)]
        else:
            # assume "First Middle Last"
            tokens = p.split()
            if len(tokens) > 1:
                last = tokens[-1]
                first = " ".join(tokens[:-1])
            else:
                last = tokens[0]
                first = ""

        # compress first/middle names into initials
        initials = " ".join(f"{t[0]}." for t in first.split() if t)
        name = f"{initials} {last}".strip()
        formatted.append(name)

    # join authors with commas, last one with 'and'
    if len(formatted) > 2:
        return ", ".join(formatted[:-1]) + ", and " + formatted[-1]
    elif len(formatted) == 2:
        return f"{formatted[0]} and {formatted[1]}"
    elif formatted:
        return formatted[0]
    else:
        return ""



def format_entry(e: Dict, number: int, show_notes: bool = False) -> str:
    """Format one entry for MediaWiki output (APS style)."""
    authors = latex_to_unicode(format_authors(e.get("author", "")))
    title = latex_to_unicode(e.get("title", ""))
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

    # ---- replacement collection loop using separate seen_published / seen_arxiv dicts ----
    seen_published: Dict[str, Dict] = {}
    seen_arxiv: Dict[str, Dict] = {}
    published_entries = []
    arxiv_entries = []
    #published_entries: List[Dict] = []
    #arxiv_entries: List[Dict] = []


    for e in db.entries:
    	# canonicalize some fields early so identifier comparisons work
    	if e.get("title"):
    		e["title"] = latex_to_unicode(e["title"])
    	if e.get("author"):
    		e["author"] = latex_to_unicode(e["author"])
    	# extract arXiv id if present
    	arxiv_id = get_arxiv_id(e)
    	if arxiv_id:
    		e["_arxiv_id"] = arxiv_id
    		# dedupe by arXiv id
    		if arxiv_id in seen_arxiv:
    			seen_arxiv[arxiv_id] = merge_records(seen_arxiv[arxiv_id], e)
    			logger.info("Merged duplicate arXiv:%s into existing arXiv record", arxiv_id)
    		else:
    			seen_arxiv[arxiv_id] = e
    	else:
    		# published or other entry -- dedupe by DOI (preferred) else normalized title
    		doi = e.get("doi", "").lower().strip()
    		title_norm = normalize_title(e.get("title", ""))
    		identifier = doi or title_norm
    		if not identifier:
    			# if still empty, skip with a warning
    			logger.warning("Skipping entry with no DOI and empty title: %s", e.get("title", "[untitled]"))
    			continue
    		if identifier in seen_published:
    			# merge (prefer cleaner core fields)
    			seen_published[identifier] = merge_records(seen_published[identifier], e)
    			logger.info("Merged duplicate published entry for identifier: %s", identifier)
    		else:
    			seen_published[identifier] = e

    # build lists from dicts
    published_entries = list(seen_published.values())
    arxiv_entries = list(seen_arxiv.values())

    # after building published_entries + arxiv_entries
    for arx in arxiv_entries:
        arxiv_id = arx.get("_arxiv_id")
        if arxiv_id and not arx.get("url"):
            arx["url"] = f"https://arxiv.org/abs/{arxiv_id}"
            logger.info("Auto-generated URL for arXiv:%s", arxiv_id)


    # ---- existing merge stage: attach arXiv to published if requested ----
    if merge_arxiv_with_published:
    	still_unmatched = []
    	for arx in arxiv_entries:
    		matched = False
    		arxiv_title = normalize_title(arx.get("title", ""))
    		arxiv_id = arx.get("_arxiv_id")

    		for i, pub in enumerate(published_entries):
    			pub_title_norm = normalize_title(pub.get("title", ""))
    			# match by DOI if pub has DOI and arx has doi, OR by normalized title OR by preprint field
    			if (
    				(pub.get("doi") and arx.get("doi") and pub["doi"].lower() == arx["doi"].lower())
    				or (pub_title_norm and pub_title_norm == arxiv_title)
    				or (pub.get("preprint") and arxiv_id and arxiv_id in pub.get("preprint", ""))
    			):
    				# merge arXiv into published record
    				published_entries[i] = merge_records(pub, arx)
    				logger.info("Merged arXiv:%s into published '%s'", arxiv_id, pub.get("title", "[untitled]"))
    				matched = True
    				break

    		if not matched:
    			still_unmatched.append(arx)

    	arxiv_entries = still_unmatched

    # warn about standalone arXiv
    for arx in arxiv_entries:
        logger.warning("Retaining standalone arXiv:%s (%s)", arx.get("_arxiv_id"), arx.get("title", "[untitled]"))
        if True:
            print("DEBUG arXiv retained:",
              f"_arxiv_id={arx.get('_arxiv_id')!r},",
              f"eprint={arx.get('eprint')!r},",
              f"year={arx.get('year')!r},",
              f"title={arx.get('title', '[untitled]')!r},",
              f"journal={arx.get('journal')!r}")



    all_entries = published_entries + arxiv_entries

    # fix missing year entries for arxiv records:
    for e in all_entries:
        if not e.get("year") and e.get("_arxiv_id"):
            e["_year_inferred"] = year_key(e)
        else:
            e["_year_inferred"] = int(e.get("year", 0))


    # Sort by year, normalised via year_key method
    all_entries.sort(key=lambda e: e["_year_inferred"], reverse=reverse)


    output_lines = []
    last_year = None
    for i, e in enumerate(all_entries, start=1):
    	if year_headings:
    		year = e.get("_year_inferred", "")
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
