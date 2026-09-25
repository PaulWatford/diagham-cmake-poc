#!/usr/bin/env python3
import argparse
import re
import os
import logging
from typing import List, Optional, Dict
import bibtexparser
import codecs
import requests
import feedparser
import time
from pathlib import Path


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


def format_entry(e: dict, index: Optional[int] = None, show_notes: bool = False) -> str:
    """Format a single BibTeX entry into MediaWiki markup."""
    # --- Author list ---
    authors = e.get("author", "")
    if authors:
        names = [a.strip() for a in re.split(r"\s+and\s+", authors)]
        formatted_authors = []
        for n in names:
            parts = n.replace(",", "").split()
            if len(parts) == 1:
                formatted_authors.append(parts[0])
            else:
                formatted_authors.append(f"{parts[-1]}, {' '.join(parts[:-1])}")
        if len(formatted_authors) == 1:
            author_str = formatted_authors[0]
        elif len(formatted_authors) == 2:
            author_str = " and ".join(formatted_authors)
        else:
            author_str = ", ".join(formatted_authors[:-1]) + ", and " + formatted_authors[-1]
    else:
        author_str = ""

    # --- Core metadata ---
    title = e.get("title", "[untitled]").strip().rstrip(".")
    journal = e.get("journal", "").strip()
    vol = e.get("volume", "").strip()
    pages = e.get("pages", "").strip()
    year = e.get("year", "").strip()
    doi = e.get("doi", "").strip()
    url = e.get("url", "").strip()
    arxiv_id = e.get("_arxiv_id") or e.get("eprint", "")
    note = e.get("note", "").strip()

    # --- Build main citation line ---
    line = f"{author_str}, ''{title}''"
    if journal:
        line += f", {journal}"
    if vol:
        line += f" '''{vol}'''"
    if pages:
        line += f", {pages}"
    if year:
        line += f" ({year})"

    # --- Build link list ---
    links = []
    has_publisher = bool(doi or (url and "arxiv.org" not in url.lower()))

    # DOI link
    if doi:
        links.append(f"[https://doi.org/{doi} doi]")

    # Publisher link (only if it's *not* arXiv)
    if url and "arxiv.org" not in url.lower():
        links.append(f"[{url} pub]")

    # arXiv link (avoid duplicate 'pub' pointing to arxiv)
    if arxiv_id:
        display_id = re.sub(r"(?i)^arxiv:", "", arxiv_id).strip()
        clean_id = re.sub(r"v\d+$", "", display_id)
        arxiv_link = f"https://arxiv.org/abs/{clean_id}"
        links.append(f"[{arxiv_link} arXiv:{display_id}]")

    # --- Combine into final output ---
    out = f"# {line}"
    if links:
        out += " " + " ".join(links)
    if show_notes and note:
        out += f" ({note})"

    return out


def format_entry_old(e: Dict, number: int, show_notes: bool = False) -> str:
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
    # --- Build link list ---
    links = []

    doi = e.get("doi", "").strip()
    url = e.get("url", "").strip()
    arxiv_id = e.get("_arxiv_id") or e.get("eprint", "")
    is_arxiv = e.get("journal", "").lower() in ["arxiv", "arxiv.org"]

    # Add DOI if available
    if doi:
        links.append(f"[https://doi.org/{doi} doi]")

    # Add publisher URL only if it's not the same as the arXiv link
    if url and ("arxiv.org" not in url.lower()):
        links.append(f"[{url} pub]")

    # Add arXiv link
    if arxiv_id:
        # canonicalize the ID
        arxiv_id = re.sub(r"(?i)^arxiv:", "", arxiv_id).strip()
        arxiv_id = re.sub(r"v\d+$", "", arxiv_id)  # strip version
        arxiv_link = f"https://arxiv.org/abs/{arxiv_id}"
        links.append(f"[{arxiv_link} arXiv:{arxiv_id}]")

    if links:
        ref += " " + " ".join(links)

    if show_notes and "note" in e:
        ref += f"\n: {e['note']}"

    return ref


# -------------------------------
# Processing
# -------------------------------
def process_bibtex(filename: str, merge_arxiv_with_published: bool = True,
                   show_notes: bool = False, reverse: bool = True,
                   year_headings: bool = True, extra_entries: Optional[List[Dict]] = None) -> List[str]:
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

    # after reading and categorizing db.entries:
    if extra_entries:
        for e in extra_entries:
            e["_from_arxiv_import"] = True
            arxiv_entries.append(e)

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

            for pub in published_entries:
                if (
                    pub.get("doi")
                    and normalize_title(pub.get("title", "")) == arxiv_title
                ):
                    # --- Merge arXiv info into the published record ---
                    if arxiv_id and not pub.get("_arxiv_id"):
                        pub["_arxiv_id"] = arxiv_id
                        # Preserve cross-links
                        if "note" in pub:
                            pub["note"] += f" Includes arXiv preprint {arxiv_id}."
                        else:
                            pub["note"] = f"Includes arXiv preprint {arxiv_id}."
                        # Always keep arXiv URL reference
                        if "url" not in pub and arx.get("url"):
                            pub["url"] = arx["url"]
                        if "eprint" not in pub and arx.get("eprint"):
                            pub["eprint"] = arx["eprint"]

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
    return output_lines, all_entries




# add arxiv import functionality

import requests
import bibtexparser
from xml.etree import ElementTree as ET

def fetch_arxiv_metadata(arxiv_ids: list[str]) -> dict:
    """Fetch metadata for one or more arXiv records and enrich with Crossref if DOI present."""
    import requests
    import xml.etree.ElementTree as ET

    if not arxiv_ids:
        return {}

    # Build proper comma-separated list for the arXiv API
    id_list_str = ",".join(arxiv_ids)
    url = f"https://export.arxiv.org/api/query?id_list={id_list_str}"

    r = requests.get(url, timeout=20)
    r.raise_for_status()

    root = ET.fromstring(r.text)
    ns = {"atom": "http://www.w3.org/2005/Atom", "arxiv": "http://arxiv.org/schemas/atom"}

    results = {}
    for entry in root.findall("atom:entry", ns):
        arxiv_id = entry.find("atom:id", ns).text.split("/")[-1]

        def get(tag):
            el = entry.find(f"atom:{tag}", ns)
            if el is not None and el.text:
                return el.text.strip()
            el = entry.find(f"arxiv:{tag}", ns)
            return el.text.strip() if el is not None and el.text else ""

        authors = [
            a.find("atom:name", ns).text
            for a in entry.findall("atom:author", ns)
            if a.find("atom:name", ns) is not None
        ]

        doi = get("doi")
        published = get("published")
        year = published[:4] if published else ""

        bib = {
            "ID": f"arxiv:{arxiv_id}",
            "ENTRYTYPE": "article",
            "title": get("title").replace("\n", " ").strip(),
            "author": " and ".join(authors),
            "year": year,
            "journal": "arXiv.org",
            "eprint": arxiv_id,
            "url": f"https://arxiv.org/abs/{arxiv_id}",
        }

        # --- Crossref enrichment if DOI present ---
        if doi:
            bib["doi"] = doi
            # Enrich from Crossref
            try:
                cr = requests.get(f"https://api.crossref.org/works/{doi}", timeout=10)
                cr.raise_for_status()
                data = cr.json().get("message", {})

                # Journal title
                bib["journal"] = data.get("container-title", ["arXiv.org"])[0]

                # Volume
                volume = data.get("volume")
                if not volume:
                    logger.warning(f"Attention – no 'volume' info for DOI {doi}")
                else:
                    bib["volume"] = volume

                # Pages / Article number (varies by publisher)
                pages = (
                    data.get("page")
                    or data.get("article-number")
                    or data.get("article_number")
                    or None
                )
                if not pages:
                    logger.warning(f"Attention – no 'page' or 'article-number' info for DOI {doi}")
                else:
                    bib["pages"] = pages

                # Year
                issued = data.get("issued", {}).get("date-parts", [[None]])
                if issued and issued[0][0]:
                    bib["year"] = str(issued[0][0])

            except Exception as e:
                logger.warning(f"Crossref lookup failed for DOI {doi}: {e}")
        results[arxiv_id] = bib

    return results


def fetch_crossref_metadata(doi):
    """Fetch metadata from Crossref."""
    try:
        r = requests.get(f"https://api.crossref.org/works/{doi}", timeout=10)
        if r.status_code != 200:
            return None
        data = r.json().get("message", {})
        meta = {
            "title": data.get("title", [""])[0],
            "authors": [
                f"{a.get('given','')} {a.get('family','')}".strip()
                for a in data.get("author", [])
            ],
            "journal": data.get("container-title", [""])[0],
            "volume": data.get("volume"),
            "year": data.get("issued", {}).get("date-parts", [[None]])[0][0],
            "doi": doi,
            "url": data.get("URL"),
        }
        return meta
    except Exception as e:
        print(f"Attention: Crossref fetch failed for {doi}: {e}")
        return None

def merge_metadata(arxiv_meta, publisher_meta):
    """Combine arXiv and publisher records, preferring publisher values."""
    merged = arxiv_meta.copy()
    merged.update({k: v for k, v in publisher_meta.items() if v})
    return merged

def to_bibtex(entry, arxiv_id):
    """Format metadata dict into BibTeX."""
    authors = " and ".join(entry["authors"])
    key = entry.get("doi", f"arXiv:{arxiv_id}").replace("/", "_")
    bibtex = f"""@article{{{key},
  title = {{{entry['title']}}},
  author = {{{authors}}},
  journal = {{{entry.get('journal', 'arXiv.org')}}},
  year = {{{entry.get('year', entry.get('published', '')[:4])}}},
  doi = {{{entry.get('doi', '')}}},
  url = {{{entry.get('url', entry.get('arxiv_url', ''))}}}
}}\n"""
    return bibtex

def import_arxiv_records(arxiv_in: str, arxiv_2_bib_out: Optional[str] = None) -> List[Dict]:
    """Import arXiv records from a file, fetch metadata in batch, and optionally save to BibTeX."""
    if not arxiv_in or not os.path.exists(arxiv_in):
        logger.warning("No valid arXiv input file specified: %s", arxiv_in)
        return []

    with open(arxiv_in, "r", encoding="utf-8") as f:
        arxiv_ids = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    if not arxiv_ids:
        logger.warning("No arXiv identifiers found in %s", arxiv_in)
        return []

    logger.info("Fetching metadata for %d arXiv records...", len(arxiv_ids))
    arxiv_data = fetch_arxiv_metadata(arxiv_ids)

    # Convert results to a flat list of BibTeX-style dicts
    records = list(arxiv_data.values())
    logger.info("Fetched %d valid arXiv entries.", len(records))

    # Optionally write fetched records to a separate BibTeX file
    if arxiv_2_bib_out:
        with open(arxiv_2_bib_out, "w", encoding="utf-8") as f:
            for e in records:
                f.write(f"@article{{e['ID']}},\n")
                for k, v in e.items():
                    if k in ["ID", "ENTRYTYPE"]:
                        continue
                    f.write(f"  {k} = {{{v}}},\n")
                f.write("}\n\n")
        logger.info("Wrote arXiv BibTeX data to %s", arxiv_2_bib_out)

    return records

# bulk arxiv import
def import_arxiv_records(arxiv_in: str, arxiv_2_bib_out: Optional[str] = None) -> List[Dict]:
    """Import arXiv records from a file, fetch metadata in batch, and optionally save to BibTeX."""
    if not arxiv_in or not os.path.exists(arxiv_in):
        logger.warning("No valid arXiv input file specified: %s", arxiv_in)
        return []

    with open(arxiv_in, "r", encoding="utf-8") as f:
        arxiv_ids = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    if not arxiv_ids:
        logger.warning("No arXiv identifiers found in %s", arxiv_in)
        return []

    logger.info("Fetching metadata for %d arXiv records...", len(arxiv_ids))
    arxiv_data = fetch_arxiv_metadata(arxiv_ids)

    # Convert results to a flat list of BibTeX-style dicts
    records = list(arxiv_data.values())
    logger.info("Fetched %d valid arXiv entries.", len(records))

    # Optionally write fetched records to a separate BibTeX file
    if arxiv_2_bib_out:
        with open(arxiv_2_bib_out, "w", encoding="utf-8") as f:
            for e in records:
                f.write(f"@article{{{e['ID']}}},\n")
                for k, v in e.items():
                    if k in ["ID", "ENTRYTYPE"]:
                        continue
                    f.write(f"  {k} = {{{v}}},\n")
                f.write("}\n\n")
        logger.info("Wrote arXiv BibTeX data to %s", arxiv_2_bib_out)

    return records


def import_arxiv_records_old(arxiv_ids: list[str], bib_filename: str, wiki_output: str, append: bool = False):
    """Import arXiv records, merge into BibTeX, and regenerate Wiki output."""
    print(f"Fetching {len(arxiv_ids)} arXiv records...")
    arxiv_data = fetch_arxiv_metadata(arxiv_ids)
    out_bib = Path("arxiv_imports.bib")
    mode = "a" if append else "w"
    with open(out_bib, mode, encoding="utf-8") as f:
        for aid, a_meta in arxiv_data.items():
            doi = a_meta.get("doi")
            if doi:
                p_meta = fetch_crossref_metadata(doi)
                if p_meta:
                    merged = merge_metadata(a_meta, p_meta)
                    print(f"OK {aid}: merged with {p_meta['journal']} ({p_meta['year']})")
                else:
                    print(f"Attention️ {aid}: no publisher data, keeping arXiv metadata.")
                    merged = a_meta
            else:
                print(f"Info️ {aid}: no DOI, arXiv only.")
                merged = a_meta
            bib_entry = to_bibtex(merged, aid)
            f.write(bib_entry)
            time.sleep(0.5)  # be kind to APIs

    print(f"📚 Saved fetched records to {out_bib}")

    # Merge with existing bib file for Wiki export
    print("🔧 Merging with existing BibTeX and regenerating wiki output...")
    combined_bib = Path("combined_import.bib")
    with open(combined_bib, "w", encoding="utf-8") as fout:
        fout.write(open(bib_filename, "r", encoding="utf-8").read())
        fout.write("\n")
        fout.write(open(out_bib, "r", encoding="utf-8").read())

    output_lines = process_bibtex(
        str(combined_bib),
        merge_arxiv_with_published=True,
        show_notes=False,
        reverse=True,
        year_headings=True,
    )
    with open(wiki_output, "w", encoding="utf-8") as f:
        f.write("\n".join(output_lines))

    print(f"✅ Wiki bibliography updated → {wiki_output}")

import unicodedata

def bibtex_entry(entry: dict, latex_output: bool = False) -> str:
    """
    Convert a Python dict representing a BibTeX entry into a formatted record.
    - latex_output: if True, converts UTF-8 diacritics to LaTeX escape sequences.
    """
    e = dict(entry)
    key = e.pop("ID", "unnamed")
    entry_type = e.pop("ENTRYTYPE", "article")

    # --- Helper: UTF-8 → LaTeX accent / spaces conversion ---
    # --- Character mappings for LaTeX escapes ---
    accent_map = {
        "ä": r"{\"a}", "ö": r"{\"o}", "ü": r"{\"u}",
        "Ä": r"{\"A}", "Ö": r"{\"O}", "Ü": r"{\"U}",
        "é": r"\'{e}", "è": r"\`{e}", "ê": r"\^{e}", "ë": r"\"{e}",
        "É": r"\'{E}", "È": r"\`{E}", "Ê": r"\^{E}", "Ë": r"\"{E}",
        "á": r"\'{a}", "à": r"\`{a}", "â": r"\^{a}", "å": r"\r{a}",
        "Á": r"\'{A}", "À": r"\`{A}", "Â": r"\^{A}", "Å": r"\r{A}",
        "ó": r"\'{o}", "ò": r"\`{o}", "ô": r"\^{o}", "õ": r"\~{o}",
        "Ó": r"\'{O}", "Ò": r"\`{O}", "Ô": r"\^{O}", "Õ": r"\~{O}",
        "í": r"\'{i}", "ì": r"\`{i}", "î": r"\^{i}", "ï": r"\"{i}",
        "Í": r"\'{I}", "Ì": r"\`{I}", "Î": r"\^{I}", "Ï": r"\"{I}",
        "ñ": r"\~{n}", "Ñ": r"\~{N}",
        "ç": r"\c{c}", "Ç": r"\c{C}",
        "ß": r"\ss{}",
        "ø": r"\o{}", "Ø": r"\O{}", "å": r"\r{a}", "Å": r"\r{A}"
    }

    # Smart quotes, dashes, ellipsis, spaces
    smart_punct_map = {
        "“": "``", "”": "''",  # smart double quotes
        "‘": "`",  "’": "'",   # smart single quotes
        "–": "--", "—": "---", # en/em dash
        "…": r"\ldots{}",      # ellipsis
        "\u00A0": "~",         # non-breaking space
    }

    def unicode_to_latex(text: str) -> str:
        """Convert accented & special Unicode characters to LaTeX equivalents."""
        out = []
        for c in text:
            if c in accent_map:
                out.append(accent_map[c])
            elif c in smart_punct_map:
                out.append(smart_punct_map[c])
            else:
                out.append(c)
        return "".join(out)

    def bib_escape(value: str) -> str:
        """Escape special characters for BibTeX."""
        if not isinstance(value, str):
            return str(value)
        value = value.replace("\\", "\\\\")
        value = value.replace("\"", "{\\textquotedbl}")
        value = value.replace("%", "\\%")
        value = value.replace("&", "\\&")
        value = value.replace("_", "\\_")
        value = value.replace("#", "\\#")
        value = value.replace("$", "\\$")
        if latex_output:
            value = unicode_to_latex(value)
        return value.strip()

    # --- Order fields for readability ---
    order = ["title", "author", "year", "journal", "volume", "pages", "doi", "url", "eprint", "note"]
    fields = []
    for k, v in e.items():
        if v:
            fields.append((k, v))
    fields.sort(key=lambda kv: (order.index(kv[0]) if kv[0] in order else 99, kv[0]))

    # --- Construct entry string ---
    field_lines = [f"  {k} = {{{bib_escape(v)}}}" for k, v in fields]
    bib_str = f"@{entry_type}{{{key},\n" + ",\n".join(field_lines) + "\n}"
    return bib_str

    
# -------------------------------
# CLI
# -------------------------------
def main():
    parser = argparse.ArgumentParser(description="Convert BibTeX to MediaWiki references.")
    parser.add_argument("bibfile", help="Input .bib file")
    parser.add_argument("--no-merge-arxiv", action="store_true", help="Do not merge arXiv with published")
    parser.add_argument("--show-notes", action="store_true", help="Show notes if present")
    parser.add_argument("--reverse", action="store_false", help="Sort bibliography in reverse chronological order (newest first)")
    parser.add_argument("--numbering", choices=["none", "global", "reverse"], default="none",
                            help="Add global numbering of all entries (global = ascending, reverse = descending).")
    parser.add_argument("--year-headings", action="store_false", help="Insert year headings in the output")

    parser.add_argument("--output", "-o",
                        help="Output file (default: same name as input, with .wiki extension; use '-' for stdout)",
                        default=None)

    parser.add_argument("--arxiv-in", help="Optional text file containing arXiv IDs to import and merge", default=None)
    parser.add_argument("--arxiv-2-bib-out", help="Optional output file for BibTeX generated from arXiv imports", default=None)
    parser.add_argument("--all-bib-out", help="Optional output file for full merged BibTeX (existing + arXiv imports)", default=None)
    parser.add_argument("--latex-output", action="store_true", help="Escape non-ASCII characters for LaTeX-compatible BibTeX output")

    parser.add_argument("--verbose", type=int, default=1, choices=[0, 1, 2],
                        help="Verbosity level (0=warnings only, 1=info, 2=debug)")
    args = parser.parse_args()

    # --- Logging setup ---
    level = logging.WARNING
    if args.verbose == 1:
        level = logging.INFO
    elif args.verbose >= 2:
        level = logging.DEBUG
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")

    # ---- Step 1 (Optional): import arXiv records ----
    imported_arxiv_records = []
    if args.arxiv_in:
        imported_arxiv_records = import_arxiv_records(
            args.arxiv_in,
            arxiv_2_bib_out=args.arxiv_2_bib_out
        )
        logger.info(f"Imported {len(imported_arxiv_records)} new arXiv records from {args.arxiv_in}")


    # --- Process main BibTeX, merging optional imported arXiv entries ---
    refs, all_entries = process_bibtex(
        args.bibfile,
        merge_arxiv_with_published=not args.no_merge_arxiv,
        show_notes=args.show_notes,
        reverse=args.reverse,
        year_headings=args.year_headings,
        extra_entries=imported_arxiv_records,  # <—— integrate imported records
    )
    # fix some output formatting - changing to global enumeration if desired.
    if args.numbering != "none":
        numbered_refs = []
        # Filter only actual entries (lines starting with "#")
        numbered_entries = [r for r in refs if r.strip().startswith("#")]
        n_total = len(numbered_entries)
        count = 0

        for ref in refs:
            ref_stripped = ref.strip()
            if not ref_stripped:
                # skip empty lines entirely
                continue
            if ref_stripped.startswith("#"):
                count += 1
                if args.numbering == "global":
                    num = count
                else:  # reverse
                    num = n_total - count + 1
                clean_ref = ref_stripped.lstrip("#").strip()
                numbered_refs.append(f"  {num}. {clean_ref}")
            else:
                # Year headings, keep unnumbered and unindented
                numbered_refs.append(ref_stripped)

        output_text = "\n".join(numbered_refs)
    else:
        # fallback: join non-empty refs without extra blank lines
        output_text = "\n".join([r.strip() for r in refs if r.strip()])


    # --- Step 3 (Optional): write merged BibTeX output ---
    if args.all_bib_out:
        with open(args.all_bib_out, "w", encoding="utf-8") as f:
            for entry in all_entries:
                f.write(bibtex_entry(entry, latex_output=args.latex_output) + "\n\n")
                f.write("\n\n")
        logger.info(f"Wrote merged BibTeX output to {args.all_bib_out}")
        
    # --- Step 4: Write Wiki output ---
    if args.output == "-" or (args.output is None and args.bibfile == "-"):
        print(output_text)
    else:
        outfile = args.output or re.sub(r"\.bib$", ".wiki", args.bibfile)
        with open(outfile, "w", encoding="utf-8") as f:
            f.write(output_text + "\n")
        logger.info("Wrote wiki output → %s", outfile)

if __name__ == "__main__":
    main()
