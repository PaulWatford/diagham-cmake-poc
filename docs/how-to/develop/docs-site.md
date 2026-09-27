# Preview and publish the documentation site

Purpose: build the `docs/` knowledge base as a website locally, and understand how it is published.
Source: new (2026-09-28); `mkdocs.yml`, `scripts_cmake/gen_mkdocs_nav.py`, `.github/workflows/pages.yml`.

The knowledge base lives in `docs/` as Markdown, versioned and reviewed with
the code (which is why the GitHub wiki is off). Material for MkDocs turns
it into a searchable site at <https://paulwatford.github.io/diagham-cmake-poc/>
on every push to `main`; pull requests that touch `docs/` build the site
with `--strict`, so a broken link fails the check.

## Locally

```
python3 -m venv ~/mkdocs-venv && ~/mkdocs-venv/bin/pip install "mkdocs-material==9.*"
~/mkdocs-venv/bin/mkdocs serve          # http://127.0.0.1:8000, rebuilds on save
~/mkdocs-venv/bin/mkdocs build --strict # what CI runs; output in build/site/
```

## What the site contains

- Every page under `docs/`: the navigation lists the hand-written pages
  (tutorials, how-to guides, reference, explanation, history, drafts); the
  700+ generated program pages and the parked wiki stubs are not in the
  navigation but are reachable from their index pages and through search.
- Maths written as `$…$` / `$$…$$` is rendered by MathJax.
- "Edit this page" links to the file on GitHub.
- The root files (`README.md`, `CONTRIBUTING.md`, `GOVERNANCE.md`, …) are
  not part of the site; the home page links to them on GitHub.

## When you add or rename a page

1. Give it the `Purpose:` and `Source:` lines like every other page.
2. Run `python3 scripts_cmake/gen_mkdocs_nav.py` — it rewrites the `nav:`
   block of `mkdocs.yml` from the tree, using each page's H1 as its title.
   CI checks that the block is up to date (`--check`).
3. `mkdocs build --strict` must pass (links, anchors, nav).
