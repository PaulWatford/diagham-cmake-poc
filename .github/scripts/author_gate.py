#!/usr/bin/env python3
"""author_gate.py: the author gate of the Guard workflow.

Decides, for one pull request, whether it is held for a named maintainer and who is asked to review it:

  * the author is not in .github/approved-developers.txt (or is a bot)  -> external-contribution, held
  * AI assistance is declared in the pull request body, or detected in the commits -> ai-assisted, held
  * an AI-assisted change touches DiagHam physics sources (src/, FQHE/, FTI/, Spin/, QuantumDots/, Base/,
    build files excepted) -> "for maintainer review" in the report (AGENTS.md: out of bounds for automation)
  * the "AI assistance" section of the pull request template has no box ticked -> the check fails until it has

A held pull request shows this check red (HELD_FOR_APPROVAL) until a login from .github/guard-reviewers.txt
adds the label `approved-by-maintainer` or approves the current head commit. Silence is not approval.
Pull requests by approved developers with no AI involvement pass without labels, comments or requests.

Environment (set by the workflow): GH_TOKEN, GITHUB_REPOSITORY, PR_NUMBER; optional GUARD_SENDER.
Offline test:  python3 author_gate.py --dry-run FIXTURE.json   (a JSON object with keys pr, commits, files,
reviews, label_events; nothing is sent to GitHub and the decided actions are printed).
"""
import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import guardlib as g  # noqa: E402

MARKER = "<!-- diagham-guard:author-gate -->"


def decide(pr, commits, files, approved, reviewers):
    author = (pr.get("user") or {}).get("login") or ""
    is_bot = author.endswith("[bot]") or (pr.get("user") or {}).get("type") == "Bot"
    is_approved = g.norm(author) in {g.norm(a) for a in approved} and not is_bot
    decl = g.disclosure(pr.get("body"))
    markers = g.ai_markers(pr, commits)
    ai = decl in ("assisted", "generated") or bool(markers) or is_bot
    physics = g.physics_files(files) if ai else []
    disclosure_missing = decl is None and not is_bot
    hold = (not is_approved) or ai
    labels = []
    if hold:
        labels.append(g.LABEL_HOLD)
    if not is_approved:
        labels.append(g.LABEL_EXTERNAL)
    if ai:
        labels.append(g.LABEL_AI)
    ask = [r for r in reviewers if g.norm(r) != g.norm(author)] if hold else []
    return {"author": author, "is_bot": is_bot, "is_approved": is_approved, "disclosure": decl,
            "markers": markers, "ai": ai, "physics": physics, "disclosure_missing": disclosure_missing,
            "hold": hold, "labels": labels, "ask": ask}


def report(d, released, release_note, reviewers):
    who = ("bot" if d["is_bot"] else "approved developer" if d["is_approved"]
           else f"not on `{g.APPROVED_FILE}`")
    if d["disclosure_missing"]:
        ai_line = "**not declared**: tick one box under \"AI assistance\" in the pull request text (the template has them)"
    else:
        ai_line = {"none": "none declared", "assisted": "declared: AI-assisted",
                   "generated": "declared: AI-generated", None: "bot author"}[d["disclosure"]]
    if d["markers"]:
        ai_line += "; detected: " + "; ".join(d["markers"][:6])
    if d["physics"]:
        phys = (f"yes, {len(d['physics'])} file(s), **for maintainer review** (AGENTS.md: DiagHam physics code is "
                "out of bounds for automated changes): " + ", ".join(f"`{p}`" for p in d["physics"][:10]))
    else:
        phys = "no" if d["ai"] else "not checked (no AI involvement)"
    if not d["hold"]:
        hold = "none"
    elif released:
        hold = release_note
    else:
        names = ", ".join(f"@{r}" for r in d["ask"]) or "(nobody listed in " + g.REVIEWERS_FILE + ")"
        hold = (f"**HELD_FOR_APPROVAL** until {names} adds the label `{g.LABEL_RELEASE}` or approves this head "
                "commit. Silence is not approval (GOVERNANCE.md; AGENTS.md).")
    lines = ["## Guard: author gate", "", "| Field | Value |", "|---|---|",
             f"| Author | @{d['author']} ({who}) |",
             f"| AI assistance | {ai_line} |",
             f"| DiagHam physics sources touched | {phys} |",
             f"| Hold | {hold} |"]
    if release_note and not released:
        lines.append(f"| Note | {release_note} |")
    lines += ["", "The build and the tests decide pass or fail; this gate only decides who must look. "
              "It never merges, never edits a file and never approves."]
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", metavar="FIXTURE", help="decide from a JSON fixture; send nothing")
    ap.add_argument("--repo-root", default=".", help="where the .github lists are read from (the base checkout)")
    a = ap.parse_args(argv)
    root = Path(a.repo_root)
    approved = g.read_logins(root / g.APPROVED_FILE)
    reviewers = g.read_logins(root / g.REVIEWERS_FILE)

    if a.dry_run:
        fx = json.loads(Path(a.dry_run).read_text(encoding="utf-8"))
        pr, commits, files = fx["pr"], fx.get("commits", []), fx.get("files", [])
        gh = g.GitHub(fx.get("repo", "owner/repo"), pr["number"], dry_run=True)
        label_events, reviews = fx.get("label_events", []), fx.get("reviews", [])
    else:
        repo, number = os.environ["GITHUB_REPOSITORY"], os.environ["PR_NUMBER"]
        gh = g.GitHub(repo, number)
        pr, commits, files = gh.pull(), gh.commits(), gh.files()
        label_events, reviews = None, None

    d = decide(pr, commits, files, approved, reviewers)
    current = {l["name"] for l in pr.get("labels") or []}
    released, note, bad_actor = (False, "", None)
    if d["hold"]:
        released, note, bad_actor = g.release_status(gh, pr, reviewers, current, label_events, reviews)

    gh.ensure_labels()
    gh.add_labels([l for l in d["labels"] if l not in current])
    stale = [l for l in (g.LABEL_HOLD, g.LABEL_EXTERNAL, g.LABEL_AI) if l not in d["labels"]]
    if released:
        stale.append(g.LABEL_HOLD)
    if bad_actor:
        stale.append(g.LABEL_RELEASE)
    gh.remove_labels(stale, current)
    if d["hold"] and not released:
        gh.request_reviews(d["ask"])
    body = report(d, released, note, reviewers)
    # one comment per gate, kept up to date; written when something is held or was held before
    if d["hold"] or d["disclosure_missing"] or current & {g.LABEL_HOLD, g.LABEL_EXTERNAL, g.LABEL_AI}:
        gh.upsert_comment(MARKER, body)

    g.step_summary(body)
    g.set_output("hold", "true" if (d["hold"] and not released) else "false")
    g.set_output("ai", "true" if d["ai"] else "false")
    print(body)
    if a.dry_run:
        print("\nactions:", *gh.actions, sep="\n  ")
    if d["disclosure_missing"]:
        print("\nauthor gate: the AI assistance section has no box ticked", file=sys.stderr)
        return 1
    if d["hold"] and not released:
        print("\nauthor gate: HELD_FOR_APPROVAL", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
