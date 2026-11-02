#!/usr/bin/env python3
"""guard_report.py: publish the damage scan of one pull request and apply the hold.

    python3 .github/scripts/guard_report.py damage.md damage.json

Reads the report that scripts_cmake/damage_scan.py wrote, keeps one comment on the pull request up to date,
sets or clears the `damage` label, and exits 1 (HELD_FOR_APPROVAL) while damage is flagged and no login from
.github/guard-reviewers.txt has released the hold (label `approved-by-maintainer`, or an approving review of
the current head). It never edits the change and never rewrites a golden.

Environment: GH_TOKEN, GITHUB_REPOSITORY, PR_NUMBER.  Offline: --dry-run FIXTURE.json (as author_gate.py).
"""
import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import guardlib as g  # noqa: E402

MARKER = "<!-- diagham-guard:damage-scan -->"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("markdown")
    ap.add_argument("json")
    ap.add_argument("--dry-run", metavar="FIXTURE")
    ap.add_argument("--repo-root", default=".")
    a = ap.parse_args(argv)
    md = Path(a.markdown).read_text(encoding="utf-8")
    result = json.loads(Path(a.json).read_text(encoding="utf-8"))
    reviewers = g.read_logins(Path(a.repo_root) / g.REVIEWERS_FILE)

    if a.dry_run:
        fx = json.loads(Path(a.dry_run).read_text(encoding="utf-8"))
        pr = fx["pr"]
        gh = g.GitHub(fx.get("repo", "owner/repo"), pr["number"], dry_run=True)
        label_events, reviews = fx.get("label_events", []), fx.get("reviews", [])
    else:
        gh = g.GitHub(os.environ["GITHUB_REPOSITORY"], os.environ["PR_NUMBER"])
        pr = gh.pull()
        label_events, reviews = None, None

    current = {l["name"] for l in pr.get("labels") or []}
    damage = bool(result.get("damage"))
    released, note, bad_actor = (False, "", None)
    if damage:
        released, note, bad_actor = g.release_status(gh, pr, reviewers, current, label_events, reviews)

    gh.ensure_labels()
    if damage:
        gh.add_labels([l for l in (g.LABEL_DAMAGE, g.LABEL_HOLD) if l not in current and not released])
        if released:
            gh.remove_labels([g.LABEL_HOLD], current)
    else:
        gh.remove_labels([g.LABEL_DAMAGE], current)
    if bad_actor:
        gh.remove_labels([g.LABEL_RELEASE], current)

    if damage and not released:
        who = ", ".join(f"@{r}" for r in reviewers) or f"(nobody listed in {g.REVIEWERS_FILE})"
        hold = (f"\n**Hold:** HELD_FOR_APPROVAL until {who} adds the label `{g.LABEL_RELEASE}` or approves this head "
                "commit. Silence is not approval.\n")
    elif damage:
        hold = f"\n**Hold:** {note}\n"
    else:
        hold = "\n**Hold:** none from this scan.\n"
    body = md + hold + (f"\n**Note:** {note}\n" if note and not released else "")
    if damage or current & {g.LABEL_DAMAGE}:
        gh.upsert_comment(MARKER, body)
    g.step_summary(body)
    g.set_output("damage", "true" if damage else "false")
    g.set_output("hold", "true" if (damage and not released) else "false")
    if a.dry_run:
        print(body)
        print("actions:", *gh.actions, sep="\n  ")
    if damage and not released:
        print("damage scan: HELD_FOR_APPROVAL", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
