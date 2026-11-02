#!/usr/bin/env python3
"""discord_post.py: announce a pull request to the volunteers' Discord channel (Guard workflow, job notify).

    python3 .github/scripts/discord_post.py                 # in the workflow: reads the environment, posts
    python3 .github/scripts/discord_post.py --dry-run FIXTURE.json   # prints the JSON it would send

Environment: DISCORD_WEBHOOK_URL, PR_NUMBER, PR_TITLE, PR_AUTHOR, PR_URL, GUARD_HOLD, GUARD_AI, GUARD_PHYSICS,
GUARD_DAMAGE, GUARD_EVENT_ACTION, GUARD_LABEL, GUARD_LABEL_ACTOR.

One message per new pull request (and one when a maintainer releases a hold with the label
approved-by-maintainer). Only the pull request's number, title, author and link are taken from it; the body is
untrusted text and is never posted. Mentions in the title are neutralised and the message allows no mentions.
A failure to post is printed and the script exits 0: a Discord outage must not turn the Guard red. Python
standard library only.
"""
import argparse
import json
import os
import sys
import urllib.error
import urllib.request

RELEASE_LABEL = "approved-by-maintainer"
TITLE_MAX = 200
CHECKLIST = ("Volunteer checklist:\n"
             "1. Does the description say what changed and why?\n"
             "2. Are the G1 to G7 lines present with their output?\n"
             "3. Does the diff touch anything in CLAUDE.md's stop-and-ask list?\n"
             "Record a check with a pull request comment starting \"Volunteer check:\".")
ZWSP = "​"


def neutralise(text):
    """Break @everyone, @here, user and role mentions: a zero-width space after every @ and every <."""
    return (text or "").replace("@", "@" + ZWSP).replace("<", "<" + ZWSP)


def build(info):
    number = info.get("number") or "?"
    title = neutralise((info.get("title") or "")[:TITLE_MAX])
    release = info.get("event_action") == "labeled" and info.get("label") == RELEASE_LABEL
    embed = {
        "title": f"#{number}  {title}"[:256],
        "url": info.get("url") or "",
        "fields": [
            {"name": "Author", "value": neutralise(info.get("author") or "not read"), "inline": True},
            {"name": "Hold", "value": info.get("hold") or "not read", "inline": True},
            {"name": "AI assistance", "value": info.get("ai") or "not read", "inline": True},
            {"name": "DiagHam physics sources touched", "value": info.get("physics") or "not read", "inline": True},
            {"name": "Damage", "value": info.get("damage") or "not read", "inline": True},
        ],
    }
    if release:
        embed["description"] = f"Released by {neutralise(info.get('label_actor') or 'not read')}"
    else:
        embed["footer"] = {"text": CHECKLIST}
    return {"username": "DiagHam Guard", "embeds": [embed], "allowed_mentions": {"parse": []}}


def from_env():
    e = os.environ
    return {"number": e.get("PR_NUMBER"), "title": e.get("PR_TITLE"), "author": e.get("PR_AUTHOR"),
            "url": e.get("PR_URL"), "hold": e.get("GUARD_HOLD"), "ai": e.get("GUARD_AI"),
            "physics": e.get("GUARD_PHYSICS"), "damage": e.get("GUARD_DAMAGE"),
            "event_action": e.get("GUARD_EVENT_ACTION"), "label": e.get("GUARD_LABEL"),
            "label_actor": e.get("GUARD_LABEL_ACTOR")}


def from_fixture(path):
    fx = json.loads(open(path, encoding="utf-8").read())
    pr = fx.get("pr") or {}
    return {"number": pr.get("number"), "title": pr.get("title"), "author": (pr.get("user") or {}).get("login"),
            "url": pr.get("html_url"), "hold": fx.get("hold"), "ai": fx.get("ai"), "physics": fx.get("physics"),
            "damage": fx.get("damage"), "event_action": fx.get("event_action"), "label": fx.get("label"),
            "label_actor": fx.get("label_actor")}


def post(url, payload):
    req = urllib.request.Request(url, data=json.dumps(payload).encode("utf-8"), method="POST",
                                 headers={"Content-Type": "application/json",
                                          "User-Agent": "diagham-guard (github.com/PaulWatford/diagham-cmake-poc)"})
    try:
        with urllib.request.urlopen(req, timeout=10) as r:
            print(f"discord: posted (HTTP {r.status})")
    except urllib.error.HTTPError as err:
        print(f"discord: not posted (HTTP {err.code}); the Guard is unaffected")
    except (urllib.error.URLError, OSError) as err:
        print(f"discord: not posted ({err}); the Guard is unaffected")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", metavar="FIXTURE", help="print the message built from a JSON fixture; send nothing")
    a = ap.parse_args(argv)
    if a.dry_run:
        print(json.dumps(build(from_fixture(a.dry_run)), indent=2, ensure_ascii=False))
        return 0
    url = os.environ.get("DISCORD_WEBHOOK_URL")
    if not url:
        print("discord: DISCORD_WEBHOOK_URL is not set; nothing posted")
        return 0
    post(url, build(from_env()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
