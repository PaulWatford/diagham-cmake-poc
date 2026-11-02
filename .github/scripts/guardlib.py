"""guardlib.py: shared pieces of the Guard workflow (.github/workflows/guard.yml).

Everything that talks to GitHub goes through the `gh` CLI with the workflow's GITHUB_TOKEN. Nothing here
merges, approves, edits a file in the repository or changes a setting: the scripts add labels, request
reviews and keep one comment per gate up to date. The maintainer lists are files on the base branch:

  .github/approved-developers.txt   logins whose pull requests are not held by the author gate
  .github/guard-reviewers.txt       logins who receive the review request and may release a hold

A hold is released when a login from guard-reviewers.txt adds the label `approved-by-maintainer` to the
pull request (the label event re-runs the workflow), or has submitted an approving review of the current
head commit. A label added by anyone else is removed again and reported.
"""
import json
import os
import re
import subprocess
import sys
from pathlib import Path

APPROVED_FILE = ".github/approved-developers.txt"
REVIEWERS_FILE = ".github/guard-reviewers.txt"
VOLUNTEERS_FILE = ".github/volunteers.txt"
LABEL_HOLD = "held-for-approval"
LABEL_EXTERNAL = "external-contribution"
LABEL_AI = "ai-assisted"
LABEL_DAMAGE = "damage"
LABEL_RELEASE = "approved-by-maintainer"
LABELS = {
    LABEL_HOLD: ("B60205", "HELD_FOR_APPROVAL: waits for a named maintainer (GOVERNANCE.md, AGENTS.md)"),
    LABEL_EXTERNAL: ("0E8A16", "author is not on .github/approved-developers.txt"),
    LABEL_AI: ("5319E7", "AI assistance declared or detected; reviewed under AGENTS.md"),
    LABEL_DAMAGE: ("D93F0B", "the damage scan flagged a locked artefact (scripts_cmake/damage_scan.py)"),
    LABEL_RELEASE: ("0052CC", "a maintainer from .github/guard-reviewers.txt releases the hold"),
}


def read_logins(path):
    p = Path(path)
    if not p.exists():
        return []
    out = []
    for line in p.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip().lstrip("@")
        if line:
            out.append(line)
    return out


def norm(login):
    return (login or "").lower()


class GitHub:
    """Thin wrapper over `gh api`. With dry_run=True nothing is sent; calls are recorded in .actions."""

    def __init__(self, repo, pr_number, dry_run=False):
        self.repo = repo
        self.pr = int(pr_number)
        self.dry_run = dry_run
        self.actions = []

    def api(self, path, *args, method=None, paginate=False, fields=None):
        cmd = ["gh", "api", path]
        if method:
            cmd += ["-X", method]
        if paginate:
            cmd += ["--paginate", "--slurp"]
        for k, v in (fields or {}).items():
            cmd += ["-f", f"{k}={v}"]
        cmd += list(args)
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            raise RuntimeError(f"{' '.join(cmd[:3])}: {p.stderr.strip()[:300]}")
        if not p.stdout.strip():
            return None
        data = json.loads(p.stdout)
        if paginate:  # --slurp returns a list of pages
            flat = []
            for page in data:
                flat.extend(page if isinstance(page, list) else [page])
            return flat
        return data

    # reads
    def pull(self):
        return self.api(f"repos/{self.repo}/pulls/{self.pr}")

    def commits(self):
        return self.api(f"repos/{self.repo}/pulls/{self.pr}/commits", paginate=True)

    def files(self):
        return self.api(f"repos/{self.repo}/pulls/{self.pr}/files", paginate=True)

    def reviews(self):
        return self.api(f"repos/{self.repo}/pulls/{self.pr}/reviews", paginate=True)

    def label_events(self):
        return self.api(f"repos/{self.repo}/issues/{self.pr}/events", paginate=True)

    def comments(self):
        return self.api(f"repos/{self.repo}/issues/{self.pr}/comments", paginate=True)

    # writes (all idempotent)
    def ensure_labels(self):
        for name, (colour, desc) in LABELS.items():
            self._do(f"label {name}", ["gh", "label", "create", name, "--color", colour, "--description", desc,
                                       "--force", "--repo", self.repo])

    def add_labels(self, names):
        if names:
            self._do(f"add labels {names}", ["gh", "pr", "edit", str(self.pr), "--repo", self.repo,
                                            "--add-label", ",".join(names)])

    def remove_labels(self, names, current):
        names = [n for n in names if n in current]
        if names:
            self._do(f"remove labels {names}", ["gh", "pr", "edit", str(self.pr), "--repo", self.repo,
                                               "--remove-label", ",".join(names)])

    def request_reviews(self, logins):
        if logins:
            self._do(f"request review {logins}", ["gh", "pr", "edit", str(self.pr), "--repo", self.repo,
                                                 "--add-reviewer", ",".join(logins)], tolerate=True)

    def upsert_comment(self, marker, body):
        body = f"{marker}\n{body}"
        if self.dry_run:
            self.actions.append(f"upsert comment {marker}")
            return
        existing = [c for c in (self.comments() or []) if marker in (c.get("body") or "")]
        if existing:
            self.api(f"repos/{self.repo}/issues/comments/{existing[0]['id']}", method="PATCH", fields={"body": body})
        else:
            self.api(f"repos/{self.repo}/issues/{self.pr}/comments", method="POST", fields={"body": body})

    def _do(self, what, cmd, tolerate=False):
        self.actions.append(what)
        if self.dry_run:
            return
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            msg = f"guard: {what} failed: {p.stderr.strip()[:300]}"
            if tolerate:
                print(msg, file=sys.stderr)
            else:
                raise RuntimeError(msg)


def release_status(gh, pr, reviewers, label_names, label_events=None, reviews=None):
    """Return (released, note, bad_label_actor).

    released: a reviewer added `approved-by-maintainer`, or approved the current head commit.
    bad_label_actor: the release label is present but its last `labeled` event came from a non-reviewer.
    """
    rev = {norm(r) for r in reviewers}
    note, bad_actor = "", None
    if LABEL_RELEASE in label_names:
        events = label_events if label_events is not None else (gh.label_events() or [])
        actor = None
        for e in events:
            if e.get("event") == "labeled" and (e.get("label") or {}).get("name") == LABEL_RELEASE:
                actor = (e.get("actor") or {}).get("login")
        if actor and norm(actor) in rev:
            return True, f"released: @{actor} added `{LABEL_RELEASE}`", None
        bad_actor = actor or "unknown"
        note = f"`{LABEL_RELEASE}` was added by @{bad_actor}, who is not in {REVIEWERS_FILE}; the label is removed"
    head_sha = ((pr.get("head") or {}).get("sha")) or ""
    for r in (reviews if reviews is not None else (gh.reviews() or [])):
        if r.get("state") == "APPROVED" and norm((r.get("user") or {}).get("login")) in rev \
                and r.get("commit_id") == head_sha:
            return True, f"released: approving review of `{head_sha[:7]}` by @{r['user']['login']}", bad_actor
    return False, note, bad_actor


def step_summary(text):
    path = os.environ.get("GITHUB_STEP_SUMMARY")
    if path:
        with open(path, "a", encoding="utf-8") as f:
            f.write(text + "\n")


def set_output(name, value):
    path = os.environ.get("GITHUB_OUTPUT")
    if path:
        with open(path, "a", encoding="utf-8") as f:
            f.write(f"{name}={value}\n")


AI_TRAILER_RE = re.compile(r"(?im)^\s*(co-authored-by|generated-by|assisted-by|signed-off-by)\s*:.*"
                           r"(claude|anthropic|copilot|chatgpt|openai|gemini|cursor|codex|devin|aider|sweep|"
                           r"jules|amazon q|codewhisperer)")
AI_TEXT_RE = re.compile(r"(?i)(generated with \[?claude code|🤖 generated|co-authored-by: github copilot|"
                        r"this (commit|change|pull request) was (generated|written) by (an )?(ai|llm|language model))")
AI_EMAIL_RE = re.compile(r"(?i)(noreply@anthropic\.com|copilot@github\.com|@users\.noreply\.github\.com$)")
AI_EMAIL_LOGIN_RE = re.compile(r"(?i)(copilot|claude|anthropic|openai|chatgpt|devin|sweep-ai|dependabot|renovate)")
PHYSICS_RE = re.compile(r"^(src|FQHE|FTI|Spin|QuantumDots|Base)/.+")
BUILD_FILE_RE = re.compile(r"(^|/)(CMakeLists\.txt|.*\.cmake|Makefile\.am)$")


def ai_markers(pr, commits):
    """Return a list of human-readable markers of AI involvement found in the body and the commits."""
    found = []
    body = pr.get("body") or ""
    if AI_TEXT_RE.search(body):
        found.append("the pull request text says it was AI-generated")
    for c in commits or []:
        sha = (c.get("sha") or "")[:7]
        commit = c.get("commit") or {}
        msg = commit.get("message") or ""
        if AI_TRAILER_RE.search(msg) or AI_TEXT_RE.search(msg):
            found.append(f"commit `{sha}` carries an AI trailer or marker")
        for who in ("author", "committer"):
            email = ((commit.get(who) or {}).get("email")) or ""
            login = ((c.get(who) or {}) or {}).get("login") or ""
            if AI_EMAIL_RE.search(email) and AI_EMAIL_LOGIN_RE.search(email):
                found.append(f"commit `{sha}` {who} e-mail `{email}`")
            if login.endswith("[bot]") or ((c.get(who) or {}).get("type") == "Bot"):
                found.append(f"commit `{sha}` {who} is the bot `{login}`")
    # de-duplicate, keep order
    seen, out = set(), []
    for f in found:
        if f not in seen:
            seen.add(f)
            out.append(f)
    return out


def disclosure(body):
    """Read the 'AI assistance' section of the pull request body.

    Returns 'none', 'assisted', 'generated', or None when no box is ticked (or the section is missing).
    """
    body = body or ""
    m = re.search(r"(?is)##\s*AI assistance(.*?)(\n##\s|\Z)", body)
    section = m.group(1) if m else ""
    ticked = re.findall(r"(?im)^\s*[-*]\s*\[[xX]\]\s*(.+)$", section)
    kinds = []
    for t in ticked:
        t = t.lower()
        if t.startswith("no ai") or "no ai tools" in t:
            kinds.append("none")
        elif "ai-generated" in t or "ai generated" in t:
            kinds.append("generated")
        elif "ai-assisted" in t or "ai assisted" in t:
            kinds.append("assisted")
    if not kinds:
        return None
    for k in ("generated", "assisted", "none"):
        if k in kinds:
            return k
    return None


def physics_files(files):
    """Paths under the DiagHam source trees that are not build files (AGENTS.md: out of bounds for automated changes)."""
    out = []
    for f in files or []:
        p = f.get("filename") or ""
        if PHYSICS_RE.match(p) and not BUILD_FILE_RE.search(p):
            out.append(p)
    return out
