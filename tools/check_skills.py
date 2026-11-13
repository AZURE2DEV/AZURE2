#!/usr/bin/env python3
"""Check the Claude Code skills in .claude/skills/ against docs/upgrade/RULES.md §3.

Rules enforced (numbers as in RULES.md §3):

  R4  SKILL.md at most 500 lines (warning; error above the hard limit of 800 or with
      --strict); every references/*.md at most 300 lines (error).
  R5  no commit hashes, job IDs, session narratives or "fixed on" notes; dates only in
      version notes ("builds before vX.Y" / "builds before YYYY-MM-DD", the date
      on the same line as "expires: YYYY-MM-DD") no more than 90 days after --today; an expired note is an error.
  R6  one owner per task: each skill directory has a SKILL.md whose `name` equals the
      directory name; names and descriptions are unique among the repository skills and
      do not repeat a user-level skill (directory $CLAUDE_USER_SKILLS_DIR, default
      ~/.claude/skills when it exists); a skill named in backticks (`azure2-...`) exists.
  R7  description at most 600 characters.
  R8  no personal absolute paths (/groups/, /Users/, /users/, /home/, /afs/, ~/...; ~/.azure2 and ~/.claude are allowed).
  Links: every docs/... path and every relative Markdown link in a skill resolves, and so
      does every references/... or scripts/... file a skill names (in the skill or the repository).
  R10 no skill (Markdown or script) uses a removed pyazr name: the names of the column "earlier name" of the
      migration table in docs/source/user_guide/pyazr.rst, read from the page, as a call (`name(`) or an
      attribute (`obj.name`); a name that also appears in the column "pyazr 3.0" is not checked, and a line of a
      version note ("builds before ...") may name it.

A line can opt out of one rule with an HTML comment `<!-- check_skills: allow R5 -->`;
use it only for a quoted example of what the rule forbids.

    python3 tools/check_skills.py [--root REPO] [--today YYYY-MM-DD] [--strict]

Exits 0 when no rule is broken, 1 otherwise; prints one line per finding as
`path:line: [rule] message`.  Pure standard library.
"""

import argparse
import datetime as _dt
import os
import re
import sys

SKILL_MAX = 500
SKILL_HARD_MAX = 800
REFERENCE_MAX = 300
DESCRIPTION_MAX = 600
EXPIRY_HORIZON_DAYS = 90

HASH_RE = re.compile(r"(?<![\w./-])([0-9a-f]{7,40})(?![\w-])")
DATE_RE = re.compile(r"(?<!\d)(?:19|20)\d\d-\d\d-\d\d(?!\d)")
EXPIRES_RE = re.compile(r"expires:\s*((?:19|20)\d\d-\d\d-\d\d)")
VERSION_NOTE_RE = re.compile(r"\bbefore\s+(?:v\d+(?:\.\d+)+|(?:19|20)\d\d-\d\d-\d\d)", re.I)
JOB_RE = re.compile(r"\b(?:jobs?|cluster|condor|qsub|slurm|pid)\s*(?:id\s*)?#?\s*\d{4,}", re.I)
NARRATIVE_RES = [
    re.compile(p, re.I)
    for p in (
        r"\bfixed (?:on|in) (?:dev|main|qt5|\d)",
        r"\b(?:on|since) dev\b",
        r"\b(?:previous|earlier) session\b",
        r"\byesterday\b",
        r"\bafter a week of\b",
    )
]
PERSONAL_PATH_RE = re.compile(
    r"(?:(?<![\w.$})])(?:/groups/|/Users/|/users/|/home/|/afs/|/scratch\d*/)"
    r"|(?:^|(?<=[\s`'\"(=:]))~/(?!\.azure2\b|\.claude\b))"
)
DOCS_PATH_RE = re.compile(r"(?<![\w./-])(docs/[\w./-]+\.(?:rst|md|py|html|txt|png|svg))")
LOCAL_FILE_RE = re.compile(r"(?<![\w./-])((?:references|scripts)/[\w.-]+\.\w+)")
MD_LINK_RE = re.compile(r"\]\(([^)\s]+)\)")
SKILL_NAME_RE = re.compile(r"`(azure2-[a-z][a-z0-9-]*|rmfit-[a-z][a-z0-9-]*)`")
ALLOW_RE = re.compile(r"<!--\s*check_skills:\s*allow\s+([R0-9a-z, ]+)-->")
PYAZR_PAGE = os.path.join("docs", "source", "user_guide", "pyazr.rst")
TABLE_HEADER_RE = re.compile(r"^\s*\* - earlier name\s*$")
LITERAL_RE = re.compile(r"``([^`]+)``")
IDENT_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def removed_pyazr_names(page_text):
    """{name: compiled pattern} of the removed names in the migration table of the pyazr page, or None if the page
    has no such table.  A row is "   * - <earlier>" then "     - <3.0>", each cell possibly continued on more
    indented lines.  A literal's name is its last dotted component before a "(" (``pyazr.azure2(file)`` gives
    azure2, ``s.sess`` gives sess); a name used anywhere in the 3.0 column is not a removed name."""
    lines = page_text.splitlines()
    start = next((i for i, l in enumerate(lines) if TABLE_HEADER_RE.match(l)), None)
    if start is None:
        return None
    earlier, current, cell = [], [], None
    for line in lines[start + 2:]:
        if line.strip() and not line.startswith("   "):
            break
        if re.match(r"^   \* - ", line):
            cell = earlier
            line = line[6:]
        elif re.match(r"^     - ", line):
            cell = current
            line = line[6:]
        if cell is not None:
            cell.append(line)
    current_names = {n for lit in LITERAL_RE.findall(" ".join(current)) for n in IDENT_RE.findall(lit)}
    patterns = {}
    for lit in LITERAL_RE.findall(" ".join(earlier)):
        head = lit.split("(", 1)[0].strip()
        if not head:
            continue
        name = head.split(".")[-1]
        if not IDENT_RE.fullmatch(name) or name in current_names:
            continue
        call_or_attr = (r"(?<![\w.])" + name + r"\s*\(") + "|" + (r"(?<=\w)\." + name + r"\b")
        patterns[name] = re.compile(call_or_attr)
    return patterns


class Finding:
    def __init__(self, path, line, rule, message, warning=False):
        self.path, self.line, self.rule, self.message, self.warning = path, line, rule, message, warning

    def __str__(self):
        kind = "warning" if self.warning else "error"
        return f"{self.path}:{self.line}: [{self.rule}] {kind}: {self.message}"


def parse_frontmatter(text):
    """Return ({key: value}, number of lines of the front matter) for a SKILL.md."""
    lines = text.splitlines()
    if not lines or lines[0].strip() != "---":
        return None, 0
    fields = {}
    for i, line in enumerate(lines[1:], start=1):
        if line.strip() == "---":
            return fields, i + 1
        m = re.match(r"^([A-Za-z_-]+):\s*(.*)$", line)
        if m:
            key, value = m.group(1), m.group(2).strip()
            if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
                value = value[1:-1]
            fields[key] = value
        elif fields and line.startswith((" ", "\t")):
            last = list(fields)[-1]
            fields[last] = (fields[last] + " " + line.strip()).strip()
    return None, 0


def allowed(line, rule):
    m = ALLOW_RE.search(line)
    return bool(m) and rule in re.split(r"[ ,]+", m.group(1).strip())


def is_hash(token):
    return any(c.isdigit() for c in token) and any(c.isalpha() for c in token)


class Checker:
    def __init__(self, root, today, strict=False, user_skills_dir=None):
        self.root = os.path.abspath(root)
        self.skills_dir = os.path.join(self.root, ".claude", "skills")
        self.today = today
        self.strict = strict
        self.user_skills_dir = user_skills_dir
        self.findings = []
        self.removed_names = {}

    def add(self, path, line, rule, message, warning=False):
        self.findings.append(Finding(os.path.relpath(path, self.root), line, rule, message, warning))

    # -- one Markdown or script file ------------------------------------------------
    def check_text(self, path, skill_dir, skill_names, markdown=True):
        with open(path, encoding="utf-8", errors="replace") as f:
            text = f.read()
        horizon = self.today + _dt.timedelta(days=EXPIRY_HORIZON_DAYS)
        for n, line in enumerate(text.splitlines(), start=1):
            if not allowed(line, "R10") and not VERSION_NOTE_RE.search(line):
                for name, rx in self.removed_names.items():
                    if rx.search(line):
                        self.add(path, n, "R10", f"removed pyazr name {name!r} (see the migration table of "
                                                 f"{PYAZR_PAGE})")
            if not allowed(line, "R8"):
                m = PERSONAL_PATH_RE.search(line)
                if m:
                    self.add(path, n, "R8", f"personal absolute path {m.group(0)!r}; use an environment variable")
            if not markdown:
                continue
            if not allowed(line, "R5"):
                self.check_history(path, n, line, horizon)
            if not allowed(line, "links"):
                self.check_links(path, n, line, skill_dir)
            if not allowed(line, "R6"):
                for name in SKILL_NAME_RE.findall(line):
                    if name not in skill_names:
                        self.add(path, n, "R6", f"links skill `{name}`, which does not exist in .claude/skills/")

    def check_history(self, path, n, line, horizon):
        for m in HASH_RE.finditer(line):
            if is_hash(m.group(1)):
                self.add(path, n, "R5", f"looks like a commit hash: {m.group(1)}")
        m = JOB_RE.search(line)
        if m:
            self.add(path, n, "R5", f"job or process ID: {m.group(0)!r}")
        for rx in NARRATIVE_RES:
            m = rx.search(line)
            if m:
                self.add(path, n, "R5", f"history or session narrative: {m.group(0)!r}")
        dates = DATE_RE.findall(line)
        if not dates:
            return
        expiry = EXPIRES_RE.search(line)
        if not (VERSION_NOTE_RE.search(line) and expiry):
            self.add(path, n, "R5", "a date outside a version note ('builds before ...' with 'expires: YYYY-MM-DD')")
            return
        try:
            expires = _dt.date.fromisoformat(expiry.group(1))
        except ValueError:
            self.add(path, n, "R5", f"invalid expiry date {expiry.group(1)}")
            return
        if expires < self.today:
            self.add(path, n, "R5", f"version note expired on {expires}; remove or renew it")
        elif expires > horizon:
            self.add(path, n, "R5", f"expiry {expires} is more than {EXPIRY_HORIZON_DAYS} days ahead")

    def check_links(self, path, n, line, skill_dir):
        here = os.path.dirname(path)
        for target in DOCS_PATH_RE.findall(line):
            if not os.path.exists(os.path.join(self.root, target)):
                self.add(path, n, "links", f"{target} does not exist")
        for target in LOCAL_FILE_RE.findall(line):
            if not any(os.path.exists(os.path.join(base, target)) for base in (skill_dir, here, self.root)):
                self.add(path, n, "links", f"{target} does not exist in {os.path.relpath(skill_dir, self.root)}")
        for target in MD_LINK_RE.findall(line):
            if re.match(r"^[a-z]+:", target) or target.startswith("#"):
                continue
            target = target.split("#", 1)[0]
            if target and not os.path.exists(os.path.normpath(os.path.join(here, target))):
                self.add(path, n, "links", f"Markdown link {target} does not resolve")

    # -- the whole skill set ----------------------------------------------------------
    def run(self):
        if not os.path.isdir(self.skills_dir):
            self.add(self.skills_dir, 0, "R6", "no .claude/skills directory")
            return self.findings
        skill_dirs = sorted(d for d in os.listdir(self.skills_dir) if os.path.isdir(os.path.join(self.skills_dir, d)))
        self.load_removed_names()
        skill_names = set(skill_dirs)
        names, descriptions = {}, {}
        for d in skill_dirs:
            sdir = os.path.join(self.skills_dir, d)
            skill_md = os.path.join(sdir, "SKILL.md")
            if not os.path.isfile(skill_md):
                self.add(sdir, 0, "R6", "skill directory without SKILL.md")
                continue
            with open(skill_md, encoding="utf-8") as f:
                text = f.read()
            nlines = len(text.splitlines())
            if nlines > SKILL_HARD_MAX:
                self.add(skill_md, nlines, "R4", f"{nlines} lines; hard limit {SKILL_HARD_MAX}")
            elif nlines > SKILL_MAX:
                self.add(skill_md, nlines, "R4", f"{nlines} lines; limit {SKILL_MAX}", warning=not self.strict)
            fields, _ = parse_frontmatter(text)
            if fields is None:
                self.add(skill_md, 1, "R6", "missing YAML front matter (--- name/description ---)")
                fields = {}
            name, desc = fields.get("name", ""), fields.get("description", "")
            if name != d:
                self.add(skill_md, 2, "R6", f"name {name!r} differs from the directory name {d!r}")
            if not desc:
                self.add(skill_md, 3, "R7", "no description")
            elif len(desc) > DESCRIPTION_MAX:
                self.add(skill_md, 3, "R7", f"description has {len(desc)} characters; limit {DESCRIPTION_MAX}")
            if name in names:
                self.add(skill_md, 2, "R6", f"name {name!r} also used by {names[name]}")
            names.setdefault(name, d)
            key = " ".join(desc.lower().split())
            if key and key in descriptions:
                self.add(skill_md, 3, "R6", f"same description as {descriptions[key]}")
            descriptions.setdefault(key, d)
            for dirpath, _, files in os.walk(sdir):
                for fn in sorted(files):
                    p = os.path.join(dirpath, fn)
                    rel = os.path.relpath(p, sdir)
                    if fn.endswith(".md"):
                        if rel.startswith("references" + os.sep):
                            with open(p, encoding="utf-8") as f:
                                k = len(f.read().splitlines())
                            if k > REFERENCE_MAX:
                                self.add(p, k, "R4", f"{k} lines; reference limit {REFERENCE_MAX}")
                        self.check_text(p, sdir, skill_names, markdown=True)
                    elif fn.endswith((".py", ".sh", ".txt", ".json", ".yaml", ".yml")):
                        self.check_text(p, sdir, skill_names, markdown=False)
        self.check_user_level(names, descriptions)
        return self.findings

    def load_removed_names(self):
        page = os.path.join(self.root, PYAZR_PAGE)
        if not os.path.isfile(page):
            self.add(page, 0, "R10", "no pyazr page: removed names not checked", warning=True)
            return
        with open(page, encoding="utf-8") as f:
            names = removed_pyazr_names(f.read())
        if names is None:
            self.add(page, 0, "R10", "the page has no migration table (a row '* - earlier name')")
            return
        self.removed_names = names

    def check_user_level(self, names, descriptions):
        udir = self.user_skills_dir
        if not udir or not os.path.isdir(udir) or os.path.realpath(udir) == os.path.realpath(self.skills_dir):
            return
        for d in sorted(os.listdir(udir)):
            skill_md = os.path.join(udir, d, "SKILL.md")
            if not os.path.isfile(skill_md):
                continue
            with open(skill_md, encoding="utf-8", errors="replace") as f:
                fields, _ = parse_frontmatter(f.read())
            if not fields:
                continue
            uname = fields.get("name", d)
            udesc = " ".join(fields.get("description", "").lower().split())
            if uname in names:
                self.add(os.path.join(self.skills_dir, names[uname], "SKILL.md"), 2, "R6",
                         f"name {uname!r} is also a user-level skill ({skill_md})")
            if udesc and udesc in descriptions:
                self.add(os.path.join(self.skills_dir, descriptions[udesc], "SKILL.md"), 3, "R6",
                         f"description repeats the user-level skill {uname!r}")


def main(argv=None):
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", default=os.path.dirname(here), help="repository root (default: parent of tools/)")
    ap.add_argument("--today", default=None, help="date for version-note expiry, YYYY-MM-DD (default: today)")
    ap.add_argument("--strict", action="store_true", help="treat a SKILL.md over 500 lines as an error")
    ap.add_argument("--user-skills", default=None,
                    help="user-level skills directory (default: $CLAUDE_USER_SKILLS_DIR, else ~/.claude/skills; "
                         "'none' to skip)")
    args = ap.parse_args(argv)
    today = _dt.date.fromisoformat(args.today) if args.today else _dt.date.today()
    user = args.user_skills or os.environ.get("CLAUDE_USER_SKILLS_DIR") or os.path.expanduser("~/.claude/skills")
    if user == "none":
        user = None
    findings = Checker(args.root, today, strict=args.strict, user_skills_dir=user).run()
    for f in findings:
        print(f)
    errors = [f for f in findings if not f.warning]
    if errors:
        print(f"check_skills: {len(errors)} error(s), {len(findings) - len(errors)} warning(s)", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
