"""DECISIONS.md: the queue of questions a campaign cannot answer from its spec (rmfit.spec).

Format (one entry per question; the header and the ``key: value`` lines are machine-read, prose may follow)::

    ## D-003 [pending] 2026-11-07 13:40 adopt-level-6+-15.23
    question: adopt the 6+ level at 15.233 MeV (gain +7,778, J not separated from 2+/3+)?
    evidence: figures/leveladd_6+_E15.2330.png; REPORT.md
    options: adopt 6+ | adopt 2+ | adopt 3+ | reject
    default_if_unanswered: wait
    answer:
    answered_by:
    answered:
    note:

``rmfit decide D-003 "adopt 6+" --by evaluator --note ...`` fills the answer; ``rmfit decisions`` lists what is pending.
A driver in WAITING state polls ``pending()``; ``default_if_unanswered`` other than ``wait`` lets a budgeted run go on
with the stated default (recorded as the answer with answered_by "default").
"""
import os
import re
import time
from typing import List, Optional

HEAD = re.compile(r"^## (D-\d+) \[(pending|answered)\] (\S+ \S+) (.*)$")
FIELDS = ("question", "evidence", "options", "default_if_unanswered", "answer", "answered_by", "answered", "note")


class Decision(dict):
    @property
    def id(self):
        return self["id"]

    @property
    def pending(self):
        return self["status"] == "pending"

    def render(self):
        L = [f"## {self['id']} [{self['status']}] {self['when']} {self['slug']}"]
        for k in FIELDS:
            L.append(f"{k}: {self.get(k, '')}")
        if self.get("prose"):
            L += ["", self["prose"].rstrip()]
        return "\n".join(L) + "\n"


class DecisionQueue:
    def __init__(self, path):
        self.path = path

    # -- io --------------------------------------------------------------------------------------
    def _read(self) -> List[Decision]:
        if not os.path.exists(self.path):
            return []
        out, cur, prose = [], None, []
        for line in open(self.path):
            m = HEAD.match(line.rstrip("\n"))
            if m:
                if cur is not None:
                    cur["prose"] = "\n".join(prose).strip()
                    out.append(cur)
                cur = Decision(id=m.group(1), status=m.group(2), when=m.group(3), slug=m.group(4))
                prose = []
                continue
            if cur is None:
                continue
            km = re.match(r"^(%s): ?(.*)$" % "|".join(FIELDS), line.rstrip("\n"))
            if km and km.group(1) not in cur:
                cur[km.group(1)] = km.group(2)
            else:
                prose.append(line.rstrip("\n"))
        if cur is not None:
            cur["prose"] = "\n".join(prose).strip()
            out.append(cur)
        return out

    def _write(self, items: List[Decision]):
        tmp = self.path + ".tmp"
        with open(tmp, "w") as fh:
            fh.write("# Decisions queue (rmfit.decisions)\n\nQuestions the campaign cannot answer from its spec.  "
                     "Answer with `rmfit <dir> decide <id> \"<choice>\" [--by who] [--note ...]`; "
                     "`rmfit <dir> decisions` lists the pending ones.\n\n")
            for d in items:
                fh.write(d.render() + "\n")
        os.replace(tmp, self.path)

    # -- api --------------------------------------------------------------------------------------
    def all(self) -> List[Decision]:
        return self._read()

    def pending(self) -> List[Decision]:
        return [d for d in self._read() if d.pending]

    def get(self, did) -> Optional[Decision]:
        for d in self._read():
            if d.id == did:
                return d
        return None

    def add(self, slug, question, evidence="", options=(), default="wait", prose="") -> str:
        items = self._read()
        n = 1 + max([int(d.id.split("-")[1]) for d in items] or [0])
        d = Decision(id=f"D-{n:03d}", status="pending", when=time.strftime("%Y-%m-%d %H:%M"), slug=slug.replace(" ", "-"),
                     question=question, evidence=evidence, options=" | ".join(options) if not isinstance(options, str) else options,
                     default_if_unanswered=default, answer="", answered_by="", answered="", note="", prose=prose)
        items.append(d)
        self._write(items)
        return d.id

    def answer(self, did, choice, by="", note="") -> Decision:
        items = self._read()
        for d in items:
            if d.id == did:
                d["answer"] = choice
                d["answered_by"] = by
                d["answered"] = time.strftime("%Y-%m-%d %H:%M")
                d["note"] = note
                d["status"] = "answered"
                self._write(items)
                return d
        raise KeyError(did)

    def apply_defaults(self, by="default") -> List[Decision]:
        """Answer every pending entry whose default is not ``wait`` with that default (a budgeted driver run)."""
        done = []
        for d in self.pending():
            dflt = (d.get("default_if_unanswered") or "wait").strip()
            if dflt and dflt.lower() != "wait":
                done.append(self.answer(d.id, dflt, by=by, note="default_if_unanswered applied"))
        return done
