"""The campaign's memory: every candidate, every score, every decision.

SQLite for the index and one ``.npz`` per candidate vector.  Nothing is ever
deleted; a candidate's history is its rows in ``scores``.  ``STATUS.md`` is
rendered from here so a fresh session (human or agent) can pick a campaign up
from the file system alone.
"""
import dataclasses
import json
import os
import sqlite3
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .keys import ParamKey, canonical_signs

_SCHEMA = """
CREATE TABLE IF NOT EXISTS structures (
    id INTEGER PRIMARY KEY, name TEXT UNIQUE, azr TEXT, keys TEXT, parent INTEGER,
    note TEXT, created REAL, fixed TEXT DEFAULT '[]');
CREATE TABLE IF NOT EXISTS candidates (
    id INTEGER PRIMARY KEY, structure_id INTEGER, parent_id INTEGER, move TEXT,
    signs TEXT, npz TEXT, status TEXT, created REAL);
CREATE TABLE IF NOT EXISTS scores (
    id INTEGER PRIMARY KEY, cid INTEGER, stage TEXT, objective REAL, chi2 REAL,
    penalty REAL, prior REAL, npts INTEGER, nfev INTEGER, wall REAL, evaluator TEXT,
    status TEXT, message TEXT, created REAL);
CREATE TABLE IF NOT EXISTS segment_scores (
    cid INTEGER, stage TEXT, seg INTEGER, chi2 REAL);
CREATE TABLE IF NOT EXISTS rounds (
    id INTEGER PRIMARY KEY, structure_id INTEGER, kind TEXT, note TEXT, stats TEXT,
    started REAL, finished REAL);
CREATE TABLE IF NOT EXISTS flips (
    structure_id INTEGER, key TEXT, outcome TEXT, delta REAL, cid INTEGER,
    round_id INTEGER, created REAL);
CREATE TABLE IF NOT EXISTS log (
    id INTEGER PRIMARY KEY, ts REAL, event TEXT, detail TEXT);
CREATE INDEX IF NOT EXISTS scores_cid ON scores(cid);
CREATE INDEX IF NOT EXISTS cand_struct ON candidates(structure_id);
"""


def keys_to_json(keys: Sequence[ParamKey]) -> str:
    return json.dumps([dataclasses.asdict(k) for k in keys])


def keys_from_json(text: str) -> List[ParamKey]:
    return [ParamKey(**d) for d in json.loads(text)]


class Ledger:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        self.dir = os.path.dirname(self.path)
        self.cdir = os.path.join(self.dir, "candidates")
        os.makedirs(self.cdir, exist_ok=True)
        self.db = sqlite3.connect(self.path, timeout=60)
        self.db.executescript(_SCHEMA)
        self.db.commit()

    def close(self):
        self.db.close()

    # -- structures -------------------------------------------------------------

    def add_structure(self, name, azr, keys, parent=None, note="", fixed=()):
        cur = self.db.execute("SELECT id FROM structures WHERE name=?", (name,))
        row = cur.fetchone()
        if row:
            return row[0]
        cur = self.db.execute(
            "INSERT INTO structures(name, azr, keys, parent, note, created, fixed) VALUES (?,?,?,?,?,?,?)",
            (name, os.path.abspath(azr), keys_to_json(keys), parent, note, time.time(), keys_to_json(fixed)))
        self.db.commit()
        self.log("structure", f"{name}: {azr} ({len(keys)} free) {note}")
        return cur.lastrowid

    def structure(self, sid) -> dict:
        r = self.db.execute("SELECT id, name, azr, keys, parent, note, fixed FROM structures WHERE id=?", (sid,)).fetchone()
        if not r:
            raise KeyError(sid)
        return dict(id=r[0], name=r[1], azr=r[2], keys=keys_from_json(r[3]), parent=r[4], note=r[5],
                    fixed=keys_from_json(r[6] or "[]"))

    def fixed_of(self, structure_id):
        r = self.db.execute("SELECT fixed FROM structures WHERE id=?", (structure_id,)).fetchone()
        return keys_from_json(r[0] or "[]") if r else []

    def structures(self) -> List[dict]:
        return [dict(id=r[0], name=r[1], azr=r[2], parent=r[3], note=r[4])
                for r in self.db.execute("SELECT id, name, azr, parent, note FROM structures ORDER BY id")]

    # -- candidates -------------------------------------------------------------

    def add_candidate(self, structure_id, x, keys, move=None, parent=None, status="new") -> int:
        x = np.asarray(x, float)
        signs = canonical_signs(x, keys, self.fixed_of(structure_id))
        cur = self.db.execute(
            "INSERT INTO candidates(structure_id, parent_id, move, signs, npz, status, created) "
            "VALUES (?,?,?,?,?,?,?)",
            (structure_id, parent, json.dumps(move) if move is not None else None, signs, "", status, time.time()))
        cid = cur.lastrowid
        npz = os.path.join(self.cdir, f"c{cid:06d}.npz")
        np.savez(npz, x=x, keys=keys_to_json(keys))
        self.db.execute("UPDATE candidates SET npz=? WHERE id=?", (npz, cid))
        self.db.commit()
        return cid

    def load(self, cid) -> Tuple[np.ndarray, List[ParamKey]]:
        r = self.db.execute("SELECT npz FROM candidates WHERE id=?", (cid,)).fetchone()
        if not r:
            raise KeyError(cid)
        d = np.load(r[0], allow_pickle=False)
        return np.asarray(d["x"], float), keys_from_json(str(d["keys"]))

    def candidate(self, cid) -> dict:
        r = self.db.execute("SELECT id, structure_id, parent_id, move, signs, status, created "
                            "FROM candidates WHERE id=?", (cid,)).fetchone()
        if not r:
            raise KeyError(cid)
        return dict(id=r[0], structure_id=r[1], parent_id=r[2],
                    move=json.loads(r[3]) if r[3] else None, signs=r[4], status=r[5], created=r[6])

    def set_status(self, cid, status):
        self.db.execute("UPDATE candidates SET status=? WHERE id=?", (status, cid))
        self.db.commit()

    def record(self, cid, result, stage, evaluator="sharded"):
        """``result`` is an engine.Score or a local.FitResult."""
        seg = getattr(result, "segment_chi2", None) or {}
        self.db.execute(
            "INSERT INTO scores(cid, stage, objective, chi2, penalty, prior, npts, nfev, wall, evaluator, "
            "status, message, created) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (cid, stage, float(result.objective), float(getattr(result, "chi2", np.nan)),
             float(getattr(result, "penalty", np.nan)), float(getattr(result, "prior", 0.0)),
             int(getattr(result, "npts", 0) or 0), int(getattr(result, "nfev", 0) or 0),
             float(getattr(result, "wall", 0.0) or 0.0), evaluator,
             getattr(result, "status", None), getattr(result, "message", None), time.time()))
        if seg:
            self.db.executemany("INSERT INTO segment_scores(cid, stage, seg, chi2) VALUES (?,?,?,?)",
                                [(cid, stage, int(k), float(v)) for k, v in seg.items()])
        self.db.commit()

    def scores(self, cid) -> List[dict]:
        return [dict(stage=r[0], objective=r[1], chi2=r[2], penalty=r[3], evaluator=r[4], status=r[5], created=r[6])
                for r in self.db.execute("SELECT stage, objective, chi2, penalty, evaluator, status, created "
                                         "FROM scores WHERE cid=? ORDER BY id", (cid,))]

    def segment_scores(self, cid, stage) -> Dict[int, float]:
        return {r[0]: r[1] for r in self.db.execute(
            "SELECT seg, chi2 FROM segment_scores WHERE cid=? AND stage=?", (cid, stage))}

    def best(self, structure_id, stage="verified") -> Optional[Tuple[int, float]]:
        r = self.db.execute(
            "SELECT s.cid, MIN(s.objective) FROM scores s JOIN candidates c ON c.id = s.cid "
            "WHERE c.structure_id=? AND s.stage=? AND c.status != 'rejected'", (structure_id, stage)).fetchone()
        if not r or r[0] is None:
            return None
        return int(r[0]), float(r[1])

    def signs_of(self, cid, structure_id=None) -> str:
        """The candidate's canonical sign string recomputed from its saved vector.  The ``signs`` column of a
        ledger written before 2026-11-12 may hold the old, not gauge-invariant form (keys.canonical_signs), and a
        ledger can mix both while older jobs still write to it, so comparisons use this, not the column.  Falls
        back to the column when the .npz is missing."""
        if structure_id is None:
            r = self.db.execute("SELECT structure_id FROM candidates WHERE id=?", (cid,)).fetchone()
            structure_id = r[0] if r else None
        try:
            x, keys = self.load(cid)
        except Exception:
            r = self.db.execute("SELECT signs FROM candidates WHERE id=?", (cid,)).fetchone()
            return r[0] if r else ""
        return canonical_signs(x, keys, self.fixed_of(structure_id))

    def population(self, structure_id, k=4, stage="verified") -> List[Tuple[int, float, str]]:
        """The ``k`` best verified candidates with distinct canonical sign patterns."""
        rows = self.db.execute(
            "SELECT c.id, MIN(s.objective) FROM scores s JOIN candidates c ON c.id = s.cid "
            "WHERE c.structure_id=? AND s.stage=? AND c.status != 'rejected' "
            "GROUP BY c.id ORDER BY MIN(s.objective)", (structure_id, stage)).fetchall()
        out, seen = [], set()
        for cid, obj in rows:
            signs = self.signs_of(cid, structure_id)
            if signs in seen:
                continue
            seen.add(signs)
            out.append((int(cid), float(obj), signs))
            if len(out) >= k:
                break
        return out

    def duplicate_of(self, structure_id, x, keys, objective, rel=5e-4, stage="verified") -> Optional[int]:
        signs = canonical_signs(x, keys, self.fixed_of(structure_id))
        for cid, obj in self.db.execute(
                "SELECT c.id, s.objective FROM scores s JOIN candidates c ON c.id = s.cid "
                "WHERE c.structure_id=? AND s.stage=?", (structure_id, stage)).fetchall():
            if abs(obj - objective) <= rel * abs(objective) and self.signs_of(cid, structure_id) == signs:
                return int(cid)
        return None

    # -- flips and rounds ---------------------------------------------------------

    def record_flip(self, structure_id, key, outcome, delta, cid=None, round_id=None):
        self.db.execute("INSERT INTO flips(structure_id, key, outcome, delta, cid, round_id, created) "
                        "VALUES (?,?,?,?,?,?,?)",
                        (structure_id, str(key), outcome, None if delta is None else float(delta),
                         cid, round_id, time.time()))
        self.db.commit()

    def flip_outcomes(self, structure_id) -> Dict[str, List[Tuple[str, float]]]:
        out: Dict[str, List[Tuple[str, float]]] = {}
        for key, outcome, delta in self.db.execute(
                "SELECT key, outcome, delta FROM flips WHERE structure_id=? ORDER BY rowid", (structure_id,)):
            out.setdefault(key, []).append((outcome, delta))
        return out

    def start_round(self, structure_id, kind, note="") -> int:
        cur = self.db.execute("INSERT INTO rounds(structure_id, kind, note, stats, started) VALUES (?,?,?,?,?)",
                              (structure_id, kind, note, "{}", time.time()))
        self.db.commit()
        self.log("round", f"start {kind} #{cur.lastrowid} on structure {structure_id} {note}")
        return cur.lastrowid

    def finish_round(self, rid, **stats):
        self.db.execute("UPDATE rounds SET stats=?, finished=? WHERE id=?",
                        (json.dumps(stats, default=float), time.time(), rid))
        self.db.commit()
        self.log("round", f"finish #{rid}: " + json.dumps(stats, default=float))

    def rounds(self, structure_id=None) -> List[dict]:
        q = "SELECT id, structure_id, kind, note, stats, started, finished FROM rounds"
        args = ()
        if structure_id is not None:
            q += " WHERE structure_id=?"
            args = (structure_id,)
        return [dict(id=r[0], structure_id=r[1], kind=r[2], note=r[3], stats=json.loads(r[4] or "{}"),
                     started=r[5], finished=r[6]) for r in self.db.execute(q + " ORDER BY id", args)]

    def log(self, event, detail=""):
        self.db.execute("INSERT INTO log(ts, event, detail) VALUES (?,?,?)", (time.time(), event, detail))
        self.db.commit()
        with open(os.path.join(self.dir, "LOG.jsonl"), "a") as fh:
            fh.write(json.dumps(dict(ts=time.strftime("%Y-%m-%d %H:%M:%S"), event=event, detail=detail)) + "\n")

    # -- status -------------------------------------------------------------------

    def status(self) -> dict:
        out = dict(structures=[], rounds=[], last_log=[])
        for s in self.structures():
            b = self.best(s["id"])
            bs = self.best(s["id"], stage="global") or self.best(s["id"], stage="polish")
            n = self.db.execute("SELECT COUNT(*) FROM candidates WHERE structure_id=?", (s["id"],)).fetchone()[0]
            out["structures"].append(dict(id=s["id"], name=s["name"], azr=s["azr"], candidates=n,
                                          best_verified=b, best_unverified=bs, note=s["note"]))
        out["rounds"] = self.rounds()[-8:]
        out["last_log"] = [dict(ts=r[0], event=r[1], detail=r[2]) for r in
                           self.db.execute("SELECT ts, event, detail FROM log ORDER BY id DESC LIMIT 12")][::-1]
        return out

    def write_status_md(self, path=None, extra: Optional[Dict[str, str]] = None):
        path = path or os.path.join(self.dir, "STATUS.md")
        st = self.status()
        lines = [f"# rmfit campaign status ({time.strftime('%Y-%m-%d %H:%M')})", ""]
        for s in st["structures"]:
            bv = s["best_verified"]
            bu = s["best_unverified"]
            lines.append(f"## structure {s['id']}: {s['name']}")
            lines.append(f"- file: `{s['azr']}`  candidates: {s['candidates']}")
            lines.append(f"- best verified: " + (f"candidate {bv[0]} objective {bv[1]:.3f}" if bv else "none"))
            if bu:
                lines.append(f"- best unverified (sharded): candidate {bu[0]} objective {bu[1]:.3f}")
            if s["note"]:
                lines.append(f"- note: {s['note']}")
            lines.append("")
        if st["rounds"]:
            lines.append("## recent rounds")
            for r in st["rounds"]:
                fin = "running" if r["finished"] is None else time.strftime("%m-%d %H:%M", time.localtime(r["finished"]))
                lines.append(f"- #{r['id']} {r['kind']} (structure {r['structure_id']}, {fin}): "
                             + json.dumps({k: v for k, v in r["stats"].items() if k != "screen"}, default=float) + (f" {r['note']}" if r["note"] else ""))
            lines.append("")
        if extra:
            for k, v in extra.items():
                lines.append(f"## {k}")
                lines.append(str(v))
                lines.append("")
        lines.append("## recent log")
        for e in st["last_log"]:
            lines.append(f"- {time.strftime('%m-%d %H:%M', time.localtime(e['ts']))} {e['event']}: {e['detail']}")
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        return path
