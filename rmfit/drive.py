"""``rmfit drive``: the closed loop (DESIGN.md Phase 2.1; retrospective finding 2.2 "nobody closes the loop").

A campaign is driven by a resumable state machine whose every state is one checkpoint in ``drive_state.json``::

    build phase (needs --base <azr> [--sav] and <dir>/spec.json; no campaign.json yet)
      PIN_CHECK -> AUDIT -> REPRODUCE -> SPEC_CHECK -> SCALE_ANCHOR -> SEED_CHANNELS -> PAIR_COMPLETENESS -> LADDER
      -> CANDIDATES -> INIT
    loop phase (campaign.json exists)
      TIGHTEN -> SIGN -> SHIFTS_FREE -> STRUCTURE_ADD -> SIGN_POST -> DATASET_OFF -> REVIEW -> EXPORT -> CYCLE
      CYCLE -> TIGHTEN (another pass) | DONE

What each state does is the skill's process (rmfit-campaign SKILL.md section 2 "First fits", section 3 steps 1-7, section 4
stopping rules) with the thresholds read from ``spec.json`` (``budget``, ``stop_rules``, ``pre_authorized``).  Anything the
spec does not answer becomes a ``DECISIONS.md`` entry; a pending entry whose default is ``wait`` ends the leg with status
WAITING, and the next leg applies the answers (``adopt-level-*`` answered "adopt" adopts the structure).  Nothing here
decides data treatment: dataset-off tests, missing pairs and review flags are recorded and queued, never acted on.

One leg = one ``rmfit <dir> drive`` call (one SGE job of ``run_crc_drive``): states run until the state machine stops
(DONE / WAITING / STOP), the wall cap ``--max-minutes`` is reached (status CONTINUE: the job script resubmits itself), or
``<dir>/STOP_DRIVE`` exists.  Engine states refuse the login node unless ``--local`` (fixture tests).  ``--dry-run`` prints
the plan from the current state without touching anything.  ``DRIVE.md`` is rewritten after every state for the supervisor.

Build-phase scope (honest limits, 2026-11-08): SCALE_ANCHOR applies the spec's data-set roles (anchor = norm fixed, free =
norm free with the prior as its systematic, excluded = segments off) to a working copy ``<stem>_drive.azr``; SEED_CHANNELS
lifts every FREE width that sits at exactly zero to ``seed_frac`` x gamma_W in the engine and bakes (zero amplitudes never
move, skill 2026-09-25); it does not add channel lines.  PAIR_COMPLETENESS lists the ENSDF excited states of the model's
residual nuclei whose partition opens below the top of the fit window and has no pair -- a finding and a non-blocking
decision for the supervisor, who writes the pair builder (an archived campaign pattern).
LADDER runs only when the spec gives ``model.ladder_Ex`` (a from-scratch device).  CANDIDATES writes ``candidates.json``
from ENSDF (broad first) when none exists; every candidate needs a ``hypothesis``.
"""
import dataclasses
import json
import os
import re
import shutil
import socket
import time
from typing import List, Optional

from . import AZURE2_ROOT  # noqa: F401
from .backend import PYAZR3
from .spec import FitSpec, Finding, write_findings_md

BUILD_STATES = ["PIN_CHECK", "AUDIT", "REPRODUCE", "SPEC_CHECK", "SCALE_ANCHOR", "SEED_CHANNELS", "PAIR_COMPLETENESS",
                "LADDER", "CANDIDATES", "INIT", "RECALL"]
LOOP_STATES = ["TIGHTEN", "SIGN", "POPULATION", "SHIFTS_FREE", "LEVELSCAN", "STRUCTURE_ADD", "SIGN_POST", "DATASET_OFF", "REVIEW", "EXPORT",
               "CYCLE", "UQ"]
STATES = BUILD_STATES + LOOP_STATES
ENGINE_STATES = {"REPRODUCE", "SEED_CHANNELS", "LADDER", "INIT", "TIGHTEN", "SIGN", "POPULATION", "SHIFTS_FREE", "LEVELSCAN", "STRUCTURE_ADD",
                 "SIGN_POST", "DATASET_OFF", "REVIEW", "EXPORT", "UQ"}
STATUS_WORDS = ("CONTINUE", "WAITING", "DONE", "STOP", "ERROR")
STATE_FILE = "drive_state.json"
STOP_FILE = "STOP_DRIVE"

# spec defaults (budget / stop_rules), SKILL.md sections 3-4 and templates/SPEC_TEMPLATE.json
DEFAULT_BUDGET = dict(polish_nfev=150, polish_rounds_max=4, sign_rounds_max=4, level_tests_max=10, level_nfev=150,
                      pinned_nfev=60, dataset_off_max=2, dataset_off_nfev=100, ladder_nfev=100, cycles_max=3,
                      node_hours_max=0, grid_check=True, figures=6,
                      # POPULATION (2026-11-08): restarts + crossover per cycle; restart_nfev 0 = polish_nfev
                      restarts_per_cycle=4, restart_nfev=0, crossover_passes=2, population_rounds_max=3)
DEFAULT_STRUCTURE = dict(seed_omitted=True, broad_first_kev=10.0, narrow_phase=True, levelscan_per_cycle=1, levelscan_top=3,
                         levelscan_gamma_kev=None, local_release_mev=0.5, local_nfev=80, recall_levels=2, recall_min_pass=0.5)
DEFAULT_STOP = dict(polish_gain_rel_min=5e-4, sign_rounds_without_accept=2, level_adds_marginal_in_a_row=2,
                    cycle_gain_rel_min=1e-3, population_gain_rel_min=1e-3, dataset_off_extra_pct=5.0, reproduce_rtol=1e-3,
                    # sign rounds only after the polish converged (skill section 3, step 1 "when done" before step 2): a flip
                    # is credited with what it gains beyond the control polish, and an unconverged incumbent moves by more
                    # than any flip (in one benchmark 18 sign rounds accepted nothing while the polishes
                    # still gained 1-2 % per round).  0 disables the gate.
                    sign_after_converged_polish=1)

SYMBOLS = ["n", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K",
           "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
           "Y", "Zr"]


class DriveStop(Exception):
    """Ends the leg with the given status word (STOP or ERROR) and message."""

    def __init__(self, status, message):
        super().__init__(message)
        self.status = status


def _now():
    return time.strftime("%Y-%m-%d %H:%M")


def _submit():
    """How to start the drive's job script on this site (qsub / sbatch / nohup ./run_crc_drive &)."""
    from . import site
    return site.submit_hint("run_crc_drive")


def in_sge_job():
    """Inside a scheduler job (SGE or Slurm; the name is historical)."""
    from . import site
    return site.in_batch_job()


class Driver:
    def __init__(self, directory, base=None, sav=None, local=False, log=None):
        self.dir = os.path.abspath(directory)
        os.makedirs(self.dir, exist_ok=True)
        self.local = local
        self._print = log or print
        self.spec_path = os.path.join(self.dir, "spec.json")
        if not os.path.exists(self.spec_path):
            raise DriveStop("ERROR", f"no spec.json in {self.dir}: the spec is the objective; copy rmfit/templates/SPEC_TEMPLATE.json, "
                                     "fill it in (or start from a spec_draft.json) and run `rmfit <dir> spec check` first")
        self.spec = FitSpec.load(self.spec_path)
        self.budget = dict(DEFAULT_BUDGET, **{k: v for k, v in (self.spec.get("budget") or {}).items() if _isnum(v) or isinstance(v, bool)})
        self.stop = dict(DEFAULT_STOP, **{k: v for k, v in (self.spec.get("stop_rules") or {}).items() if _isnum(v)})
        self.structure = dict(DEFAULT_STRUCTURE, **{k: v for k, v in (self.spec.get("structure") or {}).items()
                                                   if _isnum(v) or isinstance(v, bool) or v is None})
        self.state_path = os.path.join(self.dir, STATE_FILE)
        self.state = self._load_state()
        if base:
            base = os.path.abspath(base)
            if self.state.get("base") and os.path.abspath(self.state["base"]) != base:
                self.log(f"WARNING --base {base} differs from the recorded base {self.state['base']}; the recorded one is kept")
            else:
                self.state["base"] = base
        if sav:
            self.state["sav"] = os.path.abspath(sav)
        self._campaign = None
        self.t_leg = time.time()

    # -- state -------------------------------------------------------------------------------------
    def _load_state(self) -> dict:
        if os.path.exists(self.state_path):
            with open(self.state_path) as fh:
                st = json.load(fh)
        else:
            st = {}
        st.setdefault("state", "PIN_CHECK")
        st.setdefault("status", "CONTINUE")
        st.setdefault("cycle", 1)
        st.setdefault("history", [])
        st.setdefault("tested_candidates", [])
        st.setdefault("applied_decisions", [])
        st.setdefault("consumed_decisions", [])
        st.setdefault("dataset_off_tested", [])
        st.setdefault("counters", {})
        st.setdefault("wall_hours", 0.0)
        st.setdefault("findings", {})
        st.setdefault("created", _now())
        return st

    def save(self):
        self.state["updated"] = _now()
        tmp = self.state_path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(self.state, fh, indent=1, default=str)
        os.replace(tmp, self.state_path)

    def log(self, *a):
        msg = " ".join(str(x) for x in a)
        line = f"{time.strftime('%H:%M:%S')} drive {msg}"
        self._print(line)
        with open(os.path.join(self.dir, "campaign.log"), "a") as fh:
            fh.write(line + "\n")

    def note(self, text):
        """A dated readme entry (the campaign's lab notebook; same format as Campaign.note)."""
        if self.campaign is not None:
            self.campaign.note("drive: " + text)
        else:
            with open(os.path.join(self.dir, "readme"), "a") as fh:
                fh.write(f"\n{_now()} (rmfit drive): {text}\n")

    @property
    def campaign(self):
        if self._campaign is None and os.path.exists(os.path.join(self.dir, "campaign.json")):
            from .campaign import Campaign
            self._campaign = Campaign(self.dir, enforce_engine=True)
        return self._campaign

    @property
    def decisions(self):
        if self.campaign is not None:
            return self.campaign.decisions
        from .decisions import DecisionQueue
        return DecisionQueue(os.path.join(self.dir, "DECISIONS.md"))

    @property
    def model(self) -> str:
        """The model the build phase works on (the working copy once SCALE_ANCHOR wrote it), or the campaign's base."""
        if self.campaign is not None:
            return self.campaign.base_azr
        m = self.state.get("model") or self.state.get("base")
        if not m:
            raise DriveStop("ERROR", "no base model: run `rmfit <dir> drive --base <file.azr> [--sav <file.sav>]` the first time")
        return m

    @property
    def counters(self) -> dict:
        return self.state["counters"]

    def _stem(self):
        return os.path.splitext(os.path.basename(self.state.get("base") or self.model))[0]

    # -- decisions ---------------------------------------------------------------------------------
    def ask(self, slug, question, options, default="wait", evidence="") -> Optional[str]:
        """The answer of the latest decision with this slug if it has one and was not consumed yet; otherwise queue it
        (unless one is pending) and return None.  A default other than ``wait`` lets the leg go on with that answer."""
        slug = slug.replace(" ", "-")
        for d in reversed(self.decisions.all()):
            if d["slug"] != slug:
                continue
            if d.pending:
                return None
            if d.id in self.state["consumed_decisions"]:
                break
            self.state["consumed_decisions"].append(d.id)
            self.save()
            return d.get("answer", "")
        did = self.decisions.add(slug, question, evidence=evidence, options=options, default=default)
        self.log(f"decision {did} queued ({slug}): {question}" + (f"  [default if unanswered: {default}]" if default != "wait" else ""))
        if default and default.lower() != "wait":
            self.decisions.answer(did, default, by="default", note="default_if_unanswered applied by rmfit drive")
            self.state["consumed_decisions"].append(did)
            self.save()
            return default
        return None

    def waiting_on(self) -> List[str]:
        return [f"{d.id} {d['slug']}" for d in self.decisions.pending()
                if (d.get("default_if_unanswered") or "wait").strip().lower() == "wait"]

    def apply_decisions(self):
        """Answered ``adopt-level-*`` decisions (queued by Campaign._run_level_add as accept-pending) are applied here."""
        c = self.campaign
        if c is None:
            return
        for d in self.decisions.all():
            if d.pending or d.id in self.state["applied_decisions"] or not d["slug"].startswith("adopt-level"):
                continue
            ans = (d.get("answer") or "").strip().lower()
            m = re.search(r"structure (\d+)", d.get("question", "") + " " + d.get("evidence", ""))
            self.state["applied_decisions"].append(d.id)
            if ans.startswith("adopt") and m:
                sid = int(m.group(1))
                c.adopt_structure(sid, note=f"(decision {d.id} answered '{d.get('answer')}' by {d.get('answered_by', '')})")
                self.counters["adopted_since_sign"] = 1
                self.counters["population_exhausted"] = 0      # a new structure: the population stage may pay again
                self.note(f"decision {d.id} applied: structure {sid} adopted as the base model ({d['slug']})")
            else:
                self.log(f"decision {d.id} ({d['slug']}) answered '{d.get('answer')}': nothing to apply"
                         + ("" if m else " (no structure id in the question)"))
            self.save()

    # -- the leg -----------------------------------------------------------------------------------
    def run(self, max_minutes: float = 0.0, dry_run: bool = False) -> str:
        if dry_run:
            self.print_plan()
            return self.state["status"]
        status = "CONTINUE"
        self.state["status"] = "RUNNING"
        self.state["leg_started"] = _now()
        self.state["leg_host"] = socket.gethostname()
        self.save()
        self.log(f"leg starts at {self.state['state']} (cycle {self.state['cycle']}) on {socket.gethostname()}"
                 + (f", wall cap {max_minutes:g} min" if max_minutes else ""))
        try:
            self.apply_decisions()
            while True:
                if os.path.exists(os.path.join(self.dir, STOP_FILE)):
                    os.remove(os.path.join(self.dir, STOP_FILE))
                    self.note("STOP_DRIVE found: the drive stopped after the current state")
                    status = "STOP"
                    break
                w = self.waiting_on()
                if w:
                    self.log("waiting on " + ", ".join(w))
                    status = "WAITING"
                    break
                if self.state["state"] == "DONE":
                    status = "DONE"
                    break
                nh = float(self.budget.get("node_hours_max") or 0)
                if nh > 0 and self.state["wall_hours"] + (time.time() - self.t_leg) / 3600 > nh:
                    self.note(f"budget.node_hours_max {nh:g} h reached: the drive stopped")
                    status = "STOP"
                    break
                if max_minutes and (time.time() - self.t_leg) / 60 > max_minutes:
                    self.log(f"wall cap {max_minutes:g} min reached: the leg ends, the next one resumes at {self.state['state']}")
                    status = "CONTINUE"
                    break
                st = self.state["state"]
                if st in ENGINE_STATES and not in_sge_job() and not self.local:
                    raise DriveStop("STOP", f"state {st} needs the engine: refusing the login node ({_submit()}, or --local "
                                            "for a fixture test)")
                self.step(st)
        except DriveStop as e:
            status = e.status
            self.log(f"{status}: {e}")
            self.state["message"] = str(e)
        except Exception as err:  # noqa: BLE001 -- the leg must leave a readable state behind
            import traceback
            status = "ERROR"
            self.state["message"] = f"{type(err).__name__}: {err}"
            self.log("ERROR in state", self.state["state"], ":", traceback.format_exc())
        finally:
            self.state["wall_hours"] = float(self.state["wall_hours"]) + (time.time() - self.t_leg) / 3600
            self.state["status"] = status
            self.state["leg_finished"] = _now()
            self.save()
            self.write_md()
            if self._campaign is not None:
                self._campaign.close()
                self._campaign.ledger.write_status_md()
            self.log(f"DRIVE_STATUS: {status} at {self.state['state']} (leg {(time.time()-self.t_leg)/60:.1f} min, "
                     f"total {self.state['wall_hours']:.2f} h)")
        return status

    def step(self, st):
        handler = getattr(self, "_s_" + st.lower())
        t0 = time.time()
        self.log(f"-- {st} --")
        result = handler() or "next"
        rec = dict(state=st, cycle=self.state["cycle"], when=_now(), minutes=round((time.time() - t0) / 60, 2),
                   result=result if isinstance(result, str) else result[0],
                   note=(result[1] if isinstance(result, tuple) and len(result) > 1 else ""))
        self.state["history"].append(rec)
        if rec["result"] == "next":
            self.state["state"] = self.next_state(st)
        elif rec["result"] == "stay":
            pass
        elif rec["result"] in STATES or rec["result"] == "DONE":
            self.state["state"] = rec["result"]
        self.save()
        self.write_md()

    @staticmethod
    def next_state(st):
        i = STATES.index(st)
        return STATES[i + 1] if i + 1 < len(STATES) else "DONE"

    def print_plan(self):
        cur = self.state["state"]
        self._print(f"rmfit drive plan for {self.dir}  (status {self.state['status']}, cycle {self.state['cycle']}, "
                    f"wall {self.state['wall_hours']:.2f} h)")
        self._print(f"  spec {self.spec_path} sha {str(self.spec.sha)[:10]}; campaign.json {'present' if self.campaign is not None else 'absent'}")
        try:
            self._print(f"  model {self.model}" + (f"  sav {self.state['sav']}" if self.state.get("sav") else ""))
        except DriveStop as e:
            self._print(f"  model: {e}")
        self._print(f"  budget {json.dumps(self.budget)}")
        self._print(f"  stop_rules {json.dumps(self.stop)}")
        done = {h["state"] for h in self.state["history"] if h["result"] != "stay"}
        for s in STATES:
            mark = "->" if s == cur else ("ok" if s in done else "  ")
            need = "engine" if s in ENGINE_STATES else "      "
            self._print(f"  {mark} {s:18s} {need}  {PLAN_TEXT.get(s, '')}")
        pend = self.decisions.pending()
        self._print(f"  pending decisions: {len(pend)}" + ("".join(f"\n     {d.id} {d['slug']}: {d.get('question', '')[:100]}" for d in pend)))
        if self.state.get("message"):
            self._print(f"  last message: {self.state['message']}")

    # -- DRIVE.md ----------------------------------------------------------------------------------
    def write_md(self):
        L = [f"# rmfit drive -- {os.path.basename(self.dir)} ({_now()})", "",
             f"status **{self.state['status']}** at state **{self.state['state']}**, cycle {self.state['cycle']}, "
             f"{self.state['wall_hours']:.2f} node-hours so far; state file `{STATE_FILE}`; stop cleanly with `touch {STOP_FILE}`.", ""]
        if self.state.get("message"):
            L += [f"last message: {self.state['message']}", ""]
        try:
            if self.campaign is not None:
                cid, x, obj, keys = self.campaign.incumbent()
                L += [f"incumbent: candidate {cid}, verified objective {obj:.3f}; base model `{os.path.relpath(self.campaign.base_azr, self.dir)}`", ""]
        except Exception as err:  # noqa: BLE001
            L += [f"incumbent: n/a ({err})", ""]
        pend = self.decisions.pending()
        L += ["## Decisions pending", ""]
        L += [f"- {d.id} `{d['slug']}` (default {d.get('default_if_unanswered', 'wait')}): {d.get('question', '')}" for d in pend] or ["- none"]
        L += ["", "## States", "", "| state | status | last run | min | result | note |", "|---|---|---|---:|---|---|"]
        last = {}
        for h in self.state["history"]:
            last[h["state"]] = h
        for s in STATES:
            h = last.get(s)
            stat = "current" if s == self.state["state"] else ("done" if h and h["result"] != "stay" else ("running" if h else "pending"))
            L.append(f"| {s} | {stat} | {h['when'] if h else ''} | {h['minutes'] if h else ''} | {h['result'] if h else ''} | "
                     f"{(h['note'] or '').replace('|', '/')[:160] if h else ''} |")
        if self.state["findings"]:
            L += ["", "## Build findings", ""]
            for k, rows in self.state["findings"].items():
                L.append(f"### {k}")
                L += [f"- {r}" for r in rows] or ["- none"]
                L.append("")
        L += ["", "## History (last 15)", ""]
        L += [f"- {h['when']} {h['state']} (cycle {h['cycle']}, {h['minutes']} min): {h['result']} {h['note']}" for h in self.state["history"][-15:]]
        with open(os.path.join(self.dir, "DRIVE.md"), "w") as fh:
            fh.write("\n".join(L) + "\n")

    # -- helpers -----------------------------------------------------------------------------------
    def _finding_rows(self, key, rows):
        self.state["findings"][key] = list(rows)[:60]

    def _ex_top(self) -> Optional[float]:
        """Top of the fit window in Ex (MeV): spec model.fit_window.Ex_MeV, else the highest data point of any active segment."""
        fw = (self.spec.get("model") or {}).get("fit_window") or {}
        if fw.get("Ex_MeV"):
            return float(fw["Ex_MeV"][1])
        try:
            from .backend import SegmentSet
            from .model import rows_in_segment
            if PYAZR3:
                from .windows import pair_separation, cm_conversion
                text = self.model                  # pyazr 3: both read the project file
            else:
                from .model import cm_conversion
                from .windows import pair_separation
                text = open(self.model).read()
            segs = SegmentSet.from_file(self.model)
            base = os.path.dirname(self.model)
            top = None
            for s in segs:
                if not s.active:
                    continue
                r = rows_in_segment(s, base)
                if r is None or r.shape[0] == 0:
                    continue
                e = float(r[:, 0].max())
                e = min(e, s.energy_max) if s.energy_max > 0 else e
                ex = pair_separation(text, s.entrance_key) + cm_conversion(text, s.entrance_key) * e
                top = ex if top is None else max(top, ex)
            return top
        except Exception as err:  # noqa: BLE001
            self.log(f"fit-window top from the segments failed ({err})")
            return None

    def _seed_vector(self, ev):
        from .campaign import Campaign
        return Campaign._seed_vector(self, ev, self.state.get("sav"))

    def _job_script(self, overwrite=False) -> str:
        return write_job_script(self.dir, self.spec, self.campaign.cfg if self.campaign is not None else None, overwrite=overwrite)

    # ==================================================================================================
    # build phase
    # ==================================================================================================
    def _s_pin_check(self):
        from . import enginepin
        cur = enginepin.fingerprint()
        self.state["engine"] = cur
        if self.campaign is not None:
            rec = self.campaign.cfg.get("engine")
            if rec and rec.get("so_md5") and rec["so_md5"] != cur["so_md5"]:
                raise DriveStop("STOP", f"engine pin mismatch: campaign {enginepin.short(rec)} vs imported {enginepin.short(cur)}")
            return "next", f"campaign pin {enginepin.short(rec) if rec else '(none)'} matches"
        eng = self.spec.get("engine") or {}
        if eng.get("so_md5") and eng["so_md5"] != cur["so_md5"]:
            raise DriveStop("STOP", f"spec.engine.so_md5 {eng['so_md5'][:8]} != imported engine {enginepin.short(cur)}: run from the "
                                    f"spec's azure2_root ({eng.get('azure2_root')})")
        if os.path.realpath(cur.get("azure2_root", "")) == os.path.realpath(os.path.expanduser("~/AZURE2")) and not self.local:
            self.log("WARNING the imported engine is the SHARED checkout ~/AZURE2 (weekly releases replace it): a campaign should "
                     "pin a frozen copy; INIT will pin what this job imports")
        self._job_script()
        return "next", f"imported engine {enginepin.short(cur)}"

    def _s_audit(self):
        from .audit import audit_model, write_audit
        cfg = self.campaign.cfg if self.campaign is not None else None
        fs = audit_model(self.model, campaign_dir=self.dir, campaign_cfg=cfg, run_target_energies=True)
        # 2026-11-09 (from-scratch demo): a data set the spec EXCLUDES is not fitted, so its data findings cannot block the build;
        # they are kept in AUDIT.md at INFO with the reason (SCALE_ANCHOR switches the set off after this state)
        excluded = [ds["file"] for ds in self.spec.data_sets if ds.get("file") and ds.get("role") in ("excluded", "excluded_in_window")]
        if excluded and self.campaign is None:
            fs2 = []
            for f in fs:
                hit = next((x for x in excluded if x in f.message), None)
                if hit and f.severity in ("HIGH", "MEDIUM") and f.code.startswith(("data-", "norm-")):
                    f = Finding("INFO", f.code, f.message + f" [data set excluded by spec.json: not fitted]", f.evidence)
                fs2.append(f)
            fs = fs2
        counts = write_audit(fs, os.path.join(self.dir, "AUDIT.md"), self.model)
        highs = [f for f in fs if f.severity == "HIGH"]
        self._finding_rows("AUDIT.md", [f.line() for f in sorted(fs, key=lambda f: {"HIGH": 0, "MEDIUM": 1}.get(f.severity, 2))][:30])
        if highs:
            ans = self.ask("audit-high", f"AUDIT.md has {len(highs)} HIGH finding(s): " + "; ".join(f.message[:90] for f in highs[:3])
                           + ". Fix the input and answer 'fixed' (the audit re-runs), or 'accepted' to go on with them recorded.",
                           ("fixed", "accepted"), default="wait", evidence="AUDIT.md")
            if ans is None:
                return "stay", f"{len(highs)} HIGH -> waiting"
            if ans.lower().startswith("fixed"):
                return "stay", "re-running the audit after 'fixed'"
            self.note(f"AUDIT.md HIGH findings accepted by decision ({', '.join(f.code for f in highs)})")
        if self.campaign is None and not self.spec.pre_authorized.get("skip_subagent_audits", False):
            r = self._subagent_audit_gate(fs)
            if r is not None:
                return r
        return "next", ", ".join(f"{k} {v}" for k, v in counts.items())

    def _subagent_audit_gate(self, machine_findings):
        """The three idle-time reviews (skill section 3b: input, data, literature) BEFORE the first engine state.  In a from-scratch
        demo (2026-11-09) the lab-frame error the data audit found had been fitted for an hour, because the audits ran
        in the session's idle time after the job was submitted.  The drive cannot start a subagent: it writes the prompts
        (audits/<kind>/PROMPT.md, rmfit.audit_subagent) and waits for the supervising session to run them, verify every
        FINDINGS.md row (`rmfit <dir> audit --status KIND N confirmed|rejected|decision`) and answer 'done'; 'skip' builds
        without them (recorded in the readme; `pre_authorized.skip_subagent_audits` for fixtures and benchmarks).  A row
        left in status 'decision' is the evaluator's: one wait decision each, the answer recorded on the row."""
        from . import audit_subagent as asub
        written = []
        for kind in asub.KINDS:
            if not os.path.exists(os.path.join(self.dir, "audits", kind, "PROMPT.md")):
                try:
                    asub.write_subagent_audit(kind, self.model, self.dir, self.spec, machine_findings)
                    written.append(kind)
                except Exception as err:
                    self.log(f"audit prompt {kind} not written: {err}")
        if written:
            self.note(f"AUDIT: subagent prompts written under audits/ ({', '.join(written)}); the drive waits for their verified FINDINGS.md")
        summ = asub.findings_summary(self.dir)
        nrows = sum(x["n"] for x in summ)
        cand = sum(len(x["unverified"]) for x in summ)
        ans = self.ask("audit-subagents", f"three audit prompts under audits/{{input,data,literature}}/PROMPT.md: hand each to a background "
                       f"subagent (a cheaper model), verify every FINDINGS.md row against the files (`rmfit {self.dir} audit --status KIND N "
                       f"confirmed|rejected|decision`), then answer 'done' ({nrows} rows, {cand} unverified now); 'skip' builds without them.",
                       ("done", "skip"), default="wait", evidence="audits/*/FINDINGS.md")
        if ans is None:
            return "stay", f"subagent audits pending ({nrows} rows, {cand} unverified) -> waiting"
        if ans.lower().startswith("skip"):
            self.note("AUDIT: subagent audits skipped by decision")
            return None
        summ = asub.findings_summary(self.dir)
        cand = [(x["kind"], r["n"]) for x in summ for r in x["unverified"]]
        if cand:
            self.note(f"AUDIT: 'done' answered with {len(cand)} FINDINGS.md row(s) still 'candidate': {cand[:6]}")
            self.ask("audit-subagents", f"'done' was answered but {len(cand)} FINDINGS.md row(s) are still 'candidate' ({cand[:6]}): "
                     f"verify them (`rmfit {self.dir} audit --status KIND N ...`) and answer 'done' again, or 'skip'.",
                     ("done", "skip"), default="wait", evidence="audits/*/FINDINGS.md")
            return "stay", f"{len(cand)} audit row(s) unverified after 'done' -> waiting"
        pending, rows_out = [], []
        for x in summ:
            for r in asub.parse_findings(x["path"]):
                if r["status"] == "decision":
                    a = self.ask(f"audit-decision-{x['kind']}-{r['n']}", f"{x['kind']} audit finding {r['n']} ({r['severity']}) is the evaluator's: "
                                 f"{r['where']}: {r['what']} -- {r['why']}.  Record the treatment in spec.json (standing_decisions / data_sets) and "
                                 f"answer with it.", ("decided",), default="wait", evidence=x["path"])
                    if a is None:
                        pending.append(f"{x['kind']} {r['n']}")
                    else:
                        asub.set_status(x["path"], r["n"], "confirmed", f"evaluator: {a[:80]}")
                        self.note(f"audit {x['kind']} finding {r['n']} decided by the evaluator: {a}")
                elif r["status"] == "confirmed":
                    rows_out.append(f"{x['kind']} {r['n']} {r['severity']} confirmed: {r['where']}: {r['what'][:100]}")
        if pending:
            return "stay", f"{len(pending)} audit finding(s) await the evaluator ({', '.join(pending[:5])}) -> waiting"
        self._finding_rows("audits/*/FINDINGS.md", rows_out[:30])
        self.note(f"AUDIT: subagent audits verified ({nrows} rows: {sum(x['by_status'].get('confirmed', 0) for x in summ)} confirmed, "
                  f"{sum(x['by_status'].get('rejected', 0) for x in summ)} rejected)")
        return None

    def _s_reproduce(self):
        """Mode-1 reproduction before building (skill section 2 item 2): the base model + .sav must give the archived
        chiSquared.out total.  Plain chi2 is compared (the norm penalty's convention differs between the CLI and rmfit)."""
        from .review import output_dir_of, read_chisq_out
        from .engine import SingleSession
        from .model import single_variant
        if self.campaign is not None:
            return "next", "campaign exists: the baseline was verified at init"
        odir = output_dir_of(self.model)
        cs = os.path.join(odir, "chiSquared.out") if odir else None
        if not cs or not os.path.exists(cs):
            return "next", "no archived chiSquared.out next to the model: nothing to reproduce"
        tab = read_chisq_out(cs)
        ref = sum(d["chi2"] for d in tab.values())
        nref = sum(d["n"] for d in tab.values())
        flags = self.spec.flags()
        vdir = os.path.join(self.dir, "variants")
        os.makedirs(vdir, exist_ok=True)
        src = single_variant(self.model, "drive_reproduce", vdir)
        with SingleSession(src, cwd=vdir, flags=flags, name="drive_reproduce") as ev:
            x = self._seed_vector(ev)
            sc = ev.score(x)
        rel = abs(sc.chi2 - ref) / max(abs(ref), 1e-300)
        msg = f"plain chi2 {sc.chi2:.3f} (N {sc.npts}) vs chiSquared.out {ref:.3f} (N {nref}): {100*rel:.3f} %"
        self.state["reproduce"] = dict(chi2=sc.chi2, ref=ref, n=sc.npts, nref=nref, rel=rel)
        if rel > float(self.stop["reproduce_rtol"]) or sc.npts != nref:
            ans = self.ask("reproduce-mismatch", f"the base model does not reproduce its archived chiSquared.out ({msg}); flags "
                           f"{flags}. 'accept' to build on it anyway, 'stop' to fix the input first.", ("accept", "stop"),
                           default="wait", evidence="campaign.log REPRODUCE")
            if ans is None:
                return "stay", msg + " -> waiting"
            if ans.lower().startswith("stop"):
                raise DriveStop("STOP", "reproduction mismatch not accepted: " + msg)
            self.note(f"REPRODUCE mismatch accepted by decision: {msg}")
        else:
            self.note(f"REPRODUCE: the base model reproduces its archived chiSquared.out ({msg})")
        return "next", msg

    def _s_spec_check(self):
        cfg = self.campaign.cfg if self.campaign is not None else None
        fs = self.spec.check(self.model, campaign_cfg=cfg)
        out = os.path.join(self.dir, "SPEC_CHECK.md")
        counts = write_findings_md(out, "Spec check (rmfit drive)", fs, [f"spec `{self.spec_path}` vs model `{self.model}`", ""])
        highs = [f for f in fs if f.severity == "HIGH"]
        fixable = {"spec-excluded-active", "spec-norm-fixed"}        # SCALE_ANCHOR switches excluded sets off / fixes anchor norms
        blocking = [f for f in highs if f.code not in fixable or self.campaign is not None]
        self._finding_rows("SPEC_CHECK.md", [f.line() for f in fs][:30])
        if blocking:
            ans = self.ask("spec-high", f"spec check: {len(blocking)} HIGH disagreement(s) between spec.json and the model: "
                           + "; ".join(f.message[:90] for f in blocking[:3]) + ". Fix one of them and answer 'fixed', or 'accepted'.",
                           ("fixed", "accepted"), default="wait", evidence="SPEC_CHECK.md")
            if ans is None:
                return "stay", f"{len(blocking)} HIGH -> waiting"
            if ans.lower().startswith("fixed"):
                return "stay", "re-checking after 'fixed'"
            self.note(f"SPEC_CHECK HIGH findings accepted by decision ({', '.join(f.code for f in blocking)})")
        return "next", ", ".join(f"{k} {v}" for k, v in counts.items())

    def _s_scale_anchor(self):
        """Apply the spec's data-set roles to a working copy (First fits item 1): anchor -> norm fixed; free -> norm free,
        prior_sigma as the systematic (percent); excluded -> segments off."""
        if self.campaign is not None:
            return "next", "campaign exists: roles were applied before init"
        from .backend import AzrModel
        mdl = AzrModel.from_file(self.model)
        changed, rows = 0, []
        for ds in self.spec.data_sets:
            f = ds.get("file")
            if not f:
                continue
            role = ds.get("role", "free")
            norm = ds.get("norm") or {}
            try:
                if role == "excluded" or role == "excluded_in_window":
                    mdl.set_segment_active(f, False)
                    rows.append(f"{f}: {role} -> segments off")
                elif role in ("anchor", "anchor_geometry") or norm.get("fixed") is True:
                    n = mdl.set_segment_norm(f, vary=False)
                    rows.append(f"{f}: {role} -> norm fixed on {n} segment(s)")
                elif role == "free":
                    sig = norm.get("prior_sigma")
                    n = mdl.set_segment_norm(f, vary=(not norm.get("fixed", False)),
                                             sys_error=(100.0 * float(sig) if sig is not None else None))
                    rows.append(f"{f}: free -> norm {'free' if not norm.get('fixed') else 'fixed'}" + (f", systematic {100*float(sig):g} %" if sig is not None else "") + f" on {n} segment(s)")
                else:
                    continue
                changed += 1
            except KeyError as err:
                rows.append(f"{f}: no matching segment ({err})")
        if self.spec.pre_authorized.get("shifts_free", False) is True:
            # 2026-11-10: SHIFTS_FREE once polished in 7 evaluations without any change because the base varied no
            # shift, so no shift key existed; a shift is a parameter only when its segment line varies it (BoundsPolicy fixes
            # it until SHIFTS_FREE).  Every active segment with a shift uncertainty gets varyShift on here.
            nsh = _vary_shifts_with_error(mdl)
            rows.append(f"energy shifts: varyShift set on {nsh} segment(s) with a shift uncertainty (fixed by the policy until SHIFTS_FREE)")
        out = os.path.join(self.dir, self._stem() + "_drive.azr")
        os.makedirs(os.path.join(self.dir, "output"), exist_ok=True)
        os.makedirs(os.path.join(self.dir, "checks"), exist_ok=True)
        mdl.set_output_dir(os.path.join(self.dir, "output"))
        mdl.write(out)
        sav = self.state.get("sav")
        if sav and os.path.exists(sav):
            shutil.copy2(sav, os.path.splitext(out)[0] + ".sav")
            self.state["sav"] = os.path.splitext(out)[0] + ".sav"
        self.state["model"] = out
        self._finding_rows("SCALE_ANCHOR", rows)
        self.note(f"SCALE_ANCHOR: spec roles applied to {os.path.basename(out)} ({changed} data sets): " + "; ".join(rows))
        return "next", f"{changed} data sets -> {os.path.basename(out)}"

    def _s_seed_channels(self):
        """First fits item 2 for the channels that already exist: every FREE width at exactly zero is lifted to seed_frac x
        gamma_W (in the engine, skill 2026-09-18 / 09-25: never in keV, zero never moves) and the model is baked."""
        if self.campaign is not None:
            return "next", "campaign exists: seeds belong before init"
        from .engine import SingleSession, bake
        from .model import single_variant
        from .local import BoundsPolicy
        from .backend import AzrModel
        frac = float((self.spec.get("model") or {}).get("seed_frac", 0.01))
        eps = float((self.spec.get("bounds") or {}).get("sign_eps", BoundsPolicy().sign_eps))
        flags = self.spec.flags()
        vdir = os.path.join(self.dir, "variants")
        os.makedirs(vdir, exist_ok=True)
        # 2026-11-10 (from-scratch demo): a width of exactly 0 is FIXED whatever its flag (Golden rule 3), so the engine never
        # lists it and the search over ev.keys found "nothing to seed" on a from-scratch model whose 20 background channels were
        # all free zeros.  Every active channel flagged free with gamma == 0 first gets a placeholder (a tiny physical width, so
        # that the key exists), then the seed in units of gamma_W as before; photon channels, which have no Wigner limit, keep
        # the placeholder PHOTON_SEED_EV as their start.
        PHOTON_SEED_EV = float((self.spec.get("model") or {}).get("photon_seed_ev", 1.0))
        mdl0 = AzrModel.from_file(self.model)
        placed = []
        for lv in mdl0.active_levels:
            for ch in lv.channels:
                if not bool(ch.channel_fixed) and float(ch.gamma) == 0.0:
                    ch.gamma = PHOTON_SEED_EV if ch.is_photon else 1e-3
                    placed.append(f"{lv.jpi}@{float(lv.energy):.4f} p{ch.pair} L{ch.L} S{ch.S}" + (" (photon)" if ch.is_photon else ""))
        seed_src = self.model
        if placed:
            # next to the base, not in variants/: the model's relative data paths resolve from its own directory, and
            # single_variant links the data directories of the model it is given (demo 3: "Could Not Fill Segment #1")
            seed_src = os.path.join(os.path.dirname(self.model), self._stem() + "_placeholders.azr")
            mdl0.set_output_dir(os.path.join(self.dir, "output"))
            mdl0.write(seed_src)
            self.log(f"  {len(placed)} free channels at exactly zero given a placeholder so the engine lists them: " + ", ".join(placed[:12]))
        src = single_variant(seed_src, "drive_seed", vdir)
        with SingleSession(src, cwd=vdir, flags=flags, name="drive_seed") as ev:
            x = self._seed_vector(ev)
            gW = ev.wigner_amplitudes()
            zero = [j for j, k in enumerate(ev.keys) if k.kind == "width" and gW.get(j, 0) > 0 and abs(x[j]) < eps * gW[j]]
            if not zero and not placed:
                return "next", "no free width at zero: nothing to seed"
            sc0 = ev.score(x)
            x1 = x.copy()
            for j in zero:
                x1[j] = frac * gW[j]
            sc1 = ev.score(x1)
            if sc0.objective > 0 and sc1.objective > 1.2 * sc0.objective:
                self.log(f"seed {frac:g} gamma_W raises the objective {sc0.objective:.1f} -> {sc1.objective:.1f}: halving the seed")
                frac *= 0.5
                for j in zero:
                    x1[j] = frac * gW[j]
                sc1 = ev.score(x1)
            keys = list(ev.keys)
            labels = [str(keys[j]) for j in zero]
        out = os.path.join(self.dir, self._stem() + "_seeded.azr")
        azr, sav, scb = bake(seed_src, x1, keys, out, flags, workdir=vdir, output_dir=os.path.join(self.dir, "output"),
                             name="drive_seed_bake", separate=False)
        # non-seeded physical widths must survive the bake (seed_stage1.py width_check: a failed Brune transform at load
        # rescales a whole level)
        A, B = AzrModel.from_file(self.model), AzrModel.from_file(azr)      # the pre-placeholder model: placeholders are skipped as zeros
        bad = []
        for la, lb in zip(A.active_levels, B.active_levels):
            for ca, cb in zip(la.channels, lb.channels):
                if abs(ca.gamma) <= 1e-6:
                    continue
                if abs(cb.gamma / ca.gamma - 1) > 0.01:
                    bad.append(f"{la.jpi} {la.energy:.4f} pair {ca.pair} L{ca.L}: {ca.gamma:.5g} -> {cb.gamma:.5g}")
        if bad:
            raise DriveStop("ERROR", f"seed bake changed {len(bad)} non-seeded widths by > 1 %: " + "; ".join(bad[:4]))
        self.state["model"], self.state["sav"] = azr, sav
        self.state["seeds"] = dict(frac=frac, n=len(zero), before=sc0.objective, after=sc1.objective, baked=scb.objective)
        self._finding_rows("SEED_CHANNELS", [f"{len(zero)} free widths at zero seeded at {frac:g} gamma_W: objective {sc0.objective:.1f} -> "
                                             f"{sc1.objective:.1f} (baked {scb.objective:.1f}); {len(placed)} placeholders"] + labels[:40])
        self.note(f"SEED_CHANNELS: {len(zero)} free widths at exactly zero lifted to {frac:g} gamma_W and baked into "
                  f"{os.path.basename(azr)} (objective {sc0.objective:.1f} -> {sc1.objective:.1f}; non-seeded widths unchanged)")
        return "next", f"{len(zero)} seeds at {frac:g} gamma_W ({len(placed)} channels given a placeholder) -> {os.path.basename(azr)}"

    def _s_pair_completeness(self):
        """First fits item 3 as a check: excited states (ENSDF) of the model's residual nuclei whose partition opens below the
        top of the fit window and has no pair.  Reported and queued (default 'continue'); the pair builder is the
        supervisor's (reaction-specific radii, channel spins, summed segments)."""
        from .backend import AzrModel
        top = self._ex_top()
        if top is None:
            return "next", "fit-window top unknown: skipped"
        mdl = AzrModel.from_file(self.model)
        pairs = mdl.pair_lines()
        have = {}          # (z_light, a_light, z_heavy, a_heavy) -> {e2 values}, sepE
        rows, missing = [], []
        for key, c in pairs.items():
            if bool(getattr(c, "is_photon", False)):
                continue
            try:
                if PYAZR3:
                    z1, m1, z2, m2 = c.Z1, c.M1, c.Z2, c.M2
                    sep, e2 = c.sep_energy, c.heavy_excitation
                else:
                    z1, m1, z2, m2 = c._get("z1", int), c._get("m1", float), c._get("z2", int), c._get("m2", float)
                    sep, e2 = c._get("sepE", float), c._get("e2", float)
            except Exception:  # noqa: BLE001
                continue
            part = (z1, int(round(m1)), z2, int(round(m2)))
            d = have.setdefault(part, dict(e2=set(), sep=sep, keys=[]))
            d["e2"].add(round(e2, 3))
            d["keys"].append(key)
            rows.append(f"pair {key}: {_nuc(z1, m1)} + {_nuc(z2, m2)}(E* {e2:.4f}) threshold Ex {sep + e2:.4f}")
        try:
            from pyazr import fetch_levels                      # pyazr 2 and 3 alike
        except Exception as err:  # noqa: BLE001
            self._finding_rows("PAIR_COMPLETENESS", rows + [f"ENSDF client unavailable ({err})"])
            return "next", "pyazr.nds unavailable: pairs listed only"
        for part, d in have.items():
            z, a = part[2], part[3]
            try:
                levels = fetch_levels(f"{a}{SYMBOLS[z].lower()}" if z < len(SYMBOLS) else f"{a}z{z}")
            except Exception as err:  # noqa: BLE001
                rows.append(f"{_nuc(z, a)}: ENSDF fetch failed ({err})")
                continue
            for lv in levels:
                ex_state = float(lv.energy_mev)
                thr = d["sep"] + ex_state
                if ex_state <= 0.0 or thr >= top:
                    continue
                if any(abs(ex_state - e) < 0.005 for e in d["e2"]):
                    continue
                missing.append(f"{_nuc(part[0], part[1])} + {_nuc(z, a)}({ex_state:.4f} {lv.jp}) opens at Ex {thr:.4f} < {top:.3f}: no pair")
        self._finding_rows("PAIR_COMPLETENESS", rows + (missing or ["every ENSDF state of the residual nuclei below the window top has a pair"]))
        self.state["pairs_missing"] = missing
        if missing:
            self.note(f"PAIR_COMPLETENESS: {len(missing)} partition(s) open below Ex {top:.3f} MeV without a pair: " + "; ".join(missing[:6])
                      + ". First-fits rule 3: add a pair for each before fitting (builder pattern an archived campaign); "
                        "a window cap instead is a stop-gap to log.")
            self.ask("missing-pairs", f"{len(missing)} partition(s) open below the window top have no pair (DRIVE.md). Add them (the drive "
                     f"re-runs from SCALE_ANCHOR on the new base) or 'continue' without.", ("added", "continue"), default="continue",
                     evidence="DRIVE.md PAIR_COMPLETENESS")
        return "next", f"{len(pairs)} pairs, {len(missing)} missing partitions below Ex {top:.3f}"

    def _s_ladder(self):
        """First fits item 4 when the spec asks for it (model.ladder_Ex = [cutoffs]): data up to each cutoff, levels above
        fixed, polished and carried forward by key; the final vector is baked as the base for INIT."""
        if self.campaign is not None:
            return "next", "campaign exists"
        cuts = (self.spec.get("model") or {}).get("ladder_Ex")
        if not cuts:
            return "next", "no model.ladder_Ex in the spec: skipped"
        from .windows import ladder_variants
        from .engine import Sharded, SingleSession, bake
        from .local import polish, default_bounds, BoundsPolicy
        from .keys import take
        from .model import single_variant
        from .campaign import _dc_from
        flags = self.spec.flags()
        policy = _dc_from(BoundsPolicy, self.spec.get("bounds"))
        nshards = int(self.spec.get("nshards") or 8)
        nfev = int(self.budget["ladder_nfev"])
        vdir = os.path.join(self.dir, "variants", "ladder")
        steps = ladder_variants(self.model, vdir, [float(c) for c in cuts])
        src = single_variant(self.model, "drive_ladder_base", os.path.join(self.dir, "variants"))
        with SingleSession(src, cwd=os.path.join(self.dir, "variants"), flags=flags, name="drive_ladder_base") as ev0:
            x_full, keys_full = self._seed_vector(ev0), list(ev0.keys)
            obj0 = ev0.score(x_full).objective
        rows = []
        for i, info in enumerate(steps):
            path = os.path.join(vdir, f"step{i+1}_Ex{info['cutoff']:g}.azr")
            ev = Sharded(path, os.path.join(vdir, f"shards_step{i+1}"), nshards=nshards, flags=flags)
            try:
                xs = take(x_full, keys_full, ev.keys, ev.x0)
                s0 = ev.score(xs).objective
                fr = polish(ev, xs, bounds=default_bounds(ev, xs, policy), policy=policy, max_nfev=nfev, log=self.log)
                x_full = take(fr.x, ev.keys, keys_full, x_full)
                rows.append(f"step {i+1} Ex <= {info['cutoff']:g}: {s0:.1f} -> {fr.objective:.1f} ({fr.status}, {fr.nfev} evaluations, {len(ev.keys)} free)")
                self.log(rows[-1])
            finally:
                ev.close()
        out = os.path.join(self.dir, self._stem() + "_ladder.azr")
        azr, sav, sc = bake(self.model, x_full, keys_full, out, flags, workdir=os.path.join(self.dir, "variants"),
                            output_dir=os.path.join(self.dir, "output"), name="drive_ladder_bake")
        rows.append(f"full model: {obj0:.1f} -> {sc.objective:.1f} baked into {os.path.basename(azr)}")
        self.state["model"], self.state["sav"] = azr, sav
        self._finding_rows("LADDER", rows)
        self.note("LADDER: " + "; ".join(rows))
        return "next", rows[-1]

    def _s_candidates(self):
        """candidates.json: kept if present (every entry needs a hypothesis); else built from the ENSDF adopted levels of
        spec.compound_nucleus inside the fit window that have no active level within 30 keV, broad first."""
        path = os.path.join(self.dir, "candidates.json")
        if os.path.exists(path):
            bank = json.load(open(path)).get("candidates", [])
            lacking = [i for i, c in enumerate(bank) if not (c.get("hypothesis") or c.get("note"))]
            self._finding_rows("CANDIDATES", [f"{len(bank)} candidates in candidates.json; {len(lacking)} without a hypothesis (skipped by STRUCTURE_ADD): {lacking}"])
            return "next", f"{len(bank)} candidates kept ({len(lacking)} lack a hypothesis)"
        cn = self.spec.get("compound_nucleus")
        if not cn:
            json.dump(dict(candidates=[], note="no compound_nucleus in the spec: empty bank"), open(path, "w"), indent=1)
            return "next", "no compound_nucleus in the spec: empty candidates.json"
        from .backend import AzrModel
        top = self._ex_top()
        fw = (self.spec.get("model") or {}).get("fit_window") or {}
        lo = float(fw["Ex_MeV"][0]) if fw.get("Ex_MeV") else 0.0
        mdl = AzrModel.from_file(self.model)
        have = [(lv.jpi, float(lv.energy)) for lv in mdl.active_levels]
        try:
            from pyazr import fetch_levels                      # pyazr 2 and 3 alike
            levels = fetch_levels(cn.lower())
        except Exception as err:  # noqa: BLE001
            json.dump(dict(candidates=[], note=f"ENSDF fetch failed: {err}"), open(path, "w"), indent=1)
            return "next", f"ENSDF fetch failed ({err}): empty candidates.json"
        cands = []
        for lv in levels:
            ex = float(lv.energy_mev)
            if ex < lo or (top is not None and ex > top):
                continue
            if any(abs(ex - e) < 0.030 for _, e in have):
                continue
            jp = (lv.jp or "").strip()
            # tentative assignments as firm ones: "(0+)", "4(-)", "(3)-" (the hypothesis keeps ENSDF's text)
            # (spaces are kept: pyazr's NNDC source spells a J list "(1 2)+", which must not read as 12+)
            clean = re.fullmatch(r"(\d+(/2)?)([+-])", jp.replace("(", "").replace(")", "").strip())
            jpi = f"{clean.group(1)}{clean.group(3)}" if clean else None
            w = _ensdf_width_kev(lv)
            cands.append(dict(energy=round(ex, 4), jpi=jpi, energy_fixed=True, channels=[], source="ensdf",
                              width_kev=w, hypothesis=f"ENSDF adopted level {ex:.4f} MeV {jp or '(J unknown)'}"
                              + (f", Gamma {w:g} keV" if w else "") + " has no level in the model within 30 keV"))
        cands.sort(key=lambda c: (-(c["width_kev"] or 0.0), c["energy"]))     # broad first (skill: Gamma > 10 keV first)
        nmax = int(self.budget["level_tests_max"])
        json.dump(dict(candidates=cands[:nmax], dropped_by_budget=len(cands) - min(len(cands), nmax),
                       note=f"ENSDF {cn} levels in Ex {lo:g}-{top if top is not None else 'top'} MeV without a model level within 30 keV; "
                            f"broad first; written by rmfit drive {_now()}"), open(path, "w"), indent=1)
        self._finding_rows("CANDIDATES", [f"{c['energy']:.4f} {c['jpi'] or 'scan'} {c['hypothesis']}" for c in cands[:nmax]])
        self.note(f"CANDIDATES: candidates.json written from ENSDF {cn}: {min(len(cands), nmax)} levels "
                  f"({len(cands) - min(len(cands), nmax)} more dropped by budget.level_tests_max {nmax}), broad first")
        return "next", f"{min(len(cands), nmax)} ENSDF candidates"

    def _s_init(self):
        if self.campaign is not None:
            return "next", "campaign exists"
        from .campaign import Campaign
        c = Campaign.init(self.dir, self.model, sav=self.state.get("sav"), spec=self.spec_path,
                          note=f"rmfit drive on {os.path.basename(self.model)}",
                          sav_fixed=str((self.spec.get("engine") or {}).get("sav_fixed", "ignore")) if self.spec else "ignore")
        c.write_job_script()
        self._campaign = c
        cid, x, obj, keys = c.incumbent()
        self.state["cycle_start_objective"] = obj
        self.state["baseline_objective"] = obj
        self._job_script(overwrite=True)
        return "next", f"campaign initialised: baseline {obj:.3f} (candidate {cid})"

    def _s_recall(self):
        """Validate the level-add procedure on levels the model already has before STRUCTURE_ADD is trusted (the evaluator
        2026-10-15; an archived campaign): the structure.recall_levels broadest physical levels are removed one at a time and
        searched for with no J^pi given.  Below structure.recall_min_pass the adds of this run are queued as decisions."""
        n = int(self.structure.get("recall_levels") or 0)
        if n <= 0:
            return "next", "structure.recall_levels 0: skipped"
        if int(self.budget["level_tests_max"]) <= 0:
            return "next", "budget.level_tests_max 0: no level adds to validate"
        c = self.campaign
        res = c.run_recall(n=n, min_width_kev=float(self.structure.get("broad_first_kev") or 0.0),
                           global_nfev=int(self.budget["level_nfev"]), pinned_nfev=int(self.budget["pinned_nfev"]),
                           local_release=float(self.structure.get("local_release_mev") or 0.0),
                           local_nfev=int(self.structure.get("local_nfev") or 80), note="rmfit drive RECALL")
        rows = [f"{r['level']} {r['energy']:.4f} MeV Gamma {r['gamma_kev']:.1f} keV -> {r.get('outcome')} as {r.get('jpi_found')}"
                f" @ {r.get('energy_found', float('nan')):.4f}: "
                + ("re-found" if r["passed"] else
                   ("LOCATED but rejected by the criteria" if r.get("located") else "NOT re-found") + f" ({r['why']})")
                for r in res["results"]]
        self._finding_rows("RECALL", rows or ["no physical level broad enough inside the window: nothing removed"])
        self.state["recall"] = dict(n=res["n"], n_pass=res["n_pass"], n_located=res.get("n_located", res["n_pass"]), fraction=res["fraction"])
        if res["n"] == 0:
            return "next", "no known level broad enough below bg_energy to remove: recall skipped"
        minp = float(self.structure.get("recall_min_pass") or 0.0)
        if res["fraction"] < minp:
            self.state["adds_untrusted"] = True
            nloc = int(res.get("n_located", res["n_pass"]))
            cause = (f"the search located {nloc} of {res['n']} (right J^pi, energy, width) but the criteria rejected them: "
                     f"see the RECALL table (a theta2 cap the incumbent's own level sits at is a spec bounds setting, not a search failure)"
                     if nloc > res["n_pass"] else "the search did not locate the removed level(s)")
            did = self.decisions.add("recall-low", f"the level-add procedure re-found {res['n_pass']} of {res['n']} known levels "
                                     f"(below structure.recall_min_pass {minp:g}): {cause}; its adds are queued as accept-pending "
                                     f"decisions instead of adopted.  Answer 'trust' to pre-authorize them anyway.",
                                     evidence="DRIVE.md RECALL table; recall/<level>/ child campaigns", options=("adopt-manually", "trust"),
                                     default="adopt-manually")
            self.note(f"RECALL: {res['n_pass']} of {res['n']} known levels re-found, {nloc} located (< {minp:g}): level adds untrusted, "
                      f"decision {did} (non-blocking)")
            return "next", (f"recall {res['n_pass']}/{res['n']} below {minp:g} ({nloc} located, rejected by the criteria): "
                            f"level adds become decisions until '{did}' is answered 'trust'")
        return "next", f"recall {res['n_pass']}/{res['n']} known levels re-found: level adds trusted"

    # ==================================================================================================
    # loop phase
    # ==================================================================================================
    def _incumbent_obj(self) -> float:
        return float(self.campaign.incumbent()[2])

    def _s_tighten(self):
        c = self.campaign
        if c is None:
            raise DriveStop("ERROR", "no campaign.json: the build phase did not reach INIT")
        n = int(self.counters.get("polish_rounds", 0))
        if self.state.get("cycle_start_objective") is None:
            # a campaign that existed before the drive (perturb sibling, retrofit) never passed INIT: the cycle's
            # start is its incumbent when the loop phase begins (a benchmark run's CYCLE saw 0 % and stopped)
            self.state["cycle_start_objective"] = self._incumbent_obj()
            self.state.setdefault("baseline_objective", self.state["cycle_start_objective"])
            self.save()
        if n >= int(self.budget["polish_rounds_max"]):
            self.counters["polish_rounds"] = 0
            self.counters["polish_converged"] = 0
            return "next", f"budget.polish_rounds_max {self.budget['polish_rounds_max']} reached before convergence (last gain {100*float(self.counters.get('polish_last_gain', 0)):.3f} %)"
        before = self._incumbent_obj()
        fr = c.run_polish(max_nfev=int(self.budget["polish_nfev"]), tighten=(n == 0), note="rmfit drive TIGHTEN")
        after = self._incumbent_obj()
        gain = (before - after) / max(abs(before), 1e-300)
        self.counters["polish_rounds"] = n + 1
        self.counters["polish_last_gain"] = gain
        msg = f"round {n+1}: {before:.3f} -> {after:.3f} ({100*gain:.3f} %, {fr.status}, {fr.nfev} evaluations)"
        if gain < float(self.stop["polish_gain_rel_min"]) or fr.status == "crashed":
            self.counters["polish_rounds"] = 0
            self.counters["polish_converged"] = 1 if fr.status != "crashed" else 0
            return "next", msg + " -> converged" if fr.status != "crashed" else msg
        return "stay", msg

    def _polish_gate(self, what):
        """stop_rules.sign_after_converged_polish: SIGN and POPULATION wait for a converged TIGHTEN.  Returns the deferral note,
        or None when the stage runs: the polish converged, the gate is off, this is the last cycle (budget.cycles_max -- there is
        no next cycle to defer to), or the stage was already deferred once in this campaign (a polish that did not converge in
        polish_rounds_max rounds twice will not in a third: a from-scratch demo with 162 free parameters, 6 x 300 evaluations
        still gaining 10 % per round, every cycle deferred both stages and the drive would have ended with no sign round; the
        13N end-to-end test with its seeded channels and cycles_max 1 never reached POPULATION).  The stage then runs on the
        unconverged incumbent and the readme says so: the cure is a larger budget.polish_nfev / polish_rounds_max."""
        if not int(self.stop.get("sign_after_converged_polish", 1)) or self.counters.get("polish_converged", 1) != 0:
            return None
        gain = 100 * float(self.counters.get("polish_last_gain", 0))
        last = int(self.state["cycle"]) >= int(self.budget["cycles_max"])
        nd = int(self.counters.get(f"{what}_deferrals", 0))
        if last or nd >= 1:
            why = "the last cycle (budget.cycles_max)" if last else f"already deferred {nd} cycle(s)"
            self.note(f"{what.upper()} runs on an unconverged polish (last gain {gain:.3f} %): {why}; budget.polish_nfev "
                      f"{self.budget['polish_nfev']} x polish_rounds_max {self.budget['polish_rounds_max']} is too small for this model")
            return None
        self.counters[f"{what}_deferrals"] = nd + 1
        return (f"polish not converged (last gain {gain:.3f} %): {what} deferred to the next cycle "
                f"(stop_rules.sign_after_converged_polish)")

    def _sign_state(self, label):
        c = self.campaign
        if not self.counters.get(f"{label}_rounds"):
            gate = self._polish_gate(label)
            if gate:
                return "next", gate
        n = int(self.counters.get(f"{label}_rounds", 0))
        noacc = int(self.counters.get(f"{label}_noaccept", 0))
        if n >= int(self.budget["sign_rounds_max"]):
            self.counters[f"{label}_rounds"] = self.counters[f"{label}_noaccept"] = 0
            return "next", f"budget.sign_rounds_max {self.budget['sign_rounds_max']} reached"
        st = c.run_sign_round(note=f"rmfit drive {label}")
        acc = int(st.get("accepted", 0) or 0)
        n += 1
        noacc = 0 if acc else noacc + 1
        self.counters[f"{label}_rounds"], self.counters[f"{label}_noaccept"] = n, noacc
        msg = (f"round {n}: {st.get('proposed')} proposals, {st.get('improving_screen')} improving on screen, {acc} accepted, "
               f"{st.get('before', 0):.3f} -> {st.get('verified', st.get('after', 0)):.3f}")
        if n >= 2 and noacc >= int(self.stop["sign_rounds_without_accept"]):
            self.counters[f"{label}_rounds"] = self.counters[f"{label}_noaccept"] = 0
            self.counters["adopted_since_sign"] = 0
            return "next", msg + f" -> {noacc} rounds without an acceptance: converged for this structure"
        return "stay", msg

    def _s_sign(self):
        return self._sign_state("sign")

    def _s_sign_post(self):
        if not self.counters.get("adopted_since_sign"):
            return "next", "no structure change since the last sign rounds: skipped"
        return self._sign_state("sign_post")

    def _s_population(self):
        """Restarts + crossover (Campaign.run_population), the basin-crossing stage.  In a benchmark the drive's
        sign rounds left a basin 13.5 % above the reference (no single/pair/channel/group flip pays, even with a long
        polish), six polished restarts + one crossover pass came within 0.7 % of it.  Runs once per cycle after the sign rounds, on a
        converged polish only, until budget.population_rounds_max or a round that gains < stop_rules.population_gain_rel_min
        (then skipped until a structure change).  An improvement is a new basin: SIGN_POST finishes its last signs."""
        c = self.campaign
        k = int(self.budget.get("restarts_per_cycle", 0) or 0)
        if k <= 0:
            return "next", "budget.restarts_per_cycle 0: skipped"
        gate = self._polish_gate("population")
        if gate:
            return "next", gate
        n = int(self.counters.get("population_rounds", 0))
        if n >= int(self.budget.get("population_rounds_max", 3)):
            return "next", f"budget.population_rounds_max {self.budget.get('population_rounds_max', 3)} reached"
        if self.counters.get("population_exhausted"):
            return "next", "the last population round gained < stop_rules.population_gain_rel_min: skipped until a structure change"
        before = self._incumbent_obj()
        nfev = int(self.budget.get("restart_nfev") or 0) or int(self.budget["polish_nfev"])
        st = c.run_population(restarts=k, nfev=nfev, passes=int(self.budget.get("crossover_passes", 2)),
                              note="rmfit drive POPULATION")
        after = self._incumbent_obj()
        gain = (before - after) / max(abs(before), 1e-300)
        self.counters["population_rounds"] = n + 1
        self.counters["population_last_gain"] = gain
        msg = (f"round {n+1}: {k} restarts, {st.get('crossover_moves', 0)} crossover patterns, {st.get('polishes', 0)} polishes of "
               f"{nfev}; {before:.3f} -> {after:.3f} ({100*gain:.3f} %) via {st.get('best')}")
        if gain >= float(self.stop.get("population_gain_rel_min", 1e-3)):
            self.counters["adopted_since_sign"] = 1        # a new basin: SIGN_POST finishes the last signs
            self.counters["polish_converged"] = 1
            return "next", msg + " -> new basin adopted"
        self.counters["population_exhausted"] = 1
        return "next", msg + " -> below population_gain_rel_min: no further population rounds for this structure"

    def _s_shifts_free(self):
        c = self.campaign
        if not self.spec.pre_authorized.get("shifts_free", False) is True:
            return "next", "pre_authorized.shifts_free is not true: skipped"
        before = self._incumbent_obj()
        fr = c.run_polish(max_nfev=int(self.budget["polish_nfev"]), shifts_free=True, note="rmfit drive SHIFTS_FREE")
        after = self._incumbent_obj()
        return "next", f"{before:.3f} -> {after:.3f} ({fr.status}, {fr.nfev} evaluations)"

    def _bank(self):
        path = os.path.join(self.dir, "candidates.json")
        data = json.load(open(path)) if os.path.exists(path) else dict(candidates=[])
        return path, data, data.get("candidates", [])

    def _is_narrow(self, cand) -> bool:
        """Known width below structure.broad_first_kev (unknown widths count as broad: tested in the broad phase)."""
        bf = float(self.structure.get("broad_first_kev") or 0.0)
        w = _cand_width_kev_known(cand)
        return bf > 0 and w is not None and w < bf

    def _s_levelscan(self):
        """Residual-driven candidates when the bank's broad queue is exhausted (the evaluator 2026-09-14/10-12: one level of every
        J^pi at the width the residual shows, compared inside the model; a broad-level scan found the campaign's
        largest gain).  Campaign.run_levelscan; the best (J^pi, pattern) per residual cluster is appended to candidates.json
        with its screen as hypothesis and the next two J^pi as alternatives; STRUCTURE_ADD tests it like any candidate."""
        per = int(self.structure.get("levelscan_per_cycle") or 0)
        if per <= 0:
            return "next", "structure.levelscan_per_cycle 0: skipped"
        if int(self.counters.get("level_tests", 0)) >= int(self.budget["level_tests_max"]):
            return "next", "no level-test budget left this cycle: no scan"
        path, data, bank = self._bank()
        tested = set(int(i) for i in self.state["tested_candidates"])
        broad = [i for i in range(len(bank)) if i not in tested and not self._is_narrow(bank[i])
                 and (bank[i].get("hypothesis") or bank[i].get("note"))]
        if broad:
            return "next", f"{len(broad)} broad candidate(s) untested: no scan needed yet"
        done = int(self.counters.get("levelscan_rounds", 0))
        if done >= per:
            return "next", f"structure.levelscan_per_cycle {per} reached this cycle"
        c = self.campaign
        top = self._ex_top()
        fw = (self.spec.get("model") or {}).get("fit_window") or {}
        lo = float(fw["Ex_MeV"][0]) if fw.get("Ex_MeV") else None
        from .levelscan import to_level_candidate
        res = c.run_levelscan(top=int(self.structure.get("levelscan_top") or 3), gamma_kev=self.structure.get("levelscan_gamma_kev"),
                              ex_range=((lo if lo is not None else 0.0, top) if top is not None else None), refit=0,
                              note="rmfit drive LEVELSCAN")
        added = []
        for r in res.get("results", []):
            bg = [row for row in (r.get("best_global") or []) if row.get("gain", 0) > 0]
            if not bg:
                continue
            row = bg[0]
            ex = float(row.get("ex", r["ex"]))
            if any(abs(ex - float(b.get("energy", 1e9))) < 0.030 and b.get("jpi") in (None, row["jpi"]) for b in bank):
                continue
            lc = to_level_candidate(row, ex)
            bank.append(dict(energy=lc.energy, jpi=lc.jpi, jpis=[q["jpi"] for q in bg[:3]], channels=lc.channels, energy_fixed=True,
                             source="levelscan", width_kev=float(row["gamma_total_ev"]) / 1e3, cycle=int(self.state["cycle"]),
                             hypothesis=f"residual-run cluster at Ex {r['ex']:.3f} MeV (chi2 {float(r.get('score') or 0):.0f} in runs, "
                                        f"run width {float(r.get('run_width') or 0) * 1e3:.0f} keV); template screen of every J^pi at Gamma "
                                        f"{r['gamma_kev']:g} keV: best {row['jpi']} {row['desc']} (screen gain {row['gain']:+.0f})"
                                        + (", then " + ", ".join(f"{q['jpi']} {q['gain']:+.0f}" for q in bg[1:3]) if len(bg) > 1 else "")))
            added.append(f"{lc.jpi} @ {lc.energy:.3f}")
        data["candidates"] = bank
        json.dump(data, open(path, "w"), indent=1)
        self.counters["levelscan_rounds"] = done + 1
        self.counters["levelscan_last_added"] = len(added)
        msg = f"{len(res.get('results', []))} residual cluster(s) screened, {len(added)} candidate(s) appended: " + ", ".join(added)
        self.note("LEVELSCAN: " + msg)
        return "next", msg

    def _s_structure_add(self):
        c = self.campaign
        from .structure import LevelCandidate
        path, data, bank = self._bank()
        tested = set(int(i) for i in self.state["tested_candidates"])
        marg = int(self.counters.get("marginal_in_a_row", 0))
        ntests = int(self.counters.get("level_tests", 0))
        if ntests >= int(self.budget["level_tests_max"]):
            return "next", f"budget.level_tests_max {self.budget['level_tests_max']} reached ({len(bank) - len(tested)} untested)"
        if marg >= int(self.stop["level_adds_marginal_in_a_row"]):
            self.counters["marginal_in_a_row"] = 0
            return "next", f"{marg} marginal adds in a row: the queue ends here ({len(bank) - len(tested)} untested)"
        # a candidate without a hypothesis is never tested (every structure candidate needs one)
        for i in range(len(bank)):
            if i not in tested and not (bank[i].get("hypothesis") or bank[i].get("note")):
                self.log(f"candidate {i} ({bank[i].get('energy')} {bank[i].get('jpi')}) has no hypothesis: skipped (every structure candidate needs one)")
                self.state["tested_candidates"].append(i)
                tested.add(i)
        untested = [i for i in range(len(bank)) if i not in tested]
        # broad first (the evaluator 2026-09-14: Gamma > 10 keV first, narrow alignment later; DESIGN.md 2.2): known widths below
        # structure.broad_first_kev wait until the broad queue AND the residual scan are exhausted (structure.narrow_phase)
        broad = [i for i in untested if not self._is_narrow(bank[i])]
        narrow = [i for i in untested if self._is_narrow(bank[i])]
        if broad:
            order = sorted(broad, key=lambda i: (-_cand_width_kev(bank[i]), i))
        else:
            if not self.state.get("broad_done"):
                per = int(self.structure.get("levelscan_per_cycle") or 0)
                done = int(self.counters.get("levelscan_rounds", 0))
                if per > 0 and done < per and (done == 0 or int(self.counters.get("levelscan_last_added", 0)) > 0):
                    return "LEVELSCAN", "broad queue empty: residual-driven scan before the narrow phase"
                self.state["broad_done"] = True
                if narrow:
                    self.note(f"STRUCTURE_ADD: broad structure queue exhausted; {len(narrow)} narrow candidate(s) (Gamma < "
                              f"{self.structure.get('broad_first_kev')} keV) " + ("follow (narrow alignment)" if self.structure.get("narrow_phase", True) else "stay deferred (structure.narrow_phase false)"))
            if not narrow:
                return "next", f"every candidate tested ({len(bank)})"
            if not self.structure.get("narrow_phase", True):
                return "next", f"broad queue done; {len(narrow)} narrow candidate(s) deferred (structure.narrow_phase false)"
            order = sorted(narrow, key=lambda i: (-_cand_width_kev(bank[i]), i))
        i = order[0]
        spec = bank[i]
        cand = LevelCandidate(energy=float(spec["energy"]), jpi=spec.get("jpi"), channels=spec.get("channels", []),
                              source=spec.get("source", "drive"), note=(spec.get("hypothesis") or spec.get("note") or ""),
                              energy_fixed=spec.get("energy_fixed", True))
        jpis = list(spec.get("jpis") or ([cand.jpi] if cand.jpi else [])) or None
        self.log(f"candidate {i}: {cand.jpi or 'scan'} @ {cand.energy} MeV -- hypothesis: {cand.note}"
                 + (f" (J^pi tested: {', '.join(jpis)})" if jpis and len(jpis) > 1 else ""))
        if self.state.get("adds_untrusted"):
            for d in self.decisions.all():
                if d["slug"] == "recall-low" and not d.pending and (d.get("answer") or "").strip().lower().startswith("trust"):
                    self.state["adds_untrusted"] = False
                    self.note(f"decision {d.id} answered 'trust': level adds pre-authorized again")
        if self.state.get("adds_untrusted") and getattr(c, "spec", None) is not None:
            c.spec.data.setdefault("pre_authorized", {})["adopt_level"] = False      # this process only: adds become decisions
            self.log("  level adds untrusted after RECALL: an accept is recorded as accept-pending (decision)")
        st = c.run_level_add(cand, jpis=jpis, global_nfev=int(self.budget["level_nfev"]), pinned_nfev=int(self.budget["pinned_nfev"]),
                             local_release=float(self.structure.get("local_release_mev") or 0.0),
                             local_nfev=int(self.structure.get("local_nfev") or 80), note="rmfit drive",
                             # the J^pi group's other open channels FREE (tiny seed), not gamma 0 FIXED: an ENSDF candidate
                             # names no channels (structure.seed_omitted, default True since 2026-11-22)
                             seed_omitted=bool(self.structure.get("seed_omitted", True)))
        self.state["tested_candidates"].append(i)
        self.counters["level_tests"] = ntests + 1
        out = st.get("outcome", "?")
        if out == "accept":
            self.counters["marginal_in_a_row"] = 0
            self.counters["adopted_since_sign"] = 1
            self.counters["population_exhausted"] = 0          # a new structure: the population stage may pay again
            self.state["broad_done"] = False                   # new residuals: the broad queue / scan may have more
        elif out == "accept-pending":
            self.counters["marginal_in_a_row"] = 0
        elif out == "scan_failed":
            # not a measured marginal level: nothing could be tested (e.g. no channel open at that energy); it must not
            # end the queue for the candidates behind it
            self.log(f"  candidate {i}: untestable ({out}); not counted toward level_adds_marginal_in_a_row")
        else:
            self.counters["marginal_in_a_row"] = marg + 1
        return "stay", f"candidate {i} {st.get('jpi', cand.jpi) or 'scan'} @ {cand.energy}: {out} (delta {st.get('delta', float('nan')):+.1f})"

    def _s_dataset_off(self):
        """Step 5: the Botek test on the worst chi2/N sets; the result is a table and (if the extra recovery is large) a
        decision with default 'keep' -- whether a data set stays is the evaluator's."""
        c = self.campaign
        import statistics
        cid = c.incumbent()[0]
        rows = c.dataset_table(cid)
        if len(rows) < 2:
            return "next", "fewer than two data sets"
        med = statistics.median(r["per_pt"] for r in rows)
        done = set(self.state["dataset_off_tested"])
        worst = [r for r in rows if r["per_pt"] > max(3.0 * med, 3.0) and r["name"] not in done][: int(self.budget["dataset_off_max"])]
        if not worst or len(done) >= int(self.budget["dataset_off_max"]):
            return "next", f"{len(done)} sets tested this cycle; median chi2/N {med:.2f}"
        r = worst[0]
        st = c.run_dataset_off(r["name"], nfev=int(self.budget["dataset_off_nfev"]), note="rmfit drive")
        self.state["dataset_off_tested"].append(r["name"])
        pct = float(st.get("extra_pct", 0.0))
        msg = f"{r['name']} (chi2/N {r['per_pt']:.1f} vs median {med:.1f}): extra recovery {st.get('extra_recovery', 0):+.0f} ({pct:.2f} %)"
        slug = f"dataset-off-{r['name']}"
        if pct > float(self.stop["dataset_off_extra_pct"]) and any(d.get("slug") == slug for d in self.decisions.all()):
            msg += " (decision already queued in an earlier cycle)"
        elif pct > float(self.stop["dataset_off_extra_pct"]):
            self.ask(slug, f"{r['name']} distorts the rest of the fit: {msg}. Keep it, exclude it (role excluded in spec.json, "
                     f"then re-init), or change its treatment?", ("keep", "exclude", "retreat"), default="keep", evidence="readme dataset-off entry; REPORT.md")
        return "stay", msg

    def _baseline_review_flags(self):
        """HIGH flags of the init baseline (candidate 1, reviewed once into REVIEW_baseline.md): a flag the base model
        already carried is the evaluator's (e.g. a broad pole spanning several
        unfitted peaks in the evaluator's final model, with 'no extra levels' a standing decision), so the drive reports it
        but waits only on a NEW one.  Keyed as the bench does: (code, message up to the first colon)."""
        c = self.campaign
        if self.state.get("baseline_review") is None:
            from .review import review_campaign
            flags = []
            try:
                bcid = None
                try:
                    row = c.ledger.db.execute("SELECT id FROM candidates WHERE move LIKE '%baseline%' ORDER BY id LIMIT 1").fetchone()
                    bcid = int(row[0]) if row else None
                except Exception:
                    bcid = None
                rb = review_campaign(c, grid_check=False, figures=0, cid=bcid or 1, path=os.path.join(c.dir, "REVIEW_baseline.md"))
                flags = [[f.code, f.message.split(":")[0].strip()] for f in rb["findings"] if f.severity == "HIGH"]
                self.log(f"baseline review (candidate {bcid or 1}) -> REVIEW_baseline.md, {len(flags)} HIGH flag(s) inherited from the base model")
            except Exception as err:
                self.log(f"baseline review failed ({err}): every HIGH flag counts as new")
            self.state["baseline_review"] = flags
            self.save()
        return {(a, b) for a, b in self.state["baseline_review"]}

    def _s_review(self):
        c = self.campaign
        from .review import review_campaign
        kw = {} if self.budget.get("close_fraction") is None else {"close_fraction": float(self.budget["close_fraction"])}   # review.CLOSE_GAP_FRACTION otherwise
        r = review_campaign(c, grid_check=bool(self.budget.get("grid_check", True)), figures=int(self.budget.get("figures", 6)), **kw)
        highs = [f for f in r["findings"] if f.severity == "HIGH"]
        counts = r["counts"]
        note = ", ".join(f"{k} {v}" for k, v in counts.items())
        key = lambda f: (f.code, f.message.split(":")[0].strip())
        from .review import same_flag_subject
        has = lambda f, keys: any(f.code == k[0] and same_flag_subject(key(f)[1], k[1]) for k in keys)   # a level may move within its width
        if highs:
            inherited = self._baseline_review_flags()
            old = [f for f in highs if has(f, inherited)]
            if old:
                note += f"; {len(old)} HIGH inherited from the base model (the evaluator's, not waited on): " + "; ".join(
                    f"{f.code}: {key(f)[1]}" for f in old[:3])
            # a flag the evaluator already answered 'continue' is not asked again in a later cycle (a benchmark,
            # 2026-11-10: the same broad-pole flag was queued in cycles 1 and 2)
            cont = {tuple(k) for k in self.state.get("review_continued", [])}
            seen = [f for f in highs if f not in old and has(f, cont)]
            if seen:
                note += f"; {len(seen)} HIGH already answered 'continue' (not asked again): " + "; ".join(f"{f.code}: {key(f)[1]}" for f in seen[:3])
            highs = [f for f in highs if f not in old and f not in seen]
        if highs:
            ans = self.ask("review-high", f"REVIEW.md has {len(highs)} new HIGH flag(s): " + "; ".join(f"{f.code}: {f.message[:80]}" for f in highs[:3])
                           + ". Look at figures/review_*.png; 'continue' to export as is, 'stop' to intervene.", ("continue", "stop"),
                           default="wait", evidence="REVIEW.md; figures/review_*.png")
            if ans is None:
                return "stay", f"{len(highs)} new HIGH review flags -> waiting"
            if ans.lower().startswith("stop"):
                raise DriveStop("STOP", "review HIGH flags: evaluator asked to stop")
            self.state["review_continued"] = sorted({tuple(k) for k in self.state.get("review_continued", [])} | {key(f) for f in highs})
            self.save()
        return "next", note

    def _s_export(self):
        c = self.campaign
        azr, sav, sc = c.export(note="rmfit drive export")
        c.report()
        return "next", f"{os.path.basename(azr)} objective {sc.objective:.3f}"

    def _s_cycle(self):
        obj = self._incumbent_obj()
        start = float(self.state.get("cycle_start_objective") or obj)
        gain = (start - obj) / max(abs(start), 1e-300)
        cyc = int(self.state["cycle"])
        path = os.path.join(self.dir, "candidates.json")
        bank = json.load(open(path)).get("candidates", []) if os.path.exists(path) else []
        untested = len(bank) - len(set(self.state["tested_candidates"]))
        msg = f"cycle {cyc}: {start:.3f} -> {obj:.3f} ({100*gain:.3f} %); {untested} candidates untested"
        self.note(f"CYCLE {cyc} complete: objective {start:.3f} -> {obj:.3f} ({100*gain:.3f} %), {untested} candidates untested; "
                  + ("another pass" if gain >= float(self.stop["cycle_gain_rel_min"]) and cyc < int(self.budget["cycles_max"]) else "DONE"))
        if gain >= float(self.stop["cycle_gain_rel_min"]) and cyc < int(self.budget["cycles_max"]):
            self.state["cycle"] = cyc + 1
            self.state["cycle_start_objective"] = obj
            self.state["dataset_off_tested"] = []
            self.counters["level_tests"] = 0
            self.counters["marginal_in_a_row"] = 0
            self.counters["levelscan_rounds"] = 0
            self.counters["levelscan_last_added"] = 0
            self.state["broad_done"] = False
            return "TIGHTEN", msg + " -> another pass"
        if (self.spec.get("uq") or {}).get("enabled"):
            return "UQ", msg + " -> UQ (spec uq.enabled)"
        return "DONE", msg + " -> DONE"

    def _s_uq(self):
        """Phase 4 (2026-11-10): rmfit.uq on the incumbent after the last cycle -- export, saturation check, Laplace, NUTS chains,
        gates, linear vs sampled bands, UQ.md.  Chains share the leg's slots; a failed gate after uq.max_attempts is a decision."""
        from . import uq as _uq
        su = dict(self.spec.get("uq") or {})
        slots = int(self.spec.get("nshards") or self.campaign.nshards or 4)
        chains = max(2, int(su.get("chains", 4)))
        nprocs = max(1, min(chains, slots))
        threads = max(1, slots // nprocs)
        ans = (self.ask("uq-gates", "", [], default="wait") if self.state.get("uq_asked") and any(d["slug"] == "uq-gates" for d in self.decisions.all())
               else None)
        if ans == "rerun":
            sp = os.path.join(self.dir, "uq", "uq_state.json")
            if os.path.exists(sp):
                us = json.load(open(sp)); us["attempt"] = int(us.get("attempt", 0)) + 1
                us["stages"].pop("SAMPLE", None); us["stages"].pop("GATES", None)
                json.dump(us, open(sp, "w"), indent=1)
            su["max_attempts"] = int(su.get("max_attempts", 2)) + 1
        elif ans == "stop":
            return "DONE", "UQ gates failed; the evaluator stopped the stage"
        st = _uq.run(campaign=self.campaign, spec_uq=su, log=self.log, nprocs=nprocs, threads=threads)
        g = st["stages"].get("GATES", {})
        ba = st["stages"].get("BANDS", {})
        summary = (f"R-hat max {g.get('rhat_max', float('nan')):.3f}, min ESS {g.get('ess_min', float('nan')):.0f}, divergent {100*g.get('divergent_frac', 0):.2f} %; "
                   f"sd ratio median {g.get('sd_ratio_median', float('nan')):.2f}; " + "; ".join(
                       f"{r['label']}: {r['best0']:.4g} +- {100*r['rel_lin0']:.1f} % linear / {100*r['rel_sampled0']:.1f} % sampled at {r['e0']:.3g} MeV"
                       for r in ba.get("grids", [])[:3]))
        if g.get("passed") or ans == "accept":
            self.note("UQ " + ("PASS" if g.get("passed") else "accepted by the evaluator with failed gates") + ": " + summary + " (UQ.md)")
            return "DONE", summary
        self.state["uq_asked"] = True
        self.save()
        a2 = self.ask("uq-gates", f"the UQ chains failed their gates after {st.get('attempt', 0) + 1} attempt(s): {'; '.join(g.get('fails', []))}"
                                  + (f" ({g['no_rerun']})" if g.get("no_rerun") else "") + ". "
                                  f"accept the result as is (UQ.md says FAIL), rerun with longer chains, or stop the stage?",
                      ["accept", "rerun", "stop"], default="wait", evidence=summary)
        if a2 is None:
            return "stay", "waiting on uq-gates: " + "; ".join(g.get("fails", []))
        return "stay", f"uq-gates answered {a2}"


PLAN_TEXT = {
    "PIN_CHECK": "imported engine vs the campaign pin / spec.engine; writes run_crc_drive",
    "AUDIT": "rmfit.audit -> AUDIT.md; HIGH -> decision audit-high (wait); then the three subagent audit prompts (audits/<kind>/PROMPT.md) -> decision audit-subagents (wait) unless pre_authorized.skip_subagent_audits",
    "REPRODUCE": "base model + .sav vs the archived chiSquared.out (plain chi2); mismatch -> decision (wait)",
    "SPEC_CHECK": "spec.json vs the model -> SPEC_CHECK.md; HIGH -> decision spec-high (wait)",
    "SCALE_ANCHOR": "spec data-set roles applied to <stem>_drive.azr (anchor fixed, free with prior, excluded off)",
    "SEED_CHANNELS": "free widths at exactly zero -> seed_frac x gamma_W, baked to <stem>_seeded.azr",
    "PAIR_COMPLETENESS": "ENSDF states of the residual nuclei open below the window top without a pair -> finding + decision (continue)",
    "LADDER": "model.ladder_Ex steps polished and carried forward (skipped without it)",
    "CANDIDATES": "candidates.json kept, or built from ENSDF <compound_nucleus> (broad first, hypothesis each)",
    "INIT": "Campaign.init --spec (baseline verified, engine pinned)",
    "RECALL": "structure.recall_levels broadest known levels removed one at a time and searched for with no J^pi (validates level-add); below recall_min_pass adds become decisions",
    "TIGHTEN": "polish --tighten rounds until gain < polish_gain_rel_min (budget polish_rounds_max)",
    "SIGN": "sign rounds until sign_rounds_without_accept (>= 2 rounds; budget sign_rounds_max)",
    "POPULATION": "restarts_per_cycle sign-perturbed restarts polished to convergence + crossover_passes between them (basin crossing); gain < population_gain_rel_min -> off until a structure change",
    "SHIFTS_FREE": "polish --shifts-free if pre_authorized.shifts_free",
    "LEVELSCAN": "when the broad candidate queue is empty: residual-run clusters -> template screen of every J^pi -> candidates appended (structure.levelscan_per_cycle)",
    "STRUCTURE_ADD": "level add per candidates.json entry, broad first (structure.broad_first_kev), narrow alignment after the broad queue + scan; hypothesis required; adoption per pre_authorized.adopt_level",
    "SIGN_POST": "sign rounds again if a structure was adopted",
    "DATASET_OFF": "Botek test on the worst chi2/N sets (budget dataset_off_max); large extra recovery -> decision (keep)",
    "REVIEW": "rmfit review (+ grid check) -> REVIEW.md; a NEW HIGH flag (not in the base model's REVIEW_baseline.md) -> decision review-high (wait)",
    "EXPORT": "export + report",
    "CYCLE": "another pass from TIGHTEN if the cycle gained >= cycle_gain_rel_min and cycles < cycles_max, else UQ (spec uq.enabled) or DONE",
    "UQ": "rmfit.uq on the incumbent: export, saturated poles held, Laplace, NUTS chains, R-hat/ESS/divergence gates, linear vs sampled bands -> UQ.md; failed gates after uq.max_attempts -> decision uq-gates (wait)",
}


def _cand_width_kev_known(cand):
    """Width in keV when the candidate states one (width_kev, or non-zero seed channel widths), else None."""
    w = cand.get("width_kev")
    if _isnum(w) and float(w) > 0:
        return float(w)
    tot = sum(abs(float(ch.get("gamma", 0.0))) for ch in (cand.get("channels") or []) if _isnum(ch.get("gamma")))
    return tot / 1e3 if tot > 0 else None


HBAR_LN2_KEV_S = 6.582119569e-19 * 0.6931471805599453     # hbar ln 2 (keV s): Gamma = hbar ln 2 / T_1/2


def _ensdf_width_kev(lv) -> Optional[float]:
    """Total width (keV) of an ENSDF level from pyazr.fetch_levels: ``width_kev`` when set, else from the half-life.
    fetch_levels leaves width_kev None for levels whose width ENSDF quotes in the half-life field (all of 16O above
    S_alpha: 9.585 1- "420" keV) and gives half_life_sec = hbar ln 2 / Gamma; without the width the broad-first sort fell
    back to energy order and the budget cut dropped the top of the window."""
    w = getattr(lv, "width_kev", None)
    if _isnum(w) and float(w) > 0:
        return float(w)
    t = getattr(lv, "half_life_sec", None)
    if not _isnum(t) or not (0.0 < float(t) < float("inf")):
        return None
    return float(f"{HBAR_LN2_KEV_S / float(t):.4g}")


def _cand_width_kev(cand) -> float:
    w = cand.get("width_kev")
    if _isnum(w):
        return float(w)
    g = sum(abs(float(ch.get("gamma", 0.0))) for ch in (cand.get("channels") or []) if _isnum(ch.get("gamma", 0.0)))
    return g / 1e3 if g > 0 else -1.0


def _isnum(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _nuc(z, m):
    a = int(round(float(m)))
    if z == 0 and a == 1:
        return "n"
    sym = SYMBOLS[z] if 0 <= z < len(SYMBOLS) else f"Z{z}"
    return {(1, 1): "p", (1, 2): "d", (1, 3): "t", (2, 4): "a"}.get((z, a), f"{a}{sym}")


def _vary_shifts_with_error(mdl) -> int:
    """Set varyShift on every active ``<segmentsData>`` line whose shift uncertainty is > 0 (the fields after dataNorm:
    varyNorm, normError, shift, shiftError, varyShift; dataNorm sits two fields later for phase-shift lines).  Returns the count."""
    if PYAZR3:
        return mdl.vary_shifts_with_error()        # rmfit.azr3.Model: the shift of every data set fixed with an error
    out, inside, n = [], False, 0
    for line in mdl._suffix.splitlines():
        s = line.strip()
        if s == "<segmentsData>":
            inside = True
        elif s == "</segmentsData>":
            inside = False
        elif inside and s and not s.startswith("#"):
            t = line.split()
            try:
                i = 8 + (2 if int(t[7]) % 10 == 2 else 0)
                if t[0] == "1" and float(t[i + 4]) > 0 and t[i + 5] == "0":
                    t[i + 5] = "1"; n += 1
                    line = " ".join(t)
            except (ValueError, IndexError):
                pass
        out.append(line)
    mdl._suffix = "\n".join(out)
    return n


def write_job_script(directory, spec: FitSpec, cfg: Optional[dict], name="run_crc_drive", overwrite=False, max_minutes=600) -> str:
    """``<dir>/run_crc_drive`` from ``templates/run_crc_drive.csh``: engine guard from campaign.json (or spec.engine before init),
    one leg per job, self-resubmission on CONTINUE."""
    from . import enginepin, site
    path = os.path.join(directory, name)
    if os.path.exists(path) and not overwrite:
        return path
    sched = site.scheduler()
    slurm = sched == "slurm"
    tpl = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates",
                       {"sge": "run_crc_drive.csh", "slurm": "run_slurm_drive.sh", "shell": "run_shell_drive.sh"}[sched])
    rec = (cfg or {}).get("engine")
    if not rec:
        eng = spec.get("engine") or {}
        rec = dict(azure2_root=eng["azure2_root"], so_md5=eng["so_md5"]) if eng.get("so_md5") and eng.get("azure2_root") else None
    parent = os.path.basename(os.path.dirname(os.path.abspath(directory)))
    reaction = parent if "+" in parent else os.path.basename(os.path.abspath(directory))
    slots = int(spec.get("nshards") or (cfg or {}).get("nshards") or 8)
    text = open(tpl).read()
    if slurm:
        leg, tmin = site.drive_leg_minutes(max_minutes)
        note = "" if leg == int(max_minutes) else (f"# DRIVE_MAX_MIN {int(max_minutes)} shortened to {leg} to fit the site's slurm_time "
                                                    f"({tmin} min incl. 30 for the last state and the resubmission)")
        text = (text.replace("{slurm_header}", site.slurm_header("drive_" + reaction.replace("+", ""), slots, "drive_slurm.log",
                                                               "a", time_minutes=tmin))
                    .replace("{leg_note}", note))
        max_minutes = leg
    if sched == "shell":
        text = text.replace("{lock_lines}", site.shell_lock_lines(os.path.join(os.path.abspath(directory), site.SHELL_LOCK)))
    guard = enginepin.guard_lines_csh if sched == "sge" else enginepin.guard_lines_sh
    text = text.replace("{setup_lines}", site.setup_lines())
    text = (text
            .replace("{mail_line}", site.mail_line())
            .replace("{queue_line}", site.queue_line())
            .replace("{pe_line}", site.pe_line(slots))
            .replace("{jobname}", "drive_" + reaction.replace("+", ""))
            .replace("{campaign_dir}", os.path.abspath(directory))
            .replace("{slots}", str(slots))
            .replace("{max_minutes}", str(int(max_minutes)))
            .replace("{engine_guard}", guard(rec))
            .replace("{pythonpath}", os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, 0o755)
    return path


def main_cli(a) -> int:
    """Dispatch of ``rmfit <dir> drive ...`` (called by rmfit.cli before a Campaign is opened: the build phase has none)."""
    if getattr(a, "reset", False) and os.path.exists(os.path.join(a.campaign, STATE_FILE)):
        bak = os.path.join(a.campaign, STATE_FILE + "." + time.strftime("%Y%m%d%H%M%S"))
        os.replace(os.path.join(a.campaign, STATE_FILE), bak)
        print("drive state reset; previous state kept as", bak)
    try:
        d = Driver(a.campaign, base=getattr(a, "base", None), sav=getattr(a, "sav", None), local=getattr(a, "local", False))
    except DriveStop as e:
        print(f"{e.status}: {e}")
        return 2
    if getattr(a, "state", None):
        if a.state not in STATES:
            print("unknown state", a.state, "; states:", " ".join(STATES))
            return 2
        d.state["state"] = a.state
        d.state["status"] = "CONTINUE"
        d.save()
        d.log(f"state set to {a.state} by --state")
    if getattr(a, "write_job", False):
        print("job script:", write_job_script(d.dir, d.spec, d.campaign.cfg if d.campaign is not None else None, overwrite=True,
                                              max_minutes=getattr(a, "max_minutes", 0) or 600))
        return 0
    status = d.run(max_minutes=float(getattr(a, "max_minutes", 0) or 0), dry_run=getattr(a, "dry_run", False))
    print(f"DRIVE_STATUS: {status} state {d.state['state']}")
    return 0 if status in ("CONTINUE", "DONE", "WAITING", "STOP") else 1
