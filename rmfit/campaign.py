"""A campaign directory: configuration, evaluators, ledger, incumbent.

Layout of a campaign directory (created by ``rmfit init``):

    campaign.json     base model, engine flags, shard count, bounds policy, search config
    <name>.azr/.sav   the base model and its parameter companion (GUI-reviewable)
    data/             data files (the .azr references them relatively)
    variants/         every derived .azr (shards, single-session copies, bakes)
    ledger.sqlite, candidates/, LOG.jsonl, STATUS.md
    readme            free-text lab notebook (rmfit appends; humans append)
"""
import dataclasses
import json
import math
import os
import time
from typing import Dict, Optional, Sequence

import numpy as np

from . import ARCHIVE_FLAGS
from .engine import SingleSession, Sharded, Score, ObjectiveSpec, bake, verify
from .keys import ParamKey, take, level_group
from .ledger import Ledger
from .local import BoundsPolicy, FitResult, default_bounds, polish, polish_refreeze, sanity, free_mask
from .model import single_variant, base_segments, rows_in_segment
from .search import SearchConfig, sign_round, crossover_round, tau_flip
from .structure import LevelCandidate, add_level_variant, scan_jpi_on_window, classify, JPI_GROUPS, resolve_jpis  # noqa: F401
from . import enginepin
from .backend import PYAZR3
from .spec import FitSpec
from .decisions import DecisionQueue


def _dc_from(cls, d):
    names = {f.name for f in dataclasses.fields(cls)}
    return cls(**{k: v for k, v in (d or {}).items() if k in names})


class Campaign:
    def __init__(self, directory, enforce_engine=True):
        """``enforce_engine``: refuse to open on an engine other than the pinned one (``campaign.json["engine"]``,
        2026-11-07); the CLI passes False for status/report/log, which only warn."""
        self.dir = os.path.abspath(directory)
        with open(os.path.join(self.dir, "campaign.json")) as f:
            self.cfg = json.load(f)
        self.base_azr = os.path.normpath(os.path.join(self.dir, self.cfg["base_azr"]))   # a base outside the dir ("../x.azr") must match the ledger's absolute path
        self.flags = dict(self.cfg.get("flags", ARCHIVE_FLAGS))
        self.nshards = int(self.cfg.get("nshards", 24))
        self.policy = _dc_from(BoundsPolicy, self.cfg.get("policy"))
        self.objective = _dc_from(ObjectiveSpec, self.cfg.get("objective"))
        if self.objective.cov_reference and not os.path.isabs(self.objective.cov_reference):
            self.objective = dataclasses.replace(self.objective, cov_reference=os.path.join(
                self.dir, self.objective.cov_reference))
        self.search = _dc_from(SearchConfig, self.cfg.get("search"))
        self.vdir = os.path.join(self.dir, "variants")
        os.makedirs(self.vdir, exist_ok=True)
        self.ledger = Ledger(os.path.join(self.dir, "ledger.sqlite"))
        self._sharded: Optional[Sharded] = None
        self._log_fh = None
        # engine pin (enginepin.py): a campaign is one engine; objectives across engines are not comparable
        self.engine, mismatch = enginepin.check(self.cfg.get("engine"), strict=enforce_engine, log=self.log)
        if mismatch and enforce_engine:
            self.note("ENGINE CHANGE allowed by RMFIT_ALLOW_ENGINE_CHANGE: " + mismatch)
        # fit specification (spec.py) and the decisions queue (decisions.py), 2026-11-07
        self.spec: Optional[FitSpec] = None
        sp = os.path.join(self.dir, "spec.json")
        if os.path.exists(sp):
            try:
                self.spec = FitSpec.load(sp)
                if self.cfg.get("spec_sha") and self.spec.sha != self.cfg["spec_sha"]:
                    self.log("WARNING spec.json changed since init (sha differs from campaign.json spec_sha): "
                             "run `rmfit spec check`; re-init with --spec to re-derive the config")
            except Exception as err:
                self.log(f"WARNING spec.json unreadable: {err}")
        self.decisions = DecisionQueue(os.path.join(self.dir, "DECISIONS.md"))

    # -- logging ------------------------------------------------------------------

    def log(self, *a):
        msg = " ".join(str(x) for x in a)
        line = f"{time.strftime('%H:%M:%S')} {msg}"
        print(line, flush=True)
        with open(os.path.join(self.dir, "campaign.log"), "a") as fh:
            fh.write(line + "\n")

    def note(self, text):
        """Append a dated entry to the campaign's readme (the lab notebook)."""
        with open(os.path.join(self.dir, "readme"), "a") as fh:
            fh.write(f"\n{time.strftime('%Y-%m-%d %H:%M')} (rmfit): {text}\n")
        self.ledger.log("readme", text)

    # -- evaluators ---------------------------------------------------------------

    @property
    def sharded(self) -> Sharded:
        if self._sharded is None:
            t0 = time.time()
            self._sharded = Sharded(self.base_azr, os.path.join(self.vdir, "shards"), objective=self.objective,
                                    nshards=self.nshards, flags=self.flags,
                                    weights=self.cfg.get("weights"))
            self.log(f"sharded evaluator up: {len(self._sharded.shards)} shards, "
                     f"{self._sharded.npts} points, {len(self._sharded.keys)} free, {time.time()-t0:.0f}s")
        return self._sharded

    def close(self):
        if self._sharded is not None:
            self._sharded.close()
            self._sharded = None

    # -- structure / incumbent ---------------------------------------------------------

    @property
    def structure_id(self) -> int:
        for s in self.ledger.structures():
            if os.path.normpath(os.path.abspath(s["azr"])) == self.base_azr:
                return s["id"]
        raise RuntimeError("campaign not initialized: run `rmfit init`")

    def incumbent(self):
        """``(cid, x, objective, keys)`` of the best verified candidate."""
        b = self.ledger.best(self.structure_id, stage="verified")
        if b is None:
            raise RuntimeError("no verified candidate yet: run `rmfit init`")
        cid, obj = b
        x, keys = self.ledger.load(cid)
        return cid, x, obj, keys

    def verify_candidate(self, x, keys, move=None, parent=None, note="") -> tuple:
        """Fresh single-session score of ``x``; records it as *verified*."""
        base_keys = self._base_keys()
        missing = set(base_keys) - set(keys)
        if missing:
            raise ValueError(f"vector lacks {len(missing)} base parameters: " + ", ".join(str(k) for k in list(missing)[:5]))
        xb = take(x, keys, base_keys, self._base_x0())
        sc, _ = verify(self.base_azr, xb, self.flags, workdir=self.vdir, name="verify", objective=self.objective,
                       keys=base_keys)
        cid = self.ledger.add_candidate(self.structure_id, xb, base_keys, move=move, parent=parent, status="verified")
        self.ledger.record(cid, sc, "verified", evaluator="single")
        self.log(f"verified candidate {cid}: objective {sc.objective:.3f} (chi2 {sc.chi2:.3f} + pen {sc.penalty:.3f}) {note}"
                 + ("" if sc.chi2_self is None else f" [self-consistent chi2 {sc.chi2_self:.3f}, "
                    f"{sc.chi2_self - sc.chi2:+.3f} vs the fixed reference]"))
        self._log_target_distance(xb, base_keys, sc.objective)
        return cid, sc

    def _log_target_distance(self, x, keys, objective):
        """For a perturbation benchmark (``campaign.json`` has a ``perturbation``
        block naming the parent campaign and candidate): log the gauge-invariant
        sign distance of a verified vector to that target and the remaining gap."""
        pert = self.cfg.get("perturbation")
        if not pert:
            return
        try:
            from .keys import sign_distance
            lt = Ledger(os.path.join(pert["parent"], "ledger.sqlite"))
            xt, kt = lt.load(int(pert["parent_candidate"]))
            if [str(k) for k in kt] != [str(k) for k in keys]:
                self.log("  (target key set differs; no sign distance)")
                return
            dist, diff = sign_distance(x, xt, keys, self.ledger.fixed_of(self.structure_id))
            gap = objective / float(pert["parent_objective"]) - 1.0
            self.log(f"  sign distance to the target {dist} ({gap:+.2%} above it): "
                     + "; ".join(str(k) for k in diff))
        except Exception as err:            # diagnostics must never stop a campaign
            self.log(f"  (sign distance to the target unavailable: {err})")

    def reaction_name(self) -> str:
        """The archive's reaction directory name (``<A><El>+<proj>``) when the campaign
        sits in one, else the base model's stem."""
        parent = os.path.basename(os.path.dirname(self.dir))
        if "+" in parent:
            return parent
        return os.path.splitext(os.path.basename(self.cfg.get("base_azr", "fit.azr")))[0]

    def write_job_script(self, name="run_crc_rmfit", overwrite=False) -> str:
        """The job script for this campaign from ``rmfit/templates`` (SGE, Slurm or plain shell: site.scheduler())."""
        from . import site
        path = os.path.join(self.dir, name)
        if os.path.exists(path) and not overwrite:
            return path
        sched = site.scheduler()
        slurm = sched == "slurm"
        tpl = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates",
                           {"sge": "run_crc_rmfit.csh", "slurm": "run_slurm_rmfit.sh", "shell": "run_shell_rmfit.sh"}[sched])
        with open(tpl) as f:
            text = f.read()
        jobname = "rmfit_" + self.reaction_name().replace("+", "")
        if slurm:
            text = text.replace("{slurm_header}", site.slurm_header(jobname, int(self.cfg.get("nshards") or 24), "rmfit_slurm.log", "abe"))
        if sched == "shell":
            text = (text.replace("{lock_lines}", site.shell_lock_lines(os.path.join(self.dir, site.SHELL_LOCK)))
                        .replace("{nshards}", str(int(self.cfg.get("nshards") or 24))))
        text = text.replace("{setup_lines}", site.setup_lines())
        guard = enginepin.guard_lines_csh if sched == "sge" else enginepin.guard_lines_sh
        text = (text.replace("{mail_line}", site.mail_line())
                    .replace("{queue_line}", site.queue_line())
                    .replace("{pe_line}", site.pe_line(24))         # the SGE template's 24 slots (2026-11-07 on)
                    .replace("{jobname}", jobname)
                    .replace("{campaign_dir}", self.dir)
                    .replace("{engine_guard}", guard(self.cfg.get("engine")))
                    .replace("{pythonpath}", os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
        with open(path, "w") as f:
            f.write(text)
        os.chmod(path, 0o755)
        return path

    def _base_keys(self):
        return self.ledger.structure(self.structure_id)["keys"]

    def _base_x0(self):
        return np.zeros(len(self._base_keys()))

    # -- initialisation -------------------------------------------------------------------

    def _seed_vector(self, ev, sav=None):
        """The start vector of a campaign (sav: a parameter file seeding the free values); pyazr 3 / pyazr 2 versions below."""
        # through the class: rmfit.drive.Driver calls Campaign._seed_vector with itself as ``self``
        return (Campaign._seed_vector_p3(self, ev, sav) if PYAZR3 else Campaign._seed_vector_p2(self, ev, sav))

    def _seed_vector_p2(self, ev, sav=None):
        """The base model's free vector, with the values of a ``.sav`` companion applied."""
        x = ev.x0.copy()
        if sav:
            vals = {}
            for line in open(sav):
                t = line.split()
                if len(t) >= 2:
                    try:
                        vals[t[0]] = float(t[1])
                    except ValueError:
                        pass
            snap = ev.worker.call("snapshot", None)
            # A .sav names level parameters POSITIONALLY (energy_N, width_N_c).  After any
            # level add/remove those names point at the wrong levels, and seeding from them
            # produces a garbage vector (an objective about twenty times too large).
            # The exported <levels> block already carries the fitted R-matrix values, so if
            # the .sav's level count disagrees with the model's, seed norms/shifts only.
            sav_levels = {int(k.split("_")[1]) for k in vals if k.startswith("energy_")}
            model_levels = {int(p.name.split("_")[1]) for p in snap["parameters"]
                            if p.name.startswith("energy_")}
            skip_levels = bool(sav_levels) and (max(sav_levels) != max(model_levels | {0}))
            if skip_levels:
                self.log(f"init: WARNING {os.path.basename(sav)} indexes {max(sav_levels)} levels but the "
                         f"model has {max(model_levels | {0})}: a positional .sav is invalid after a level "
                         f"add/remove, so its level parameters are IGNORED and only norms/shifts are seeded")
            n = 0
            for p in snap["parameters"]:
                if p.fixed or p.free_index is None or p.name not in vals:
                    continue
                if skip_levels and not p.name.startswith("segment_"):
                    continue
                x[p.free_index] = vals[p.name]; n += 1
            self.log(f"init: {n} free parameters taken from {os.path.basename(sav)}"
                     + (" (norms/shifts only)" if skip_levels else ""))
        return x

    def _seed_vector_p3(self, ev, sav=None):
        """The base model's free vector, with the values of a ``.sav`` companion applied."""
        x = ev.x0.copy()
        if sav:
            from pyazr.parameters import read_parameter_file
            file_basis, vals = read_parameter_file(sav)
            snap = ev.worker.call("snapshot", None)
            if file_basis is not None and str(file_basis) != str(snap["basis"]):
                raise ValueError(f"init: {os.path.basename(sav)} holds amplitudes of the basis {file_basis!r}, the "
                                 f"campaign's sessions fit in {snap['basis']!r}: its values are not this vector's")
            # A parameter file names level parameters POSITIONALLY (energy_N, width_N_c).  After any
            # level add/remove those names point at the wrong levels, and seeding from them
            # produces a garbage vector (an objective about twenty times too large).
            # The exported <levels> block already carries the fitted R-matrix values, so if
            # the file's level count disagrees with the model's, seed the nuisances only.
            sav_levels = {int(k.split("_")[1]) for k in vals if k.startswith("energy_")}
            model_levels = {int(p["name"].split("_")[1]) for p in snap["parameters"] if p["kind"] == "energy"}
            skip_levels = bool(sav_levels) and (max(sav_levels) != max(model_levels | {0}))
            if skip_levels:
                self.log(f"init: WARNING {os.path.basename(sav)} indexes {max(sav_levels)} levels but the "
                         f"model has {max(model_levels | {0})}: a positional parameter file is invalid after a level "
                         f"add/remove, so its level parameters are IGNORED and only norms/shifts are seeded")
            n = 0
            for p in snap["parameters"]:
                if p["fixed"] or p["free_index"] is None or p["name"] not in vals:
                    continue
                if skip_levels and p["kind"] in ("energy", "width"):
                    continue
                x[p["free_index"]] = vals[p["name"]]; n += 1
            self.log(f"init: {n} free parameters taken from {os.path.basename(sav)}"
                     + (" (norms/shifts only)" if skip_levels else ""))
        return x

    @staticmethod
    def apply_sav_fixed(base_azr, sav, flags, log=print, mode="ignore"):
        """Reconcile a parameter file's FIXED values with the model (see the backend versions below: _p2 pyazr 2, _p3 pyazr 3)."""
        return (Campaign.apply_sav_fixed_p3(base_azr, sav, flags, log, mode) if PYAZR3 else Campaign.apply_sav_fixed_p2(base_azr, sav, flags, log, mode))

    @staticmethod
    def apply_sav_fixed_p2(base_azr, sav, flags, log=print, mode="ignore"):
        """What the FIXED parameters of a ``.azr`` + ``.sav`` pair mean (2026-11-10).

        AZURE2's parameter-file reader sets every parameter it names, FIXED level energies and widths included
        (``AZUREParams::ReadUserParameters``; pyazr ``azure2.load_sav``), but ``_seed_vector`` takes only the FREE
        values from the file and the sharded evaluator takes the fixed ones from the ``.azr``.  When the two disagree
        there are two models, and the archive uses BOTH conventions: a CLI fit's ``param.sav`` holds the fixed values
        the fit ran with (in one archived fit the two differed in two dozen fixed background widths by up
        to 10x, and a benchmark initialised at more than twice the archived objective on the wrong model), while a
        campaign built by a script (the ``.azr`` carries the frozen block and the added level, the
        ``.sav`` is the previous export, 125-237 fixed values apart) means the ``.azr`` and uses the file as a seed of
        the free values only.  So the disagreement is never silent, and which file wins is the caller's:
        ``mode="ignore"`` (default, the behaviour since the first campaign) keeps the ``.azr``'s fixed values and logs a
        WARNING with the objective both ways; ``mode="apply"`` (``init --sav-fixed apply``: the pair is a CLI fit) writes
        ``<stem>_sav.azr`` (+ ``.sav``) beside the base with the file applied (``pyazr.azure2.save_fit``, fresh-session
        verified) and returns that pair.  Returns ``(azr, sav, note_or_None)``.  A positional ``.sav`` whose level count
        disagrees with the model is left alone (its level parameters are already ignored by ``_seed_vector``)."""
        import numpy as np
        from pyazr import azure2
        base_azr = os.path.abspath(base_azr); sav = os.path.abspath(sav)
        vals = {}
        for line in open(sav):
            t = line.split()
            if len(t) >= 2:
                try:
                    vals[t[0]] = float(t[1])
                except ValueError:
                    pass
        sav_levels = {int(k.split("_")[1]) for k in vals if k.startswith("energy_")}
        kw = {k: v for k, v in dict(flags or {}).items() if k in ("use_brune", "use_park", "ignore_externals", "use_long_wavelength")}
        with azure2(base_azr, data_mode=True, **kw) as m:
            model_levels = {int(p.name.split("_")[1]) for p in m.parameters if p.name.startswith("energy_")}
            if sav_levels and max(sav_levels) != max(model_levels | {0}):
                return base_azr, sav, None
            cur = np.asarray(m.sess.params_all_rwa(), float).copy()
            full = m.full_rwa_from_sav(sav)
            diff = [(p.name, float(cur[p.index]), float(full[p.index])) for p in m.parameters
                    if p.fixed and 0 <= p.index < full.size
                    and abs(cur[p.index] - full[p.index]) > 1e-6 * max(1e-6, abs(cur[p.index]), abs(full[p.index]))]   # a .sav holds 8 digits
            if not diff:
                return base_azr, sav, None
            worst = max(diff, key=lambda d: abs(d[2] - d[1]) / max(abs(d[1]), abs(d[2]), 1e-300))
            free_idx = [p.index for p in m.parameters if not p.fixed and p.index >= 0]
            obj_azr = float(m.objective(full[free_idx]))         # the .azr's fixed values, the file's free ones (what init seeds)
            x = m.load_sav(sav)
            obj_sav = float(m.objective(x))                       # the pair as the CLI reads it
            head = (f"{len(diff)} FIXED parameter(s) of {os.path.basename(sav)} differ from {os.path.basename(base_azr)} "
                    f"(largest {worst[0]}: .azr {worst[1]:.6g}, .sav {worst[2]:.6g}); plain objective {obj_azr:.3f} with the .azr's "
                    f"fixed values, {obj_sav:.3f} with the .sav's (the CLI's reading)")
            if mode != "apply":
                note = (head + " -- the .azr's are USED (the .sav seeds the free values only); if this pair is a CLI fit and its "
                        "param.sav holds the fixed values the fit ran with, re-init with --sav-fixed apply")
                log("init: WARNING " + note)
                return base_azr, sav, note
            out = os.path.splitext(base_azr)[0] + "_sav.azr"
            azr2, sav2 = m.save_fit(out, x, param_sav=True, verify=True)
            note = head + f" -- --sav-fixed apply: the campaign's base is {os.path.basename(azr2)} (+ .sav) = the .azr with the .sav applied, fresh-session verified"
            log("init: WARNING " + note)
        return azr2, sav2, note

    @staticmethod
    def apply_sav_fixed_p3(base_azr, sav, flags, log=print, mode="ignore"):
        """What the FIXED parameters of a ``.azr`` + ``.sav`` pair mean (2026-11-10).

        AZURE2's parameter-file reader sets every parameter it names, FIXED level energies and widths included
        (``pyazr.Session(file, parameter_file=)``), but ``_seed_vector`` takes only the FREE
        values from the file and the sharded evaluator takes the fixed ones from the ``.azr``.  When the two disagree
        there are two models, and the archive uses BOTH conventions: a CLI fit's ``param.sav`` holds the fixed values
        the fit ran with (in one archived fit the two differed in two dozen fixed background widths by up
        to 10x, and a benchmark initialised at more than twice the archived objective on the wrong model), while a
        campaign built by a script (the ``.azr`` carries the frozen block and the added level, the
        ``.sav`` is the previous export, 125-237 fixed values apart) means the ``.azr`` and uses the file as a seed of
        the free values only.  So the disagreement is never silent, and which file wins is the caller's:
        ``mode="ignore"`` (default, the behaviour since the first campaign) keeps the ``.azr``'s fixed values and logs a
        WARNING with the objective both ways; ``mode="apply"`` (``init --sav-fixed apply``: the pair is a CLI fit) writes
        ``<stem>_sav.azr`` (+ ``.sav``) beside the base with the file applied (``pyazr.Session.save_fit``, verified by
        reopening) and returns that pair.  Returns ``(azr, sav, note_or_None)``.  A positional ``.sav`` whose level count
        disagrees with the model is left alone (its level parameters are already ignored by ``_seed_vector``)."""
        import numpy as np
        from pyazr import Session
        from pyazr.parameters import read_parameter_file
        base_azr = os.path.abspath(base_azr); sav = os.path.abspath(sav)
        _, vals = read_parameter_file(sav)
        sav_levels = {int(k.split("_")[1]) for k in vals if k.startswith("energy_")}
        kw = {k: v for k, v in dict(flags or {}).items() if k in ("use_brune", "use_park", "ignore_externals", "use_long_wavelength")}
        # two sessions of one project: as written, and with the file taken as the program takes it (every value
        # of the file, the fixed ones included: Session(file, parameter_file=))
        with Session(base_azr, **kw) as m, Session(base_azr, parameter_file=sav, **kw) as ms:
            model_levels = {int(p.name.split("_")[1]) for p in m.parameters if p.kind == "energy"}
            if sav_levels and max(sav_levels) != max(model_levels | {0}):
                return base_azr, sav, None
            diff = [(p.name, float(p.start), float(q.start)) for p, q in zip(m.parameters, ms.parameters)
                    if p.fixed and abs(p.start - q.start) > 1e-6 * max(1e-6, abs(p.start), abs(q.start))]   # a file holds 8 digits
            if not diff:
                return base_azr, sav, None
            worst = max(diff, key=lambda d: abs(d[2] - d[1]) / max(abs(d[1]), abs(d[2]), 1e-300))
            obj_azr = float(m.objective()(m.parameters.read(sav)))   # the .azr's fixed values, the file's free ones (what init seeds)
            x = ms.parameters.x0
            obj_sav = float(ms.objective()(x))                       # the pair as the CLI reads it
            head = (f"{len(diff)} FIXED parameter(s) of {os.path.basename(sav)} differ from {os.path.basename(base_azr)} "
                    f"(largest {worst[0]}: .azr {worst[1]:.6g}, .sav {worst[2]:.6g}); plain objective {obj_azr:.3f} with the .azr's "
                    f"fixed values, {obj_sav:.3f} with the .sav's (the CLI's reading)")
            if mode != "apply":
                note = (head + " -- the .azr's are USED (the .sav seeds the free values only); if this pair is a CLI fit and its "
                        "param.sav holds the fixed values the fit ran with, re-init with --sav-fixed apply")
                log("init: WARNING " + note)
                return base_azr, sav, note
            out = os.path.splitext(base_azr)[0] + "_sav.azr"
            azr2, sav2 = ms.save_fit(out, x, param_sav=True, verify=True)
            note = head + f" -- --sav-fixed apply: the campaign's base is {os.path.basename(azr2)} (+ .sav) = the .azr with the .sav applied, fresh-session verified"
            log("init: WARNING " + note)
        return azr2, sav2, note

    @classmethod
    def init(cls, directory, base_azr, sav=None, flags=None, nshards=24, policy=None, search=None, note="",
             objective=None, cov_reference=None, spec=None, sav_fixed="ignore"):
        """Create campaign.json, register the base model and verify the baseline.

        Data covariance with model-scaled (``relvector``) terms (2026-10-25, the evaluator: "use a fixed reference model"):
        the campaign fixes those terms at ONE reference model for all its evaluations -- polish, verify, export --
        so every number in its ledger is the same objective.  By default the reference is the seed's own model
        (frozen at the baseline evaluation and saved as ``cov_reference.npz``); seeds are built from the last
        adopted fit, so that is the best finished fit.  ``cov_reference`` = path of another fit's .npz (e.g. to give
        a control and a variant the same reference).  The reference is updated only by starting the next campaign.
        An existing campaign's reference is kept on re-init."""
        directory = os.path.abspath(directory)
        base_azr = os.path.abspath(base_azr)
        os.makedirs(directory, exist_ok=True)
        # spec (2026-11-07): copied to <dir>/spec.json, its flags/policy/search/objective/nshards used where the explicit
        # arguments are absent, its sha recorded; `rmfit spec check` reports disagreements with the .azr before the baseline
        spec_sha = None
        if spec:
            import shutil
            spec_src = os.path.abspath(spec)
            spec_dst = os.path.join(directory, "spec.json")
            if os.path.abspath(spec_src) != os.path.abspath(spec_dst):
                shutil.copy2(spec_src, spec_dst)
            fs = FitSpec.load(spec_dst)
            derived = fs.derive_campaign_cfg()
            flags = flags if flags is not None else derived.get("flags")
            policy = policy if policy else derived.get("policy")
            search = search if search else derived.get("search")
            objective = objective if objective else derived.get("objective")
            if "nshards" in derived and nshards == 24:
                nshards = derived["nshards"]
            spec_sha = fs.sha
            probs = fs.check(base_azr)
            highs = [f for f in probs if f.severity == "HIGH"]
            for f in probs:
                print("spec check:", f.line())
            if highs:
                raise ValueError(f"spec check: {len(highs)} HIGH disagreement(s) between spec.json and the .azr -- fix one of them first: "
                                 + "; ".join(f.message for f in highs[:3]))
        # an existing campaign.json keeps its policy / search / weights unless
        # new ones are given (re-running init must not silently reset them)
        prev = {}
        cj = os.path.join(directory, "campaign.json")
        if os.path.exists(cj):
            with open(cj) as f:
                prev = json.load(f)
        cfg = dict(prev)
        cfg.update(base_azr=os.path.relpath(base_azr, directory),
                   flags=(prev.get("flags", ARCHIVE_FLAGS) if flags is None else flags),
                   nshards=nshards, policy=(policy if policy else prev.get("policy", {})),
                   search=(search if search else prev.get("search", {})),
                   objective=(objective if objective else prev.get("objective", {})),
                   created=prev.get("created", time.strftime("%Y-%m-%d %H:%M")),
                   reinit=time.strftime("%Y-%m-%d %H:%M"))
        if spec_sha:
            cfg["spec_sha"] = spec_sha
        # engine pin: recorded once at the first init; re-pinned only with RMFIT_ALLOW_ENGINE_CHANGE=1
        old_pin = prev.get("engine")
        if "engine" not in prev or os.environ.get("RMFIT_ALLOW_ENGINE_CHANGE", "") not in ("", "0"):
            cfg["engine"] = enginepin.fingerprint()
        repinned = bool(old_pin and old_pin.get("so_md5") and old_pin.get("so_md5") != cfg["engine"].get("so_md5"))
        if cov_reference:
            import shutil
            shutil.copy2(os.path.abspath(cov_reference), os.path.join(directory, "cov_reference.npz"))
            cfg["objective"] = dict(cfg["objective"], cov_reference="cov_reference.npz")
        with open(cj, "w") as f:
            json.dump(cfg, f, indent=2)
        c = cls(directory)
        c.log(f"init: engine pinned {enginepin.short(c.cfg['engine'])}" +
              (" (new pin)" if "engine" not in prev else
               f" (re-pinned from {enginepin.short(old_pin)})" if repinned else " (kept)"))
        if repinned:
            c.note(f"init: engine re-pinned (RMFIT_ALLOW_ENGINE_CHANGE): {enginepin.short(old_pin)} -> "
                   f"{enginepin.short(c.cfg['engine'])}; objectives before and after are not comparable")
            enginepin.repin_job_scripts(directory, old_pin, c.cfg["engine"], log=c.log)
        if cov_reference:
            c.log(f"init: relvector covariance model fixed at the reference {os.path.abspath(cov_reference)} "
                  f"(copied to cov_reference.npz)")
        if sav:
            base0 = base_azr
            base_azr, sav, sav_note = cls.apply_sav_fixed(base_azr, sav, c.flags, log=c.log, mode=sav_fixed)
            if base_azr != base0:
                c.cfg["base_azr"] = cfg["base_azr"] = os.path.relpath(base_azr, directory)
                with open(cj, "w") as f:
                    json.dump(cfg, f, indent=2)
                c.base_azr = os.path.normpath(base_azr)
                c.note(f"init: {os.path.basename(sav)} carries fixed level parameters that differ from {os.path.basename(base0)}; "
                       f"the campaign's base is {os.path.basename(base_azr)} (the .azr with the .sav applied, as the CLI reads the pair)")
            elif sav_note:
                c.note("WARNING at init: " + sav_note)
        src = single_variant(base_azr, "init", c.vdir)
        if ((c.objective.loss == "sivia" and c.objective.loss_scale != "none") or c.objective.dataset_weight == "scaled") \
                and not c.objective.scales:
            # reference reduced chi-squared of every data set (or one global
            # value) from the baseline under the plain objective
            with SingleSession(src, cwd=c.vdir, flags=c.flags, name="init_plain") as ev0:
                x0 = c._seed_vector(ev0, sav)
                sc0 = ev0.score(x0)
                seg_ds = ev0.dataset_of_key()
                n_of, chi_of = {}, {}
                for key, (a, b) in ev0.row_ranges:
                    ds = str(seg_ds.get(int(key), int(key)))
                    n_of[ds] = n_of.get(ds, 0) + (b - a)
                    chi_of[ds] = chi_of.get(ds, 0.0) + sc0.segment_chi2.get(key, 0.0)
                if c.objective.loss_scale == "global":
                    scales = {"*": sc0.chi2 / max(ev0.npts, 1)}
                else:
                    scales = {ds: max(chi_of[ds] / max(n_of[ds], 1), 1.0) for ds in n_of}
            c.cfg["objective"]["scales"] = scales
            c.objective = dataclasses.replace(c.objective, scales=scales)
            with open(cj, "w") as f:
                json.dump(c.cfg, f, indent=2)
            c.log(f"objective: {c.objective.label}; reference chi2/N per data set from the plain baseline "
                  f"({sc0.chi2:.1f} / {ev0.npts}): " + ", ".join(f"{k} {v:.1f}" for k, v in sorted(scales.items())))
        with SingleSession(src, cwd=c.vdir, flags=c.flags, name="init", objective=c.objective) as ev:
            keys = ev.keys
            free_rel = ev.free_norms_on_relvector()
            bad = {f: s for f, s in free_rel.items() if f not in set(c.objective.cov_free_norm_ok or [])}
            if bad:
                raise ValueError(
                    "free data normalization on relvector-covariance rows: " + "; ".join(
                        f"{f}: segment(s) {', '.join(map(str, s))}" for f, s in sorted(bad.items())) +
                    ".  A free norm there is biased whichever way the relvector reference follows it (in one large model: "
                    "a norm ran up fivefold with a data-unit reference and down threefold norm-following).  t0 method (NNPDF, JHEP 05 (2010) "
                    "075; an archived campaign): fix the norm (varyNorm 0) and put its prior in the covariance file as a "
                    "relvector term (all ones x the prior over the rows it normalizes).  To keep the norm free anyway, "
                    "name the file in the objective's cov_free_norm_ok list.")
            for f, s in sorted(free_rel.items()):
                c.log(f"init: WARNING free norm on relvector rows of {f} (segments {', '.join(map(str, s))}), allowed by "
                      f"cov_free_norm_ok -- biased in principle (10-26-26_t0_cov); watch these norms")
            x = c._seed_vector(ev, sav)
            sc = ev.score(x)
            if ev.cov is not None and ev.cov.relative and not c.objective.cov_reference:
                # the first evaluation froze the relvector terms at the seed: keep that model for the campaign
                from .covariance import save_reference
                save_reference(os.path.join(c.dir, "cov_reference.npz"), ev.cov.frozen_models(),
                               note=f"seed {os.path.basename(base_azr)} at init {time.strftime('%Y-%m-%d %H:%M')}")
                c.cfg["objective"] = dict(c.cfg.get("objective") or {}, cov_reference="cov_reference.npz")
                with open(cj, "w") as f:
                    json.dump(c.cfg, f, indent=2)
                c.objective = dataclasses.replace(c.objective, cov_reference=os.path.join(c.dir, "cov_reference.npz"))
                ev.spec = c.objective
                c.log("init: relvector covariance model fixed at the seed's model for this campaign "
                      "(cov_reference.npz; polish, verify and export all use it)")
            elif c.objective.cov_reference:
                c.log(f"init: relvector covariance model fixed at {c.objective.cov_reference}; self-consistent chi2 "
                      f"at the seed {sc.chi2_self:.3f}" if sc.chi2_self is not None else
                      f"init: covariance reference {c.objective.cov_reference} (no relvector terms in this model)")
            sid = c.ledger.add_structure(os.path.basename(base_azr), base_azr, keys, note=note,
                                         fixed=ev.fixed_amplitudes)
            cid = c.ledger.add_candidate(sid, x, keys, move=dict(kind="baseline"), status="verified")
            c.ledger.record(cid, sc, "verified", evaluator="single")
            rep = sanity(ev, x, c.policy, default_bounds(ev, x, c.policy))
        # the archived chiSquared.out in the base model's output directory (skill section 2 item 2: reproduce the last
        # finished fit before building on it).  A seed that does not reproduce it is on the wrong flags, the wrong or a
        # positional .sav, or a stale output/ (a benchmark in which ignore_externals False gave a chi2 a third above the
        # archived one; the campaign and its perturbation benchmark were built on the wrong objective).  Plain chi2 is
        # compared (the per-segment scores stay plain whatever the objective); a warning, not a refusal.
        try:
            from .review import output_dir_of, read_chisq_out
            odir = output_dir_of(base_azr)
            cs = os.path.join(odir, "chiSquared.out") if odir else ""
            if cs and os.path.exists(cs):
                tab = read_chisq_out(cs)
                ref_chi2 = sum(d["chi2"] for d in tab.values())
                ref_n = sum(d["n"] for d in tab.values())
                plain = sum(sc.segment_chi2.values()) if getattr(sc, "segment_chi2", None) else sc.chi2
                rel = abs(plain - ref_chi2) / max(abs(ref_chi2), 1e-300)
                msg = (f"init: archived {os.path.relpath(cs, directory)} chi2 {ref_chi2:.3f} (N {ref_n}) vs the seed's plain chi2 "
                       f"{plain:.3f} (N {sc.npts}): {100 * rel:.4f} %")
                if rel > 1e-3 or ref_n != sc.npts:
                    c.log("WARNING " + msg + " -- the seed does NOT reproduce the archived fit: check the flags (ignore_externals, "
                          "long wavelength, Brune), the .sav (positional names after a level change) and whether output/ is stale, "
                          "before building on it (skill section 2 item 2)")
                    c.note("WARNING at init: " + msg + "; the seed does not reproduce the archived fit")
                else:
                    c.log(msg + " -- reproduces")
        except Exception as err:              # a diagnostic must never stop an init
            c.log(f"init: archived chiSquared.out comparison skipped ({err})")
        c.log(f"init: baseline candidate {cid} objective {sc.objective:.3f} (chi2 {sc.chi2:.3f} + pen {sc.penalty:.3f}), "
              f"N={sc.npts}, sanity {rep.summary()}")
        c.note(f"campaign initialised on {os.path.basename(base_azr)}: baseline objective {sc.objective:.3f} "
               f"(chi2 {sc.chi2:.3f} + penalties {sc.penalty:.3f}, N={sc.npts}, {len(keys)} free); sanity {rep.summary()}")
        c.ledger.write_status_md()
        return c

    # -- operations ---------------------------------------------------------------------------

    def run_polish(self, max_nfev=200, shifts_free=None, note="warm-start continuation",
                   tighten=False, tighten_factor=10.0, cid=None, refreeze_every=0) -> FitResult:
        """Global polish of the incumbent on the sharded evaluator.

        ``cid`` (2026-10-05): polish THAT candidate instead of the best verified one -- for a start vector injected
        into the ledger (Campaign.verify_candidate) that is deliberately worse than the incumbent at the start, e.g. a
        solution transplanted from a data-only window fit (an archived campaign).

        With ``tighten``, amplitudes outside the policy's theta^2 caps are
        walked into the box in stages: every stage relaxes the caps by a
        factor ``tighten_factor`` less than the previous one (starting from
        the largest current excess) and polishes, so the fit follows the
        shrinking box instead of being clipped in one step.

        ``refreeze_every`` (2026-10-24; default 0 since 2026-10-25): with a data covariance that has model-scaled
        (``relvector``) terms and NO fixed reference, each stage runs as passes of at most this many evaluations with
        the covariance re-frozen between them (local.polish_refreeze).  On a large model that
        iteration oscillated, so it is off by default; campaigns now fix the
        relvector model at a reference (Campaign.init, ObjectiveSpec.cov_reference) and this option is ignored then.
        """
        if cid is None:
            cid, x, obj, keys = self.incumbent()
        else:
            x, keys = self.ledger.load(cid)
            sc = [r for r in self.ledger.scores(cid) if r.get("stage") == "verified"]
            obj = float(sc[-1]["objective"]) if sc else float("nan")
            self.log(f"polish starts from candidate {cid} (objective {obj:.3f}), not the incumbent")
        ev = self.sharded
        xs = take(x, keys, ev.keys, ev.x0)
        # None = leave the campaign policy alone.  A hard-coded False here silently
        # overrode a campaign whose policy said shifts_free=True, so a "free the
        # BandH shift" test polished with the shift excluded (twice).
        policy = self.policy if shifts_free is None else dataclasses.replace(self.policy, shifts_free=shifts_free)
        self.log(f"polish policy: shifts_free={policy.shifts_free}")
        if getattr(ev, "cov", None) is not None and ev.cov.relative:
            if ev.cov.has_reference:
                if refreeze_every and refreeze_every > 0:
                    self.log(f"polish: --refreeze-every {refreeze_every} ignored (the relvector model is fixed)")
                refreeze_every = 0
                self.log(f"polish: relvector covariance -> model fixed at the campaign reference "
                         f"({os.path.basename(self.objective.cov_reference)})")
            else:
                self.log(f"polish: relvector covariance -> " + (f"re-frozen every <= {refreeze_every} evaluations"
                                                                if refreeze_every and refreeze_every > 0
                                                                else "frozen at the start point throughout"))
        rid = self.ledger.start_round(self.structure_id, "polish", note + (" (tighten)" if tighten else ""))
        t0 = time.time()
        scales = [1.0]
        if tighten:
            gW = ev.wigner_amplitudes()
            excess = 1.0
            for j, k in enumerate(ev.keys):
                if k.kind == "width" and j in gW and gW[j] > 0:
                    e_seed = ev.level_energy.get(k.level, 0.0)
                    cap = policy.cap_for(e_seed)
                    excess = max(excess, (xs[j] / gW[j]) ** 2 / cap)
            n = int(np.ceil(np.log(excess) / np.log(tighten_factor))) if excess > 1 else 0
            scales = [s for s in (excess / tighten_factor ** i for i in range(1, n + 1)) if s > 1.0] + [1.0]
            self.log(f"tighten: largest theta^2 excess {excess:.3g} -> {len(scales)} stages, cap scales "
                     + ", ".join(f"{s:.3g}" for s in scales))
        fr = None
        x_cur = xs
        for stage, scale in enumerate(scales):
            b = default_bounds(ev, x_cur, policy, cap_scale=scale)
            fr = polish_refreeze(ev, x_cur, bounds=b, policy=policy, max_nfev=max_nfev, every=refreeze_every,
                                 log=self.log)
            self.log(f"polish stage {stage+1}/{len(scales)} (cap scale {scale:.3g}): {fr.start_objective:.3f} -> "
                     f"{fr.objective:.3f} ({fr.status}; nfev {fr.nfev}; at bound {len(fr.bounds_hit)}; "
                     f"{(time.time()-t0)/60:.1f} min)")
            if fr.bounds_hit:
                self.log("  at bound: " + ", ".join(f"{ev.keys[j]}={fr.x[j]:.4g}" for j in fr.bounds_hit[:10]))
            # every stage is recorded (sharded score) so it can be inspected
            scs = ev.score(fr.x)
            scid = self.ledger.add_candidate(self.structure_id, fr.x, ev.keys,
                                             move=dict(kind="polish_stage", stage=stage + 1, cap_scale=scale, parent=cid),
                                             parent=cid, status="stage")
            self.ledger.record(scid, scs, "global")
            from .diagnose import dataset_table, format_dataset_table
            self.log("  worst datasets:\n" + format_dataset_table(dataset_table(ev, scs), top=6))
            x_cur = fr.x
            if fr.status == "crashed":
                break
        fr.start_objective = float(obj)
        stats = dict(start=obj, end=fr.objective, status=fr.status, nfev=fr.nfev, stages=len(scales))
        if fr.objective < obj:
            ncid, sc = self.verify_candidate(fr.x, ev.keys, move=dict(kind="polish", parent=cid), parent=cid)
            stats["verified"] = sc.objective
            self.note(f"warm-start polish of candidate {cid}: {obj:.3f} -> {sc.objective:.3f} (verified, candidate {ncid}); "
                      f"{fr.status}, {fr.nfev} evaluations; parameters at a bound: "
                      + (", ".join(str(ev.keys[j]) for j in fr.bounds_hit[:8]) or "none"))
        self.ledger.finish_round(rid, **stats)
        self.ledger.write_status_md()
        return fr

    def run_sign_round(self, jpi: Optional[Sequence[str]] = None, note="") -> dict:
        cid, x, obj, keys = self.incumbent()
        ev = self.sharded
        xs = take(x, keys, ev.keys, ev.x0)
        cfg = dataclasses.replace(self.search, jpi=list(jpi) if jpi else self.search.jpi)
        # flips known to revert or to be worse on this structure are skipped
        skip = [k for k, outs in self.ledger.flip_outcomes(self.structure_id).items()
                if outs and outs[-1][0] in ("reverts", "worse", "unphysical", "crashed")]
        rid = self.ledger.start_round(self.structure_id, "sign", note or (",".join(cfg.jpi) if cfg.jpi else "all"))
        x_best, obj_best, stat, trials = sign_round(ev, xs, obj, cfg, self.policy, skip=skip, log=self.log)
        return self._record_round(rid, cid, ev, obj, x_best, obj_best, stat, trials, kind="sign round")

    def _record_round(self, rid, cid, ev, obj, x_best, obj_best, stat, trials, kind="sign round") -> dict:
        """Ledger bookkeeping shared by sign and crossover rounds."""
        for tr in trials:
            tcid = self.ledger.add_candidate(self.structure_id, tr.x, ev.keys, move=tr.move.to_json(),
                                            parent=cid, status=tr.outcome)
            if tr.stage2 is not None:
                self.ledger.record(tcid, tr.stage2, "global")
            for k in tr.move.targets:
                self.ledger.record_flip(self.structure_id, k, tr.outcome, tr.gain, tcid, rid)
        stats = dict(proposed=stat.n_proposed, screened=stat.n_screened, improving_screen=stat.n_improving_screen,
                     polished=stat.n_polished, accepted=stat.n_accepted, reverted=stat.n_reverted,
                     kept=stat.n_union_kept, before=stat.best_before, after=stat.best_after, wall_min=stat.wall / 60,
                     control=stat.control_after, screen=stat.screen)
        if obj_best < obj:
            ncid, sc = self.verify_candidate(x_best, ev.keys, move=dict(kind=kind.replace(" ", "_"), round=rid), parent=cid)
            stats["verified"] = sc.objective
            stat.verified_after = sc.objective
        self.ledger.finish_round(rid, **stats)
        self.note(f"{kind} #{rid} ({stats.get('proposed')} proposals, {stats['improving_screen']} improving on screen, "
                  f"{stats['accepted']} accepted, {stats['reverted']} reverted, {stats['kept']} kept): "
                  f"{obj:.3f} -> {stats.get('verified', obj_best):.3f}"
                  + (" (verified)" if "verified" in stats else " (no improvement)")
                  + "; trials: " + "; ".join(f"{t['move']}: {t['outcome']} {t['gain']:+.0f}" for t in stat.trials))
        self.ledger.write_status_md()
        return stats

    def run_crossover(self, donor_cid, donor_dir=None, note="") -> dict:
        """Basin crossover with a donor candidate (of this campaign, or of another
        campaign directory with the same base structure): the donor's sign
        pattern applied per channel / group / all, every proposal polished."""
        cid, x, obj, keys = self.incumbent()
        ev = self.sharded
        src = self.ledger if not donor_dir else Ledger(os.path.join(os.path.abspath(donor_dir), "ledger.sqlite"))
        xd, kd = src.load(int(donor_cid))
        missing = set(keys) - set(kd)
        bad = [k for k in missing if k.kind in ("energy", "width")]
        if bad:
            raise ValueError(f"donor lacks {len(bad)} R-matrix parameters of this model: "
                             + ", ".join(str(k) for k in bad[:5]))
        if missing:
            self.log(f"  donor lacks {len(missing)} norm/shift parameters of this model (e.g. "
                     f"{', '.join(str(k) for k in list(missing)[:3])}): taken from the incumbent")
        xs = take(x, keys, ev.keys, ev.x0)
        xds = take(take(xd, kd, keys, x), keys, ev.keys, ev.x0)   # donor pattern, incumbent's norms where absent
        label = f"candidate {donor_cid}" + (f" of {os.path.basename(os.path.abspath(donor_dir))}" if donor_dir else "")
        rid = self.ledger.start_round(self.structure_id, "crossover", note or label)
        x_best, obj_best, stat, trials = crossover_round(ev, xs, obj, xds, self.search, self.policy,
                                                         log=self.log, label=label)
        return self._record_round(rid, cid, ev, obj, x_best, obj_best, stat, trials, kind="crossover round")

    def run_population(self, restarts=4, n_flips=(6, 8, 10, 12), per_group=3, top_moves=40, jitter=0.2,
                       nfev=600, rounds_max=8, gain_tol=5e-4, passes=2, seed=None, note="") -> dict:
        """Restarts + crossover: the basin-crossing stage (2026-11-08, from a benchmark post-mortem).

        The sign rounds judge single and paired flips on short screens; a benchmark post-mortem
        showed that the way from a basin 13 % above the reference to it was a coordinated 17-flip pattern whose unpolished
        objective was three times the incumbent's -- no single, pair, channel or group flip gained more than a few units even with a
        600-evaluation polish, while a population of six sign-perturbed restarts of the incumbent, each polished to convergence,
        plus crossover between them reached 5326 in one pass (pop3).  So: ``restarts`` sign-perturbed copies of the incumbent
        (a random set of ``n_flips`` flips among the ``top_moves`` highest-exposure amplitudes, at most ``per_group`` per J^pi,
        widths jittered by ``jitter``), each polished sign-locked once and then released in rounds of ``nfev`` until the gain is
        below ``gain_tol`` (at most ``rounds_max``); the control (the incumbent polished the same way) joins them; then up to
        ``passes`` crossover passes from the best member with every other member as donor (moves.crossover_moves: the donor's
        pattern per channel / per J^pi group / all), each pattern polished once, the winner polished to convergence.  Everything
        is judged on polished objectives.  Every member is a ledger candidate (status "population"); only the final best is
        verified and becomes the incumbent if it beats it."""
        from .moves import exposure, single_flips, random_flip_set, apply, crossover_moves
        from .keys import sign_distance
        cid, x, obj, keys = self.incumbent()
        ev = self.sharded
        xs = take(x, keys, ev.keys, ev.x0)
        fixed = ev.fixed_amplitudes
        pos = {k: j for j, k in enumerate(ev.keys)}
        rid = self.ledger.start_round(self.structure_id, "population",
                                      note or f"{restarts} restarts, {passes} crossover passes")
        rng = np.random.default_rng(rid if seed is None else seed)
        t0 = time.time()
        n_pol = [0]

        def rounds(x0, label, first_bounds=None, max_rounds=rounds_max):
            xx = np.asarray(x0, float)
            prev = float(ev.objective(xx))
            for i in range(max_rounds):
                b = first_bounds if (i == 0 and first_bounds is not None) else default_bounds(ev, xx, self.policy)
                f = polish(ev, xx, bounds=b, policy=self.policy, max_nfev=nfev)
                n_pol[0] += 1
                gain = (prev - f.objective) / prev if prev else 0.0
                self.log(f"    {label} round {i+1}: {prev:.3f} -> {f.objective:.3f} ({100*gain:.3f} %, {f.status}, "
                         f"{(time.time()-t0)/60:.0f} min)")
                xx, prev = f.x, float(f.objective)
                if f.status == "crashed" or gain < gain_tol:
                    break
            return xx, prev

        def member(name, xv, ov, kind):
            mcid = self.ledger.add_candidate(self.structure_id, xv, ev.keys, move=dict(kind=kind, round=rid, name=name),
                                             parent=cid, status="population")
            self.ledger.record(mcid, ev.score(xv), "global")
            return dict(name=name, x=xv, obj=ov, cid=mcid)

        self.log(f"population round #{rid}: incumbent {cid} {obj:.3f}; {restarts} restarts of {n_flips} flips (<= {per_group} per "
                 f"J^pi group from the {top_moves} highest-exposure amplitudes, jitter {jitter}), polish {nfev} x <= {rounds_max} "
                 f"until < {100*gain_tol:.2f} %, then <= {passes} crossover passes")
        pop = []
        xc, oc = rounds(xs, "control")
        pop.append(member("control", xc, oc, "population_control"))
        expo = exposure(ev, xs)
        moves = single_flips(xs, ev.keys, exposure_map=expo, fixed=fixed)[:top_moves]
        for i in range(int(restarts)):
            n = int(n_flips) if isinstance(n_flips, (int, np.integer)) else int(rng.choice(list(n_flips)))
            chosen = random_flip_set(moves, rng, n, per_group=per_group)
            xp = xs.copy()
            for mv in chosen:
                xp = apply(xp, ev.keys, mv)
            if jitter > 0:
                for j, k in enumerate(ev.keys):
                    if k.kind == "width" and xp[j] != 0:
                        xp[j] *= rng.uniform(1 - jitter, 1 + jitter)
            lock = {k: (1 if xp[pos[k]] > 0 else -1) for mv in chosen for k in mv.targets}
            b1 = default_bounds(ev, xs, self.policy, sign_lock=lock)
            xp = b1.clip(xp)
            self.log(f"  restart {i+1}: {len(chosen)} flips " + ", ".join(str(m.targets[0]) for m in chosen)
                     + f"; raw {ev.objective(xp):.1f}")
            x1, _ = rounds(xp, f"restart {i+1} locked", first_bounds=b1, max_rounds=1)
            x2, o2 = rounds(x1, f"restart {i+1}")
            d = sign_distance(x2, xs, ev.keys, fixed)[0]
            self.log(f"  restart {i+1}: {o2:.3f}, sign distance to the incumbent {d}" + (" (same basin)" if d == 0 else ""))
            pop.append(member(f"restart {i+1}", x2, o2, "population_restart"))
        self.log("  population: " + "; ".join(f"{m['name']} {m['obj']:.3f}" for m in sorted(pop, key=lambda m: m["obj"])))
        best = min(pop, key=lambda m: m["obj"])
        n_moves, source = 0, best["name"]
        for pas in range(1, int(passes) + 1):
            cand, seen = [], []
            for d in pop:
                if d is best:
                    continue
                # two members that polished into the same sign pattern donate the same moves: polish them once
                dup = next((s for s in seen if sign_distance(d["x"], s["x"], ev.keys, fixed)[0] == 0), None)
                if dup is not None:
                    self.log(f"  crossover pass {pas}: {d['name']} has the sign pattern of {dup['name']}, not polished again")
                    continue
                seen.append(d)
                for mv in crossover_moves(best["x"], d["x"], ev.keys, fixed=fixed, label=d["name"]):
                    xm = apply(best["x"], ev.keys, mv)
                    xm, om = rounds(xm, mv.description, max_rounds=1)
                    n_moves += 1
                    cand.append((om, xm, mv.description))
            if not cand:
                self.log(f"  crossover pass {pas}: no member differs from {best['name']} in sign pattern")
                break
            om, xm, desc = min(cand, key=lambda t: t[0])
            self.log(f"  crossover pass {pas} from {best['name']} {best['obj']:.3f}: {len(cand)} patterns polished, "
                     f"best '{desc}' {om:.3f}")
            if om >= best["obj"]:
                break
            xm, om = rounds(xm, f"crossover pass {pas} winner")
            best = member(f"crossover pass {pas}: {desc}", xm, om, "population_crossover")
            pop.append(best)
            source = best["name"]
        stats = dict(before=obj, after=best["obj"], restarts=restarts, passes=passes, polishes=n_pol[0],
                     crossover_moves=n_moves, members={m["name"]: m["obj"] for m in pop}, best=source,
                     wall_min=(time.time() - t0) / 60)
        improved = best["obj"] < obj
        if improved:
            ncid, sc = self.verify_candidate(best["x"], ev.keys, move=dict(kind="population", round=rid, source=source),
                                             parent=cid)
            stats["verified"] = sc.objective
            improved = sc.objective < obj
        self.ledger.finish_round(rid, **stats)
        self.note(f"population round #{rid} ({restarts} restarts, {n_moves} crossover patterns, {n_pol[0]} polishes of {nfev}, "
                  f"{stats['wall_min']:.0f} min): {obj:.3f} -> {stats.get('verified', best['obj']):.3f} via {source}"
                  + (" (verified)" if "verified" in stats else " (no improvement)")
                  + "; members: " + "; ".join(f"{m['name']} {m['obj']:.1f}" for m in sorted(pop, key=lambda m: m["obj"])))
        self.ledger.write_status_md()
        stats["improved"] = improved
        return stats

    def run_recall(self, n: int = 2, min_width_kev: float = 10.0, levels: Optional[Sequence[str]] = None,
                   global_nfev: int = 150, pinned_nfev: int = 60, local_release: float = 0.0, local_nfev: int = 80,
                   offset_kev: float = 10.0, note: str = "") -> dict:
        """RECALL test of the level-add procedure on levels the model already has (the evaluator 2026-10-15: "go through the
        resonances that we have in the fit ... make sure that they wouldn't be excluded by our criteria"; worked example
        an archived campaign).  The ``n`` broadest physical levels (below the policy's bg_energy, physical width >=
        ``min_width_kev``) are each removed from the incumbent -- a child campaign in ``recall/<tag>/`` on the baked
        incumbent with that level switched off -- and the SAME level-add procedure STRUCTURE_ADD uses (window scan over
        every J^pi, pinned + released polish, same-budget control, classification) is asked to find a level at the
        true energy + ``offset_kev`` with no J^pi given.  A level is re-found when the outcome is accept(-pending), the
        J^pi is the true one, the fitted energy is within max(50 keV, 3 Gamma) and the width within a factor 3; a level
        that meets the last three but is rejected by the criteria is "located" (result key) and still not re-found.
        ``levels`` = explicit key labels ("3/2-#1") instead of the automatic choice.  Returns dict(results, n, n_pass,
        n_located, fraction, passed_levels); ledger round "recall"; the child campaigns are kept for inspection."""
        from .structure import LevelCandidate, deactivate_level_variant
        from .backend import AzrModel
        cid, x, obj, keys = self.incumbent()
        ev = self.sharded
        xs = take(x, keys, ev.keys, ev.x0)
        rid = self.ledger.start_round(self.structure_id, "recall", note or f"{n} broadest known levels")
        t0 = time.time()
        # physical width of every level from the incumbent: Gamma_c = (gamma_c / gamma_W,c)^2 Gamma_W,c (eV)
        gWa, gWw = ev.wigner_amplitudes(), ev.wigner_widths(xs)
        width_ev: Dict[str, float] = {}
        energy: Dict[str, float] = {}
        theta2_max: Dict[str, float] = {}    # the incumbent's largest theta^2 per level: a level the model itself holds AT
        for j, k in enumerate(ev.keys):      # the policy cap is rejected as "unphysical" by the criteria (13N 1/2+ 2.368, 2026-11-09)
            if k.kind == "width" and gWa.get(j, 0) > 0 and j in gWw:
                t2 = (xs[j] / gWa[j]) ** 2
                width_ev[k.level] = width_ev.get(k.level, 0.0) + t2 * gWw[j]
                theta2_max[k.level] = max(theta2_max.get(k.level, 0.0), float(t2))
            elif k.kind == "energy":
                energy[k.level] = float(xs[j])
        mdl = AzrModel.from_file(self.base_azr)
        count: Dict[str, int] = {}
        for lv in mdl.active_levels:               # the keys' labels: J^pi#ordinal over the ACTIVE levels in file order
            count[lv.jpi] = count.get(lv.jpi, 0) + 1
            energy.setdefault(f"{lv.jpi}#{count[lv.jpi]}", float(lv.energy))
        table = sorted(((width_ev.get(l, 0.0), l) for l in energy if energy[l] < self.policy.bg_energy), reverse=True)
        if levels:
            chosen = [(width_ev.get(l, 0.0), l) for l in levels if l in energy]
        else:
            chosen = [(w, l) for w, l in table if w >= 1e3 * float(min_width_kev)][:int(n)]
        self.log(f"recall: {len(chosen)} level(s) chosen of {len(table)} below bg_energy {self.policy.bg_energy:g} MeV: "
                 + ", ".join(f"{l} {energy[l]:.4f} MeV Gamma {w / 1e3:.1f} keV" for w, l in chosen))
        if not chosen:
            self.ledger.finish_round(rid, n=0, n_pass=0, fraction=1.0, results="[]", wall_min=(time.time() - t0) / 60)
            self.note(f"recall round #{rid}: no physical level with Gamma >= {min_width_kev:g} keV below bg_energy -- nothing to re-find")
            return dict(results=[], n=0, n_pass=0, fraction=1.0, passed_levels=[])
        from .model import link_data_dirs
        rdir = os.path.join(self.dir, "recall"); os.makedirs(rdir, exist_ok=True)
        link_data_dirs(self.base_azr, rdir)         # the baked incumbent and the children keep the model's relative data paths
        azr_inc, sav_inc, sc = self.export(path=os.path.join(rdir, "incumbent.azr"), note="recall")
        self.close()                                # the children run their own shards; the parent's come back lazily
        obj_cfg = dict(self.cfg.get("objective") or {})
        covref = obj_cfg.pop("cov_reference", None)
        if covref and not os.path.isabs(covref):
            covref = os.path.join(self.dir, covref)
        results = []
        for w, label in chosen:
            jpi, e_true = label.split("#")[0], energy[label]
            tag = label.replace("/", "").replace("#", "_") + f"_E{e_true:.4f}"
            cdir = os.path.join(rdir, tag)
            os.makedirs(os.path.join(cdir, "output"), exist_ok=True)
            link_data_dirs(self.base_azr, cdir)
            res = dict(level=label, jpi=jpi, energy=e_true, gamma_kev=w / 1e3, dir=os.path.relpath(cdir, self.dir))
            try:
                var = deactivate_level_variant(azr_inc, os.path.join(cdir, "base.azr"), jpi, e_true,
                                               output_dir=os.path.join(cdir, "output"))
                # with the incumbent's .sav: the baked .azr carries the levels but every free norm and shift at its nominal
                # value (bake writes those to the .sav only; campaigns initialised from the .azr alone started at nominal norms)
                child = Campaign.init(cdir, var["path"], sav=sav_inc, flags=self.flags, nshards=self.nshards, policy=self.cfg.get("policy"),
                                      search=self.cfg.get("search"), objective=obj_cfg or None, cov_reference=covref,
                                      note=f"recall child of {self.dir}: {label} removed")
                try:
                    c0 = child.incumbent()[2]
                    cand = LevelCandidate(energy=round(e_true + offset_kev * 1e-3, 4), jpi=None, channels=[], source="recall",
                                          note=f"recall: the known {jpi} level at {e_true:.4f} MeV (Gamma {w / 1e3:.1f} keV) removed; "
                                               f"can the procedure find it again with no J^pi given?", energy_fixed=True)
                    st = child.run_level_add(cand, jpis=None, global_nfev=global_nfev, pinned_nfev=pinned_nfev,
                                             local_release=local_release, local_nfev=local_nfev, note="recall")
                finally:
                    child.close()
                out = st.get("outcome", "?")
                jf = st.get("jpi")
                ef = st.get("energy_fit", cand.energy)
                gf = st.get("level_gamma_ev")
                why, located = [], []
                if jf != jpi:
                    located.append(f"J^pi {jf} != {jpi}")
                if abs(float(ef) - e_true) > max(0.050, 3.0 * w * 1e-6):
                    located.append(f"energy {float(ef):.4f} vs {e_true:.4f}")
                if gf and w > 0 and not (w / 3.0 <= float(gf) <= 3.0 * w):
                    located.append(f"Gamma {float(gf) / 1e3:.1f} keV vs {w / 1e3:.1f}")
                # 2026-11-09 (13N drive test 1718060): the 1/2+ at 2.368 MeV was located exactly (J^pi, energy -0.4 keV, Gamma
                # 122 vs 133 keV, gain 8.6e6) and still classified "unphysical", because its s-wave amplitude ends ON the
                # theta^2 cap -- where the incumbent's own level sits.  That is a verdict on the POLICY (cap too tight for a
                # level the model already carries), not on the search; it is reported as such and still counts as not re-found,
                # since STRUCTURE_ADD would reject the same level.
                cap_true = self.policy.cap_for(e_true)
                inc_at_cap = theta2_max.get(label, 0.0) >= cap_true * (1 - 1e-3)
                if out not in ("accept", "accept-pending"):
                    if out == "unphysical" and inc_at_cap:
                        why.append(f"outcome unphysical: amplitude at the theta2 cap {cap_true:g}, where the incumbent's own level "
                                   f"sits too (theta2 {theta2_max.get(label, 0.0):.3g}) -- the policy cap rejects a level the model "
                                   f"carries; raise theta2_cap / theta2_cap_windows in spec bounds, not the procedure")
                    else:
                        why.append(f"outcome {out}")
                why += located
                res.update(removed_objective=c0, outcome=out, jpi_found=jf, energy_found=float(ef),
                           gamma_found_kev=(float(gf) / 1e3 if gf else None), delta=st.get("delta"), scan=st.get("scan", ""),
                           at_bound=st.get("at_bound", ""), incumbent_theta2=theta2_max.get(label), cap=cap_true,
                           located=not located, passed=not why, why="; ".join(why))
            except Exception as err:  # noqa: BLE001 -- one failed child must not end the round
                res.update(outcome="error", passed=False, why=f"{type(err).__name__}: {err}")
                self.log(f"recall of {label} failed: {err}")
            results.append(res)
            self.log(f"recall {label}: " + ("RE-FOUND" if res["passed"]
                                            else ("LOCATED but rejected by the criteria (" if res.get("located") else "NOT re-found (")
                                            + res["why"] + ")")
                     + f" -- outcome {res.get('outcome')} as {res.get('jpi_found')} @ {res.get('energy_found', float('nan')):.4f}")
        n_pass = sum(1 for r in results if r["passed"])
        n_located = sum(1 for r in results if r.get("located"))
        frac = n_pass / len(results)
        self.ledger.finish_round(rid, n=len(results), n_pass=n_pass, fraction=frac, results=json.dumps(results),
                                 wall_min=(time.time() - t0) / 60)
        self.note(f"recall round #{rid} ({len(results)} known levels removed and searched for, {(time.time() - t0) / 60:.0f} min): "
                  f"{n_pass} of {len(results)} re-found ({n_located} located: right J^pi, energy and width) -- "
                  + "; ".join(f"{r['level']} {r['energy']:.3f} MeV Gamma {r['gamma_kev']:.0f} keV: "
                              + ("re-found" if r["passed"] else f"NOT ({r['why']})") for r in results)
                  + "; children in recall/")
        self.ledger.write_status_md()
        return dict(results=results, n=len(results), n_pass=n_pass, n_located=n_located, fraction=frac,
                    passed_levels=[r["level"] for r in results if r["passed"]])

    def run_levelscan(self, ex: Optional[Sequence[float]] = None, top: int = 3, gamma_kev: Optional[float] = 10.0,
                      scales: Sequence[float] = (0.3, 0.6, 1.2, 2.5, 5.0, 10.0), jpis: Optional[Sequence[str]] = None,
                      gn_top: int = 20, max_signs: int = 32, refit: int = 0,
                      ex_offsets: Sequence[float] = (-0.04, 0.0, 0.04), differential_window_kev: float = 250.0,
                      grid: Optional[Sequence[float]] = (5.0, 20.0), stage_scales: Optional[Sequence[float]] = (1.2, 5.0),
                      refine_top: int = 3, merge_kev: float = 60.0, gamma_bounds_kev: Sequence[float] = (10.0, 1000.0),
                      seg_match: Optional[str] = None, rank: str = "chi2",
                      ex_range: Optional[Sequence[float]] = None, note="") -> dict:
        """Residual-driven level search (plan Phase 7): candidate energies from
        the incumbent's residual runs (or the given ``ex`` list), a template
        screen of every J^pi x width partition x sign pattern at each, and
        optionally a refit of the best ``refit`` candidates via run_level_add.

        ``gamma_kev=None`` (CLI ``--gamma auto``, 2026-09-24): the template
        total width of each screen is the cluster's own residual-run width
        (levelscan.cluster_width) clamped into ``gamma_bounds_kev``, so the
        width grid ``scales`` x Gamma is centred on what the data show instead
        of a guess.

        2026-09-24: ``seg_match`` (regex on data-file names: only those segments'
        residuals make candidates), ``rank`` ("chi2" | "density", see
        levelscan.candidate_energies) and ``ex_range`` (lo, hi MeV: residual
        candidates outside are dropped BEFORE the ``top`` cut -- replaces the
        stage-1 / pick_window_energies / stage-2 double screen)."""
        from .levelscan import (candidate_energies, template_variant, template_screen, best_per_jpi, to_level_candidate,
                                incumbent_on_variant, auto_gamma_kev)
        from .model import single_variant
        cid, x, obj, keys = self.incumbent()
        gamma_label = "auto" if gamma_kev is None else f"{gamma_kev:g} keV"
        rid = self.ledger.start_round(self.structure_id, "levelscan", note or f"top {top}, Gamma {gamma_label}")
        t0 = time.time()
        cands = []
        if ex:
            cands = [dict(ex=float(e), score=float("nan"), npts=0, segments=[]) for e in ex]
        else:
            src = single_variant(self.base_azr, "levelscan", self.vdir)
            with SingleSession(src, cwd=self.vdir, flags=self.flags, name="levelscan", objective=self.objective) as ev1:
                xs = take(x, keys, ev1.keys, ev1.x0)
                found = candidate_energies(ev1, xs, top=max(4 * top, 40), include_match=seg_match, rank=rank)
            if ex_range:
                lo, hi = float(ex_range[0]), float(ex_range[1])
                nall = len(found)
                found = [c for c in found if lo <= c.ex <= hi]
                self.log(f"  ex_range {lo:g}-{hi:g} MeV: {len(found)} of {nall} residual clusters kept")
            if seg_match:
                self.log(f"  residual search restricted to segments whose data file matches /{seg_match}/")
            found = found[:max(4 * top, 20)]
            cands = [c.to_json() for c in found]
            self.log(f"  residual-run candidates, ranked by {rank} (Ex MeV, chi2 in runs, chi2/MeV, points, run width [differential / total] and span keV, segments): "
                     + "; ".join(f"{c['ex']:.3f}: {c['score']:.0f} ({c.get('density', 0.0):.3g}/MeV, {c['npts']} pts, width {c.get('width', 0.0) * 1e3:.0f} / "
                                 f"{c.get('width_total', 0.0) * 1e3:.0f}, span {c.get('span', 0.0) * 1e3:.0f} keV, "
                                 f"segs {c['segments']})" for c in cands))
            # 2026-09-14: three residual candidates 47-86 keV apart all screened (with the
            # +-40 keV offsets) to the same 10.41 MeV level, 40 min each.  Keep only the
            # strongest of any group closer than `merge_kev`; the offsets cover the rest.
            if merge_kev and merge_kev > 0:
                kept = []
                for c in cands:                       # already strongest-first
                    if all(abs(c["ex"] - k["ex"]) > merge_kev * 1e-3 for k in kept):
                        kept.append(c)
                if len(kept) < len(cands):
                    self.log(f"  {len(cands) - len(kept)} candidate(s) within {merge_kev:g} keV of a stronger one dropped; screening "
                             + ", ".join(f"{c['ex']:.3f}" for c in kept))
                cands = kept
        # 2026-11-04: the model's own J^pi groups unless --jpi; raises on a J^pi no pair can couple to
        jpis = resolve_jpis(self.base_azr, jpis, log=self.log)
        results = []
        for c in cands[:top]:
            tag = f"E{c['ex']:.3f}".replace(".", "p")
            wdir = os.path.join(self.vdir, f"levelscan_{tag}")
            g_kev = auto_gamma_kev(c, gamma_bounds_kev) if gamma_kev is None else float(gamma_kev)
            if gamma_kev is None:
                self.log(f"  template width at Ex {c['ex']:.4f}: differential run width {c.get('width', 0.0) * 1e3:.0f} keV "
                         f"(total-segment runs {c.get('width_total', 0.0) * 1e3:.0f}, span {c.get('span', 0.0) * 1e3:.0f} keV) "
                         f"-> Gamma_total {g_kev:g} keV, grid "
                         f"{min(scales) * g_kev:g}-{max(scales) * g_kev:g} keV (bounds {gamma_bounds_kev[0]:g}-{gamma_bounds_kev[1]:g})")
            os.makedirs(wdir, exist_ok=True)
            var = template_variant(self.base_azr, os.path.join(wdir, f"templates_{tag}.azr"), c["ex"],
                                   jpis=jpis, output_dir=os.path.join(wdir, "output"))
            if var.get("skipped"):
                self.log(f"  J^pi not screened at Ex {c['ex']:.4f}: no particle channel with Gamma_W >= "
                         f"{var['min_open_ev']:g} eV (" + ", ".join(f"{j} {g:.2g} eV" for j, g in var["skipped"].items()) + ")")
            if grid:
                nchg = self.set_grid_setting(var["path"], float(grid[0]), float(grid[1]))
                if nchg:
                    self.log(f"  screening grid setting {float(grid[0]):g}/{float(grid[1]):g} applied to {nchg} "
                             f"<targetInt> line(s) of the screen variant (the converged 20/50 stays for anything reported)")
            ev2 = Sharded(var["path"], os.path.join(wdir, "shards"), nshards=self.nshards, flags=self.flags,
                          objective=self.objective)
            try:
                x_inc = incumbent_on_variant(ev2, x, keys, var)
                # Phase 7 redesign (2026-09-12): rank the GN'd candidates on the DIFFERENTIAL
                # segments near the candidate energy -- the total objective barely
                # distinguishes J^pi (a new 11/2- once won by 181k on a model 1.5M from its
                # optimum).  The Ex ranges come from a single session of the base model.
                score_keys = None
                if differential_window_kev and differential_window_kev > 0:
                    try:
                        from .levelscan import scoring_segments
                        score_keys = scoring_segments(ev2, var["ex"], half_width=differential_window_kev * 1e-3,
                                                      differential_only=True, ex_ranges=self._segment_ex_ranges())
                        self.log(f"  local re-rank on {len(score_keys)} differential segments within "
                                 f"+-{differential_window_kev:.0f} keV of Ex {var['ex']:.4f}")
                    except Exception as err:
                        self.log(f"  differential re-rank unavailable ({err}); ranking on the total objective")
                rows = template_screen(ev2, x_inc, var, gamma_total_ev=g_kev * 1e3, scales=scales,
                                       max_signs=max_signs, gn_top=gn_top, policy=self.policy, jpis=jpis,
                                       ex_offsets=ex_offsets, score_keys=score_keys,
                                       stage_scales=stage_scales, refine_top=refine_top, log=self.log)
            finally:
                ev2.close()
            best = best_per_jpi(rows)                                               # the screen's order (local re-rank when active)
            best_global = best_per_jpi(sorted(rows, key=lambda r: -r["gain"]))     # by the global GN gain
            # 2026-09-14 (gate C'): the GLOBAL GN gain is the existence statistic -- it
            # predicted the accepted 1/2+ 9.862 to 16% while the local re-rank put that
            # level 5th (its gain is n-total and every (a,a) angle, not the +-250 keV
            # window).  The local order stays in the note for J^pi discrimination; the
            # refit takes the best by global gain.
            results.append(dict(ex=c["ex"], score=c.get("score"), run_width=c.get("width"), run_span=c.get("span"),
                                gamma_kev=g_kev, best_global=[
                dict(jpi=r["jpi"], desc=r["desc"], signs=r["signs"], gamma_total_ev=r["gamma_total_ev"], ex=r["ex"],
                     frozen=r["frozen"], gn=r["gn"], gain=r["gain"], local_gain=r.get("local_gain"),
                     widths_ev=r["widths_ev"], channels=r["channels"])
                for r in best_global[:6]], best=[
                dict(jpi=r["jpi"], desc=r["desc"], signs=r["signs"], gamma_total_ev=r["gamma_total_ev"], ex=r["ex"],
                     frozen=r["frozen"], gn=r["gn"], gain=r["gain"], widths_ev=r["widths_ev"], channels=r["channels"])
                for r in best[:6]]))
            self.note(f"levelscan at Ex {c['ex']:.4f} MeV (run chi2 {c.get('score', float('nan')):.0f}, run width "
                      f"{(c.get('width') or 0.0) * 1e3:.0f} keV, template Gamma {g_kev:g} keV): best per J^pi by GLOBAL gain -- "
                      + "; ".join(f"{r['jpi']}: {r['desc']} signs {r['signs']} @{r['ex']:.3f} gain {r['gain']:+.0f}"
                                  + (f" (local {r['local_gain']:+.0f})" if r.get('local_gain') is not None else "") for r in best_global[:6])
                      + (" | local order -- " + "; ".join(f"{r['jpi']} {r['local_gain']:+.0f}" for r in best[:4] if r.get('local_gain') is not None)
                         if any(r.get('local_gain') is not None for r in best[:4]) else ""))
        self.ledger.finish_round(rid, candidates=cands, results=results, wall_min=(time.time() - t0) / 60)
        self.ledger.write_status_md()
        if refit:
            picks = sorted([(r["best_global"][0]["gain"], r["best_global"][0]["ex"], r["best_global"][0])
                            for r in results if r.get("best_global")], reverse=True)[:refit]
            for gain, exv, row in picks:
                row = dict(row, gain=gain)
                cand = to_level_candidate(row, exv)
                self.log(f"  refit of the levelscan candidate {cand.jpi} @ {exv:.4f} MeV ({row['desc']}, screen gain {gain:+.0f})")
                self.run_level_add(cand, jpis=[cand.jpi], note="from levelscan")
        return dict(candidates=cands, results=results)

    def review_figure(self, x_before, keys_before, x_after, keys_after, ex: float, name: str, title: str = "",
                      label_before: str = "incumbent", label_after: str = "candidate", half_width: float = 0.06,
                      azr: Optional[str] = None):
        """Data with the two calculations around ``ex`` (excitation functions and
        angular distributions), written to figures/<name>.png|.pdf via a fresh
        single session on the base model."""
        from .levelscan import review_figure
        from .model import single_variant
        fdir = os.path.join(self.dir, "figures"); os.makedirs(fdir, exist_ok=True)
        path = os.path.join(fdir, f"{name}.png")
        try:
            src = single_variant(azr or self.base_azr, "review", self.vdir)
            with SingleSession(src, cwd=self.vdir, flags=self.flags, name="review", objective=self.objective) as ev:
                xb = take(x_before, keys_before, ev.keys, ev.x0); xa = take(x_after, keys_after, ev.keys, ev.x0)
                known = set(keys_before)
                for j, k in enumerate(ev.keys):          # parameters the incumbent does not have (a new level): off
                    if k not in known and k.kind == "width":
                        xb[j] = 0.0
                out = review_figure(ev, xb, xa, ex, path, half_width=half_width, title=title,
                                    label_before=label_before, label_after=label_after)
            self.log(f"  review figure: {out['path']} ({out['panels']} panels)")
            return out["path"] if out["panels"] else None
        except Exception as err:                 # a figure must never stop a campaign
            self.log(f"  (review figure failed: {err})")
            return None

    def _segment_ex_ranges(self) -> dict:
        """``{segment key: (Ex_lo, Ex_hi)}`` of the base model, from one single
        session (cached for the campaign's lifetime)."""
        if getattr(self, "_ex_ranges_cache", None) is None:
            from .levelscan import segment_ex_ranges
            azr = single_variant(self.base_azr, "exranges", self.vdir)
            with SingleSession(azr, cwd=self.vdir, flags=self.flags, name="exranges",
                               objective=self.objective) as ev:
                self._ex_ranges_cache = segment_ex_ranges(ev)
        return self._ex_ranges_cache

    @staticmethod
    def set_grid_setting(azr_path: str, mult: float, ppw: float) -> int:
        """Set the sub-point grid of every folded data set for a SCREENING variant (backend versions below).  Returns the number changed."""
        return (Campaign.set_grid_setting_p3(azr_path, mult, ppw) if PYAZR3 else Campaign.set_grid_setting_p2(azr_path, mult, ppw))

    @staticmethod
    def set_grid_setting_p2(azr_path: str, mult: float, ppw: float) -> int:
        """Rewrite the two trailing numbers of every active <targetInt> line --
        resonanceWidthMultiplier and pointsPerWidth (TargetEffect.cpp:92-102) --
        for a SCREENING variant.  Measured on a large target-integrated segment:
        5/20 is 6.7x cheaper than the converged 20/50 for a 0.03% objective
        shift; 2/x is NOT safe (-0.2..-0.7%, comparable to the width error).
        Never use on a model whose numbers will be reported.  Returns the
        number of lines changed."""
        import re
        t = open(azr_path).read()
        m = re.search(r"(<targetInt>\n)(.*?)(</targetInt>)", t, re.S)
        if not m:
            return 0
        lines = m.group(2).split("\n"); n = 0
        for i, l in enumerate(lines):
            f = l.split()
            if len(f) >= 4 and f[0] == "1":
                f[-2], f[-1] = f"{mult:g}", f"{ppw:g}"; lines[i] = " ".join(f); n += 1
        open(azr_path, "w").write(t[:m.start()] + m.group(1) + "\n".join(lines) + m.group(3) + t[m.end():])
        return n

    @staticmethod
    def set_grid_setting_p3(azr_path: str, mult: float, ppw: float) -> int:
        """Set ``width_multiplier`` and ``per_width`` of the ``quadrature`` field of every active data set that
        has one (the sub-point grid of a folded data set: how far around a narrow level the grid is refined and
        with how many points per width) -- for a SCREENING variant.  Measured on a large target-integrated
        segment: 5/20 is 6.7x cheaper than the converged 20/50 for a 0.03% objective
        shift; 2/x is NOT safe (-0.2..-0.7%, comparable to the width error).
        Never use on a model whose numbers will be reported.  Returns the
        number of data sets changed."""
        import pyazr
        project = pyazr.Project.open(os.path.abspath(azr_path))
        n = 0
        for name, rec in list(project.records().items()):
            if rec.record_kind != "dataset" or str(rec.get("active", "yes")) == "no":
                continue
            q = rec.get("quadrature")
            if not q:
                continue
            words = str(q).split()
            parts = dict(zip(words[::2], words[1::2]))
            parts["width_multiplier"], parts["per_width"] = f"{mult:g}", f"{ppw:g}"
            project.set(name, quadrature=" ".join(f"{k} {v}" for k, v in parts.items()))
            n += 1
        if n:
            project.save()
        return n

    def run_relocate(self, level: str, ex_new: float, half_window_kev: float = 5.0, max_nfev: int = 100,
                     note="") -> dict:
        """Move an existing level (key label such as ``9/2+#1``) to ``ex_new`` and
        re-polish with its energy confined to a narrow window; verify; the result
        becomes the incumbent if it beats it by tau."""
        from .levelscan import relocate
        from .search import tau_flip
        cid, x, obj, keys = self.incumbent()
        ev = self.sharded
        xs = take(x, keys, ev.keys, ev.x0)
        rid = self.ledger.start_round(self.structure_id, "relocate", f"{level} -> {ex_new} MeV {note}")
        # CONTROL first: polish the UNMOVED incumbent with the same budget, and judge the
        # relocation against that.  Without it a relocation scores whatever a warm-start
        # polish would have found anyway -- a level that walked straight back to
        # its starting energy still showed a "gain" of 309 (100 evaluations of drift).
        frC = polish(ev, xs, bounds=default_bounds(ev, xs, self.policy), policy=self.policy,
                     max_nfev=max_nfev, log=self.log)
        control = float(frC.objective)
        self.log(f"  control polish of the unmoved incumbent: {obj:.3f} -> {control:.3f} ({frC.status}); "
                 f"the relocation is judged against {control:.3f}")
        fr = relocate(ev, xs, level, ex_new, self.policy, half_window_kev=half_window_kev, max_nfev=max_nfev, log=self.log)
        stats = dict(level=level, ex_new=ex_new, before=obj, control=control, after=fr.objective,
                     gain_vs_control=control - fr.objective, status=fr.status, nfev=fr.nfev)
        tau = tau_flip(obj, ev.npts, self.search)
        if fr.status != "crashed":
            stats["figure"] = self.review_figure(xs, ev.keys, fr.x, ev.keys, ex_new,
                                                 f"relocate_{level.replace('/', '').replace('#', '_')}_E{ex_new:.4f}",
                                                 title=f"relocate {level} -> {ex_new:.4f} MeV: {obj:.0f} -> {fr.objective:.0f}",
                                                 label_before="incumbent", label_after="relocated")
        if fr.status != "crashed" and fr.objective < control - tau:
            ncid, sc = self.verify_candidate(fr.x, ev.keys, move=dict(kind="relocate", level=level, ex=ex_new), parent=cid)
            stats["verified"] = sc.objective
            self.note(f"relocation of {level} to {ex_new:.4f} MeV: {obj:.3f} -> {sc.objective:.3f} (verified, candidate {ncid}; "
                      f"control polish {control:.3f}; {fr.status}, {fr.nfev} evaluations)")
        else:
            tcid = self.ledger.add_candidate(self.structure_id, fr.x, ev.keys, move=dict(kind="relocate", level=level, ex=ex_new),
                                             parent=cid, status="worse")
            self.note(f"relocation of {level} to {ex_new:.4f} MeV rejected: {obj:.3f} -> {fr.objective:.3f} "
                      f"(sharded; tau {tau:.0f}; candidate {tcid} kept for the record)")
        self.ledger.finish_round(rid, **stats)
        self.ledger.write_status_md()
        return stats

    def export(self, path=None, note="export") -> tuple:
        """Bake the incumbent into a GUI-reviewable .azr + .sav (fresh-session verified)."""
        cid, x, obj, keys = self.incumbent()
        path = path or os.path.join(self.dir, self.reaction_name() + ".azr")
        azr, sav, sc = bake(self.base_azr, x, keys, path, self.flags, workdir=self.vdir, objective=self.objective,
                            output_dir=os.path.join(self.dir, "output"), name="export")
        self.log(f"exported candidate {cid} -> {azr} (+ .sav); objective {sc.objective:.3f}; round-trip {sc.roundtrip}"
                 + ("" if sc.chi2_self is None else f"; self-consistent chi2 {sc.chi2_self:.3f}"))
        self.note(f"exported candidate {cid} (objective {obj:.3f}) to {os.path.basename(azr)} + .sav; "
                  f"fresh-session objective {sc.objective:.3f}")
        return azr, sav, sc


    # -- reporting -------------------------------------------------------------------------

    def dataset_table(self, cid, stage="verified"):
        """Per-dataset chi2 of a ledger candidate from its stored per-segment
        scores (no engine session needed)."""
        seg_chi2 = self.ledger.segment_scores(cid, stage)
        segs = base_segments(self.base_azr)
        base_dir = os.path.dirname(self.base_azr)
        by = {}
        for s in segs:
            if not s.active or s.key not in seg_chi2:
                continue
            name = os.path.splitext(os.path.basename(s.data_file))[0]
            d = by.setdefault(name, dict(name=name, chi2=0.0, n=0, segments=[]))
            d["chi2"] += seg_chi2[s.key]
            d["n"] += int(rows_in_segment(s, base_dir).shape[0])
            d["segments"].append(s.key)
        rows = list(by.values())
        for d in rows:
            d["per_pt"] = d["chi2"] / max(d["n"], 1)
        rows.sort(key=lambda d: -d["per_pt"])
        return rows

    def report(self, path=None):
        """REPORT.md: incumbent, history of verified objectives, per-dataset table."""
        from .diagnose import format_dataset_table
        path = path or os.path.join(self.dir, "REPORT.md")
        cid, x, obj, keys = self.incumbent()
        st = self.ledger.status()
        lines = [f"# rmfit report -- {os.path.basename(self.dir)} ({time.strftime('%Y-%m-%d %H:%M')})", "",
                 f"Incumbent: candidate {cid}, verified objective {obj:.3f}; base model `{os.path.basename(self.base_azr)}`.", ""]
        hist = self.ledger.db.execute(
            "SELECT s.cid, s.objective, s.chi2, s.penalty, s.created, c.move FROM scores s JOIN candidates c ON c.id=s.cid "
            "WHERE c.structure_id=? AND s.stage='verified' ORDER BY s.id", (self.structure_id,)).fetchall()
        lines += ["## verified candidates", "", "| candidate | objective | chi2 | penalty | move | when |", "|---:|---:|---:|---:|---|---|"]
        for h in hist:
            mv = json.loads(h[5]).get("kind", "") if h[5] else ""
            lines.append(f"| {h[0]} | {h[1]:.3f} | {h[2]:.3f} | {h[3]:.3f} | {mv} | "
                         f"{time.strftime('%m-%d %H:%M', time.localtime(h[4]))} |")
        lines += ["", "## rounds", ""]
        for r in st["rounds"]:
            lines.append(f"- #{r['id']} {r['kind']}: " + json.dumps(r["stats"], default=float))
        lines += ["", f"## per-dataset chi-squared of candidate {cid}", "", format_dataset_table(self.dataset_table(cid)), ""]
        with open(path, "w") as fh:
            fh.write("\n".join(lines))
        return path


    # -- structure moves ----------------------------------------------------------------------

    def adopt_structure(self, sid, note=""):
        """Make a ledger structure the campaign's base model (its best verified
        candidate becomes the incumbent)."""
        s = self.ledger.structure(sid)
        self.cfg["base_azr"] = os.path.relpath(s["azr"], self.dir)
        with open(os.path.join(self.dir, "campaign.json"), "w") as f:
            json.dump(self.cfg, f, indent=2)
        self.base_azr = s["azr"]
        self.close()
        self.log(f"adopted structure {sid} ({s['name']}) as the base model {note}")
        self.ledger.log("adopt", f"structure {sid} {s['name']} {note}")

    def _spec_allows_adopt(self, jpi, cand, gamma_ev, delta, obj_ref, near, top_share) -> bool:
        """The spec's ``pre_authorized.adopt_level``: absent -> True (old behaviour); False -> never; a dict ->
        rules max_gamma_kev, min_gain_rel, forbid_near_same_jpi (near list non-empty), max_top_dataset_share."""
        if self.spec is None:
            return True
        rule = self.spec.pre_authorized.get("adopt_level", True)
        if rule is True or rule is None:
            return True
        if rule is False:
            return False
        if isinstance(rule, dict):
            ok = True
            if gamma_ev is not None and rule.get("max_gamma_kev") and gamma_ev > 1e3 * float(rule["max_gamma_kev"]):
                ok = False
            if rule.get("min_gain_rel") and obj_ref and delta < float(rule["min_gain_rel"]) * obj_ref:
                ok = False
            if rule.get("forbid_near_same_jpi") and near:
                ok = False
            if rule.get("max_top_dataset_share") and top_share is not None and top_share > float(rule["max_top_dataset_share"]):
                ok = False
            return ok
        return False

    def run_level_add(self, cand: LevelCandidate, *args, **kwargs) -> dict:
        """_run_level_add, with the policy's ``new_level_theta2_cap`` (2026-09-27) applied to the new level
        through a +-0.5 keV cap window at its seed energy: in force for the whole test, kept in
        campaign.json only when the level is accepted (adopt_structure writes the policy then)."""
        cap = float(getattr(self.policy, "new_level_theta2_cap", 0.0) or 0.0)
        if cap <= 0:
            return self._run_level_add(cand, *args, **kwargs)
        win = BoundsPolicy.seed_window(cand.energy, cap)
        self.policy.theta2_cap_windows = [win] + list(self.policy.theta2_cap_windows or [])
        pol = self.cfg.setdefault("policy", {})
        pol["theta2_cap_windows"] = [win] + list(pol.get("theta2_cap_windows") or [])
        self.log(f"  new-level theta2 cap {cap:g} for the level at {cand.energy} MeV (window {win[0]}-{win[1]} MeV)")
        stats = None
        try:
            stats = self._run_level_add(cand, *args, **kwargs)
            return stats
        finally:
            if not (isinstance(stats, dict) and stats.get("outcome") == "accept"):
                self.policy.theta2_cap_windows = [w for w in self.policy.theta2_cap_windows if w != win]
                pol["theta2_cap_windows"] = [w for w in pol.get("theta2_cap_windows", []) if w != win]

    def _run_level_add(self, cand: LevelCandidate, jpis=None, half_width=0.5, pinned_nfev=60,
                       global_nfev=150, note="", near_same_jpi_kev: float = 20.0,
                       dead_channel_ev: float = 1.0, local_release: float = 0.0, local_nfev: int = 80,
                       seed_omitted: Optional[bool] = None) -> dict:
        """Test one candidate level: J^pi scan on a window (pinned energy,
        sign patterns), then the best J^pi on the full model (pinned polish,
        released polish, verification), then classification.

        ``seed_omitted``: passed to ``add_level_variant`` for the full model -- the J^pi group's other open, alive
        particle channels are written FREE with a tiny seed instead of gamma 0 FIXED (None: the RMFIT_SEED_OMITTED
        environment variable, off by default).  A candidate that names no channels (an ENSDF level) otherwise gets only
        the lowest wave of the group's first two open pairs, and its other open channels stay fixed at zero in every
        later structure (2026-11-22: a recovered level in a benchmark lacked four of its open channels; freeing them
        was worth 1,305 of chi2).  Seeded channels count in tau only if the fit uses them (structure.effective_k).

        Accepting means the new structure becomes the base model.
        """
        from .diagnose import compare
        cid, x, obj, keys = self.incumbent()
        tag = f"E{cand.energy:.3f}".replace(".", "p")
        rid = self.ledger.start_round(self.structure_id, "level_add",
                                      f"{cand.jpi or 'scan'} @ {cand.energy} MeV ({cand.source}) {note}")
        t0 = time.time()
        # 1. window scan
        scan_dir = os.path.join(self.vdir, f"scan_{tag}")
        jl = [cand.jpi] if (cand.jpi and not jpis) else (jpis or None)
        results = scan_jpi_on_window(self.base_azr, scan_dir, cand, self.flags, objective=self.objective, jpis=jl,
                                     half_width=half_width, policy=self.policy, x_base=x, base_keys=keys,
                                     max_nfev=pinned_nfev, log=self.log, nshards=self.nshards)
        good = [r for r in results if r.status != "crashed" and r.x is not None]
        scan_summary = "; ".join(f"{r.jpi}:{r.signs}:{r.delta:+.0f}" for r in results)
        self.log(f"window scan: {scan_summary}")
        if not good:
            self.ledger.finish_round(rid, outcome="scan_failed", scan=scan_summary)
            self.note(f"level candidate {cand.jpi or 'scan'} @ {cand.energy} MeV: window scan failed ({scan_summary})")
            return dict(outcome="scan_failed", scan=results)
        best = good[0]
        # 2. full-model variant with the level, seeded from the window result
        path = os.path.join(self.vdir, f"add_{best.jpi.replace('/', '_')}_{tag}.azr")
        # 2026-09-28: the ledger keys structures by NAME and returns the existing one -- a re-test of a
        # candidate tested before (the run-off comparison copies) would inherit the old structure's
        # candidates and its incumbent.  Give a repeat a fresh name.
        taken = {s_["name"] for s_ in self.ledger.structures()}
        if os.path.basename(path) in taken:
            n_ = 2
            while os.path.basename(path).replace(".azr", f"_r{n_}.azr") in taken:
                n_ += 1
            path = path.replace(".azr", f"_r{n_}.azr")
            self.log(f"  structure name taken; this test is registered as {os.path.basename(path)}")
        # 2026-09-19: the level is written with a FREE energy.  Until now cand.energy_fixed (True
        # by default, and True from levelscan) went into the file as levelFix = 1, so the
        # "released" polish never released the new level's energy and the adopted structure
        # kept it fixed for every later campaign (all levels added 09-14 .. 09-19 sat at round
        # numbers; freeing four of them by hand gained clearly).  Pinning is now done
        # by the free mask (stage A) and a narrow window (stage B), not by the file.
        import dataclasses
        release_energy = bool(getattr(cand, "release_energy", True))
        cand_full = dataclasses.replace(cand, energy_fixed=not release_energy)
        info = add_level_variant(self.base_azr, path, cand_full, best.jpi,
                                 output_dir=os.path.join(self.vdir, f"output_add_{tag}"), seed_omitted=seed_omitted)
        if info.get("seeded_omitted"):
            self.log(f"  {len(info['seeded_omitted'])} omitted open channel(s) of the {best.jpi} group seeded free: "
                     + ", ".join(f"p{q['pair']} L{q['L']}" for q in info["seeded_omitted"]))
        lvl = info["level"]
        self.log(f"full model with {best.jpi} @ {cand.energy} MeV as {lvl}: {os.path.basename(path)}")
        ev2 = Sharded(path, os.path.join(self.vdir, f"shards_add_{tag}"), nshards=self.nshards, flags=self.flags,
                      objective=self.objective)
        try:
            sid2 = self.ledger.add_structure(os.path.basename(path), path, ev2.keys, parent=self.structure_id,
                                             note=f"+{best.jpi} @ {cand.energy} MeV ({cand.source})", fixed=ev2.fixed_amplitudes)
            # incumbent parameters by key; the new level's amplitudes from the window scan
            x2 = take(x, keys, ev2.keys, ev2.x0)
            # seed the new level's amplitudes from the window result: the window
            # variant appended the same level, so its keys match ev2's
            wx = best.x
            if wx is not None and getattr(best, "keys", None):
                own_keys = [k for k in best.keys if k.level == lvl]
                pos_w = {k: i for i, k in enumerate(best.keys)}
                pos_2 = {k: i for i, k in enumerate(ev2.keys)}
                for k in own_keys:
                    if k in pos_2:
                        x2[pos_2[k]] = wx[pos_w[k]]
            # 2026-10-12: channels that exist only because of RMFIT_SEED_OMITTED go back to the file's tiny seed.
            # The window model's data do not constrain them; in one window polish theta^2
            # ~1 went into three unconstrained channels and the full model started six times above its optimum (the test never recovered).
            _omit = {(q["pair"], q["L"], q["S"]) for q in (info.get("seeded_omitted") or [])}
            if _omit:
                _reset = [j for j, k in enumerate(ev2.keys)
                          if k.level == lvl and k.kind == "width" and (k.pair, k.L, k.S) in _omit]
                for j in _reset:
                    x2[j] = ev2.x0[j]
                if _reset:
                    self.log(f"  {len(_reset)} seeded-omitted channel(s) of {lvl} reset to their file seeds after the window scan")
            own = [j for j, k in enumerate(ev2.keys) if k.level == lvl]
            # 2026-09-14: channels of the new level with no physical width at its energy
            # (closed, or Gamma_W below `dead_channel_ev`: high-L / barely open) are zeroed
            # and FIXED.  Their amplitudes do nothing to the cross sections, so a polish
            # leaves them wherever the seed put them and the accepted 7/2+ 10.413 came out
            # with formal theta^2 ~0.6 in four such channels -- meaningless numbers that
            # would be reported as widths.  They are also not counted in tau's k.
            dead = []
            try:
                from .backend import AzrModel
                photon_pairs = {c.pair for lv in AzrModel.from_file(path).levels for c in lv.channels
                                if getattr(c, "is_photon", False)}
                gWw = ev2.wigner_widths(x2)          # particle channels only: photon amplitudes are linear, never "dead"
                for j in own:
                    k = ev2.keys[j]
                    if k.kind == "width" and k.pair not in photon_pairs and gWw.get(j, 0.0) < dead_channel_ev:
                        x2[j] = 0.0; dead.append(j)
            except Exception as err:
                self.log(f"  dead-channel check skipped ({err})")
            if dead:
                self.log(f"  {len(dead)} channel(s) of {lvl} with Gamma_W < {dead_channel_ev:g} eV at {cand.energy} MeV zeroed and fixed: "
                         + ", ".join(str(ev2.keys[j]).split(' ', 1)[1].rstrip(']') for j in dead))
            def _fix_dead(bnd):
                for j in dead:
                    bnd.fixed[j] = True
                return bnd
            # the evaluator 2026-09-28: while the new level is fitted, every OTHER free level wider than
            # policy.freeze_broad_during_test_ev is held fixed (energy + widths) in the released
            # stage; a final polish then frees everything before verification.  The level is
            # placed against a fixed broad background instead of running off to fill it.
            frozen_broad, gsum = [], {}
            fz = float(getattr(self.policy, "freeze_broad_during_test_ev", 0.0) or 0.0)
            if fz > 0:
                try:
                    phys = np.asarray(ev2.transform(x2), float)
                    for jj, kk in enumerate(ev2.keys[:ev2.n_rmatrix]):
                        if kk.kind == "width" and jj < phys.size:
                            gsum[kk.level] = gsum.get(kk.level, 0.0) + abs(float(phys[jj]))
                    frozen_broad = sorted(L for L, g in gsum.items() if L != lvl and g > fz)
                except Exception as err:
                    self.log(f"  freeze-broad: width evaluation failed ({err}); nothing frozen")
                self.log(f"  freeze-broad: {len(frozen_broad)} level(s) wider than {fz / 1e3:.0f} keV held fixed while "
                         f"the new level is fitted" + (": " + ", ".join(f"{L} ({gsum[L] / 1e3:.0f} keV)" for L in frozen_broad)
                                                       if frozen_broad else ""))
            fz_idx = [j for j, k in enumerate(ev2.keys) if k.level in frozen_broad]

            def _fix_broad(bnd):
                for j in fz_idx:
                    bnd.fixed[j] = True
                return bnd
            s2 = ev2.score(x2)
            self.log(f"full model at the seeded level: objective {s2.objective:.3f} (incumbent {obj:.3f})")
            # 3. pinned polish: the level's widths + norms; then released global polish
            free = free_mask(ev2, levels=[lvl], with_norms=True)
            e_own = [j for j in own if ev2.keys[j].kind == "energy"]
            for j in e_own:
                free[j] = False                                   # pinned: the energy stays at the candidate's
            b = _fix_dead(default_bounds(ev2, x2, self.policy))
            fA = polish(ev2, x2, free=free, bounds=b, policy=self.policy, max_nfev=pinned_nfev, log=self.log)
            self.log(f"pinned polish: {fA.start_objective:.3f} -> {fA.objective:.3f} ({fA.status})")
            bB = _fix_broad(_fix_dead(default_bounds(ev2, fA.x, self.policy)))
            half = None
            if e_own:
                # released: the energy may move +-max(50 keV, 3 Gamma), never more than the policy's window
                try:        # Gamma_c = theta_c^2 * Gamma_W,c over the level's own particle channels (MeV)
                    gWa, gWw2 = ev2.wigner_amplitudes(), ev2.wigner_widths(fA.x)
                    gam = sum((fA.x[jj] / gWa[jj]) ** 2 * gWw2[jj] for jj in own
                              if ev2.keys[jj].kind == "width" and gWa.get(jj, 0) > 0 and jj in gWw2) * 1e-6
                except Exception as err:
                    self.log(f"  width estimate for the energy window failed ({err}); using the candidate's seed widths")
                    gam = sum(abs(float(c.get("gamma", 0.0))) for c in (cand.channels or [])) * 1e-6
                half = min(max(0.050, 3.0 * gam), float(self.policy.e_window))
                for j in e_own:
                    bB.lb[j] = max(bB.lb[j], cand.energy - half); bB.ub[j] = min(bB.ub[j], cand.energy + half)
                self.log(f"  released stage: {lvl} energy free within +-{1000 * half:.0f} keV of {cand.energy:.4f} MeV "
                         f"(Gamma ~ {1000 * gam:.1f} keV)")
            x_rel = fA.x
            if local_release and local_release > 0:
                # 2026-10-02 local release: before the global stage, every level of ANY J^pi within +-local_release MeV
                # of the candidate, all background widths and the norms are polished together, so a new level is judged
                # after its neighbours and the background have re-partitioned (the evaluator's proposal)
                from .local import local_mask
                mL = local_mask(ev2, fA.x, [cand.energy], local_release, bg_energy=self.policy.bg_energy)
                fL = polish(ev2, fA.x, free=mL, bounds=bB, policy=self.policy, max_nfev=local_nfev, log=self.log)
                self.log(f"local-release polish (+-{1000 * local_release:.0f} keV, {int(mL.sum())} free): "
                         f"{fL.start_objective:.3f} -> {fL.objective:.3f} ({fL.status})")
                x_rel = fL.x
            fB = polish(ev2, x_rel, bounds=bB, policy=self.policy, max_nfev=global_nfev, log=self.log)
            self.log(f"released polish: {fB.start_objective:.3f} -> {fB.objective:.3f} ({fB.status})"
                     + (f" [{len(frozen_broad)} broad level(s) fixed]" if fz_idx else ""))
            b_final = bB
            if fz_idx:
                bR = _fix_dead(default_bounds(ev2, fB.x, self.policy))
                if e_own and half is not None:
                    for j in e_own:
                        bR.lb[j] = max(bR.lb[j], cand.energy - half); bR.ub[j] = min(bR.ub[j], cand.energy + half)
                fR = polish(ev2, fB.x, bounds=bR, policy=self.policy, max_nfev=global_nfev, log=self.log)
                self.log(f"release polish (broad levels freed): {fR.start_objective:.3f} -> {fR.objective:.3f} ({fR.status})")
                fB = fR
                b_final = bR
            # control: the incumbent (new level switched off = its seeds zeroed)
            # polished with the same budget on the same evaluator, so the level
            # is credited only with what it gains beyond plain continuation
            x_ctrl0 = x2.copy()
            for j in own:
                if ev2.keys[j].kind == "width":
                    x_ctrl0[j] = 0.0
            ctrl_free = np.array([not (k.level == lvl) for k in ev2.keys])
            if local_release and local_release > 0:
                # the control gets the SAME local-then-global budget, without the new level
                from .local import local_mask
                mC = local_mask(ev2, x_ctrl0, [cand.energy], local_release, bg_energy=self.policy.bg_energy,
                                exclude_levels=[lvl])
                fCL = polish(ev2, x_ctrl0, free=mC, bounds=default_bounds(ev2, x_ctrl0, self.policy),
                             policy=self.policy, max_nfev=local_nfev, log=self.log)
                self.log(f"control local-release polish (level off, {int(mC.sum())} free): "
                         f"{fCL.start_objective:.3f} -> {fCL.objective:.3f} ({fCL.status})")
                x_ctrl0 = fCL.x
            fC = polish(ev2, x_ctrl0, free=ctrl_free, bounds=default_bounds(ev2, x_ctrl0, self.policy),
                        policy=self.policy, max_nfev=global_nfev, log=self.log)
            self.log(f"control polish (level off): {fC.start_objective:.3f} -> {fC.objective:.3f} ({fC.status})")
            # 4. verify on a single session of the variant
            sc, keys2 = verify(path, fB.x, self.flags, objective=self.objective, keys=ev2.keys,
                               workdir=self.vdir, name=f"verify_add_{tag}")
            ncid = self.ledger.add_candidate(sid2, take(fB.x, ev2.keys, keys2, np.zeros(len(keys2))), keys2,
                                             move=dict(kind="level_add", jpi=best.jpi, energy=cand.energy), parent=cid,
                                             status="verified")
            self.ledger.record(ncid, sc, "verified", evaluator="single")
            # 5. classify
            gW = ev2.wigner_amplitudes()
            theta2 = {str(ev2.keys[j]): float((fB.x[j] / gW[j]) ** 2) for j in own
                      if ev2.keys[j].kind == "width" and j in gW and gW[j] > 0}
            e_moved = 0.0
            for j in own:
                if ev2.keys[j].kind == "energy":
                    e_moved = float(fB.x[j] - cand.energy)
            # k for tau: "effective parameters" (2026-09-20, the evaluator) -- the scan pattern's channels
            # always count; a channel seeded by RMFIT_SEED_OMITTED counts only if the fit used it
            # (structure.effective_k).  With nothing seeded this is the old len(own) - len(dead).
            from .structure import effective_k
            try:
                _gWw = ev2.wigner_widths(fB.x)
            except Exception:
                _gWw = {}
            _params = []
            for j in own:
                if j in dead:
                    continue
                kk = ev2.keys[j]
                if kk.kind == "width":
                    t2 = float((fB.x[j] / gW[j]) ** 2) if (j in gW and gW[j] > 0) else 0.0
                    _params.append(dict(kind="width", pair=kk.pair, L=kk.L, S=kk.S, theta2=t2,
                                        gamma_ev=t2 * float(_gWw.get(j, 0.0))))
                else:
                    _params.append(dict(kind=kk.kind))
            _seed = [(q["pair"], q["L"], q["S"]) for q in (info.get("seeded_omitted") or [])]
            _k = effective_k(_params, _seed)
            k_new = _k["k_effective"]
            s2n = obj / max(sc.npts, 1)
            _lnN = np.log(max(sc.npts, 2))
            tau = max(s2n * k_new * _lnN, 0.002 * obj)
            _kstats = dict(k_pattern=_k["k_pattern"], k_effective=_k["k_effective"], k_all=_k["k_all"],
                           tau_pattern=float(max(s2n * _k["k_pattern"] * _lnN, 0.002 * obj)),
                           tau_all=float(max(s2n * _k["k_all"] * _lnN, 0.002 * obj)))
            if _seed:
                self.log(f"  k for tau: pattern {_k['k_pattern']}, effective {_k['k_effective']}, all {_k['k_all']} "
                         f"-> tau {tau:.0f} (pattern-only {_kstats['tau_pattern']:.0f}, all-channels {_kstats['tau_all']:.0f}); "
                         f"seeded channels the fit used: "
                         + (", ".join(f"p{a} L{b} S{cc:g}" for a, b, cc in _k["used"]) or "none"))
            # redistribution against the incumbent's verified per-segment table
            inc_seg = self.ledger.segment_scores(cid, "verified")
            rows = []
            gains = losses = 0.0
            for kseg, v in sc.segment_chi2.items():
                d = v - inc_seg.get(kseg, v)
                if d < 0: gains -= d
                else: losses += d
            rho = losses / gains if gains > 0 else (np.inf if losses > 0 else 0.0)
            cap = self.policy.cap_for(cand.energy)
            # reference = the better of the verified incumbent and the control (sharded, same evaluator
            # as fB; the two differ from verified numbers by the same ~1e-6)
            obj_ref = min(obj, fC.objective) if fC.status != "crashed" else obj
            # the new level's width amplitudes the released polish left ON a bound: a level
            # pinned at the theta^2 cap is an absorber, not a resonance (see classify)
            at_bound = [str(ev2.keys[j]) for j in getattr(fB, "bounds_hit", []) or []
                        if j in own and ev2.keys[j].kind == "width"]
            if self.policy.width_bound_for(cand.energy) is not None:
                from .local import width_bound_pins
                wbp = width_bound_pins(ev2.keys, own, b_final.ub, ev2.wigner_amplitudes(), cap)
                pinned_wb = [k for k in at_bound if k in wbp]
                if pinned_wb:
                    self.log(f"  {len(pinned_wb)} width(s) on the WIDTH bound ({self.policy.width_bound_for(cand.energy) / 1e3:.0f} keV "
                             f"per channel), not counted as unphysical: " + ", ".join(pinned_wb))
                at_bound = [k for k in at_bound if k not in wbp]
            # total-width floor (policy.total_width_floor_windows): the floor HOLDING the new level (its total width at
            # the floor) is the "no broad level of this J^pi wanted" signal, as an amplitude on a per-channel floor was
            for _f in getattr(b_final, "floors", None) or []:
                if _f.level == lvl:
                    self.log(f"  total width of {lvl} after the released polish: {_f.gamma(fB.x) / 1e3:.1f} keV "
                             f"(floor {_f.gmin / 1e3:g} keV)")
                    if _f.binding(fB.x):
                        at_bound.append(f"Gamma_total[{lvl}] at the {_f.gmin / 1e3:g} keV floor")
            if at_bound:
                self.log(f"  new level's amplitudes at a bound after the released polish: {', '.join(at_bound)}")
            # the new level's summed physical width (eV), from the per-channel list built for effective_k
            level_gamma_ev = (float(sum(q.get("gamma_ev", 0.0) for q in _params if q.get("kind") == "width"))
                              if _gWw else None)
            # 2026-10-12 (the evaluator): judge a resonance on the data it touches -- chi2 with and without it inside
            # E_R +- 1/2/3 Gamma, candidate and control on ONE fresh single session (same build: a pyazr change
            # during the job moved verified-vs-sharded numbers by +53 on 10-12).  The control's vector is kept in
            # the ledger (status/stage "control", never "verified", so it can never become the incumbent).
            win, ctrl_cid, ctrl_single = None, None, None
            try:
                ctrl_cid = self.ledger.add_candidate(sid2, take(fC.x, ev2.keys, keys2, np.zeros(len(keys2))), keys2,
                                                     move=dict(kind="level_add_control", jpi=best.jpi,
                                                               energy=cand.energy), parent=cid, status="control")
                from .diagnose import window_delta, format_window_delta
                wpath = single_variant(path, f"window_{tag}", self.vdir)
                with SingleSession(wpath, cwd=self.vdir, flags=self.flags, name=f"window_{tag}",
                                   objective=self.objective) as evw:
                    xw = take(fB.x, ev2.keys, evw.keys, evw.x0)
                    xo = take(fC.x, ev2.keys, evw.keys, evw.x0)
                    scC = evw.score(xo)
                    ctrl_single = float(scC.objective)
                    self.ledger.record(ctrl_cid, scC, "control", evaluator="single")
                    win = window_delta(evw, xw, xo, cand.energy + e_moved, (level_gamma_ev or 0.0) * 1e-6,
                                       min_half_mev=float(getattr(self.policy, "window_min_half_mev", 0.0) or 0.0))
                self.log(f"  control on a fresh single session: {ctrl_single:.3f} (its polish, sharded: {fC.objective:.3f}); "
                         f"level - control on the same session {sc.objective - ctrl_single:+.1f} (objective)")
                for line in format_window_delta(win):
                    self.log("  " + line)
            except Exception as err:
                self.log(f"  +-m Gamma window table skipped ({err})")
            narrow_floor = float(getattr(self.policy, "min_level_width_ev", 0.0) or 0.0) or None
            broad_ceiling = float(getattr(self.policy, "max_level_width_ev", 0.0) or 0.0) or None
            # 2026-09-21: is the add a SPLIT of a same-J^pi neighbour?  Summed (free-channel) physical
            # widths of the other levels of this J^pi at the incumbent (x2) and after the released polish
            split = None
            try:
                from .structure import split_evidence
                _gW0 = ev2.wigner_widths(x2) if _gWw else {}
                _nb = {}
                for j, kk in enumerate(ev2.keys):
                    if kk.level is None or kk.level == lvl or kk.jpi != best.jpi:
                        continue
                    d = _nb.setdefault(kk.level, dict(label=kk.level, energy=None, gamma_before_ev=0.0, gamma_after_ev=0.0))
                    if kk.kind == "energy":
                        d["energy"] = float(fB.x[j])
                    elif kk.kind == "width" and j in gW and gW[j] > 0:
                        d["gamma_before_ev"] += float((x2[j] / gW[j]) ** 2) * float(_gW0.get(j, 0.0))
                        d["gamma_after_ev"] += float((fB.x[j] / gW[j]) ** 2) * float(_gWw.get(j, 0.0))
                # fixed-energy levels and levels at/above the background energy are poles, not partners (2026-10-16)
                _bgE = float(getattr(self.policy, "bg_energy", 1e9))
                _nbl = [d for d in _nb.values() if d["energy"] is not None and d["energy"] < _bgE]
                if level_gamma_ev:
                    split = split_evidence(cand.energy + e_moved, level_gamma_ev, _nbl)
                if split:
                    split["label"] = f"{split['label']}@{[d['energy'] for d in _nbl if d['label'] == split['label']][0]:.4f}"
            except Exception as err:
                self.log(f"  split check skipped ({err})")
            # share of the summed dataset gains carried by the single best dataset
            top_share, top_name = None, ""
            try:
                from types import SimpleNamespace
                from .diagnose import compare as _compare
                _cmp = _compare(ev2, SimpleNamespace(segment_chi2=inc_seg, objective=obj), sc)
                if _cmp["gains"] > 0:
                    top_name = _cmp["rows"][0]["name"]; top_share = float(-_cmp["rows"][0]["delta"] / _cmp["gains"])
            except Exception as err:
                self.log(f"  per-dataset share skipped ({err})")
            outcome = classify(obj_ref, sc.objective, theta2, e_moved, cap, tau, rho, at_bound=at_bound,
                               e_window=(half if e_own else 0.05),
                               level_gamma_ev=level_gamma_ev, narrow_width_ev=narrow_floor, splits=bool(split),
                               broad_width_ev=broad_ceiling)
            delta = obj_ref - sc.objective
            for hw in (0.15, 0.5):                     # broad levels need a wider view
                fig = self.review_figure(x, keys, fB.x, ev2.keys, cand.energy,
                                         f"leveladd_{best.jpi.replace('/', '')}_E{cand.energy:.4f}",
                                         title=f"level {best.jpi} @ {cand.energy} MeV: {obj_ref:.0f} -> {sc.objective:.0f} ({outcome})",
                                         label_before="incumbent", label_after="with the level", azr=path, half_width=hw)
                if fig:
                    break
            self.log(f"level {best.jpi} @ {cand.energy} MeV: verified {sc.objective:.3f} vs incumbent {obj:.3f} "
                     f"/ control {fC.objective:.3f} (delta {delta:+.1f}, tau {tau:.0f}, rho {rho:.2f}, "
                     f"max theta2 {max(theta2.values(), default=0):.3g}, "
                     f"energy {cand.energy + e_moved:.4f} = seed {1000 * e_moved:+.1f} keV, "
                     f"Gamma {('%.4g eV' % level_gamma_ev) if level_gamma_ev is not None else 'n/a'}) -> {outcome}"
                     + (f"  [narrower than the campaign's floor of {narrow_floor:g} eV: recorded, not adopted]"
                        if outcome == "narrow" else ""))
            if top_share is not None:
                self.log(f"  {100 * top_share:.0f}% of the summed dataset gains come from {top_name}")
            if split:
                self.log(f"  SPLIT: the existing {split['label']} is {split['sep_kev']:.1f} keV away (unresolved) and gave up "
                         f"{split['lost_ev']:.4g} eV = {100 * split['lost_frac']:.0f}% of the new level's width"
                         + (" -- one lineshape described by two levels: recorded, NOT adopted; the evaluator decides"
                            if outcome == "splits" else ""))
            # 2026-09-14: a candidate that lands on top of an existing level of the same J^pi
            # is a LINESHAPE hypothesis, not a level -- a level a few hundred eV wide was once accepted
            # about a keV from a few-keV level of the same J^pi because BandH carries no resolution
            # convolution.  The criteria test whether a level HELPS; they cannot tell it from
            # a resolution or single-point effect, so the verdict is marked provisional and
            # the competing test (resolution on the driving segments, level removed) is named.
            near = []
            try:
                from .backend import AzrModel
                for lv in AzrModel.from_file(self.base_azr).find(jpi=best.jpi):
                    if lv.channels and all(bool(c.channel_fixed) and float(c.gamma) == 0.0 for c in lv.channels):
                        continue       # a switched-off level (deactivate_level_variant: recall children) is not a neighbour
                    if abs(float(lv.energy) - float(cand.energy)) <= near_same_jpi_kev * 1e-3:
                        near.append(f"{best.jpi}@{float(lv.energy):.4f}")
            except Exception as err:
                self.log(f"  same-J^pi proximity check skipped ({err})")
            if near:
                self.log(f"  PROVISIONAL: within {near_same_jpi_kev:g} keV of the existing {', '.join(near)} -- a lineshape "
                         f"hypothesis (doublet vs. resolution/target convolution on the driving segments vs. one "
                         f"point); test the alternatives before adopting it as a level")
            # The spec gate runs BEFORE the ledger row, the readme note and the returned dict are
            # built (fixture a queue job, 2026-11-07: the gate logged accept-pending but the
            # result still said accept).
            if outcome == "accept" and not self._spec_allows_adopt(best.jpi, cand, level_gamma_ev, delta, obj_ref, near, top_share):
                outcome = "accept-pending"
                did = self.decisions.add(f"adopt-level-{best.jpi}-{cand.energy:.3f}",
                                         f"adopt the {best.jpi} level at {cand.energy + e_moved:.4f} MeV (gain {delta:+.1f} of {obj_ref:.0f}, "
                                         f"Gamma {('%.4g eV' % level_gamma_ev) if level_gamma_ev is not None else 'n/a'}, structure {sid2})?",
                                         evidence=f"ledger structure {sid2} candidate {ncid}; figures/leveladd_*E{cand.energy:.4f}*; REPORT.md",
                                         options=("adopt", "reject", "test alternatives"), default="wait")
                self.log(f"  accept NOT pre-authorized by spec.json -> recorded as accept-pending, decision {did} queued in DECISIONS.md")
                self.ledger.log("decision", f"{did} level {best.jpi} @ {cand.energy} structure {sid2}")
            stats = dict(**_kstats, jpi=best.jpi, signs=best.signs, energy=cand.energy, delta=delta, tau=tau, rho=rho,
                         theta2=theta2, outcome=outcome, verified=sc.objective, control=fC.objective,
                         structure=sid2, candidate=ncid, scan=scan_summary, wall_min=(time.time() - t0) / 60,
                         near_same_jpi=", ".join(near), at_bound=", ".join(at_bound),
                         e_moved=float(e_moved), energy_fit=float(cand.energy + e_moved),
                         level_gamma_ev=level_gamma_ev, split=json.dumps(split) if split else "",
                         top_dataset=top_name, top_dataset_share=top_share,
                         control_cid=ctrl_cid, control_single=ctrl_single,
                         window=json.dumps(dict(e_r=win["e_r"], gamma_mev=win["gamma_mev"], windows=win["windows"],
                                                datasets=win["datasets"][:8] + win["datasets"][-8:])) if win else "")
            self.ledger.finish_round(rid, **{k: (v if not isinstance(v, dict) else json.dumps(v)) for k, v in stats.items()})
            self.note(f"level candidate {best.jpi} @ {cand.energy} MeV ({cand.source}): window scan {scan_summary}; "
                      f"full model verified {sc.objective:.3f} vs incumbent {obj:.3f} / control {fC.objective:.3f} "
                      f"(delta {delta:+.1f}, tau {tau:.0f}, rho {rho:.2f}); "
                      f"theta2 " + ", ".join(f"{k.split(' ')[1] if ' ' in k else k}={v:.2g}" for k, v in theta2.items())
                      + (f"; at a bound: {', '.join(at_bound)}" if at_bound else "")
                      + f" -> {outcome}" + (" -- adopted as the base model" if outcome == "accept" else "")
                      + (f"; PROVISIONAL (within {near_same_jpi_kev:g} keV of {', '.join(near)}: lineshape hypothesis)" if near else ""))
            if outcome == "accept":
                self.adopt_structure(sid2, note=f"level {best.jpi} @ {cand.energy} MeV accepted")
            self.ledger.write_status_md()
            return stats
        finally:
            ev2.close()


    def run_dataset_off(self, file_substr: str, nfev: int = 100, note: str = "") -> dict:
        """Switch every segment whose data file matches ``file_substr`` off,
        polish, verify, and report the *extra* recovery: how much the rest of
        the data improve beyond simply dropping that dataset's own chi2 (the
        Botek diagnostic from an earlier campaign: a dataset that distorts the
        fit shows a large extra recovery; a merely noisy one shows none)."""
        from .structure import dataset_off_variant, carry_moved_norms
        cid, x, obj, keys = self.incumbent()
        tag = "".join(ch if ch.isalnum() else "_" for ch in file_substr)[:30]
        rid = self.ledger.start_round(self.structure_id, "dataset_off", f"{file_substr} {note}")
        t0 = time.time()
        path = os.path.join(self.vdir, f"off_{tag}.azr")
        info = dataset_off_variant(self.base_azr, path, file_substr,
                                   output_dir=os.path.join(self.vdir, f"output_off_{tag}"))
        inc_seg = self.ledger.segment_scores(cid, "verified")
        segs_off = [s.key for s in base_segments(self.base_azr) if s.active and file_substr in s.data_file]
        own_chi2 = sum(inc_seg.get(k, 0.0) for k in segs_off)
        self.log(f"dataset off: {file_substr} -> {info['segments']} segments deactivated; their chi2 in the incumbent {own_chi2:.1f}")
        moved = info.get("moved") or {}
        if moved:
            self.log("dataset off: normgroup leader(s) switched off, the group's norm fields and fitted norm carried "
                     "to the next active member: " + ", ".join(f"seg {o} -> seg {n}" for o, n in moved.items()))
        ev2 = Sharded(path, os.path.join(self.vdir, f"shards_off_{tag}"), nshards=self.nshards, flags=self.flags,
                      objective=self.objective)
        try:
            # a switched-off leader's group: its fitted norm now lives under the new leader's key
            x2 = carry_moved_norms(x, keys, take(x, keys, ev2.keys, ev2.x0), ev2.keys, moved)
            s0 = ev2.score(x2)
            fr = polish(ev2, x2, bounds=default_bounds(ev2, x2, self.policy), policy=self.policy,
                        max_nfev=nfev, log=self.log)
            sc, keys2 = verify(path, fr.x, self.flags, objective=self.objective, keys=ev2.keys,
                               workdir=self.vdir, name=f"verify_off_{tag}")
            rest_before = obj - own_chi2
            extra = rest_before - sc.objective
            # who gains when this dataset is removed: per-dataset chi2 of the REST, frozen
            # (s0) vs refitted (s1), on the same evaluator -- the pairwise tension map
            from .diagnose import compare
            s1 = ev2.score(fr.x)
            cmp = compare(ev2, s0, s1)
            by_dataset = [dict(name=r["name"], before=float(r["before"]), after=float(r["after"]),
                               delta=float(r["delta"]), n=int(r["n"])) for r in cmp["rows"]]
            self.log(f"dataset off {file_substr}: per-dataset change of the rest (frozen -> refitted; rho {cmp['rho']:.2f}):\n"
                     + "\n".join(f"      {r['name']:34s} {r['before']:11.1f} -> {r['after']:11.1f}  {r['delta']:+10.1f}  (N {r['n']})"
                                 for r in by_dataset if abs(r["delta"]) >= 1.0))
            self.log(f"dataset off {file_substr}: rest of the data {rest_before:.1f} -> {sc.objective:.1f} after refit "
                     f"(extra recovery {extra:+.1f} = {100 * extra / max(rest_before, 1):.2f}% of the rest; "
                     f"the dataset's own chi2 was {own_chi2:.1f})")
            stats = dict(file=file_substr, segments=info["segments"], own_chi2=own_chi2, rest_before=rest_before,
                         rest_after=sc.objective, extra_recovery=extra, extra_pct=100 * extra / max(rest_before, 1),
                         by_dataset=by_dataset, wall_min=(time.time() - t0) / 60)
            self.ledger.finish_round(rid, **stats)
            self.note(f"dataset-off test {file_substr} ({info['segments']} segments, own chi2 {own_chi2:.0f}): "
                      f"the rest of the data refit from {rest_before:.0f} to {sc.objective:.0f}, extra recovery "
                      f"{extra:+.0f} ({stats['extra_pct']:.2f}%). A large extra recovery means the dataset was "
                      f"distorting the fit; a small one that it is merely poorly described.")
            self.ledger.write_status_md()
            return stats
        finally:
            ev2.close()


    def run_remove(self, top: int = 3, nfev: int = 100, keep_last_of_jpi: bool = True,
                   note: str = "") -> dict:
        """Backward elimination, the plan's Phase 4(a).  Rank every level by the
        frozen cost of switching it off, then HONESTLY refit the cheapest
        ``top`` of them and keep a removal whose refitted objective is no worse
        than ``max(s^2 k ln N, 0.1% obj)``.

        The frozen ranking is only an upper bound on a level's importance --
        the whole point of the honest refit is that the rest of the model
        compensates, and that can turn a level's apparent cost of
        millions into a few thousand.  Never removes the last level of a J^pi
        group: that one carries the group's hard-sphere term."""
        from .structure import frozen_removal_screen, deactivate_level_variant
        cid, x, obj, keys = self.incumbent()
        rid = self.ledger.start_round(self.structure_id, "remove", note)
        t0 = time.time()
        ev = Sharded(self.base_azr, os.path.join(self.vdir, "shards_remove"), nshards=self.nshards,
                     flags=self.flags, objective=self.objective)
        rows, results = [], []
        try:
            xs = take(x, keys, ev.keys, ev.x0)
            rows = frozen_removal_screen(ev, xs, log=self.log)
            n_by_jpi = {}
            for k in ev.keys:
                if k.kind == "width" and k.level:
                    n_by_jpi.setdefault(level_group(k.level), set()).add(str(k.level))
            # the hard-sphere term is carried by ANY active level of the group, fixed ones included: counting only the
            # levels with free widths called a spurious added 6+ "the only level of 6+" on a windowed model whose other
            # 6+ levels are fixed.  Single-nucleus labels only ("J^pi#n"); a
            # multi-nucleus label ("...@23Ne") keeps the free-level count.
            from .backend import AzrModel
            n_active = {}
            for lv in AzrModel.from_file(self.base_azr).active_levels:
                n_active[lv.jpi] = n_active.get(lv.jpi, 0) + 1
            npts = max(sum(b - a for _, (a, b) in ev.row_ranges), 1)
            s2 = obj / npts
        finally:
            ev.close()
        self.log(f"removal screen: {len(rows)} levels ranked; s^2 = {s2:.2f}, N = {npts}")
        # CONTROL: polish the UNMODIFIED model from the same start with the same budget.
        # Without this a removal looks good whenever the incumbent simply is not at a local
        # minimum for that budget -- on the 13N fixture two different removals both "gained"
        # 265k, which was the polish escaping, not the levels.  Every cost below is measured
        # against the control, not against the incumbent.
        evC = Sharded(self.base_azr, os.path.join(self.vdir, "shards_rm_control"), nshards=self.nshards,
                      flags=self.flags, objective=self.objective)
        try:
            xC = take(x, keys, evC.keys, evC.x0)
            frC = polish(evC, xC, bounds=default_bounds(evC, xC, self.policy), policy=self.policy,
                         max_nfev=nfev, log=self.log)
            scC, _ = verify(self.base_azr, frC.x, self.flags, objective=self.objective, keys=evC.keys,
                            workdir=self.vdir, name="verify_rm_control")
        finally:
            evC.close()
        control = float(scC.objective)
        self.log(f"control polish of the full model, same budget: {obj:.1f} -> {control:.1f} "
                 f"({control - obj:+.1f}); removals are judged against this")
        for row in rows[:max(top, 0)]:
            lvl = str(row["level"])
            jpi = lvl.split("#")[0]
            grp = level_group(lvl)
            n_grp = len(n_by_jpi.get(grp, ()))
            if "@" not in str(lvl):
                n_grp = max(n_grp, n_active.get(str(grp), 0))
            if keep_last_of_jpi and n_grp <= 1:
                self.log(f"  {lvl}: the only level of {jpi}; kept (it carries the group's hard-sphere term)")
                continue
            energy = self._level_energy(lvl)
            if energy is None:
                self.log(f"  {lvl}: could not resolve its energy from the model; skipped")
                continue
            tag = "".join(ch if ch.isalnum() else "_" for ch in lvl)[:24]
            path = os.path.join(self.vdir, f"rm_{tag}.azr")
            try:
                info = deactivate_level_variant(self.base_azr, path, jpi, energy,
                                                nucleus=ParamKey("energy", level=lvl).nucleus_tag or None,
                                                output_dir=os.path.join(self.vdir, f"output_rm_{tag}"))
            except Exception as err:
                self.log(f"  {lvl}: could not build the removal variant ({err}); skipped")
                continue
            ev2 = Sharded(path, os.path.join(self.vdir, f"shards_rm_{tag}"), nshards=self.nshards,
                          flags=self.flags, objective=self.objective)
            try:
                x2 = take(x, keys, ev2.keys, ev2.x0)
                fr = polish(ev2, x2, bounds=default_bounds(ev2, x2, self.policy), policy=self.policy,
                            max_nfev=nfev, log=self.log)
                sc, _ = verify(path, fr.x, self.flags, objective=self.objective, keys=ev2.keys,
                               workdir=self.vdir, name=f"verify_rm_{tag}")
            finally:
                ev2.close()
            k_free = int(row["n_free"]) + 1          # its widths plus its energy
            tau = max(s2 * k_free * math.log(npts), 0.001 * obj)
            cost = sc.objective - control          # against the control, not the incumbent
            verdict = "REMOVE" if cost <= tau else "keep"
            self.log(f"  {lvl} at {energy:.4f} MeV: frozen {row['delta']:+.1f}, honest refit "
                     f"{obj:.1f} -> {sc.objective:.1f} vs control {control:.1f} "
                     f"(cost {cost:+.1f}, threshold {tau:.1f}) -> {verdict}")
            results.append(dict(level=lvl, jpi=jpi, energy=energy, frozen_delta=row["delta"],
                                refit=sc.objective, control=control, cost=cost, threshold=tau,
                                verdict=verdict, n_free=k_free, variant=path,
                                removed_with_it=info.get("removed", [])))
        stats = dict(screened=len(rows), refitted=len(results), control=control,
                     control_gain=control - obj,
                     removable=[r["level"] for r in results if r["verdict"] == "REMOVE"],
                     results=results, s2=s2, wall_min=(time.time() - t0) / 60)
        self.ledger.finish_round(rid, **stats)
        if results:
            self.note(f"backward elimination ({len(rows)} levels screened, {len(results)} refitted; "
                      f"control polish of the full model {obj:.0f} -> {control:.0f}): "
                      + "; ".join(f"{r['level']} at {r['energy']:.3f} MeV cost {r['cost']:+.0f} "
                                  f"vs threshold {r['threshold']:.0f} -> {r['verdict']}" for r in results)
                      + ".  A frozen ranking is an upper bound only; the honest refit is what decides.")
        else:
            self.note(f"backward elimination: {len(rows)} levels screened, none refitted "
                      f"(all candidates were the last level of their J^pi group, or unresolvable).")
        self.ledger.write_status_md()
        return stats

    def _level_energy(self, label: str):
        """The file energy of a level label like ``5/2+#9`` (``5/2+#9@23Ne`` with several compound nuclei), for the
        variant builders.  The ordinal counts the ACTIVE levels of that J^pi (and nucleus) in file order, as the
        engine does (pyazr ``AzrModel.engine_level_keys``)."""
        from .backend import AzrModel
        k = ParamKey("energy", level=str(label))
        try:
            jpi, idx, nuc = k.jpi, k.ordinal, k.nucleus_tag
        except (IndexError, ValueError):
            return None
        m = AzrModel.from_file(self.base_azr)
        lv = [l for l in m.find(jpi=jpi) if l.active]
        if nuc:
            lv = [l for l in lv if l.channels and str(m.nucleus_of_pair(l.channels[0].pair)) == nuc]
        if 1 <= idx <= len(lv):
            return float(lv[idx - 1].energy)
        return None

    # -- the recovery benchmark ------------------------------------------------------------

    def perturb(self, out_dir, n_flips=8, seed=0, jitter=0.3, e_kick=0.02, polish_nfev=60, note="") -> dict:
        """Create a sibling campaign whose baseline is the incumbent with
        ``n_flips`` high-exposure relative signs flipped (spread over J^pi
        groups), amplitudes jittered by x U(1-jitter, 1+jitter), free level
        energies kicked by +-e_kick MeV, then polished into the nearest local
        minimum -- a genuinely worse basin the sign rounds must find their way
        out of.  The plan's Validation A."""
        from .moves import exposure, single_flips, apply
        from .model import base_segments as _bs
        import shutil
        cid, x, obj, keys = self.incumbent()
        out_dir = os.path.abspath(out_dir)
        os.makedirs(out_dir, exist_ok=True)
        from .model import link_data_dirs
        # from the base model's own directory first (a base outside the campaign directory
        # has its data/ there, not in the campaign), then from the campaign directory as before
        link_data_dirs(self.base_azr, out_dir)
        link_data_dirs(self.base_azr, out_dir, src_dir=self.dir)
        base_name = os.path.basename(self.base_azr)
        base_copy = os.path.join(out_dir, base_name)
        if PYAZR3:
            from .azr3 import write_variant
            os.makedirs(os.path.join(out_dir, "output"), exist_ok=True)
            write_variant(self.base_azr, base_copy, os.path.join(out_dir, "output"))
        else:
            text = open(self.base_azr).read()
            from .model import set_output_dir_text
            with open(base_copy, "w") as fh:
                fh.write(set_output_dir_text(text, os.path.join(out_dir, "output")))
        ev = self.sharded
        xs = take(x, keys, ev.keys, ev.x0)
        rng = np.random.default_rng(seed)
        expo = exposure(ev, xs)
        moves = single_flips(xs, ev.keys, exposure_map=expo, fixed=ev.fixed_amplitudes)
        chosen, groups = [], {}
        for mv in moves:                      # highest exposure first, at most 2 per J^pi group
            g = mv.jpi
            if groups.get(g, 0) >= 2:
                continue
            chosen.append(mv); groups[g] = groups.get(g, 0) + 1
            if len(chosen) >= n_flips:
                break
        xp = xs.copy()
        for mv in chosen:
            xp = apply(xp, ev.keys, mv)
        for j, k in enumerate(ev.keys):
            if k.kind == "width" and xp[j] != 0:
                xp[j] *= rng.uniform(1 - jitter, 1 + jitter)
            elif k.kind == "energy":
                xp[j] += rng.uniform(-e_kick, e_kick)
        # settle into the new basin with the flipped signs locked first (an
        # unlocked polish from a raw flip walks straight back through zero),
        # then release everything
        pos = {k: j for j, k in enumerate(ev.keys)}
        lock = {k: (1 if xp[pos[k]] > 0 else -1) for mv in chosen for k in mv.targets}
        b1 = default_bounds(ev, xs, self.policy, sign_lock=lock)
        xp = b1.clip(xp)
        o_raw = ev.objective(xp)
        f1 = polish(ev, xp, bounds=b1, policy=self.policy, max_nfev=polish_nfev, log=self.log)
        b = default_bounds(ev, f1.x, self.policy)
        fr = polish(ev, f1.x, bounds=b, policy=self.policy, max_nfev=polish_nfev, log=self.log)
        from .keys import canonical_signs as _cs
        same_basin = _cs(fr.x, ev.keys, ev.fixed_amplitudes) == _cs(xs, ev.keys, ev.fixed_amplitudes)
        self.log(f"perturb: {len(chosen)} flips " + ", ".join(m.description for m in chosen)
                 + f"; raw {o_raw:.1f} -> locked {f1.objective:.1f} -> released {fr.objective:.1f} "
                 f"(incumbent {obj:.1f}); " + ("signs REVERTED to the incumbent's pattern" if same_basin
                                               else "a different sign basin"))
        # the new campaign: same config, its baseline = the perturbed, polished vector (verified there)
        cfg = dict(self.cfg); cfg["base_azr"] = base_name
        # a fixed relvector covariance model travels with the objective: the sibling must score with the SAME reference
        # (campaign.json named cov_reference.npz, the file stayed behind and the
        # sibling's first session stopped with FileNotFoundError)
        if self.objective.cov_reference:
            shutil.copy2(self.objective.cov_reference, os.path.join(out_dir, "cov_reference.npz"))
            cfg["objective"] = dict(cfg.get("objective") or {}, cov_reference="cov_reference.npz")
        cfg["perturbation"] = dict(parent=self.dir, parent_candidate=cid, parent_objective=obj,
                                   flips=[m.description for m in chosen], seed=seed, jitter=jitter, e_kick=e_kick)
        with open(os.path.join(out_dir, "campaign.json"), "w") as f:
            json.dump(cfg, f, indent=2)
        c2 = Campaign(out_dir)
        src = single_variant(base_copy, "init", c2.vdir)
        with SingleSession(src, cwd=c2.vdir, flags=c2.flags, name="init", objective=c2.objective) as ev1:
            x1 = take(fr.x, ev.keys, ev1.keys, ev1.x0)
            sc = ev1.score(x1)
            sid = c2.ledger.add_structure(base_name, base_copy, ev1.keys, note=f"perturbed from {self.dir} candidate {cid}",
                                          fixed=ev1.fixed_amplitudes)
            pcid = c2.ledger.add_candidate(sid, x1, ev1.keys, move=dict(kind="perturbed", parent=cid,
                                                                       flips=[m.description for m in chosen]), status="verified")
            c2.ledger.record(pcid, sc, "verified", evaluator="single")
        c2.log(f"perturbed baseline candidate {pcid}: verified {sc.objective:.3f} (parent incumbent {obj:.3f}, "
               f"target to recover: {obj:.3f})")
        c2.note(f"recovery benchmark: baseline = {os.path.basename(self.dir)} candidate {cid} ({obj:.1f}) with "
                f"{len(chosen)} relative signs flipped (" + "; ".join(m.description for m in chosen)
                + f"), widths jittered x U({1-jitter:.2f},{1+jitter:.2f}), energies kicked +-{e_kick} MeV, "
                f"polished {polish_nfev} iterations (sign-locked, then released) -> verified {sc.objective:.1f}"
                + (" -- NOTE: the released polish walked back to the incumbent's sign pattern" if same_basin else "")
                + f". Success = sign rounds bring it back to <= {obj:.1f} x (1 + 1e-4). {note}")
        c2.ledger.write_status_md()
        c2.close()
        return dict(dir=out_dir, candidate=pcid, objective=sc.objective, target=obj, flips=[m.description for m in chosen],
                    same_basin=same_basin)
