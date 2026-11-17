"""Per-campaign fit specification: ``spec.json`` (authoritative, machine-read) + ``SPEC.md`` (rendered for the evaluator).

Why (DESIGN.md, finding 2.1): in every campaign of Aug-Nov 2026 the decisions that define
the objective -- which data sets, which scale anchors, systematics, excluded rows, fit window, frozen block, caps,
parametrization -- were made in chat, late, often after a campaign had run.  Where they were written down with an
owner (an MCMC study's "Standing decisions") it worked.  The spec is that document in a form rmfit can read:

* ``rmfit spec check``  compares the .azr (segment table, pairs, flags) and campaign.json with the spec and lists
  every disagreement (nothing is changed);
* ``rmfit init --spec spec.json`` derives flags / policy / search / objective / nshards from it and records the
  spec's sha256 in campaign.json; ``Campaign`` warns when the spec changed since;
* ``rmfit spec render`` writes SPEC.md (standing decisions and open questions at the end);
* the driver and the level-add adoption read ``pre_authorized`` and ``stop_rules``; anything outside the spec goes to
  DECISIONS.md (``rmfit.decisions``), never to a chat gate.

Schema (every section optional; unknown keys are kept and reported, never refused)::

    reaction, compound_nucleus, campaign, owners[]
    engine: {azure2_root, so_md5, rmfit, flags{...}}         flags -> campaign flags
    parametrization: brune | park | standard                   -> use_brune / use_park flags
    model: {fit_window{E_lab_MeV|Ex_MeV: [lo, hi]}, frozen_block, pairs_required[{key,label}], radius_fm,
            bg_poles{energy_MeV, cap}, candidates_file, ...}
    data_sets[]: {file (substring of the data file name) | segments [keys], role anchor|free|excluded|screen,
                  norm {nominal, fixed, prior_sigma (fraction), group}, window [lo, hi] (lab MeV), excluded_rows,
                  systematic_source, covariance_file, target_int {active, energy_kind}, note}
    bounds:    BoundsPolicy fields        objective: ObjectiveSpec fields        search: SearchConfig fields
    nshards, budget{...}, pre_authorized{adopt_level: bool | {max_gamma_kev, min_gain_rel, require_jpi_in_literature,
    forbid_near_same_jpi_kev}, flip_adopt, dataset_off_adopt, drop_rows, change_norm_prior, extend_window,
    add_pairs_on_window_extension, shifts_free}, stop_rules{...}
    standing_decisions[]: {date, who, text}      open_decisions[]: strings
    uq: {enabled, chains, warmup, draws, target_accept, j_saturated, sample_norms, priors{name: [mean, sigma]}, bands[], gates...}
        (rmfit.uq: the drive runs it after its last cycle when enabled)
"""
import dataclasses
import hashlib
import json
import os
import time
from typing import List, Optional

KNOWN_TOP = {"_status", "reaction", "compound_nucleus", "campaign", "owners", "engine", "parametrization", "model",
             "data_sets", "bounds", "objective", "search", "nshards", "budget", "pre_authorized", "stop_rules", "structure",
             "standing_decisions", "open_decisions", "notes", "uq"}
ROLES = {"anchor", "anchor_geometry", "free", "excluded", "excluded_in_window", "screen"}


@dataclasses.dataclass
class Finding:
    severity: str      # HIGH | MEDIUM | LOW | INFO
    code: str
    message: str
    evidence: str = ""

    def line(self):
        return f"- **{self.severity}** `{self.code}`: {self.message}" + (f"  ({self.evidence})" if self.evidence else "")


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


class FitSpec:
    """A loaded spec.json.  ``spec[key]`` / ``spec.get(key, default)`` read the raw dict."""

    def __init__(self, data: dict, path: Optional[str] = None):
        self.data = data
        self.path = path
        self.sha = sha256_of(path) if path and os.path.exists(path) else None

    @classmethod
    def load(cls, path):
        with open(path) as fh:
            return cls(json.load(fh), os.path.abspath(path))

    def __getitem__(self, k):
        return self.data[k]

    def get(self, k, default=None):
        return self.data.get(k, default)

    # -- sections --------------------------------------------------------------------------------
    @property
    def data_sets(self) -> List[dict]:
        return list(self.data.get("data_sets") or [])

    @property
    def pre_authorized(self) -> dict:
        return dict(self.data.get("pre_authorized") or {})

    @property
    def stop_rules(self) -> dict:
        return dict(self.data.get("stop_rules") or {})

    def flags(self) -> Optional[dict]:
        """Engine flags from ``engine.flags``, with ``parametrization`` applied on top."""
        eng = self.data.get("engine") or {}
        flags = dict(eng.get("flags") or {}) if eng.get("flags") else None
        par = (self.data.get("parametrization") or "").lower()
        if par:
            flags = dict(flags or {})
            flags.update({"brune": dict(use_brune=True, use_park=False),
                          "park": dict(use_brune=False, use_park=True),
                          "standard": dict(use_brune=False, use_park=False)}.get(par, {}))
        return flags

    def validate(self) -> List[Finding]:
        """Structural problems of the spec itself (unknown sections, bad roles, non-numeric windows)."""
        out = []
        for k in self.data:
            if k not in KNOWN_TOP:
                out.append(Finding("LOW", "spec-unknown-key", f"unknown top-level key {k!r} (kept, not read by rmfit)"))
        for i, ds in enumerate(self.data_sets):
            if not (ds.get("file") or ds.get("segments")):
                out.append(Finding("MEDIUM", "spec-dataset-id", f"data_sets[{i}] names neither file nor segments"))
            role = ds.get("role", "free")
            if role not in ROLES:
                out.append(Finding("MEDIUM", "spec-role", f"data_sets[{i}] role {role!r} not in {sorted(ROLES)}"))
            w = ds.get("window")
            if w is not None and not (isinstance(w, (list, tuple)) and len(w) == 2 and all(_isnum(v) for v in w)):
                out.append(Finding("LOW", "spec-window", f"data_sets[{i}] window {w!r} is not [lo, hi] (ignored by check)"))
        par = (self.data.get("parametrization") or "").lower()
        if par and par not in ("brune", "park", "standard"):
            out.append(Finding("MEDIUM", "spec-parametrization", f"parametrization {par!r} unknown"))
        for sec, cls in (("bounds", "rmfit.local.BoundsPolicy"), ("objective", "rmfit.engine.ObjectiveSpec"),
                         ("search", "rmfit.search.SearchConfig")):
            d = self.data.get(sec)
            if d:
                names = _field_names(cls)
                bad = [k for k in d if k not in names]
                if bad:
                    out.append(Finding("LOW", f"spec-{sec}-keys", f"{sec} keys not in {cls.split('.')[-1]}: {bad}"))
        return out

    # -- derived campaign configuration ---------------------------------------------------------------
    def derive_campaign_cfg(self) -> dict:
        """What ``rmfit init`` takes from the spec: flags, policy, search, objective, nshards (only the sections given)."""
        out = {}
        fl = self.flags()
        if fl:
            out["flags"] = fl
        for sec, cls in (("bounds", "rmfit.local.BoundsPolicy"), ("search", "rmfit.search.SearchConfig"),
                         ("objective", "rmfit.engine.ObjectiveSpec")):
            d = self.data.get(sec)
            if d:
                names = _field_names(cls)
                out[{"bounds": "policy"}.get(sec, sec)] = {k: v for k, v in d.items() if k in names}
        if self.data.get("nshards"):
            out["nshards"] = int(self.data["nshards"])
        return out

    # -- rendering -------------------------------------------------------------------------------------
    def render_md(self) -> str:
        d = self.data
        L = [f"# Fit specification: {d.get('reaction', '?')} -- {d.get('campaign', '')}", ""]
        if d.get("_status"):
            L += [f"_{d['_status']}_", ""]
        L += [f"Owners: {', '.join(d.get('owners') or ['?'])}.  Compound nucleus {d.get('compound_nucleus', '?')}.  "
              f"Parametrization {d.get('parametrization', '?')}."]
        if self.path:
            L.append(f"Source `{os.path.basename(self.path)}` sha256 {self.sha[:12] if self.sha else '?'} (rendered {time.strftime('%Y-%m-%d %H:%M')}).")
        eng = d.get("engine") or {}
        if eng:
            L += ["", "## Engine", "", f"- AZURE2_ROOT `{eng.get('azure2_root', '?')}`, .so md5 `{eng.get('so_md5', '?')}`, rmfit {eng.get('rmfit', '?')}",
                  f"- flags `{json.dumps(eng.get('flags') or {})}`"]
        m = d.get("model") or {}
        if m:
            L += ["", "## Model", ""]
            for k, v in m.items():
                L.append(f"- **{k}**: {_md(v)}")
        if self.data_sets:
            L += ["", "## Data sets", "", "| file / segments | role | norm | window | notes |", "|---|---|---|---|---|"]
            for ds in self.data_sets:
                ident = ds.get("file") or f"segments {ds.get('segments')}"
                L.append(f"| {ident} | {ds.get('role', 'free')} | {_md(ds.get('norm', ''))} | {_md(ds.get('window', ''))} | "
                         f"{_md({k: v for k, v in ds.items() if k not in ('file', 'segments', 'role', 'norm', 'window')})} |")
        for sec in ("bounds", "objective", "search", "budget", "pre_authorized", "stop_rules", "structure"):
            if d.get(sec):
                L += ["", f"## {sec}", "", "```json", json.dumps(d[sec], indent=1), "```"]
        if d.get("standing_decisions"):
            L += ["", "## Standing decisions (do not reopen without the owners agreeing)", ""]
            for sd in d["standing_decisions"]:
                L.append(f"- {sd.get('date', '?')} ({sd.get('who', '?')}): {sd.get('text', '')}")
        if d.get("open_decisions"):
            L += ["", "## Open decisions (for the evaluator; see DECISIONS.md)", ""]
            for od in d["open_decisions"]:
                L.append(f"- {od}")
        return "\n".join(L) + "\n"

    # -- checks against a model / campaign ---------------------------------------------------------------
    def check(self, azr_path: str, campaign_cfg: Optional[dict] = None) -> List[Finding]:
        """Disagreements between the spec and the .azr (segments, pairs, flags) / campaign.json.  Nothing is changed."""
        from .backend import SegmentSet
        out = self.validate()
        segs = SegmentSet.from_file(azr_path)
        # data sets
        for i, ds in enumerate(self.data_sets):
            matched = _match_segments(segs, ds)
            ident = ds.get("file") or f"segments {ds.get('segments')}"
            if not matched:
                if ds.get("role") in ("screen",):
                    continue
                out.append(Finding("MEDIUM", "spec-dataset-missing", f"data_sets[{i}] {ident!r}: no segment of the .azr matches"))
                continue
            role = ds.get("role", "free")
            act = [s for s in matched if s.active]
            if role in ("excluded",) and act:
                out.append(Finding("HIGH", "spec-excluded-active", f"{ident!r} is EXCLUDED in the spec but segments "
                                   f"{[s.key for s in act]} are active"))
            if role in ("anchor", "free", "anchor_geometry") and not act:
                out.append(Finding("MEDIUM", "spec-role-inactive", f"{ident!r} has role {role} but every matching segment is inactive"))
            norm = ds.get("norm") or {}
            hits = {}            # (severity, code, message) -> [segment keys]: one finding per data set, not per segment

            def hit(sev, code, msg, s):
                hits.setdefault((sev, code, msg), []).append(s.key)
            for s in act:
                if norm.get("fixed") is True and s.vary_norm:
                    hit("HIGH", "spec-norm-fixed", "spec says the norm is FIXED, the .azr varies it", s)
                if norm.get("fixed") is False and not s.vary_norm:
                    hit("MEDIUM", "spec-norm-free", "spec says the norm is free, the .azr fixes it (t0 method with a covariance file?)", s)
                if _isnum(norm.get("nominal")) and abs(s.norm - float(norm["nominal"])) > 1e-6:
                    hit("MEDIUM", "spec-norm-nominal", f"nominal norm in the .azr != spec {norm['nominal']}", s)
                if _isnum(norm.get("prior_sigma")) and s.vary_norm:
                    want = 100.0 * float(norm["prior_sigma"])
                    if abs(s.norm_error - want) > 0.05 * max(want, 1e-9) + 1e-6:
                        hit("MEDIUM", "spec-norm-prior", f"normError in the .azr != spec prior_sigma {norm['prior_sigma']} = {want:.3g} %", s)
                w = ds.get("window")
                if isinstance(w, (list, tuple)) and len(w) == 2 and all(_isnum(v) for v in w):
                    lo, hi = float(w[0]), float(w[1])
                    if s.energy_min < lo - 1e-6 or s.energy_max > hi + 1e-6:
                        hit("MEDIUM", "spec-window", f".azr window outside the spec window [{lo}, {hi}]", s)
                if ds.get("covariance_file") and (s.covariance or "") and os.path.basename(s.covariance) != os.path.basename(ds["covariance_file"]):
                    hit("LOW", "spec-covariance", f"covariance file differs from the spec's {ds['covariance_file']!r}", s)
                if ds.get("covariance_file") and not s.covariance:
                    hit("MEDIUM", "spec-covariance", "spec names a covariance file, the segment has none", s)
            for (sev, code, msg), ks in hits.items():
                out.append(Finding(sev, code, f"{ident!r}: {msg}", f"segments {_ranges(ks)} ({len(ks)} of {len(act)})"))
        # pairs
        pr = (self.data.get("model") or {}).get("pairs_required")
        if pr:
            try:
                from .backend import AzrModel
                keys = set(pair_lines(AzrModel.from_file(azr_path)))
                for p in pr:
                    k = p.get("key") if isinstance(p, dict) else (p if isinstance(p, int) else None)
                    if k is not None and int(k) not in keys:
                        out.append(Finding("HIGH", "spec-pair-missing", f"pair {k} ({p.get('label', '') if isinstance(p, dict) else ''}) required by the spec is not in the .azr"))
            except Exception as err:   # pyazr unavailable: say so, do not fail
                out.append(Finding("INFO", "spec-pairs-unchecked", f"pairs_required not checked ({err})"))
        # engine / flags / config vs campaign.json
        # 2026-11-09 (from-scratch demo): model.radius_fm {pair: fm} vs the channel radius the .azr carries -- the spec said 6.0 fm
        # (copied from a readme), the file 3.5 fm; nothing compared them
        rad = (self.data.get("model") or {}).get("radius_fm") or {}
        if rad:
            try:
                from .backend import AzrModel
                pl = pair_lines(AzrModel.from_file(azr_path))
                for key, want in rad.items():
                    c = pl.get(int(key))
                    if c is None:
                        out.append(Finding("MEDIUM", "spec-radius", f"model.radius_fm names pair {key}, which the .azr does not have"))
                    elif not getattr(c, "is_photon", False) and abs(float(c.channel_radius) - float(want)) > 1e-6:
                        out.append(Finding("HIGH", "spec-radius", f"pair {key}: spec radius_fm {float(want):g} fm, the .azr has {float(c.channel_radius):g} fm"))
            except Exception as err:  # noqa: BLE001
                out.append(Finding("LOW", "spec-radius", f"radius check skipped ({err})"))
        if campaign_cfg:
            eng = self.data.get("engine") or {}
            ce = campaign_cfg.get("engine") or {}
            if eng.get("so_md5") and ce.get("so_md5") and eng["so_md5"] != ce["so_md5"]:
                out.append(Finding("HIGH", "spec-engine", f"spec engine .so md5 {eng['so_md5'][:8]} != campaign pin {ce['so_md5'][:8]}"))
            fl = self.flags()
            if fl:
                cf = campaign_cfg.get("flags") or {}
                for k, v in fl.items():
                    if k in cf and cf[k] != v:
                        out.append(Finding("HIGH", "spec-flags", f"flag {k}: spec {v}, campaign {cf[k]}"))
            for sec, ckey in (("bounds", "policy"), ("search", "search"), ("objective", "objective")):
                d = self.data.get(sec) or {}
                c = campaign_cfg.get(ckey) or {}
                for k, v in d.items():
                    if k in c and c[k] != v and not (isinstance(v, (list, tuple)) and list(v) == list(c[k])):
                        out.append(Finding("MEDIUM", f"spec-{sec}", f"{sec}.{k}: spec {v!r}, campaign.json {c[k]!r}"))
            if campaign_cfg.get("spec_sha") and self.sha and campaign_cfg["spec_sha"] != self.sha:
                out.append(Finding("LOW", "spec-changed", "spec.json changed since init (sha differs from campaign.json spec_sha)"))
        return out


def _field_names(dotted):
    mod, _, cls = dotted.rpartition(".")
    import importlib
    c = getattr(importlib.import_module(mod), cls)
    return {f.name for f in dataclasses.fields(c)}


def _isnum(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _ranges(keys):
    """[23, 24, 25, 30] -> '23-25, 30'."""
    ks = sorted(set(int(k) for k in keys))
    out, start, prev = [], None, None
    for k in ks + [None]:
        if start is None:
            start = prev = k
        elif k is not None and k == prev + 1:
            prev = k
        else:
            out.append(f"{start}-{prev}" if prev != start else f"{start}")
            start = prev = k
    return ", ".join(out)


def pair_lines(model) -> dict:
    """``{pair key: first channel line naming it}`` -- ``AzrModel.pair_lines`` where it exists (pyazr >= 2026-10-29),
    else built from the channels (older frozen snapshots predate it)."""
    if hasattr(model, "pair_lines"):
        return model.pair_lines()
    out = {}
    for lv in model.levels:
        for c in lv.channels:
            out.setdefault(c.pair, c)
    return dict(sorted(out.items()))


def _md(v):
    if isinstance(v, (dict, list)):
        return "`" + json.dumps(v) + "`"
    return str(v)


def _match_segments(segs, ds):
    if ds.get("segments"):
        want = {int(k) for k in ds["segments"]}
        return [s for s in segs if s.key in want]
    sub = ds.get("file")
    if not sub:
        return []
    return [s for s in segs if sub in os.path.basename(s.data_file) or sub in s.data_file]


def write_findings_md(path, title, findings: List[Finding], header_lines=()):
    sev_order = {"HIGH": 0, "MEDIUM": 1, "LOW": 2, "INFO": 3}
    fs = sorted(findings, key=lambda f: sev_order.get(f.severity, 9))
    counts = {s: sum(1 for f in fs if f.severity == s) for s in sev_order}
    L = [f"# {title}", "", f"{time.strftime('%Y-%m-%d %H:%M')}  " + ", ".join(f"{k} {v}" for k, v in counts.items()), ""]
    L += list(header_lines)
    L += [f.line() for f in fs] or ["- no findings"]
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    return counts
