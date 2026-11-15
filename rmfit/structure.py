"""Structure moves: add a level, remove a level, switch a dataset off.

Each one is a new model file (a new *structure* in the ledger) whose free
vector shares every key with the parent except the level being added or
removed, so parameters move between the two by key.  Levels are always
appended at the end of the file and removed by deactivation (widths zeroed
and fixed, line kept), which keeps the engine's level numbering -- and hence
every existing key -- stable.

The add-level test follows the archive's hard-won rules: pin the candidate's
energy and magnitudes while the interference sign is scanned (a free level
relocates itself into a background-pole basin and the comparison is void),
judge on a local window *and* the global total, and classify the outcome
(accept / collapse / blowup / live-harmful / relocated / reverts) instead of
reading a single chi-squared.
"""
import dataclasses
import json
import os
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

from . import AZURE2_ROOT  # noqa: F401
from .engine import SingleSession, Sharded, EvalError
from .keys import ParamKey, take, jpi_string
from .local import BoundsPolicy, default_bounds, free_mask, polish, sanity, gn_step
from .model import single_variant
from .windows import window_variant
from pyazr import AzrModel

# The half-integer-J list (e.g. a 17O compound nucleus).  Until 2026-11-04 it was the DEFAULT J^pi list of levelscan,
# structure add and scan_jpi_on_window; on an integer-J compound nucleus (e.g. 14N) that silently built unphysical dummy levels
# (J = 1/2 with alpha0 l = 3, s = 3).  The default is now the model's own groups (default_jpis), and every
# requested J^pi is checked against the model's particle pairs (resolve_jpis).  Kept for explicit use only.
JPI_GROUPS = ["1/2+", "1/2-", "3/2+", "3/2-", "5/2+", "5/2-", "7/2+", "7/2-", "9/2+", "9/2-", "11/2+", "11/2-"]


@dataclass
class LevelCandidate:
    """A level to try: J^pi (or None to scan), energy (MeV), seed channels
    ``[{pair, L, S, gamma}]`` (physical widths in eV, as the <levels> field wants
    them; omitted channels are seeded from the group's structure at a small
    fraction of the Wigner limit), and provenance."""
    energy: float
    jpi: Optional[str] = None
    channels: List[dict] = field(default_factory=list)
    source: str = ""
    note: str = ""
    energy_fixed: bool = True

    def to_json(self):
        return dataclasses.asdict(self)


def jpi_parse(jpi: str):
    j, par = jpi[:-1], 1 if jpi.endswith("+") else -1
    J = float(j.split("/")[0]) / 2 if "/" in j else float(j)
    return J, par


def lowest_L(J: float, parity: int, j1: float, pi1: int, j2: float, pi2: int):
    """The lowest partial waves (L, S) allowed for (J, parity) with particle
    spins/parities j1, j2 -- the natural seed channels of a new level.

    Coupling J = L + S needs J - S to be an INTEGER (2026-11-04): without that check
    J = 1/2 with S = 3 passed the triangle test at L = 3, so a half-integer level
    received channels in an integer-spin system.  Empty when J cannot couple."""
    out = []
    smin, smax = abs(j1 - j2), j1 + j2
    S = smin
    while S <= smax + 1e-9:
        if abs((J - S) - round(J - S)) > 1e-6:
            S += 1
            continue
        for L in range(0, 12):
            if abs(J - S) - 1e-9 <= L <= J + S + 1e-9 and pi1 * pi2 * (-1) ** L == parity:
                out.append((L, S))
                break
        S += 1
    return out


def _as_model(model):
    return model if isinstance(model, AzrModel) else AzrModel.from_file(model)


def model_pairs(model) -> dict:
    """``{pair: (j1, pi1, j2, pi2)}`` of every PARTICLE pair a channel of the model uses
    (photon channels skipped).  ``model`` is an AzrModel or a path."""
    mdl = _as_model(model)
    out = {}
    for lv in mdl.levels:
        for c in lv.channels:
            if c.pair in out or bool(getattr(c, "is_photon", False)):
                continue
            out[c.pair] = (float(c.J1), int(c.parity1), float(c.J2), int(c.parity2))
    return out


def jpi_allowed(model, jpi: str) -> bool:
    """Can a level of this J^pi couple to at least one particle pair of the model?"""
    J, par = jpi_parse(jpi)
    return any(lowest_L(J, par, *sp) for sp in model_pairs(model).values())


def default_jpis(model) -> List[str]:
    """The J^pi groups that already hold an ACTIVE level of the model (a background
    pole counts), ordered by J then parity (+ first).  The default J^pi list of
    levelscan / structure add / scan_jpi_on_window since 2026-11-04: a dummy of an
    existing group reuses that group's channels, so the template variant adds no
    partial wave the model does not already have (a new group adds channels whose
    hard-sphere terms change the elastic channels even at zero width: a clear change
    in one model)."""
    mdl = _as_model(model)
    seen = {}
    for lv in mdl.levels:
        if lv.active:
            seen[lv.jpi] = (float(lv.J), -int(lv.parity))
    return sorted(seen, key=lambda j: seen[j])


def resolve_jpis(model, jpis: Optional[Sequence[str]] = None, log=None) -> List[str]:
    """The J^pi list to use: ``jpis`` as given, or ``default_jpis`` when empty/None.
    Raises ValueError if any J^pi cannot couple to any particle pair of the model
    (e.g. half-integer J in an integer-spin compound nucleus), instead of skipping it
    or building an unphysical level.  ``log`` (optional) gets one line naming the list
    and any requested group the model does not have yet (it adds partial waves)."""
    mdl = _as_model(model)
    given = bool(jpis)
    out = list(jpis) if given else default_jpis(mdl)
    if not out:
        raise ValueError("no J^pi to scan: the model has no active level and none was given (--jpi)")
    bad = [j for j in out if not jpi_allowed(mdl, j)]
    if bad:
        spins = ", ".join(f"pair {p}: {sp[0]:g}{'+' if sp[1] > 0 else '-'} + {sp[2]:g}{'+' if sp[3] > 0 else '-'}"
                          for p, sp in sorted(model_pairs(mdl).items()))
        raise ValueError(f"J^pi {bad} cannot couple to any particle pair of the model ({spins}); "
                         f"the model's own groups are {default_jpis(mdl)}")
    if log is not None:
        new = [j for j in out if not mdl.find(jpi=j)]
        log(f"  J^pi groups ({'given' if given else 'the model own, default'}): {' '.join(out)}"
            + (f"; NEW groups {new} add partial waves the model does not have (hard-sphere terms change "
               f"the elastic channels even at zero width)" if new else ""))
    return out


def omitted_open_channels(existing, named, energy: float, seed_ev: float = 1.0,
                          dead_channel_ev: float = 1.0, margin: float = 10.0) -> List[dict]:
    """The group's particle channels NOT in ``named`` that are open and comfortably alive
    at ``energy``: Wigner-limit width ``Gamma_W = 2 P_L * 3 hbar^2 / (2 mu a^2)`` at least
    ``margin * dead_channel_ev`` (computed OFFLINE with pyazr's penetrability, so nothing
    returned here is anywhere near the ``Gamma_W < dead_channel_ev`` rule by which
    ``run_level_add`` zeroes a channel through its bounds -- an exactly-zero FREE amplitude
    does not survive bake's write/re-read; two keys were once lost).  Each comes
    back free with a tiny seed ``min(seed_ev, 1% of Gamma_W)`` (theta^2 <= 0.01): enough to
    leave the zero stationary point, not a guess at the strength."""
    from pyazr.transform import penetrability, HBARC, UCONV
    have = {(int(c["pair"]), int(c["L"]), round(float(c["S"]), 3)) for c in named}
    out, seen = [], set()
    for lv in existing:
        for c in lv.channels:
            key = (int(c.pair), int(round(c.L)), round(float(c.S), 3))
            if key in have or key in seen:
                continue
            seen.add(key)
            if bool(getattr(c, "is_photon", False)):
                continue
            m1, m2 = c._get("m1", float), c._get("m2", float)
            e_rel = float(energy) - (c._get("sepE", float) + c._get("e2", float))
            if e_rel <= 0.0 or m1 <= 0.0 or m2 <= 0.0:
                continue
            mu, a = m1 * m2 / (m1 + m2), c._get("chRad", float)
            try:
                P = penetrability(key[1], c._get("z1", int), c._get("z2", int), mu, a, e_rel)
            except Exception:
                continue
            gw_ev = 2.0 * P * 3.0 * HBARC ** 2 / (2.0 * mu * UCONV * a * a) * 1e6
            if not (gw_ev >= margin * dead_channel_ev):
                continue
            out.append(dict(pair=key[0], L=key[1], S=float(c.S), gamma=min(seed_ev, 0.01 * gw_ev),
                            fixed=False, wigner_ev=gw_ev))
    return out


def effective_k(params: Sequence[dict], seeded, frac: float = 0.01, theta2_floor: float = 0.05) -> dict:
    """How many parameters a new level is CHARGED for in tau = s^2 k ln N (2026-09-20, the evaluator).

    ``params``: the level's free, non-dead parameters, each ``{kind, pair, L, S, theta2,
    gamma_ev}`` with ``gamma_ev ~ theta2 * Gamma_W`` (energies: ``{kind: "energy"}``).
    ``seeded``: the ``(pair, L, S)`` that ``omitted_open_channels`` added.  Channels of the scan
    pattern always count.  A SEEDED channel counts only if the fit actually used it: its partial
    width is >= ``frac`` of the level's summed width, or its theta^2 >= ``theta2_floor``.  (A
    theta^2 floor alone would be wrong: the seed is min(1 eV, 1% Gamma_W), so on a channel with
    Gamma_W = 6 MeV a final 5 keV is theta^2 ~ 1e-3 and still dominates a narrow level.)
    Counting all ~22 seeded channels put tau above every level accepted in
    the low-energy campaign that month, although most of them end where they were seeded.
    With nothing seeded, k_effective == k_all == the old k_new."""
    seeded = {(int(a), int(b), round(float(c), 3)) for a, b, c in seeded}
    total = sum(float(q.get("gamma_ev", 0.0)) for q in params if q.get("kind") == "width")
    k_pattern, used = 0, []
    for q in params:
        key = None
        if q.get("kind") == "width" and q.get("pair") is not None:
            key = (int(q["pair"]), int(q["L"]), round(float(q["S"]), 3))
        if key is None or key not in seeded:
            k_pattern += 1
            continue
        g, t2 = float(q.get("gamma_ev", 0.0)), float(q.get("theta2", 0.0))
        if (total > 0 and g >= frac * total) or t2 >= theta2_floor:
            used.append(key)
    return dict(k_pattern=k_pattern, k_effective=k_pattern + len(used), k_all=len(params), used=used)


def add_level_variant(base_azr: str, path: str, cand: LevelCandidate, jpi: str,
                      output_dir: Optional[str] = None, seed_fraction: float = 0.05,
                      seed_omitted: Optional[bool] = None) -> dict:
    """Write ``path`` = ``base_azr`` with the candidate appended as a level of
    ``jpi``.  Channels named in the candidate get their gamma; the group's
    other channels are inert.  If ``jpi`` has no level yet, the lowest
    partial waves into pairs 1 and 2 define the group with small seeds.

    ``seed_omitted`` (2026-09-20; default None = read the environment variable
    ``RMFIT_SEED_OMITTED``, unset/0 = the behaviour above): also add the group's OTHER
    open, non-dead particle channels, free, with a tiny seed (``omitted_open_channels``).
    Without it a level found by a sparse scan pattern ("channel 7 only") has that one
    channel free for ever -- pyazr add_level writes the rest gamma = 0, FIXED -- so it can
    never carry inelastic strength and its accept/reject verdict is reached with those
    channels frozen (in one window, 5 of 7 added levels had an inelastic channel fixed at zero).  The
    extra channels count in k_new, i.e. they raise tau: the conservative reading.
    """
    if seed_omitted is None:
        seed_omitted = os.environ.get("RMFIT_SEED_OMITTED", "0").strip().lower() not in ("", "0", "false", "no")
    mdl = AzrModel.from_file(base_azr)
    resolve_jpis(mdl, [jpi])                      # raises for a J^pi no particle pair can couple to (2026-11-04)
    J, par = jpi_parse(jpi)
    existing = mdl.find(jpi=jpi)
    if existing:
        allowed = [(c.pair, c.L, c.S) for c in existing[0].channels]
        chans = []
        for spec in cand.channels:
            key = (int(spec["pair"]), int(spec["L"]), float(spec["S"]))
            match = [a for a in allowed if a[0] == key[0] and a[1] == key[1] and abs(a[2] - key[2]) < 1e-6]
            if match:
                chans.append(dict(pair=key[0], L=key[1], S=key[2], gamma=float(spec.get("gamma", spec.get("gamma_eV", spec.get("gamma_ev", 0.0)))), fixed=False))
        if not chans:
            # seed the lowest allowed wave of pairs 1 and 2 present in the group
            for pair in (1, 2):
                ch = sorted([a for a in allowed if a[0] == pair], key=lambda a: a[1])
                if ch:
                    chans.append(dict(pair=ch[0][0], L=ch[0][1], S=ch[0][2], gamma=1000.0, fixed=False))
    else:
        tmpl = mdl.levels[0].channels[0]
        chans = []
        for pair in (1, 2):
            c = mdl._pair_template(pair)
            for L, S in lowest_L(J, par, c.J1, c.parity1, c.J2, c.parity2)[:1]:
                g = next((float(s.get("gamma", s.get("gamma_eV", s.get("gamma_ev", 0.0)))) for s in cand.channels
                          if int(s["pair"]) == pair and int(s["L"]) == L), 1000.0)
                chans.append(dict(pair=pair, L=L, S=S, gamma=g, fixed=False))
    seeded = []
    if seed_omitted and existing:
        seeded = omitted_open_channels(existing, chans, cand.energy)
        chans = chans + [{k: v for k, v in c.items() if k != "wigner_ev"} for c in seeded]
    lv = mdl.add_level(J=J, parity=par, energy=cand.energy, channels=chans, level_fixed=cand.energy_fixed)
    if output_dir is not None:
        mdl.set_output_dir(output_dir)
    mdl.write(path)
    return dict(path=path, jpi=jpi, energy=cand.energy, channels=chans,
                level=f"{jpi}#{len(existing) + 1}", seeded_omitted=seeded)


def deactivate_level_variant(base_azr: str, path: str, jpi: str, energy: float,
                             output_dir: Optional[str] = None, tol=1e-3, nucleus: Optional[str] = None) -> dict:
    """``nucleus`` (``"23Ne"``, from a multi-nucleus level label ``"1/2+#3@23Ne"``) limits the match to that compound
    nucleus: two nuclei can have a level of one J^pi at the same energy."""
    mdl = AzrModel.from_file(base_azr)
    victims = mdl.find(jpi=jpi, energy=energy, tol=tol)
    if nucleus:
        victims = [lv for lv in victims if lv.channels and str(mdl.nucleus_of_pair(lv.channels[0].pair)) == nucleus]
    if not victims:
        raise KeyError(f"no level {jpi}{'@' + nucleus if nucleus else ''} at {energy} MeV")
    for lv in victims:
        lv.set_fixed(True)                 # the energy leaves the free vector too
        for c in lv.channels:
            c.gamma = 0.0
            c.channel_fixed = True
    if output_dir is not None:
        mdl.set_output_dir(output_dir)
    mdl.write(path)
    return dict(path=path, jpi=jpi, energy=energy, removed=[str(v) for v in victims])


def dataset_off_variant(base_azr: str, path: str, file_substr: str, output_dir=None) -> dict:
    mdl = AzrModel.from_file(base_azr)
    n = mdl.set_segment_active(file_substr, False)
    if output_dir is not None:
        mdl.set_output_dir(output_dir)
    mdl.write(path)
    return dict(path=path, file=file_substr, segments=n)


# -- the pinned J^pi scan on a window ------------------------------------------------

@dataclass
class ScanResult:
    jpi: str
    signs: str
    objective: float
    delta: float                 # vs the window baseline
    x: Optional[np.ndarray]
    status: str
    note: str = ""
    keys: Optional[list] = None


def split_evidence(new_energy: float, new_gamma_ev: float, neighbours: Sequence[dict],
                   transfer_frac: float = 0.5, own_frac: float = 0.05,
                   max_sep_kev: float = 1000.0) -> Optional[dict]:
    """Is a level add really a SPLIT of an existing level of the same J^pi?

    ``neighbours``: ``[{label, energy (MeV), gamma_before_ev, gamma_after_ev}]`` -- the other levels
    of the new level's J^pi with their summed physical widths at the incumbent and after the
    released polish.  A neighbour is a split partner when BOTH hold:

    * the pair is unresolved: separation <= max(Gamma_new, Gamma_neighbour before);
    * width moved from it to the new level: it lost >= ``transfer_frac`` x Gamma_new, AND that loss
      is >= ``own_frac`` of its own width.  (Without the second test a 2 MeV-wide background level
      2.2 MeV away, which shed 16 keV = 0.7% of itself in the same refit, was named the partner of
      an 11 keV level -- found by replaying the real case before the rule was installed.)

    * and it is within ``max_sep_kev`` (default 1 MeV) of the new level.  2026-10-16: a 50 MeV background
      pole has a "physical width" of tens of MeV at the reference energy, so the overlap test called it
      unresolved from 40 MeV away and its 7.7 MeV of re-shuffled width passed the transfer test: the
      recall test labelled a re-added broad level a
      "split" of a background pole at 50 MeV, and another add likewise named its group's background pole.
      The caller also leaves out levels at or above the campaign's background energy.

    Overlap alone is not enough -- broad same-J^pi levels legitimately overlap (a 360 keV
    level 54 keV from a same-J^pi neighbour was clearly real).  2026-09-21: a
    level added 5.6 keV from a 44 keV level of the same J^pi, which gave up most of its width; 96%
    of the gain came from one 1-keV-step dataset.  One lineshape described by two levels.
    Returns the NEAREST partner as ``{label, sep_kev, lost_ev, lost_frac}`` or ``None``."""
    if not new_gamma_ev or new_gamma_ev <= 0:
        return None
    best = None
    for nb in neighbours:
        sep_ev = abs(float(nb["energy"]) - float(new_energy)) * 1e6
        g0, g1 = float(nb["gamma_before_ev"]), float(nb["gamma_after_ev"])
        lost = g0 - g1
        if (sep_ev <= max(new_gamma_ev, g0) and sep_ev <= max_sep_kev * 1e3
                and lost >= transfer_frac * new_gamma_ev and lost >= own_frac * g0):
            row = dict(label=str(nb["label"]), sep_kev=float(sep_ev / 1e3), lost_ev=float(lost),
                       lost_frac=float(lost / new_gamma_ev))
            if best is None or row["sep_kev"] < best["sep_kev"]:
                best = row
    return best


def classify(before_obj, after_obj, level_theta2: Dict[str, float], moved_energy: float,
             cap: float, tau: float, rho: Optional[float],
             at_bound: Sequence[str] = (), e_window: float = 0.05,
             level_gamma_ev: Optional[float] = None, collapse_width_ev: float = 1.0,
             narrow_width_ev: Optional[float] = None, splits: bool = False,
             broad_width_ev: Optional[float] = None) -> str:
    """Outcome of a level test.  ``at_bound``: the new level's width keys that the
    released polish left ON a bound.  2026-09-16: two levels
    were accepted with large gains while three or four channels sat
    exactly at the theta^2 cap (1.0, then 0.5) and the rest had collapsed to ~1e-35 --
    broad absorbers whose Brune transformation fails ("Denominator less than zero"),
    so the model could not even be baked.  The blowup test only fires ABOVE the cap;
    a level pinned AT the cap, or with theta^2 within 0.1% of it, is now
    "unphysical" and is not adopted.  Tightening the cap does not help (the polish
    drives the amplitudes to whatever cap exists), the gate does."""
    if not np.isfinite(after_obj):
        return "crashed"
    # 2026-09-20: theta^2 alone mistakes a NARROW level for a vanished one when the channel's
    # Wigner limit is large (a level with max theta^2 5.6e-5 was ~35 eV wide, Gamma_W 0.3-1
    # MeV, and was called "collapse").  When the level's summed physical width is
    # known, "collapse" additionally requires it to be below collapse_width_ev.
    if max(level_theta2.values(), default=0.0) < 1e-4 and (level_gamma_ev is None or level_gamma_ev < collapse_width_ev):
        return "collapse"
    if max(level_theta2.values(), default=0.0) > cap * (1 + 1e-6):
        return "blowup"
    if at_bound or max(level_theta2.values(), default=0.0) >= cap * (1 - 1e-3):
        return "unphysical"
    # 2026-09-19: the released stage now frees the new level's energy inside +-e_window (MeV);
    # a level that ends AT the edge of that window wants to be somewhere else = relocated.
    # (Before, the energy could not move at all and this test never fired; the old fixed
    # 50 keV would have rejected a legitimate 60 keV move of a 300 keV-wide level.)
    if abs(moved_energy) >= 0.98 * e_window:
        return "relocated"
    if after_obj > before_obj - tau:
        return "live-harmful" if after_obj > before_obj else "marginal"
    if rho is not None and rho > 0.5:
        return "redistributes"
    # 2026-09-21: passes everything, but it took its width from an unresolved neighbour of the
    # same J^pi (split_evidence): one lineshape described by two levels -- the evaluator's call
    if splits:
        return "splits"
    # passes everything, but the campaign only adopts levels at least narrow_width_ev wide
    # (BoundsPolicy.min_level_width_ev; broad structure first, narrow alignment at the end stage)
    if narrow_width_ev and level_gamma_ev is not None and level_gamma_ev < narrow_width_ev:
        return "narrow"
    # 2026-09-22 (the evaluator): the mirror image -- passes everything but ends WIDER than
    # broad_width_ev (BoundsPolicy.max_level_width_ev): a level fitting the average cross
    # section, not a resonance shape.  Not adopted; the window is rebuilt with narrow levels.
    if broad_width_ev and level_gamma_ev is not None and level_gamma_ev > broad_width_ev:
        return "broad"
    return "accept"


def _open(azr, workdir, flags, nshards, name, objective=None):
    """A sharded evaluator when the window is big (the dense 17O data sit at
    the top of the range, so a +-0.5 MeV window there still holds ~15k points
    and a single-thread Jacobian costs minutes), else one session."""
    if nshards and nshards > 1:
        return Sharded(azr, os.path.join(workdir, f"shards_{name}"), nshards=nshards, flags=flags, objective=objective)
    return SingleSession(azr, cwd=workdir, flags=flags, name=name, objective=objective)


def scan_jpi_on_window(base_azr: str, workdir: str, cand: LevelCandidate, flags: dict,
                       jpis: Optional[Sequence[str]] = None, half_width: float = 0.5,
                       policy: BoundsPolicy = BoundsPolicy(), x_base=None, base_keys=None,
                       max_nfev=60, log=print, nshards: int = 1, polish_top: int = 2,
                       objective=None) -> List[ScanResult]:
    """For every J^pi: window model around the candidate energy, level added
    with pinned energy; every sign pattern of its seeded channels is screened
    with one damped Gauss-Newton step on the level's columns (+ norms) and
    the best ``polish_top`` patterns are polished sign-locked.  Returns the
    best result per J^pi, ordered by objective.  ``x_base``/``base_keys``
    carry the incumbent's parameters into each variant by key.
    ``jpis`` None/empty = the model's own groups (resolve_jpis, 2026-11-04)."""
    jpis = resolve_jpis(base_azr, jpis, log=log)
    os.makedirs(workdir, exist_ok=True)
    base_dir = os.path.dirname(os.path.abspath(base_azr))
    from .model import link_data_dirs
    link_data_dirs(base_azr, workdir)
    ex_lo, ex_hi = cand.energy - half_width, cand.energy + half_width
    win = os.path.join(workdir, "window_base.azr")
    info = window_variant(base_azr, win, ex_lo, ex_hi, level_margin=0.3,
                          output_dir=os.path.join(workdir, "output_window_base"))
    log(f"  window [{ex_lo:.2f}, {ex_hi:.2f}] MeV: {info['segments_kept']} segments kept, "
        f"{info['segments_dropped']} dropped, {len(info['levels_fixed'])} levels fixed")
    with _open(win, workdir, flags, nshards, "window_base", objective=objective) as ev:
        x0 = ev.x0 if x_base is None else take(x_base, base_keys, ev.keys, ev.x0)
        b = default_bounds(ev, x0, policy)
        f0 = polish(ev, x0, bounds=b, policy=policy, max_nfev=max_nfev, log=log)
        obj0 = f0.objective
        x_win, keys_win = f0.x, ev.keys
    log(f"  window baseline objective {obj0:.3f}")
    results = []
    for jpi in jpis:
        tag = jpi.replace("/", "_")
        path = os.path.join(workdir, f"window_{tag}.azr")
        try:
            info = add_level_variant(win, path, cand, jpi, output_dir=os.path.join(workdir, f"output_{tag}"))
        except Exception as err:
            results.append(ScanResult(jpi, "", np.inf, np.inf, None, "crashed", f"variant: {err}"))
            continue
        try:
            with _open(path, workdir, flags, nshards, f"win_{tag}", objective=objective) as ev:
                x = take(x_win, keys_win, ev.keys, ev.x0)
                lvl = info["level"]
                own = [j for j, k in enumerate(ev.keys) if k.kind == "width" and k.level == lvl and x[j] != 0]
                if not own:
                    results.append(ScanResult(jpi, "", np.inf, np.inf, None, "crashed", "no free amplitude"))
                    continue
                n = len(own)
                # only the first 16 patterns are screened (below): generate those, not all 2^(n-1).  Building the full list
                # first ran out of memory on a model (n ~ 40 seeded amplitudes for one level, where nearly every
                # channel is open): three drive legs SIGKILLed in this screen, 2026-11-13.  Same patterns, same order.
                patterns = [[1] * n] if n == 1 else [[1] + [1 if (m >> i) & 1 else -1 for i in range(n - 1)]
                                                    for m in range(min(2 ** (n - 1), 16))]
                cols = np.where(free_mask(ev, levels=[lvl], with_norms=True))[0]
                screened = []
                for pat in patterns[:16]:
                    xs = x.copy()
                    for j, s_ in zip(own, pat):
                        xs[j] = abs(xs[j]) * s_
                    lock = {ev.keys[j]: s_ for j, s_ in zip(own, pat)}
                    b = default_bounds(ev, xs, policy, sign_lock=lock)
                    try:
                        x1, o1, o_raw = gn_step(ev, b.clip(xs), cols, bounds=b)
                    except EvalError:
                        continue
                    screened.append((min(o1, o_raw), pat, x1 if o1 <= o_raw else b.clip(xs)))
                screened.sort(key=lambda t: t[0])
                log(f"  {jpi}: {len(screened)} sign patterns screened; best "
                    + ", ".join("".join("+" if s_ > 0 else "-" for s_ in p) + f":{o:.0f}" for o, p, _ in screened[:4]))
                best = None
                for o, pat, xstart in screened[:polish_top]:
                    lock = {ev.keys[j]: s_ for j, s_ in zip(own, pat)}
                    b = default_bounds(ev, xstart, policy, sign_lock=lock)
                    fr = polish(ev, xstart, free=free_mask(ev, levels=[lvl], with_norms=True), bounds=b,
                                policy=policy, max_nfev=max_nfev)
                    if fr.status == "crashed":
                        continue
                    if best is None or fr.objective < best.objective:
                        best = ScanResult(jpi, "".join("+" if s_ > 0 else "-" for s_ in pat), fr.objective,
                                          fr.objective - obj0, fr.x, fr.status, keys=list(ev.keys))
                if best is None:
                    best = ScanResult(jpi, "", np.inf, np.inf, None, "crashed")
                log(f"  {jpi}: best pattern {best.signs} objective {best.objective:.3f} (delta {best.delta:+.1f}, {best.status})")
                results.append(best)
        except Exception as err:
            results.append(ScanResult(jpi, "", np.inf, np.inf, None, "crashed", repr(err)[:200]))
            log(f"  {jpi}: failed ({repr(err)[:120]})")
    results.sort(key=lambda r: r.objective)
    return results


def frozen_removal_screen(ev, x, log=print) -> List[dict]:
    """Delta chi2 of switching each level off at fixed parameters (upper bound
    on its importance): zero its free widths and re-evaluate."""
    x = np.asarray(x, float)
    o0 = ev.objective(x)
    levels = {}
    for j, k in enumerate(ev.keys):
        if k.kind == "width":
            levels.setdefault(k.level, []).append(j)
    out = []
    for lvl, cols in levels.items():
        if not any(x[cols] != 0):
            continue
        xz = x.copy()
        xz[cols] = 0.0
        try:
            o = ev.objective(xz)
        except EvalError:
            o = np.inf
        out.append(dict(level=lvl, delta=o - o0, n_free=len(cols)))
    out.sort(key=lambda d: d["delta"])
    for d in out[:10]:
        log(f"  {d['level']}: frozen delta chi2 {d['delta']:+.1f} ({d['n_free']} free widths)")
    return out
