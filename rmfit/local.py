"""Bounded local optimization in reduced-width-amplitude space.

``scipy.optimize.least_squares(method="trf")`` on the residual vector an
:class:`rmfit.engine.Evaluator` returns (data rows + penalty rows + optional
prior rows), with the analytic Jacobian, box bounds, and the guards the
azure2-eval skill prescribes: dead (all-zero) Jacobian columns are dropped
before the solver sees them, an engine failure inside a trial step returns a
large residual instead of aborting the fit, and every result is run through
:func:`sanity` before anyone believes it.

Why bounds rather than a "revert guard": the archive's Levenberg-Marquardt
runaways (E -> 1e21 MeV, NaN widths) came from unbounded steps along
near-singular Jacobian directions; a trust-region step inside a box cannot
leave the box, and a converged fit sitting *on* a bound is reported as such
so the modeler can decide whether the bound is physics or a nuisance.
"""
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import inspect
import numpy as np
from scipy.optimize import least_squares

from .engine import Evaluator, EvalError
from .keys import ParamKey, level_group

BIG_RESIDUAL = 1.0e3


@dataclass
class BoundsPolicy:
    theta2_cap: float = 1.0          # |gamma|^2 <= theta2_cap * gamma_W^2 for physical levels
    bg_theta2_cap: float = 3.0       # ... for background poles
    theta2_cap_windows: list = field(default_factory=list)   # [[e_min, e_max, cap], ...] (MeV): levels whose seed energy
                                     # lies in a window get that cap instead (the evaluator 2026-09-25: relax the limit only near
                                     # the n1 threshold); first matching window wins
    bg_energy: float = 11.5          # levels at or above this Ex (MeV) count as background
    e_window: float = 0.3            # MeV either side of the current energy (physical levels)
    bg_energy_free: bool = False     # free the energies of background poles?
    norm_range: tuple = (0.2, 5.0)   # x nominal
    shift_nsigma: float = 5.0
    shifts_free: bool = False        # energy shifts are finite-differenced (2 evals each): fix them in exploration
    sign_eps: float = 1e-6           # sign lock: |gamma| >= eps * gamma_W
    min_level_width_ev: float = 0.0  # level adds: a new level that passes every test but ends NARROWER than this (eV)
    max_level_width_ev: float = 0.0  # level adds: a new level that passes every test but ends WIDER than this (eV) is "broad", not adopted (0 = off)
                                     # is classified "narrow" and not adopted (broad-structure-first campaigns); 0 = off
    threshold_margin_mev: float = 0.003   # a free level energy never comes closer than this (MeV) to a particle-pair
                                          # threshold; 0 disables.  2026-09-23: a 3/2- polished onto the n2 threshold
                                          # (0.4 keV above it) and the Brune transform of the export failed there
    freeze_broad_during_test_ev: float = 0.0   # level adds: every OTHER free level whose summed |Gamma| exceeds this (eV) is
                                          # held fixed while the new level is fitted, then everything is released in a
                                          # final polish before verification (the evaluator 2026-09-28; 0 = off)
    width_bound_windows: list = field(default_factory=list)   # [[e_min, e_max, gamma_max_ev(, pair)], ...] (MeV, eV): every free
                                          # partial width of a level whose seed energy lies in a window is bounded to
                                          # |Gamma_c| <= gamma_max_ev, i.e. theta_c^2 <= gamma_max_ev / Gamma_W,c on top of
                                          # the theta^2 cap (the evaluator 2026-09-28: keep a level added for a narrow feature
                                          # from running off while it is fitted).  First matching window wins.
    width_floor_windows: list = field(default_factory=list)   # [[e_min, e_max, gamma_min_ev, pair], ...] (MeV, eV, pair
                                          # key): a level whose seed energy lies in a window keeps its partial width in
                                          # channel `pair` at Gamma_c >= gamma_min_ev (theta_c^2 >= gamma_min_ev / Gamma_W,c,
                                          # at most 0.95 of the cap), on the side of zero its amplitude is on at x_ref
                                          # (a one-sided sign lock; flipping a whole level is a symmetry, so one side is
                                          # enough while its other channels are free).  the evaluator 2026-10-05: test whether a
                                          # BROAD level at Ex 9.4 is wanted, without letting it collapse.  First match wins.
    total_width_floor_windows: list = field(default_factory=list)   # [[e_min, e_max, gamma_min_ev], ...] (MeV, eV): a level
                                          # whose seed energy lies in a window keeps its TOTAL width sum_c Gamma_W,c theta_c^2
                                          # (open particle channels, Gamma_W at x_ref) >= gamma_min_ev.  Any channel may carry
                                          # it; single amplitudes are free, including through zero.  A one-sided penalty row
                                          # ((gamma_min - Gamma) / (total_floor_sigma * gamma_min) when Gamma < gamma_min) in
                                          # polish, not a box bound.  the evaluator 2026-10-12: "the floor should apply to the total
                                          # width, not just the n0 width ... as long as at least one of them is large enough".
    total_floor_sigma: float = 1.0e-4     # relative width of the floor penalty: 0.01 % below the floor costs 1.  Stiff on
                                          # purpose: a data pull of P chi2 per unit relative width sags the floor by ~P sigma^2/2
                                          # (a 1e5 pull -> 0.05 %); at 5e-3 the same pull sagged it by a third (unit test)
    photon_width_max_ev: float = 1.0e4    # |Gamma_gamma| <= this (eV) for every free photon channel of a physical level, and
    bg_photon_width_max_ev: float = 1.0e5  # ... of a background pole (0 = unbounded).  Photon amplitudes have no Wigner limit, so
                                          # the theta^2 cap never touched them: a from-scratch demo (2026-11-10) polished a
                                          # 5/2- background pole to Gamma_gamma 1e12 eV.  Generous stops against a runaway, not
                                          # Weisskopf limits (Gamma_W(E1) at 20 MeV in A=13 is ~3 keV, E2 ~5 eV); the bound is
                                          # |x| <= sqrt(Gamma_max / c) with Gamma_gamma = c x^2 at x_ref (the transform's local
                                          # factor), rebuilt at every polish.  A photon width at the bound shows in the review.
    new_level_theta2_cap: float = 0.0     # level adds: the NEW level's theta^2 cap (0 = the ordinary cap).  Applied as a
                                          # +-0.5 keV theta2_cap_windows entry at the candidate's seed energy, kept in the
                                          # campaign policy when the level is accepted (the evaluator 2026-09-27: "increase the cap
                                          # for the narrow resonances to 3")
    narrow_window: Optional[dict] = None   # {gamma_kev, factor, min_kev}: a level whose total width is below
                                           # gamma_kev keeps its energy within +-max(min_kev, factor*Gamma) instead of e_window

    @staticmethod
    def seed_window(energy: float, cap: float, half_kev: float = 0.5) -> list:
        """A theta2_cap_windows entry covering one level's seed energy only."""
        return [round(float(energy) - half_kev * 1e-3, 7), round(float(energy) + half_kev * 1e-3, 7), float(cap)]

    def width_bound_for(self, e_level: float, pair=None) -> Optional[float]:
        """gamma_max_ev of the first width_bound_windows entry containing e_level (MeV), else None.  An entry may carry a
        4th element, a pair key (2026-10-05, the evaluator: bound only the n2 widths of three levels): it then applies to that
        channel only and is skipped for the others.  ``pair=None`` (a level-wide question) matches any entry."""
        for w in self.width_bound_windows or []:
            if float(w[0]) <= e_level <= float(w[1]):
                if len(w) >= 4 and w[3] is not None and pair is not None and int(w[3]) != int(pair):
                    continue
                return float(w[2])
        return None

    def total_width_floor_for(self, e_level: float) -> Optional[float]:
        """gamma_min_ev of the first total_width_floor_windows entry containing e_level (MeV), else None."""
        for w in self.total_width_floor_windows or []:
            if float(w[0]) <= e_level <= float(w[1]):
                return float(w[2])
        return None

    def width_floor_for(self, e_level: float, pair) -> Optional[float]:
        """gamma_min_ev of the first width_floor_windows entry containing e_level (MeV) for this pair, else None."""
        for w in self.width_floor_windows or []:
            if float(w[0]) <= e_level <= float(w[1]) and (len(w) < 4 or w[3] is None or int(w[3]) == int(pair)):
                return float(w[2])
        return None

    def cap_for(self, e_level: float) -> float:
        """theta^2 cap for a level at excitation energy e_level (MeV): an energy window's cap if one matches, else the
        background-pole cap above bg_energy, else the physical-level cap."""
        for w in self.theta2_cap_windows or []:
            lo, hi, cap = float(w[0]), float(w[1]), float(w[2])
            if lo <= e_level <= hi:
                return cap
        return self.bg_theta2_cap if e_level >= self.bg_energy else self.theta2_cap


@dataclass
class Bounds:
    lb: np.ndarray
    ub: np.ndarray
    fixed: np.ndarray                # bool mask of parameters the policy fixes
    floors: list = field(default_factory=list)   # TotalWidthFloor rows (policy.total_width_floor_windows)

    def clip(self, x, margin=1e-9):
        x = np.array(x, float, copy=True)
        span = np.where(np.isfinite(self.ub - self.lb), self.ub - self.lb, 1.0)
        lo = np.where(np.isfinite(self.lb), self.lb + margin * span, -np.inf)
        hi = np.where(np.isfinite(self.ub), self.ub - margin * span, np.inf)
        return np.clip(x, lo, hi)


@dataclass
class TotalWidthFloor:
    """One level's total-width floor: Gamma(x) = sum_c gW_width[c] * (x[c] / gW_amp[c])^2 >= gmin (eV), as the
    one-sided residual row r = max(0, gmin - Gamma) / (sigma * gmin).  Gamma_W fixed at the reference vector."""
    level: str
    cols: np.ndarray
    gW_amp: np.ndarray
    gW_width: np.ndarray
    gmin: float
    sigma: float

    def gamma(self, x) -> float:
        x = np.asarray(x, float)
        return float(np.sum(self.gW_width * (x[self.cols] / self.gW_amp) ** 2))

    def residual(self, x) -> float:
        return max(0.0, self.gmin - self.gamma(x)) / (self.sigma * self.gmin)

    def jac_row(self, x, n) -> np.ndarray:
        row = np.zeros(n)
        if self.gamma(x) < self.gmin:
            x = np.asarray(x, float)
            row[self.cols] = -2.0 * self.gW_width * x[self.cols] / self.gW_amp ** 2 / (self.sigma * self.gmin)
        return row

    def binding(self, x, tol=0.01) -> bool:
        """The floor holds the level: its total width is within ``tol`` (relative) of the floor or below it."""
        return self.gamma(x) <= self.gmin * (1.0 + tol)


def floor_rows(floors, x, n, jac=False):
    """Residuals (and Jacobian rows) of the total-width floors at x."""
    if not floors:
        return np.zeros(0), (np.zeros((0, n)) if jac else None)
    r = np.array([f.residual(x) for f in floors])
    return r, (np.vstack([f.jac_row(x, n) for f in floors]) if jac else None)


@dataclass
class SanityReport:
    finite: bool = True
    transform_finite: bool = True
    max_theta2: float = 0.0
    theta2_violations: List[str] = field(default_factory=list)
    energy_violations: List[str] = field(default_factory=list)
    runaway: List[str] = field(default_factory=list)
    bounds_hit: List[str] = field(default_factory=list)
    # 2026-09-19: AZURE2 builds before dev e239a12 (2026-10-24) merge same-J^pi levels closer than 1 keV when they READ a
    # model, so such a fit cannot be written and read back on them.  merge_pairs fail the check; close_pairs (< 1.5 keV) are
    # a warning.  2026-11-08 (the evaluator): on a newer engine levels with different IDs stay two levels at any gap, so a close
    # pair is only a warning there; merge_pairs is then just a pair the reader would still take as one level.
    merge_pairs: List[str] = field(default_factory=list)
    close_pairs: List[str] = field(default_factory=list)
    # 2026-10-10 (the evaluator): a pair < 1 keV apart whose broader member can be moved to SEPARATE_GAP_KEV for less
    # than SEPARATE_MAX_REL_SHIFT of its width is not a failure: bake separates it before writing (see
    # separate_close_levels).  Narrow pairs stay in merge_pairs -- a lineshape question for the evaluator.
    separable_pairs: List[str] = field(default_factory=list)

    @property
    def ok(self):
        return (self.finite and self.transform_finite and not self.theta2_violations
                and not self.energy_violations and not self.runaway and not self.merge_pairs)

    def summary(self):
        if self.ok:
            return ("ok" + (f" ({len(self.bounds_hit)} at bound)" if self.bounds_hit else "")
                    + (f"; WARNING same-J^pi levels within 1.5 keV: {'; '.join(self.close_pairs)}" if self.close_pairs else "")
                    + (f"; same-J^pi pair(s) < 1 keV, broad enough to be separated at export: {'; '.join(self.separable_pairs)}"
                       if self.separable_pairs else ""))
        bits = []
        if not self.finite:
            bits.append("non-finite")
        if not self.transform_finite:
            bits.append("transform non-finite")
        if self.theta2_violations:
            bits.append(f"{len(self.theta2_violations)} amplitude(s) above their theta2 cap")
        if self.energy_violations:
            bits.append(f"E out of window x{len(self.energy_violations)}")
        if self.runaway:
            bits.append(f"runaway x{len(self.runaway)}")
        if self.merge_pairs:
            bits.append(f"same-J^pi levels AZURE2 would read as one level: {'; '.join(self.merge_pairs)}")
        return ", ".join(bits)


@dataclass
class FitResult:
    x: np.ndarray
    objective: float
    chi2: float
    penalty: float
    prior: float
    status: str                      # converged | max_nfev | crashed | rejected
    message: str
    nfev: int
    njev: int
    free: np.ndarray                 # bool mask of the parameters that moved
    dead: List[int]
    bounds_hit: List[int]
    wall: float
    start_objective: float
    sanity: Optional[SanityReport] = None
    passes: Optional[list] = None     # polish_refreeze: one dict per pass (2026-10-24)

    @property
    def gain(self):
        return self.start_objective - self.objective


# -- bounds ---------------------------------------------------------------------

_THRESHOLD_CACHE: Dict[str, tuple] = {}


_PHOTON_CACHE: Dict[str, frozenset] = {}


def photon_pairs(azr_path: str) -> frozenset:
    """File pair keys of the model's photon (capture) pairs, from the .azr's channel lines.  Cached per path."""
    if azr_path in _PHOTON_CACHE:
        return _PHOTON_CACHE[azr_path]
    try:
        from pyazr import AzrModel
        m = AzrModel.from_file(azr_path)
        out = frozenset(int(c.pair) for lv in m.levels for c in lv.channels if c.is_photon)
    except Exception:
        out = frozenset()
    _PHOTON_CACHE[azr_path] = out
    return out


def particle_thresholds(azr_path: str) -> tuple:
    """Excitation energies (MeV) at which the model's particle pairs open, from the .azr's channel
    lines (separation energy + residual excitation; photon pairs excluded).  Cached per path."""
    if azr_path in _THRESHOLD_CACHE:
        return _THRESHOLD_CACHE[azr_path]
    try:
        from pyazr import AzrModel
        m = AzrModel.from_file(azr_path)
        thr = tuple(sorted({round(float(c.sep_energy) + float(c.excitation), 6)
                            for lv in m.levels for c in lv.channels if not c.is_photon}))
    except Exception:
        thr = ()
    _THRESHOLD_CACHE[azr_path] = thr
    return thr


def _park_value(ev, r) -> float:
    """The penalty rows the engine keeps in its scalar chi2: Park's J > 0 wall (engine.ParkWall) and, since
    2026-11-07, the measured scattering lengths (engine.ScatLenRows); 0 when the model has neither."""
    f = getattr(ev, "extra_penalty_value", None) or getattr(ev, "park_value", None)
    return float(f(r)) if f is not None else 0.0


def keep_off_thresholds(lb: float, ub: float, x: float, thresholds, margin: float):
    """Cut an energy window [lb, ub] so it stays on x's side of every threshold, at least ``margin``
    away.  A level sitting inside the margin is pushed to the margin on its own side (the polish
    clips x into the box).  Pure function; returns (lb, ub)."""
    if margin <= 0:
        return lb, ub
    for T in thresholds:
        if lb < T + margin and ub > T - margin:          # the window touches the forbidden band
            if x >= T:
                lb = max(lb, T + margin)
            else:
                ub = min(ub, T - margin)
    return lb, ub


def default_bounds(ev: Evaluator, x_ref, policy: BoundsPolicy = BoundsPolicy(),
                   sign_lock: Optional[Dict[ParamKey, int]] = None,
                   cap_scale: float = 1.0) -> Bounds:
    """Box bounds from the policy, around ``x_ref``.

    Widths: ``|gamma| <= sqrt(cap * gamma_W^2)`` with the level's cap by its
    seed energy; ``sign_lock`` maps a width key to +1/-1 and keeps that
    amplitude on that side of zero (``|gamma| >= sign_eps * gamma_W``).
    Energies: a window around the current value for physical levels,
    fixed for background poles unless ``bg_energy_free``.  Norms: a range
    times nominal.  Shifts: ``shift_nsigma`` sigma around nominal, or fixed.

    ``cap_scale`` > 1 relaxes every theta^2 cap by that factor -- but never
    below the amplitude's current theta^2 if that already exceeds the policy
    cap -- so a fit that starts outside the physical box (one archived
    fit had dozens of amplitudes with theta^2 > 1) can be walked into it
    in stages instead of being clipped in one step (see
    :meth:`rmfit.campaign.Campaign.run_polish` with ``tighten``).
    """
    x_ref = np.asarray(x_ref, float)
    n = x_ref.size
    lb = np.full(n, -np.inf)
    ub = np.full(n, np.inf)
    fixed = np.zeros(n, bool)
    gW = ev.wigner_amplitudes()
    _azr = getattr(ev, "azr", None) or getattr(ev, "base_azr", None)
    _thresholds = particle_thresholds(_azr) if (_azr and policy.threshold_margin_mev > 0) else ()
    sign_lock = sign_lock or {}
    narrow_width: Dict[str, float] = {}
    if policy.narrow_window:
        # physical partial widths (eV) of the free R-matrix block at x_ref: the
        # first n_rmatrix entries of the engine's transform align with ev.keys
        try:
            phys = np.asarray(ev.transform(x_ref), float)
            for jj, kk in enumerate(ev.keys[:ev.n_rmatrix]):
                if kk.kind == "width" and jj < phys.size:
                    narrow_width[kk.level] = narrow_width.get(kk.level, 0.0) + abs(float(phys[jj])) * 1e-3
        except Exception:
            narrow_width = {}
    gWw: Dict[int, float] = {}
    if policy.width_bound_windows or policy.width_floor_windows or policy.total_width_floor_windows:
        try:
            gWw = dict(ev.wigner_widths(x_ref))
        except Exception:
            gWw = {}
    # photon channels: no Wigner limit, so an absolute width cap (policy.photon_width_max_ev / bg_photon_width_max_ev)
    # through the transform's asymptotic factor Gamma_gamma = c x^2.  With external capture the physical width is
    # c (x - x0)^2 with x0 of order 1 (the external amplitude; 13N 1/2+ E1: x0 = 1.0, so Gamma(x_ref = 0.22) / x_ref^2
    # is 12x the asymptotic c), hence c is probed at a large amplitude, one transform per photon channel, where
    # x0 does not matter; the bound is a stop against a runaway, not a Weisskopf limit
    ph_pairs = getattr(ev, "photon_pairs", None)
    if ph_pairs is None:
        ph_pairs = photon_pairs(_azr) if _azr else frozenset()
    ph_c: Dict[int, float] = {}
    if ph_pairs and (policy.photon_width_max_ev > 0 or policy.bg_photon_width_max_ev > 0):
        try:
            X = 1.0e4
            for j, k in enumerate(ev.keys[:ev.n_rmatrix]):
                if k.kind == "width" and k.pair in ph_pairs:
                    probe = x_ref.copy(); probe[j] = X
                    ph_c[j] = abs(float(np.asarray(ev.transform(probe), float)[j])) / X ** 2
        except Exception:
            ph_c = {}
    for j, k in enumerate(ev.keys):
        if k.kind == "width":
            e_seed = ev.level_energy.get(k.level, 0.0)
            cap = policy.cap_for(e_seed)
            w = gW.get(j)
            if w is None or not np.isfinite(w) or w <= 0:
                c_ph = ph_c.get(j)
                gmax_ph = policy.bg_photon_width_max_ev if e_seed >= policy.bg_energy else policy.photon_width_max_ev
                if c_ph is not None and np.isfinite(c_ph) and c_ph > 0 and gmax_ph > 0:
                    hi = float(np.sqrt(gmax_ph / c_ph))
                    if cap_scale > 1.0:
                        hi = max(hi, min(abs(float(x_ref[j])) * (1 + 1e-6), hi * np.sqrt(cap_scale)))
                    lb[j], ub[j] = -hi, hi
                    s = sign_lock.get(k)
                    if s is not None:
                        if s > 0:
                            lb[j] = policy.sign_eps * hi
                        else:
                            ub[j] = -policy.sign_eps * hi
                continue
            if cap_scale > 1.0:
                cur = (x_ref[j] / w) ** 2
                cap = max(cap, min(cur * (1 + 1e-6), cap * cap_scale))
            gmax = policy.width_bound_for(e_seed, k.pair)
            if gmax is not None and gWw.get(j, 0.0) > 0:
                cap = min(cap, gmax / float(gWw[j]))          # theta_c^2 <= Gamma_max / Gamma_W,c
            hi = np.sqrt(cap) * w
            lb[j], ub[j] = -hi, hi
            s = sign_lock.get(k)
            eps = policy.sign_eps * w
            gmin = policy.width_floor_for(e_seed, k.pair) if k.pair is not None else None
            # an amplitude that is exactly 0 is a level SWITCHED OFF (run_level_add's control zeroes the new level's
            # widths): no floor, or the clip would switch it back on (2026-10-05, 10-5-26_ex94_floor 1/2+ / 3/2- controls)
            if gmin is not None and gWw.get(j, 0.0) > 0 and x_ref[j] != 0.0:
                # one-sided floor |gamma| >= sqrt(gmin / Gamma_W,c) * gamma_W, never above 0.95 of the cap
                eps = max(eps, min(np.sqrt(gmin / float(gWw[j])) * w, 0.95 * hi))
                if s is None:
                    s = -1 if x_ref[j] < 0 else 1
            if s is not None:
                if s > 0:
                    lb[j] = eps
                else:
                    ub[j] = -eps
        elif k.kind == "energy":
            e_seed = ev.level_energy.get(k.level, x_ref[j])
            if e_seed >= policy.bg_energy and not policy.bg_energy_free:
                fixed[j] = True
            else:
                half = policy.e_window
                nw = policy.narrow_window
                if nw and k.level in narrow_width:
                    g = narrow_width[k.level]                       # keV, sum of the level's free partial widths
                    if g < float(nw.get("gamma_kev", 5.0)):
                        half = max(float(nw.get("min_kev", 1.0)), float(nw.get("factor", 3.0)) * g) * 1e-3
                lb[j] = x_ref[j] - half
                ub[j] = x_ref[j] + half
                lb[j], ub[j] = keep_off_thresholds(lb[j], ub[j], x_ref[j], _thresholds, policy.threshold_margin_mev)
        elif k.kind == "norm":
            d = ev.segments.by_key(k.seg)
            nom = d.norm if d is not None else 1.0
            lb[j], ub[j] = policy.norm_range[0] * nom, policy.norm_range[1] * nom
        elif k.kind == "shift":
            d = ev.segments.by_key(k.seg)
            if not policy.shifts_free or d is None or not d.energy_shift_error:
                fixed[j] = True
            else:
                lb[j] = d.energy_shift - policy.shift_nsigma * d.energy_shift_error
                ub[j] = d.energy_shift + policy.shift_nsigma * d.energy_shift_error
        elif k.kind in ("shift_sqrt", "shift_fn"):
            # 2026-10-22: the sqrt(E) coefficient and shift-function parameters follow the shift policy:
            # fixed unless shifts_free and the term has an uncertainty, then nominal +- shift_nsigma sigma
            d = ev.segments.by_key(k.seg)
            if d is None:
                nom, s = 0.0, 0.0
            elif k.kind == "shift_sqrt":
                nom, s = d.energy_shift_sqrt, d.energy_shift_sqrt_error
            else:
                nom, s = float(d.shift_fn_params[k.term][0]), float(d.shift_fn_params[k.term][1])
            if not policy.shifts_free or not s:
                fixed[j] = True
            else:
                lb[j] = nom - policy.shift_nsigma * s
                ub[j] = nom + policy.shift_nsigma * s
        elif k.kind == "weight":
            # a mixture member's abundance / stoichiometric weight (rmfit.mixture)
            lb[j], ub[j] = 0.0, max(1.0, 2.0 * x_ref[j])
        elif k.kind == "density":
            # a fitted target density: between a hundredth of the line's and the largest its sub-point grids span
            # (beyond that AZURE2 computes at the end of the range and warns)
            d = getattr(ev, "density_info", {}).get(k.effect)
            if d is not None:
                lb[j], ub[j] = 0.01 * d["density"], d["max"]
        elif k.kind == "doppler":
            # a fitted Doppler temperature: the engine's own range [T0/4, the layout temperature its grids span]
            d = getattr(ev, "doppler_info", {}).get(k.effect)
            if d is not None:
                lb[j], ub[j] = 0.25 * d["temperature"], d["max"]
        elif k.kind == "screening":
            # an electron-screening potential (eV): the engine's own limits, U_e >= 0, no ceiling
            lb[j], ub[j] = 0.0, np.inf
        elif k.kind == "thickness":
            # a fitted sample thickness (atoms/barn): positive, no ceiling (the engine's limit is 0)
            d = getattr(ev, "thickness_info", {}).get(k.seg)
            if d is not None:
                lb[j], ub[j] = 1.0e-6 * d["n"], np.inf
    floors = []
    if policy.total_width_floor_windows and gWw:
        by_level: Dict[str, list] = {}
        for j, k in enumerate(ev.keys):
            if k.kind == "width" and k.pair is not None:
                gmin = policy.total_width_floor_for(ev.level_energy.get(k.level, 0.0))
                if gmin is not None:
                    by_level.setdefault(k.level, []).append(j)
        for lvl, js in by_level.items():
            # open particle channels only (Gamma_W > 0 at x_ref); a level whose amplitudes are ALL exactly zero is switched
            # off (run_level_add's control): no floor, as for width_floor_windows
            js = [j for j in js if gWw.get(j, 0.0) > 0 and gW.get(j, 0.0) and np.isfinite(gW.get(j, np.nan))]
            if not js or all(x_ref[j] == 0.0 for j in js):
                continue
            gmin = policy.total_width_floor_for(ev.level_energy.get(lvl, 0.0))
            floors.append(TotalWidthFloor(lvl, np.array(js, int), np.array([gW[j] for j in js], float),
                                          np.array([gWw[j] for j in js], float), float(gmin), float(policy.total_floor_sigma)))
    return Bounds(lb, ub, fixed, floors)


# -- sanity -----------------------------------------------------------------------

# -- near-degenerate same-J^pi levels ------------------------------------------------
# 2026-10-10 (the evaluator: "If the levels are broad, then we could probably just get away with doing a check where we give the
# energy of the broader level a little adjustment in energy when this issue is going to cause a problem").  AZURE2 merges
# same-J^pi levels < 1 keV apart when it READS a model, so a fit can hold such a pair in memory but not write it.  For a
# broad level a few keV is nothing: measured on a large merged model, moving a ~6 MeV-wide level by 1 / 2 / 5 keV
# cost +1.2 / +3.9 / +21 in chi2 (frozen, on an objective of several hundred thousand).  So before a write the broader member of a pair closer than
# SEPARATE_MIN_GAP_KEV is moved away to SEPARATE_GAP_KEV, if that shift is < SEPARATE_MAX_REL_SHIFT of its width;
# narrow pairs (the 9/2+ doublet of 2026-09-19) are not touched and still fail.
SEPARATE_MIN_GAP_KEV = 1.5
SEPARATE_GAP_KEV = 2.5
SEPARATE_MAX_REL_SHIFT = 1e-2
# 2026-11-08: with an engine that keeps close levels apart (enginepin.keeps_close_levels) nothing above is applied; a pair
# closer than this (MeV-scale round-off of the printed energy, far below any fitted gap) still counts as one level.
SAME_ENERGY_KEV = 1e-6


def level_total_widths(ev, x) -> dict:
    """Gamma (eV) per level label: the sum of |physical width| over the level's FREE width keys at ``x``
    (fixed widths are not in the vector; for the levels this is used on, the broad free ones, that is the bulk)."""
    phys = np.asarray(ev.transform(np.asarray(x, float)), float)
    out = {}
    for j, k in enumerate(ev.keys):
        if k.kind == "width" and k.level and j < phys.size and np.isfinite(phys[j]):
            out[k.level] = out.get(k.level, 0.0) + abs(float(phys[j]))
    return out


def _separation_plan(l1, l2, gap_kev, gam, free_energy_levels, target_kev=SEPARATE_GAP_KEV,
                     max_rel=SEPARATE_MAX_REL_SHIFT):
    """For a same-J^pi pair (l1 below l2, ``gap_kev`` apart): ``(label, signed shift in keV)`` moving the broader member
    with a FREE energy away from the other to ``target_kev``, or None if no member qualifies (shift >= max_rel * Gamma)."""
    shift = target_kev - gap_kev
    for g, lab, sign in sorted([(gam.get(l1, 0.0), l1, -1.0), (gam.get(l2, 0.0), l2, +1.0)], reverse=True):
        if lab in free_energy_levels and g > 0 and shift * 1e3 / g < max_rel:
            return lab, sign * shift
    return None


def separate_close_levels(ev, x, min_gap_kev=SEPARATE_MIN_GAP_KEV, target_kev=SEPARATE_GAP_KEV,
                          max_rel=SEPARATE_MAX_REL_SHIFT):
    """A copy of ``x`` in which every same-J^pi pair closer than ``min_gap_kev`` has its broader member (free energy,
    shift < ``max_rel`` of its width) moved away to ``target_kev``.  Returns ``(x, moved, refused)``: ``moved`` = list of
    dicts (level, other, gap_kev, shift_kev, gamma_kev, e_from, e_to), ``refused`` = messages for pairs that could not
    be separated (narrow, or no free energy).  The caller re-scores and decides whether the cost is acceptable."""
    x = np.array(x, float, copy=True)
    ekey = {k.level: j for j, k in enumerate(ev.keys) if k.kind == "energy" and k.level}
    en = {lab: float(e) for lab, e in getattr(ev, "level_energy", {}).items()}
    for lab, j in ekey.items():
        en[lab] = float(x[j])
    gam = None
    moved, refused = [], []
    groups = {}
    for lab, e in en.items():
        groups.setdefault(level_group(lab), []).append(lab)     # J^pi and, with several nuclei, the nucleus
    for jpi, labs in groups.items():
        for _ in range(len(labs)):                       # repeat until no close pair is left (chains of three)
            lv = sorted((en[l], l) for l in labs)
            bad = [(e1, l1, e2, l2) for (e1, l1), (e2, l2) in zip(lv, lv[1:]) if 1000.0 * (e2 - e1) < min_gap_kev]
            bad = [b for b in bad if f"{b[1]}/{b[3]}" not in {r.split(':')[0] for r in refused}]
            if not bad:
                break
            e1, l1, e2, l2 = bad[0]; gap = 1000.0 * (e2 - e1)
            if gam is None:
                gam = level_total_widths(ev, x)
            plan = _separation_plan(l1, l2, gap, gam, set(ekey), target_kev, max_rel)
            if plan is None:
                refused.append(f"{l1}/{l2}: {gap:.3f} keV apart, widths {gam.get(l1, 0) / 1e3:.3g} / {gam.get(l2, 0) / 1e3:.3g} keV"
                               f" -- no member broad enough (shift < {max_rel:g} of Gamma) with a free energy")
                continue
            lab, dk = plan; j = ekey[lab]; e0 = float(x[j])
            x[j] = e0 + dk * 1e-3; en[lab] = float(x[j])
            moved.append(dict(level=lab, other=(l2 if lab == l1 else l1), gap_kev=round(gap, 4), shift_kev=round(dk, 4),
                              gamma_kev=round(gam.get(lab, 0.0) / 1e3, 3), e_from=e0, e_to=float(x[j])))
    return x, moved, refused


def sanity(ev: Evaluator, x, policy: BoundsPolicy = BoundsPolicy(),
           bounds: Optional[Bounds] = None, transform=True, x_ref=None) -> SanityReport:
    """``x_ref`` is the vector the energy window is measured from (the
    incumbent); without it the file's seed energies are used, which is wrong
    once a campaign has legitimately moved a level by more than one window."""
    x = np.asarray(x, float)
    ref = None if x_ref is None else np.asarray(x_ref, float)
    rep = SanityReport(finite=bool(np.all(np.isfinite(x))))
    gW = ev.wigner_amplitudes()
    for j, k in enumerate(ev.keys):
        if k.kind == "width":
            w = gW.get(j)
            if w:
                t2 = float((x[j] / w) ** 2)
                rep.max_theta2 = max(rep.max_theta2, t2)
                e_seed = ev.level_energy.get(k.level, 0.0)
                cap = policy.cap_for(e_seed)
                if t2 > cap * (1 + 1e-6):
                    rep.theta2_violations.append(f"{k}: theta2={t2:.3g}")
            if abs(x[j]) > 1e3:
                rep.runaway.append(f"{k}: {x[j]:.3g}")
        elif k.kind == "energy":
            if abs(x[j]) > 60.0:
                rep.runaway.append(f"{k}: E={x[j]:.3g}")
            e_seed = ev.level_energy.get(k.level, x[j])
            e_ref = float(ref[j]) if ref is not None else e_seed
            if e_seed < policy.bg_energy and abs(x[j] - e_ref) > policy.e_window * (1 + 1e-6):
                rep.energy_violations.append(f"{k}: {x[j]:.4f} vs reference {e_ref:.4f}")
        elif k.kind == "norm":
            if not (0 < x[j] < 100):
                rep.runaway.append(f"{k}: norm={x[j]:.3g}")
    if bounds is not None:
        tol = 1e-6
        for j in range(x.size):
            if bounds.fixed[j]:
                continue
            span = bounds.ub[j] - bounds.lb[j]
            if np.isfinite(span) and span > 0:
                if x[j] - bounds.lb[j] <= tol * span or bounds.ub[j] - x[j] <= tol * span:
                    rep.bounds_hit.append(str(ev.keys[j]))
    # same-J^pi spacing: energies from the evaluator's level table, free ones taken from x
    try:
        en = {lab: float(e) for lab, e in getattr(ev, "level_energy", {}).items()}
        for j, k in enumerate(ev.keys):
            if k.kind == "energy" and k.level in en:
                en[k.level] = float(x[j])
        groups = {}
        for lab, e in en.items():
            groups.setdefault(level_group(lab), []).append((e, lab))
        gam = None
        from .enginepin import keeps_close_levels
        keeps = keeps_close_levels()
        for jpi, lv in groups.items():
            lv.sort()
            for (e1, l1), (e2, l2) in zip(lv, lv[1:]):
                gap = 1000.0 * (e2 - e1)
                if gap < 1.5:
                    msg = f"{l1} {e1:.6f} / {l2} {e2:.6f} ({gap:.3f} keV)"
                    if keeps:
                        (rep.merge_pairs if gap < SAME_ENERGY_KEV else rep.close_pairs).append(msg)
                        continue
                    if gap >= 1.0:
                        rep.close_pairs.append(msg)
                        continue
                    if gam is None:
                        gam = level_total_widths(ev, x) if transform else {}
                    plan = _separation_plan(l1, l2, gap, gam, {k.level for k in ev.keys if k.kind == "energy"})
                    if plan is not None:
                        rep.separable_pairs.append(msg + f" -> {plan[0]} moves {plan[1]:+.2f} keV")
                    else:
                        rep.merge_pairs.append(msg)
    except Exception:
        pass
    if transform and rep.finite:
        try:
            phys = np.asarray(ev.transform(x), float)
            rep.transform_finite = bool(np.all(np.isfinite(phys)))
        except EvalError:
            rep.transform_finite = False
    return rep


# -- the polish -------------------------------------------------------------------


def _ev_residuals(ev, x, jac, cols):
    """``ev.residuals`` asking only for the Jacobian columns ``cols`` when the evaluator takes them (the engine's
    held finite-difference columns are then skipped, 2026-10-30); every column otherwise.  The columns outside
    ``cols`` may come back zero, so callers read ``J[:, cols]`` only."""
    if jac and cols is not None:
        try:
            takes = "cols" in inspect.signature(ev.residuals).parameters
        except (TypeError, ValueError):
            takes = False
        if takes:
            return ev.residuals(x, jac=True, cols=cols)
    return ev.residuals(x, jac=jac)

def polish(ev: Evaluator, x0, free: Optional[np.ndarray] = None,
           bounds: Optional[Bounds] = None, policy: BoundsPolicy = BoundsPolicy(),
           max_nfev=200, ftol=1e-9, xtol=1e-9, gtol=1e-10, x_scale="jac",
           verbose=0, log=None) -> FitResult:
    """Minimize the evaluator's objective over the ``free`` parameters from ``x0``.

    ``free`` is a bool mask (default: everything the bounds do not fix).  The
    Jacobian at ``x0`` decides which of those columns are dead (identically
    zero -- a parameter no active datum sees) and drops them, because a zero
    column with ``x_scale="jac"`` makes TRF report convergence after one
    evaluation.
    """
    t0 = time.time()
    x0 = np.asarray(x0, float).copy()
    n = x0.size
    if bounds is None:
        bounds = default_bounds(ev, x0, policy)
    if free is None:
        free = ~bounds.fixed
    free = np.asarray(free, bool) & ~bounds.fixed
    x0 = bounds.clip(x0)
    say = log or (print if verbose else (lambda *a: None))

    floors = list(getattr(bounds, "floors", None) or [])

    n_ev = {"rows": None}             # the evaluator's own rows; the floor rows go after them

    def _resid(x, jac, want=None):
        r, J = _ev_residuals(ev, x, jac, want)
        n_ev["rows"] = r.size
        if floors:
            rf, Jf = floor_rows(floors, x, n, jac=jac)
            r = np.concatenate([r, rf])
            if jac:
                J = np.vstack([J, Jf])
        return r, J

    try:
        r0, J0 = _resid(x0, True, np.where(free)[0])   # only the free columns are differenced
    except EvalError as err:
        return FitResult(x0, np.inf, np.inf, 0.0, 0.0, "crashed", f"start: {err}", 0, 0,
                         free, [], [], time.time() - t0, np.inf)
    if floors:
        say("  total-width floor(s): " + "; ".join(f"{f.level} >= {f.gmin / 1e3:g} keV (now {f.gamma(x0) / 1e3:.1f} keV, "
                                                  f"{f.cols.size} channels)" for f in floors))

    def _ev_rows(rr):
        # 2026-11-07: the engine-side penalty rows (Park wall, scattering lengths) sit at the END of the
        # evaluator's rows, so they must be read before the floor rows are appended
        return rr[:n_ev["rows"]] if n_ev["rows"] is not None else rr
    start_obj = ev.rows_objective(r0)
    cols = np.where(free)[0]
    alive = np.max(np.abs(J0[:, cols]), axis=0) > 0
    dead = [int(c) for c in cols[~alive]]
    cols = cols[alive]
    if dead:
        say(f"  {len(dead)} dead columns dropped: " + ", ".join(str(ev.keys[c]) for c in dead[:6]))
    if cols.size == 0:
        return FitResult(x0, start_obj, ev.rows_chi2(r0),
                         ev.penalties.value(x0) + _park_value(ev, _ev_rows(r0)), ev.prior.value(x0) if ev.prior else 0.0,
                         "converged", "nothing to fit", 1, 1, free, dead, [],
                         time.time() - t0, start_obj)

    cache = {"x": x0.copy(), "r": r0, "J": J0}
    counts = {"f": 1, "j": 1}

    def full(z):
        x = x0.copy()
        x[cols] = z
        return x

    def fun(z):
        x = full(z)
        if np.array_equal(x, cache["x"]) and cache["r"] is not None:
            return cache["r"]
        counts["f"] += 1
        try:
            r, _ = _resid(x, False)
        except EvalError as err:
            say(f"  engine failure in a trial step ({err}); returning a large residual")
            return np.full(r0.size, BIG_RESIDUAL)
        cache.update(x=x.copy(), r=r, J=None)
        return r

    def jac(z):
        x = full(z)
        if np.array_equal(x, cache["x"]) and cache["J"] is not None:
            return cache["J"][:, cols]
        counts["j"] += 1
        try:
            r, J = _resid(x, True, cols)
        except EvalError as err:
            say(f"  engine failure in a Jacobian ({err}); returning zeros")
            return np.zeros((r0.size, cols.size))
        cache.update(x=x.copy(), r=r, J=J)
        if counts["j"] % 10 == 0:
            say(f"    iteration {counts['j']}: objective {ev.rows_objective(r):.3f} ({time.time()-t0:.0f}s)")
        return J[:, cols]

    lb, ub = bounds.lb[cols], bounds.ub[cols]
    try:
        loss = ev.ls_loss
        res = least_squares(fun, x0[cols], jac=jac, bounds=(lb, ub), method="trf",
                            x_scale=x_scale, max_nfev=max_nfev, ftol=ftol, xtol=xtol,
                            gtol=gtol, verbose=2 if verbose > 1 else 0,
                            **({} if loss is None else dict(loss=loss, f_scale=1.0)))
    except Exception as err:
        return FitResult(x0, start_obj, np.nan, np.nan, np.nan, "crashed", repr(err),
                         counts["f"], counts["j"], free, dead, [], time.time() - t0, start_obj)
    x = full(res.x)
    # final consistent numbers at the solution
    try:
        r, _ = _resid(x, False)
        counts["f"] += 1
    except EvalError as err:
        return FitResult(x0, start_obj, np.nan, np.nan, np.nan, "crashed", f"end: {err}",
                         counts["f"], counts["j"], free, dead, [], time.time() - t0, start_obj)
    obj = ev.rows_objective(r)
    chi2 = ev.rows_chi2(r)
    pen = ev.penalties.value(x) + _park_value(ev, _ev_rows(r))
    pri = ev.prior.value(x) if ev.prior is not None else 0.0
    status = "converged" if res.status in (1, 2, 3, 4) else ("max_nfev" if res.status == 0 else "crashed")
    hit = []
    for c, lo, hi in zip(cols, lb, ub):
        span = hi - lo
        if np.isfinite(span) and span > 0 and (x[c] - lo <= 1e-6 * span or hi - x[c] <= 1e-6 * span):
            hit.append(int(c))
    if obj > start_obj:                      # never hand back something worse
        x, obj, chi2, pen, pri = x0, start_obj, ev.rows_chi2(r0), \
            ev.penalties.value(x0) + _park_value(ev, _ev_rows(r0)), (ev.prior.value(x0) if ev.prior else 0.0)
        status, hit = "rejected", []
    say(f"  polish: {start_obj:.3f} -> {obj:.3f} ({status}, {res.message.strip()}; "
        f"nfev={counts['f']} njev={counts['j']}, {time.time()-t0:.0f}s)")
    return FitResult(x, obj, chi2, pen, pri, status, res.message.strip(), counts["f"],
                     counts["j"], free, dead, hit, time.time() - t0, start_obj)



def polish_refreeze(ev: Evaluator, x0, bounds: Optional[Bounds] = None, policy: BoundsPolicy = BoundsPolicy(),
                    max_nfev=200, every=100, tol_rel=1e-5, tol_abs=1.0, max_passes=50, log=None,
                    **kw) -> FitResult:
    """:func:`polish` with the data covariance's model-scaled (``relvector``) terms re-frozen between passes.

    A ``relvector`` term is frozen at the model of the evaluator's first evaluation.  One long polish from a poor
    seed therefore minimises chi2 with the systematics scaled by the SEED's model, and a fresh session at the result
    disagrees (by 4 % in one case).  Here each pass is
    an ordinary least-squares run on a FIXED objective (frozen at the pass's start point) of at most ``every``
    evaluations; then the covariance is re-frozen at the pass's end point and that point is re-scored.  This is the
    usual fixed-point treatment of model-scaled systematics: the covariance never follows the model inside a
    least-squares run (which would confuse the trust region and let the fit inflate the model to grow its own
    variance), and it is brought up to date between runs.  Stops when a pass ends ``converged``/``rejected`` and the
    re-freeze moved the objective by at most ``max(tol_abs, tol_rel * objective)``, or when ``max_nfev`` (total over
    passes) or ``max_passes`` is used up.  The returned ``objective``/``chi2``/``penalty`` are the self-consistent
    ones (frozen at the returned ``x``, as a fresh session sees them) and ``start_objective`` is frozen at ``x0``.
    Without a ``relvector`` covariance, with a fixed covariance reference (``cov.has_reference``) or with
    ``every <= 0``, this IS :func:`polish` (same call, same result).

    2026-10-25: the iteration need not contract -- on a large model it oscillated (re-freeze jumps from +3k to +140k)
    and the last pass ended well above the seed.  The returned point is therefore the
    BEST self-consistent one over the passes (its pass number in ``message``), not the last.  rmfit campaigns no
    longer use this by default; they fix the relvector model at a reference (ObjectiveSpec.cov_reference)."""
    cov = getattr(ev, "cov", None)
    if cov is None or not cov.relative or getattr(cov, "has_reference", False) or every is None or every <= 0:
        return polish(ev, x0, bounds=bounds, policy=policy, max_nfev=max_nfev, log=log, **kw)
    say = log or (lambda *a: None)
    t0 = time.time()
    x = np.asarray(x0, float).copy()
    if bounds is None:
        bounds = default_bounds(ev, x, policy)
    left, nfev, njev, start, k = int(max_nfev), 0, 0, None, 0
    fr = sc = None
    passes = []
    best = None                                           # (self-consistent objective, x, score, pass)
    while True:
        k += 1
        ev.refreeze_covariance()                        # the pass's first evaluation (its start point) freezes it
        fr = polish(ev, x, bounds=bounds, policy=policy, max_nfev=max(1, min(int(every), left)), log=log, **kw)
        nfev, njev, left = nfev + fr.nfev, njev + fr.njev, left - fr.nfev
        if start is None:
            start = fr.start_objective
        if fr.status == "crashed":
            say(f"  refreeze pass {k}: polish crashed ({fr.message}); stopping")
            break
        x = fr.x
        ev.refreeze_covariance()
        sc = ev.score(x)                                  # frozen at the pass's end point = what a fresh session sees
        jump = sc.objective - fr.objective
        tol = max(tol_abs, tol_rel * abs(sc.objective))
        passes.append(dict(pass_=k, start=fr.start_objective, end=fr.objective, refrozen=sc.objective, jump=jump,
                           status=fr.status, nfev=fr.nfev))
        if best is None or sc.objective < best[0]:
            best = (sc.objective, np.array(x, float, copy=True), sc, k)
        say(f"  refreeze pass {k}: {fr.start_objective:.3f} -> {fr.objective:.3f} ({fr.status}, nfev {fr.nfev}); "
            f"re-frozen at its end {sc.objective:.3f} (jump {jump:+.3f}, tol {tol:.3g}); {left} evaluations left")
        if fr.status in ("converged", "rejected") and abs(jump) <= tol:
            status, msg = "converged", f"converged with the covariance frozen at the result ({k} passes)"
            break
        if left <= 0 or k >= max_passes:
            status, msg = "max_nfev", f"evaluation budget used ({k} passes; last re-freeze jump {jump:+.3f})"
            break
    if fr.status == "crashed":
        if sc is None:                                    # crashed in the first pass: polish's own result
            fr.nfev, fr.njev = nfev, njev
            return fr
        status, msg = "crashed", fr.message                # keep the best good pass's point
    if best is not None and best[3] != len(passes):
        say(f"  refreeze: returning pass {best[3]}'s point ({best[0]:.3f}), the best self-consistent one; the last "
            f"pass ended at {sc.objective:.3f}")
        msg = f"{msg}; best self-consistent point = pass {best[3]} of {len(passes)}"
    if best is not None:
        _, x, sc, _ = best
        ev.refreeze_covariance()
        ev.score(x)                                       # leave the evaluator frozen at the returned point
    out = FitResult(x, sc.objective, sc.chi2, sc.penalty, sc.prior, status, msg, nfev, njev, fr.free, fr.dead,
                    fr.bounds_hit, time.time() - t0, start)
    out.passes = passes
    return out

def width_bound_pins(keys, own, ub, gW, cap) -> List[str]:
    """Width keys of a level (positions ``own``) whose upper bound is BELOW the theta^2-cap bound
    sqrt(cap) * gamma_W -- i.e. set by a width_bound_windows entry, not by the cap.  An amplitude
    sitting on such a bound is the bound doing its job (the evaluator 2026-09-28: the level is held at the
    feature's width), not an absorber, so it must not count as 'unphysical'."""
    out = []
    for j in own:
        k = keys[j]
        if k.kind != "width" or j not in gW or not (gW[j] > 0):
            continue
        cap_hi = np.sqrt(cap) * gW[j]
        if np.isfinite(ub[j]) and ub[j] < cap_hi * (1 - 1e-9):
            out.append(str(k))
    return out


def free_mask(ev: Evaluator, jpi: Optional[Sequence[str]] = None,
              levels: Optional[Sequence[str]] = None, kinds=("energy", "width"),
              with_norms=True, with_shifts=False, extra: Sequence[ParamKey] = ()) -> np.ndarray:
    """A bool mask over ``ev.keys``: the R-matrix parameters of the given
    J^pi groups and/or levels (all if neither is given), plus norms/shifts."""
    jpi = set(jpi) if jpi else None
    levels = set(levels) if levels else None
    extra = set(extra)
    m = np.zeros(len(ev.keys), bool)
    for j, k in enumerate(ev.keys):
        if k in extra:
            m[j] = True
        elif k.kind in kinds and k.level is not None:
            if jpi is None and levels is None:
                m[j] = True
            elif (jpi and k.jpi in jpi) or (levels and k.level in levels):
                m[j] = True
        elif k.kind == "norm" and with_norms:
            m[j] = True
        elif k.kind in ("shift", "shift_sqrt", "shift_fn") and with_shifts:
            m[j] = True
    return m


def level_energies_now(ev: Evaluator, x) -> dict:
    """Current energy (MeV) of every level: the fitted value where the energy is a free key, else the file value."""
    out = dict(ev.level_energy)
    for j, k in enumerate(ev.keys):
        if k.kind == "energy" and k.level is not None:
            out[k.level] = float(x[j])
    return out


def local_mask(ev: Evaluator, x, centers: Sequence[float], half_width: float, bg_energy: Optional[float] = None,
               with_norms: bool = True, exclude_levels: Sequence[str] = ()) -> np.ndarray:
    """2026-10-02 "local release" (the evaluator: a flip or a new level should be judged after "a limited polish of the nearby
    other levels and background levels").  A bool mask over ``ev.keys``: the energies and widths of EVERY level, of any
    J^pi, whose current energy lies within +-``half_width`` MeV of any of ``centers``; the widths of every background
    level (current energy >= ``bg_energy``, when given); and the normalizations.  ``exclude_levels`` are left out."""
    E = level_energies_now(ev, x)
    near = {lab for lab, e in E.items() if any(abs(e - c) <= half_width for c in centers)} - set(exclude_levels)
    bg = {lab for lab, e in E.items() if bg_energy is not None and e >= bg_energy} - set(exclude_levels)
    m = np.zeros(len(ev.keys), bool)
    for j, k in enumerate(ev.keys):
        if k.level is not None and k.kind in ("energy", "width"):
            m[j] = k.level in near or (k.kind == "width" and k.level in bg)
        elif k.kind == "norm" and with_norms:
            m[j] = True
    return m


def polish_group(ev: Evaluator, x0, jpi, with_norms=True, **kw) -> FitResult:
    """Polish one J^pi group (its energies and widths) plus the normalizations."""
    return polish(ev, x0, free=free_mask(ev, jpi=[jpi], with_norms=with_norms), **kw)


def gn_step(ev: Evaluator, x, cols, damping=1e-2, bounds: Optional[Bounds] = None):
    """One damped Gauss-Newton step on ``cols`` from ``x``; returns
    ``(x1, objective_at_x1, objective_at_x)`` -- three evaluations' worth, the
    cheap screen for a proposal (flip, kick) before any real polish."""
    x = np.asarray(x, float)
    r, J = _ev_residuals(ev, x, True, cols)
    obj0 = ev.rows_objective(r)
    Jc = J[:, cols]
    alive = np.max(np.abs(Jc), axis=0) > 0
    cols = np.asarray(cols)[alive]
    Jc = Jc[:, alive]
    if cols.size == 0:
        return x.copy(), obj0, obj0
    A = Jc.T @ Jc
    A[np.diag_indices_from(A)] *= (1.0 + damping)
    g = Jc.T @ r
    try:
        step = -np.linalg.solve(A, g)
    except np.linalg.LinAlgError:
        step = -np.linalg.lstsq(A, g, rcond=None)[0]
    x1 = x.copy()
    x1[cols] += step
    if bounds is not None:
        x1 = bounds.clip(x1)
    try:
        obj1 = ev.objective(x1)
    except EvalError:
        obj1 = np.inf
    return x1, obj1, obj0
