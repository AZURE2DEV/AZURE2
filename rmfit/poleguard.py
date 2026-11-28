"""Pole guard: keep a fit from building narrow poles that sit on no level of their J^pi (pyazr 3 only).

Why (an archived campaign, 2026-11-24..27): a collective pole of one J^pi's level matrix (large sub-threshold
amplitudes and their Brune closed-channel shifts) made a narrow peak in a bare total cross section MeV away from any
level of that J^pi.  A theta^2 cap only moved it (by about a MeV, down to an eV-scale width); removing the level that
carried it let the refit rebuild it from the other levels of the group.  It is too narrow for any data point, so the chi2 never objects.
the evaluator 2026-11-27 chose "a check in the polish that rejects such steps".

Method: for each J^pi group, the phase of det A^-1 (``pyazr.Session.jpi_block``) on the real E_x axis.  A pole E_p - i G/2
near the axis turns the phase by ~pi across ~G, so the steepest slope gives G ~ 2 / max d(arg)/dE.  A coarse grid cannot step
over it (the whole turn shows between two neighbours however narrow the pole); intervals where the phase moves by more than
0.3 rad are zoomed until the slope is resolved.  A narrow pole (G < gamma_max) is SPURIOUS when no level of its own J^pi lies
within max(near, near_factor * G) of it, the level energies taken at the same parameter vector (jpi_block's ``levels``).

Use: the policy key ``pole_guard`` (BoundsPolicy.pole_guard), e.g.
    {"e_min": 4.15, "e_max": 14.0, "step": 0.002, "gamma_max": 0.001, "near": 0.002, "near_factor": 10, "jpi": "all"}
(energies in MeV of E_x).  Below e_min (from bound_min) the guard also finds bound-state poles (zeros of the real det A^-1
where every channel is closed) and calls them spurious off-level the same way.  ``polish`` rejects every trial step with
more spurious poles than its start point, counted separately in the continuum and below threshold.
"""
import os
from typing import Dict, List, Optional, Tuple

import numpy as np

DEFAULTS = dict(e_min=None, e_max=None, step=0.002, gamma_max=0.001, near=0.002, near_factor=10.0, jpi="all",
                bound_min=0.0, bound_step=0.002)


def _phase(sess, jg, E, x):
    b = sess.jpi_block(jg, E, x)
    return np.unwrap(np.angle(np.asarray(b["det_re"]) + 1j * np.asarray(b["det_im"]))), b


def narrow_poles(sess, jg, e_lo, e_hi, x=None, step=2e-3, gmax=1e-3, depth=10) -> List[Tuple[float, float]]:
    """[(E_pole, Gamma)] (MeV) of the poles of J group ``jg`` with Gamma < ``gmax`` in [e_lo, e_hi] of E_x."""
    found = []

    def zoom(a, b, d):
        if not b - a > 1e-12 * max(1.0, abs(a)):         # the turn is narrower than the float spacing: Gamma ~ 0
            return 0.5 * (a + b), 0.0
        E = np.linspace(a, b, 41)
        ph, _ = _phase(sess, jg, E, x)
        dph = np.abs(np.diff(ph))
        i = int(np.argmax(dph))
        slope = dph[i] / (E[1] - E[0])
        if dph[i] < 0.3 or d == 0:
            return 0.5 * (E[i] + E[i + 1]), (2.0 / slope if slope > 0 else np.inf)
        return zoom(E[max(i - 1, 0)], E[min(i + 2, 40)], d - 1)

    E = np.arange(e_lo, e_hi + step / 2, step)
    ph, _ = _phase(sess, jg, E, x)
    dph = np.abs(np.diff(ph))
    for i in np.flatnonzero(dph > 0.3):
        e, g = zoom(E[i], E[i + 1], depth)
        if g < gmax:
            found.append((float(e), float(g)))
    merged: List[Tuple[float, float]] = []
    for e, g in sorted(found):
        if merged and abs(e - merged[-1][0]) < 2 * step:
            if g < merged[-1][1]:
                merged[-1] = (e, g)
        else:
            merged.append((e, g))
    return merged


def bound_poles(sess, jg, e_lo, e_hi, x=None, step=5e-4) -> List[float]:
    """E_x (MeV) of the real-axis zeros of det A^-1 of J group ``jg`` in [e_lo, e_hi] where every channel is closed
    (det_im = 0): bound-state poles.  A level matrix pole pushed below the lowest threshold becomes one of these (an
    archived campaign, 2026-11-27: zeroing or scaling down a level moved the collective pole a few hundred keV, to just
    below the lowest threshold), where the phase test sees nothing."""
    if not e_hi > e_lo:
        return []
    E = np.arange(e_lo, e_hi + step / 2, step)
    b = sess.jpi_block(jg, E, x)
    d, di = np.asarray(b["det_re"]), np.asarray(b["det_im"])
    ok = (di == 0) & np.isfinite(d)
    out = []
    for i in np.flatnonzero(ok[:-1] & ok[1:] & (np.sign(d[:-1]) * np.sign(d[1:]) < 0)):
        a, c = E[i], E[i + 1]
        out.append(float(a - d[i] * (c - a) / (d[i + 1] - d[i])))     # linear interpolation of the zero
    return out


class PoleGuard:
    """``guard(x)`` -> the spurious narrow poles at the evaluator's vector ``x``: [(jpi, E_pole, Gamma, nearest level)]."""

    def __init__(self, ev, cfg: dict, log=None):
        from .backend import PYAZR3
        if not PYAZR3:
            raise NotImplementedError("pole_guard needs pyazr 3 (Session.jpi_block)")
        import pyazr
        from .keys import build_keys
        self.cfg = dict(DEFAULTS, **(cfg or {}))
        azr = getattr(ev, "base_azr", None) or ev.azr
        flags = dict(getattr(ev, "flags", {}) or {})
        self.sess = pyazr.Session(azr, cwd=os.path.dirname(os.path.abspath(azr)), **flags)
        P = self.sess.parameters
        self.x0 = np.array(P.x0, float)
        gkeys = build_keys(self.sess)
        pos = {str(k): j for j, k in enumerate(ev.keys)}
        # guard free index -> evaluator column (level parameters must all map; nuisances do not matter here)
        self.map = []
        for i, k in enumerate(gkeys):
            if k is None:
                continue
            j = pos.get(str(k))
            if j is None:
                if k.kind in ("energy", "width"):
                    raise ValueError(f"pole_guard: level parameter {k} of {os.path.basename(azr)} is not among the evaluator's keys")
                continue
            self.map.append((i, j))
        self.groups: Dict[int, str] = {}
        for p in P:
            if p.kind in ("width", "energy") and p.jgroup is not None:
                self.groups[int(p.jgroup)] = p.jpi
        want = self.cfg["jpi"]
        if want not in (None, "all"):
            want = set([want] if isinstance(want, str) else want)
            self.groups = {g: j for g, j in self.groups.items() if j in want}
        lv = [p.level_energy for p in P if p.kind in ("width", "energy")]
        self.e_min = float(self.cfg["e_min"] if self.cfg["e_min"] is not None else min(lv) - 0.5)
        self.e_max = float(self.cfg["e_max"] if self.cfg["e_max"] is not None else max(e for e in lv if e < 30) + 0.5)
        self._last = (None, None)
        if log:
            log(f"  pole guard: {len(self.groups)} J^pi group(s) {sorted(set(self.groups.values()))}, E_x {self.e_min:g}-{self.e_max:g} MeV, "
                f"step {self.cfg['step'] * 1e3:g} keV, narrow < {self.cfg['gamma_max'] * 1e3:g} keV, off-level beyond "
                f"max({self.cfg['near'] * 1e3:g} keV, {self.cfg['near_factor']:g} Gamma)")

    def vector(self, x):
        g = self.x0.copy()
        x = np.asarray(x, float)
        for i, j in self.map:
            g[i] = x[j]
        return g

    def __call__(self, x) -> List[Tuple[str, float, float, float]]:
        if self._last[0] is not None and np.array_equal(self._last[0], x):
            return self._last[1]
        g = self.vector(x)
        out = []
        for jg, jpi in sorted(self.groups.items()):
            poles = narrow_poles(self.sess, jg, self.e_min, self.e_max, g, step=self.cfg["step"], gmax=self.cfg["gamma_max"])
            bound = [(e, 0.0) for e in bound_poles(self.sess, jg, float(self.cfg["bound_min"]), self.e_min, g,
                                                   step=self.cfg["bound_step"])]
            if not poles and not bound:
                continue
            levels = [float(L["E_level"]) for L in self.sess.jpi_block(jg, [self.e_min], g)["levels"]]
            for e, gam in poles + bound:
                near = min(levels, key=lambda L: abs(L - e)) if levels else np.inf
                if abs(near - e) > max(self.cfg["near"], self.cfg["near_factor"] * gam):
                    out.append((jpi, e, gam, near))
        self._last = (np.array(x, float), out)
        return out

    def counts(self, bad) -> Tuple[int, int]:
        """(continuum, bound) numbers of spurious poles: a bound one is below e_min (where the phase test starts)."""
        nb = sum(1 for _, e, _, _ in bad if e < self.e_min)
        return len(bad) - nb, nb

    def worse(self, bad, start) -> bool:
        """More spurious poles than ``start`` in the continuum or below threshold (a bound pole climbing into the
        continuum counts as worse)."""
        c, b = self.counts(bad)
        c0, b0 = self.counts(start)
        return c > c0 or b > b0

    @staticmethod
    def text(bad) -> str:
        return "; ".join(f"{j} E_x {e:.5f} MeV " + (f"Gamma {g * 1e6:.1f} eV" if g > 0 else "BOUND") + f" (nearest level {n:.4f})"
                         for j, e, g, n in bad) or "none"

    def close(self):
        try:
            self.sess.close()
        except Exception:
            pass
