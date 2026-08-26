"""Adaptive tabulation of AZURE2 cross sections.

Produces, for any (entrance, exit) particle-pair combination of a model, a
*non-uniform* energy grid on which the cross section is tabulated such that
log-linear interpolation between the tabulated points reproduces the engine
to a requested relative tolerance.  Knots concentrate across resonances and
thin out on smooth stretches, so a table that would need thousands of
uniform points typically needs a few hundred adaptive ones.

The scheme is derivative-free interval bisection: an interval between two
tabulated points is accepted once the engine value at its midpoint agrees
with the interpolant within ``rel_tol``; otherwise the midpoint becomes a
knot and both halves are re-examined.  Every evaluated midpoint is kept in
the table (the information is already paid for).  All pairs requested are
refined together, one engine session per refinement round.

Typical use::

    from pyazr.tabulate import tabulate

    tabs = tabulate("fit.azr", pairs=[(1, 2), (1, 3)],
                    e_min=0.5, e_max=3.1, rel_tol=5e-3)
    tabs[(1, 2)].save("sigma_p0.dat")

Notes
-----
- ``e_min``/``e_max`` are centre-of-mass energies of the *entrance* pair by
  default (``frame="cm"``); the data files handed to the engine are written
  in the laboratory frame using the entrance-pair masses from the model.
- Points the engine drops (below a threshold, outside limits) simply do not
  appear in the table; the returned object records the covered range.
- The engine's own parameters are used unless an ``x`` vector (in the
  session's ``params_rwa`` space) is given.
- This tabulates central values.  The analytic parameter Jacobian could be
  used on the same knots to add uncertainty columns from a fit covariance;
  that is deliberately left to the user (pyazr does not own fitting).
"""

import os
import re
import shutil
import tempfile

import numpy as np

from .azrfile import AzrModel
from .azure2 import azure2

__all__ = ["tabulate", "TabulatedCrossSection"]

_SEG_RE = re.compile(r"<segmentsData>.*</segmentsData>", re.S)


class TabulatedCrossSection:
    """Result of an adaptive tabulation for one (entrance, exit) pair."""

    def __init__(self, entrance, exit, e_cm, e_lab, sigma, rel_tol, n_rounds):
        self.entrance = entrance
        self.exit = exit
        self.e_cm = np.asarray(e_cm, float)
        self.e_lab = np.asarray(e_lab, float)
        self.sigma = np.asarray(sigma, float)
        self.rel_tol = rel_tol
        self.n_rounds = n_rounds

    def __repr__(self):
        return (f"<TabulatedCrossSection {self.entrance}->{self.exit} "
                f"n={len(self.e_cm)} E_cm={self.e_cm.min():.4g}-"
                f"{self.e_cm.max():.4g} MeV rel_tol={self.rel_tol}>")

    def __call__(self, e_cm):
        """Log-linear interpolation at arbitrary c.m. energies."""
        e = np.atleast_1d(np.asarray(e_cm, float))
        ln = np.interp(e, self.e_cm, np.log(np.maximum(self.sigma, 1e-300)))
        out = np.exp(ln)
        return out if out.size > 1 else float(out[0])

    def save(self, path, header=True):
        """Write ``E_cm  E_lab  sigma`` columns (MeV, MeV, barn)."""
        with open(path, "w") as f:
            if header:
                f.write(f"# AZURE2 adaptive tabulation: entrance pair "
                        f"{self.entrance} -> exit pair {self.exit}\n"
                        f"# {len(self.e_cm)} points, log-linear interpolation "
                        f"good to {self.rel_tol:g} relative\n"
                        f"# E_cm[MeV]      E_lab[MeV]     sigma[b]\n")
            for ec, el, s in zip(self.e_cm, self.e_lab, self.sigma):
                f.write(f"{ec:.8e} {el:.8e} {s:.8e}\n")


def _entrance_masses(model, entrance):
    for lv in model.levels:
        for c in lv.channels:
            if c.pair == entrance:
                return float(c.M1), float(c.M2)
    raise ValueError(f"entrance pair {entrance} not found in the model")


def tabulate(azr_file, pairs, e_min, e_max, rel_tol=5e-3, frame="cm",
             x=None, min_points=33, max_points=4000, max_rounds=24,
             de_min=None, observable="angle-integrated", cwd=None,
             session_opts=None, verbose=False):
    """Adaptively tabulate cross sections for the given pair combinations.

    ``pairs`` is a list of ``(entrance, exit)`` tuples (file pair keys).
    Returns ``{(entrance, exit): TabulatedCrossSection}``.
    """
    session_opts = dict(session_opts or {})
    pairs = [tuple(p) for p in pairs]
    src = os.path.abspath(azr_file)
    base_cwd = cwd or os.path.dirname(src)

    model0 = AzrModel.from_file(src)
    m1, m2 = _entrance_masses(model0, pairs[0][0])
    lab = {}
    for ent, _ in pairs:
        a, b = _entrance_masses(model0, ent)
        lab[ent] = (a + b) / b       # E_lab = E_cm * (m1+m2)/m2

    if frame == "cm":
        span = (float(e_min), float(e_max))
    elif frame == "lab":
        f0 = lab[pairs[0][0]]
        span = (float(e_min) / f0, float(e_max) / f0)
    else:
        raise ValueError("frame must be 'cm' or 'lab'")
    if de_min is None:
        de_min = (span[1] - span[0]) / 2 ** 16

    work = tempfile.mkdtemp(prefix="azr_tab_", dir=base_cwd)
    os.makedirs(os.path.join(work, "output"), exist_ok=True)
    os.makedirs(os.path.join(work, "checks"), exist_ok=True)
    txt = open(src).read()
    txt = _SEG_RE.sub("<segmentsData>\n</segmentsData>", txt)
    # keep engine side effects inside the scratch directory
    txt = re.sub(r"^\s*\S*output\S*/?\s*(#Full Path to Output Directory)",
                 "output/                                    \\1", txt, flags=re.M)
    txt = re.sub(r"^\s*\S*checks\S*/?\s*(#Full Path to Checks Directory)",
                 "checks/                                    \\1", txt, flags=re.M)
    stripped = os.path.join(work, "_stripped.azr")
    open(stripped, "w").write(txt)

    knots = {p: None for p in pairs}     # (E_cm sorted, sigma) arrays
    active = {p: None for p in pairs}    # list of (lo, hi) c.m. intervals

    def evaluate(requests):
        """requests: {pair: sorted c.m. energy array} -> {pair: (E_cm, sigma)}"""
        m = AzrModel.from_file(stripped)
        keys = [p for p in pairs if p in requests and len(requests[p])]
        for p in keys:
            ent, ex = p
            e = np.asarray(requests[p], float)
            fn = os.path.join(work, f"tab_{ent}_{ex}.dat")
            np.savetxt(fn, np.c_[e * lab[ent], np.zeros(e.size),
                                 np.ones(e.size), np.ones(e.size)],
                       fmt="%.10g")
            m.add_data_segment(fn, entrance=ent, exit=ex,
                               observable=observable, vary_norm=False,
                               energy_min=0.0, energy_max=1e5)
        run = os.path.join(work, "_run.azr")
        m.write(run)
        out = {}
        with azure2(run, cwd=work, **session_opts) as a:
            v = np.asarray(a.params_rwa, float) if x is None else np.asarray(x, float)
            c = a.calculate_rwa(v)
            e = a.calculate_energies(v)
            for i, p in enumerate(keys):
                out[p] = (np.asarray(e[i], float), np.asarray(c[i], float))
        return out

    # round 0: uniform seed
    seed = np.linspace(span[0], span[1], max(int(min_points), 5))
    res = evaluate({p: seed for p in pairs})
    for p in pairs:
        ec, sg = res.get(p, (np.array([]), np.array([])))
        o = np.argsort(ec)
        knots[p] = [ec[o], np.maximum(sg[o], 0.0)]
        active[p] = list(zip(ec[o][:-1], ec[o][1:]))

    tol = np.log1p(rel_tol)
    n_rounds = 1
    for _ in range(max_rounds):
        req = {}
        for p in pairs:
            mids = [0.5 * (lo + hi) for lo, hi in active[p]
                    if (hi - lo) > de_min and len(knots[p][0]) < max_points]
            if mids:
                req[p] = np.array(sorted(mids))
        if not req:
            break
        res = evaluate(req)
        n_rounds += 1
        for p in pairs:
            if p not in res:
                active[p] = []
                continue
            em, sm = res[p]
            E, S = knots[p]
            lnS = np.log(np.maximum(S, 1e-300))
            new_active = []
            for lo, hi in active[p]:
                if (hi - lo) <= de_min or len(E) >= max_points:
                    continue
                mid = 0.5 * (lo + hi)
                j = np.argmin(np.abs(em - mid))
                if abs(em[j] - mid) > 0.26 * (hi - lo):
                    continue                     # engine dropped this point
                s_true = max(sm[j], 1e-300)
                s_hat = np.exp(np.interp(em[j], E, lnS))
                if abs(np.log(s_true) - np.log(s_hat)) > tol:
                    new_active += [(lo, em[j]), (em[j], hi)]
            # merge every evaluated midpoint into the table
            E = np.concatenate([E, em]); S = np.concatenate([S, np.maximum(sm, 0.0)])
            o = np.argsort(E)
            E, S = E[o], S[o]
            keep = np.concatenate([[True], np.diff(E) > 1e-12])
            knots[p] = [E[keep], S[keep]]
            active[p] = new_active
        if verbose:
            tots = {p: (len(knots[p][0]), len(active[p])) for p in pairs}
            print(f"round {n_rounds}: " +
                  " ".join(f"{p}:{n}({a})" for p, (n, a) in tots.items()),
                  flush=True)

    out = {}
    for p in pairs:
        E, S = knots[p]
        out[p] = TabulatedCrossSection(p[0], p[1], E, E * lab[p[0]], S,
                                       rel_tol, n_rounds)
    shutil.rmtree(work, ignore_errors=True)
    return out
