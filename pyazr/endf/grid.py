"""The energy grid of the sections the writer fills from the engine (ENDF-102, 2023: MF=3 and MF=6 tabulate on the
lab energy of the projectile, pp. 103, 134), and the cosine grid of a tabulated distribution.

One rule for every grid: :func:`pyazr.refine_grid` (the engine's refinement of an ``adaptive`` extrapolation grid,
the same C++ function) splits an interval while its midpoint differs from the linear interpolation of its ends by
more than ``tolerance x max(|f(m)|, 1e-3 max|f|)`` for any component.  A midpoint can lie on the chord of a curve
that departs from it elsewhere in the interval (the flank of a narrow level), so :func:`refine` checks the quarter
points of every interval with the same rule afterwards and refines again from the points that fail.  The model is
evaluated once per abscissa (a cache), so the check costs one evaluation of the quarter points.
"""

import numpy as np

#: Evenly spaced points of the start grid over the energy range (the level energies and thresholds are added).
START_POINTS = 41
#: Rounds of the quarter-point check: refine, check, refine from the failing points again.
CHECK_ROUNDS = 4


def start_grid(session, x, energy_range_MeV):
    """The start grid of every written section (lab MeV of the projectile, sorted): ``START_POINTS`` evenly spaced
    points over the range, and inside it every threshold of an exit pair and every level energy of the compound with
    the two energies half its width away (E_r -+ Gamma/2).  A refinement that starts there finds a level narrower
    than the spacing (its peak and half-maximum points are grid points), which an even grid finds only when a midpoint
    falls near it.

    Where the lists come from: the levels are pyazr's, the physical values of the session's parameters at ``x``
    (``Session.transform``): the energy of each level (observed E_x in MeV, fixed or free) and its width, the sum of
    its partial widths (eV) of the channels open at E_x (a channel below its pair's threshold holds an ANC, not a
    width); the thresholds are the session's pair records (``Session.pairs``: E_x = separation energy + excitation of
    each exit pair).  E_x becomes the lab energy of the projectile (the entrance pair's first particle) by the
    engine's row ``cm_to_lab_energy`` at E_cm = E_x - S of the entrance pair (non-relativistic: the writer refuses a
    relativistic project before)."""
    from .. import physics
    lo, hi = (float(v) for v in energy_range_MeV)
    P = session.parameters
    rows = [p.to_dict() for p in P]
    physical = session.transform(P.x0 if x is None else np.asarray(x, float))
    pairs = {p.key: p.to_dict() for p in session.pairs}
    entrance = next(p for p in pairs.values() if p["is_entrance"])
    s_ent = entrance["sep_energy"] + entrance["excitation"]
    opens = {k: p["sep_energy"] + p["excitation"] for k, p in pairs.items() if p["pair_kind"] == "particle"}
    levels = {}
    for r, v in zip(rows, physical):
        if r["kind"] == "energy":
            levels.setdefault((r["jgroup"], r["level"]), [0.0, 0.0])[0] = float(v)
    for r, v in zip(rows, physical):
        if r["kind"] != "width":
            continue
        lv = levels[(r["jgroup"], r["level"])]
        if r["radiation_type"] != "P" or lv[0] > opens.get(r["pair_key"], np.inf):
            lv[1] += abs(float(v)) * 1e-6  # eV -> MeV
    points_cm = [p["sep_energy"] + p["excitation"] - s_ent for p in pairs.values()
                 if p["pair_kind"] == "particle" and not p["is_entrance"]]
    for e, width in levels.values():
        points_cm += [e - s_ent, e - s_ent - width / 2, e - s_ent + width / 2]
    e_cm = np.array([e for e in points_cm if e > 0], float)
    lab = np.atleast_1d(physics.cm_to_lab_energy(entrance["M1"], entrance["M2"], e_cm)) if e_cm.size else e_cm
    seeds = [float(e) for e in lab if lo < e < hi]
    return sorted(set(np.linspace(lo, hi, START_POINTS).tolist()) | set(seeds))


def _deviation(values_at, chord, scale):
    """The rule of :func:`pyazr.refine_grid` per point and component: |f - chord| / max(|f|, scale)."""
    return np.abs(values_at - chord) / np.maximum(np.abs(values_at), scale)


def refine(start, model, tolerance, max_points):
    """The grid of :func:`pyazr.refine_grid` from ``start`` with the quarter-point check (the module's text).
    ``model(xs)``: one sequence of components per abscissa.  Returns the dict of :func:`pyazr.refine_grid` (``x``,
    ``values``, ``budget_reached``, ...) with ``check``: the largest deviation at the quarter points by the rule (0
    when every quarter point passes), ``check_rounds`` and ``evaluations`` (abscissae the model computed)."""
    from ..grids import refine_grid
    cache = {}

    def cached(xs):
        todo = [x for x in xs if x not in cache]
        if todo:
            for x, v in zip(todo, model(todo)):
                cache[x] = np.atleast_1d(np.asarray(v, float))
        return [cache[x] for x in xs]

    grid = sorted(float(x) for x in start)
    worst = 0.0
    for rounds in range(1, CHECK_ROUNDS + 1):
        refined = refine_grid(grid, cached, tolerance, max_points)
        if refined.get("budget_reached"):
            break
        x = np.asarray(refined["x"], float)
        v = np.asarray(refined["values"], float).reshape(len(x), -1)
        quarters = np.concatenate([x[:-1] + (x[1:] - x[:-1]) / 4, x[:-1] + 3 * (x[1:] - x[:-1]) / 4])
        fq = np.array(cached(quarters.tolist()), float).reshape(len(quarters), -1)
        chord = np.column_stack([np.interp(quarters, x, v[:, k]) for k in range(v.shape[1])])
        scale = 1e-3 * np.max(np.abs(v), axis=0)
        dev = _deviation(fq, chord, scale).max(axis=1)
        worst = float(dev.max()) if dev.size else 0.0
        failing = quarters[dev > tolerance]
        if not failing.size:
            worst = 0.0
            break
        if len(x) + failing.size > max_points:
            refined["budget_reached"] = True
            refined["largest_deviation"] = worst
            refined["largest_deviation_at"] = float(quarters[int(np.argmax(dev))])
            break
        grid = sorted(set(x.tolist()) | set(failing.tolist()))
    refined["check"] = worst
    refined["check_rounds"] = rounds
    refined["evaluations"] = len(cache)
    return refined


def thin(x, f, tolerance):
    """The subset of a refined grid ``x`` that one component ``f`` needs: the points kept so that linear
    interpolation between them gives every dropped point of ``x`` to the rule of :func:`refine`,
    |f - chord| <= tolerance x max(|f|, 1e-3 max|f|) (a recursive split at the worst point, as Douglas and
    Peucker's).  For a table whose rows are tabulated each on its own abscissae (the cosines of a LAW=5 panel), from a
    grid refined for all rows at once.  Returns the indices kept, increasing (the ends always)."""
    x = np.asarray(x, float)
    f = np.asarray(f, float)
    scale = 1e-3 * np.max(np.abs(f)) if f.size else 0.0
    keep = {0, len(x) - 1}
    stack = [(0, len(x) - 1)]
    while stack:
        i, j = stack.pop()
        if j - i < 2:
            continue
        inner = np.arange(i + 1, j)
        chord = f[i] + (f[j] - f[i]) * (x[inner] - x[i]) / (x[j] - x[i])
        dev = _deviation(f[inner], chord, scale)
        k = int(np.argmax(dev))
        if dev[k] > tolerance:
            m = int(inner[k])
            keep.add(m)
            stack += [(i, m), (m, j)]
    return np.array(sorted(keep))
