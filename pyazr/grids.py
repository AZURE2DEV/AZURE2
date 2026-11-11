"""Extrapolation grids of project format 3 (docs/upgrade/project_format3.md section 3.8), through the engine.

The engine holds the one implementation of the grid forms (``experiment/src/Grids.cpp``) and of the adaptive
refinement (``experiment/src/AdaptiveGrid.cpp``); this module hands them over, so that a grid computed here is the
grid an ``extrapolation`` record with the same text gives.
"""

__all__ = ["extrapolation_grid", "refine_grid"]


def extrapolation_grid(grid, masses=None, relativistic=False, project_dir="", angles=False):
    """The values of a grid written as in an ``extrapolation`` record's ``energies`` (or ``angles``) field.

    ``grid`` is the field's text, e.g. ``"range 20 800 keV cm points 200 log"``, ``"list 0.1 0.2 MeV lab"``,
    ``'file "grid.txt" keV cm'``.  Returns a list of lab energies in MeV (angles: deg, in the frame written).  A c.m.
    energy grid needs ``masses = (m1, m2)`` of the entrance pair in u (as ``<levels>`` holds them): the energies are
    taken to the lab frame as the engine does, ``E_lab = E_cm (m1 + m2)/m2`` (``relativistic=True``: the relativistic
    inverse).  A grid file is read relative to ``project_dir``.  An ``adaptive`` grid gives its start grid (the
    refinement needs the model: :func:`refine_grid`, or the engine's extrapolation run).  A grid that does not parse
    or cannot be built raises ``ValueError`` with the engine's message.
    """
    from .described import _engine as _engine_module
    m1, m2 = masses if masses is not None else (0.0, 0.0)
    ok, values, message = _engine_module().extrapolation_grid(str(grid), float(m1), float(m2), bool(relativistic),
                                                              str(project_dir), bool(angles))
    if not ok:
        raise ValueError(message)
    return list(values)


def refine_grid(start, model, tolerance, max_points):
    """Refine a grid as the engine refines an ``adaptive`` extrapolation grid (the same C++ function).

    ``start``: the abscissae to start from; ``model(xs)``: the values at a list of abscissae, one number or one
    sequence of numbers (all of the same length) per abscissa.  Each sweep evaluates the midpoints of the intervals
    made in the sweep before in one call and splits an interval where, for any component,
    ``|f(m) - (f(a) + f(b))/2| > tolerance * max(|f(m)|, 1e-3 max|f|)``; it stops when none is split or the grid
    would grow beyond ``max_points`` (docs: "Extrapolation Grids").  Returns a dict: ``x`` (the grid), ``values``
    (one list of components per point), ``budget_reached``, ``largest_deviation`` and ``largest_deviation_at`` (when
    stopped by ``max_points``), ``sweeps``, ``evaluations``.
    """
    from .described import _engine as _engine_module

    def batch(xs):
        out = []
        for v in model(list(xs)):
            out.append([float(c) for c in v] if hasattr(v, "__len__") else [float(v)])
        return out

    r = _engine_module().refine_grid([float(x) for x in start], batch, float(tolerance), int(max_points))
    return dict(r)
