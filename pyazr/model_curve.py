"""The model of a data set as a curve on a fine grid over the data set's range (``Session.model_curve``).

A data set's model is computed at its points (``Session.evaluate(x).data[i].model``).  For a data set that is a
function of one variable -- an excitation function (fixed angle, energy varying) or an angular distribution (fixed
energy, angle varying) -- :func:`model_curves` gives the model between the points as well: the engine's own
evaluation of an extrapolation record made from the data set's record (its reaction, quantity, unit, pair keys,
target, beam, detector, quadrature and the quantity's fields; ``kinds()["records"]["extrapolation"]["fields"]``
says which fields an extrapolation takes, and only those are copied) on a grid, in a session of the extrapolations
of the same project and options, at the data session's fit vector carried by name
(:func:`pyazr.parameters.carry_model_parameters`).

The grid.  Per branch (one fixed value of the other variable) ``points`` values evenly spaced between the branch's
first and last point, together with the data set's own points.  Energies are the beam's lab energies of the points
(``extra["energy_lab"]``, after the data set's energy shift at x), so the curve at a point's energy is the model of
that point without the shift applied twice; the fixed angle of an excitation function is the frame in which the
points' angle is one value (where it is one value in both frames, as for elastic scattering at a lab angle, the
frame of the file's angle column); an angular distribution's angles are c.m.  A narrow resonance
between two grid values is not seen: ``points`` sets the spacing (the engine's adaptive grid needs a run of the
program, ``pyazr.refine_grid``).

What the curve is at the data set's points.  For a data set without a fold the curve's value at a point is the
model of :meth:`Session.evaluate` to rounding (T1, T2, an alpha + 12C ratio to Rutherford: identical).  A folded
point (a target, a beam's energy spread, a detector's resolution) is integrated over the nodes of a mesh, which the
engine plans per compound nucleus from the folded points of the session (``docs/upgrade/energy_grid.md`` section
2.2; a data set with an energy shift, and THM, on grids of their own).  :func:`model_curves` makes the curves of
several data sets in one session, so their folds share one mesh as the data's do: on T2 an unshifted folded data set
agrees with the data session to 4e-9, a shifted one or THM to 8e-5 (the quadrature).  The curve of one data set alone
(:meth:`Session.model_curve`) is folded on a mesh of its own points: 2.4e-3 on T2.

Refused, with the three parts (:class:`pyazr.ParameterError`): a data set not in the evaluation (inactive or
unknown: :class:`pyazr.ParameterKeyError`), a session of the extrapolations, a data set whose columns give each
point something a grid point has no value for (a beam-profile window, a Legendre order, a decay energy: the roles
outside :data:`CURVE_COLUMNS`), columns read from the file's header, and points that are not one function of one
variable (scattered in energy and angle, more than :data:`MAX_BRANCHES` branches, a branch of one point).
"""

import dataclasses
import os
import re
import shutil
import tempfile

import numpy as np

from .errors import FormatError, ParameterError, ParameterKeyError, close_names, listed

#: The column roles of a data set that a curve can be made for: the point's coordinates and its datum.  Any other
#: role (``schema()["column_roles"]``) gives each point a value a grid point does not have.
CURVE_COLUMNS = ("E_lab", "E_cm", "theta_lab", "theta_cm", "value", "error", "error_percent", "_")

#: The grid values per branch between its first and last point, by default.
DEFAULT_POINTS = 200

#: The most branches (fixed angles of excitation functions, or fixed energies of angular distributions) of one curve.
MAX_BRANCHES = 64

#: The most grid points of one data set's curve (all branches).
MAX_POINTS = 20000

# two coordinates are one value when they agree to this relative tolerance (the engine's own numbers of one file
# value through the same conversion are identical; the tolerance only absorbs the last bits)
_SAME = 1e-9


@dataclasses.dataclass(frozen=True)
class CurveBranch:
    """One line of a curve: the grid entries ``start`` to ``stop`` (a slice of the curve's arrays), along which
    ``variable`` varies, at the fixed ``fixed`` (``"angle_lab"``, ``"angle_cm"`` in deg, ``"energy_lab"`` in MeV;
    ``None`` for an angle-integrated quantity) of value ``value``."""
    start: int
    stop: int
    fixed: str
    value: float


@dataclasses.dataclass
class ModelCurve:
    """The model of one data set on a grid (:func:`model_curves`).

    - ``name``, ``quantity``, ``unit`` (the engine's unit of ``model``: the unit of ``Evaluation.data[i].model``),
      ``variable`` (``"energy"`` or ``"angle"``: what varies along each branch).
    - Per grid point: ``energies`` (MeV, c.m.), ``energy_lab`` (MeV, the beam's lab energy), ``angles`` (deg, c.m.),
      ``excitation_energies`` (MeV), ``model``, ``sfactor_conversion`` (the engine's factor from a cross section in b
      to the S factor in MeV b, 0 where it has none), ``at_data`` (True where the grid value is one of the data
      set's own points).  ``extra`` holds ``energy_lab``, as ``EvaluatedData.extra`` does, so that
      ``pyazr.endf.Material.curves(curve, session=s)`` takes a curve as it takes a data set.
    - ``branches``: :class:`CurveBranch` es in grid order; ``x`` a copy of the data session's fit vector.
    """
    name: str
    quantity: str
    unit: str
    variable: str
    energies: np.ndarray
    energy_lab: np.ndarray
    angles: np.ndarray
    excitation_energies: np.ndarray
    model: np.ndarray
    sfactor_conversion: np.ndarray
    at_data: np.ndarray
    branches: tuple
    x: np.ndarray

    @property
    def extra(self):
        return {"energy_lab": self.energy_lab}

    @property
    def n(self):
        """The number of grid points."""
        return len(self.model)

    def table(self):
        """The grid as text: per point its branch, c.m. energy, lab energy, c.m. angle and the model."""
        lines = [f"{self.name}  ({self.quantity}, {self.unit}; {self.n} points along the {self.variable} in "
                 f"{len(self.branches)} branch{'es' if len(self.branches) != 1 else ''})",
                 f"{'branch':>6} {'E_cm/MeV':>14} {'E_lab/MeV':>14} {'angle/deg':>10} {'model':>14}"]
        for b_index, b in enumerate(self.branches):
            for k in range(b.start, b.stop):
                lines.append(f"{b_index:>6} {self.energies[k]:14.7g} {self.energy_lab[k]:14.7g} "
                             f"{self.angles[k]:10.5g} {self.model[k]:14.7g}")
        return "\n".join(lines)

    def __repr__(self):
        return (f"ModelCurve({self.name!r}, {self.quantity}, {self.n} points along the {self.variable}, "
                f"{len(self.branches)} branches)")


@dataclasses.dataclass(frozen=True)
class CurvePlan:
    """The extrapolation records of the curves of some data sets at one x (:func:`plan_curves`): per data set the
    records' fields (one record per branch) and the branches' layout, or the refusal; ``text`` names the plan (two
    plans with the same text need the same session)."""
    names: tuple
    records: dict      # name -> [fields of the extrapolation record of each branch]
    layouts: dict      # name -> (variable, [(fixed, value, grid (MeV lab or deg c.m.), at_data mask)])
    refused: dict      # name -> ParameterError
    quantities: dict   # name -> (quantity, unit)
    text: str


def _refuse(name, what, allowed, fix):
    return ParameterError(f"no curve of {name!r}: {what}; allowed: {allowed}; fix: {fix}", what=what,
                          allowed=allowed, fix=fix)


def _unique(values):
    """The distinct values of ``values`` (sorted), two values that agree to :data:`_SAME` taken as one."""
    v = np.sort(np.asarray(values, float))
    if v.size == 0:
        return v
    keep = [v[0]]
    for a in v[1:]:
        if abs(a - keep[-1]) > _SAME * max(abs(a), abs(keep[-1]), 1e-300):
            keep.append(a)
    return np.asarray(keep)


def _same(a, b):
    return abs(a - b) <= _SAME * max(abs(a), abs(b), 1e-300)


def _grid(values, points):
    """``points`` values evenly spaced between the first and last of ``values``, with ``values`` themselves; sorted,
    and a mask of the entries that are ``values``."""
    own = _unique(values)
    even = np.linspace(own[0], own[-1], max(int(points), 2)) if own.size > 1 else own
    merged = _unique(np.concatenate([own, even]))
    # an even value that agrees with a datum to the tolerance is the datum's value (merged keeps the smaller)
    at_data = np.array([any(_same(g, o) for o in own) for g in merged]) if own.size <= 2000 else np.isin(merged, own)
    for i, g in enumerate(merged):
        if at_data[i]:
            merged[i] = own[np.argmin(np.abs(own - g))]
    return merged, at_data


def _numbers(values):
    return " ".join(repr(float(v)) for v in values)


def _branches(name, data, points, angle_frame=None):
    """The branches of a data set's points (``data``: its ``EvaluatedData``): ``(variable, [(fixed, value, grid,
    at_data)])``, the grid in MeV lab for an energy branch, in deg c.m. for an angle branch.  ``angle_frame``: the
    frame of the data file's angle column (``"angle_lab"``, ``"angle_cm"``), which a fixed angle takes when the
    points' angle is one value in both frames."""
    lab = np.asarray(data.extra["energy_lab"], float)
    if "angle_lab" not in data.extra:  # an angle-integrated quantity: one excitation function
        if _unique(lab).size < 2:
            raise _refuse(name, "every point at one energy", "a data set with points at two energies or more",
                          "Session.evaluate(x).dataset(name) shows its points")
        grid, at = _grid(lab, points)
        return "energy", [(None, None, grid, at)]
    angle_lab, angle_cm = np.asarray(data.extra["angle_lab"], float), np.asarray(data.angles, float)
    energies = _unique(lab)
    # the frame in which the fixed angles are fewest: a detector at a lab angle, or c.m. angles from the file; on a
    # tie (elastic scattering: one c.m. angle for one lab angle) the file's own, else c.m. (the engine computes the
    # quantities that have no lab-angle form, an analyzing power, at c.m. angles)
    frames = (("angle_cm", angle_cm), ("angle_lab", angle_lab))
    if angle_frame == "angle_lab":
        frames = frames[::-1]
    frame, angles = min(frames, key=lambda f: _unique(f[1]).size)
    fixed_angles = _unique(angles)
    by_angle = [(a, np.array([_same(v, a) for v in angles])) for a in fixed_angles]
    by_energy = [(e, np.array([_same(v, e) for v in lab])) for e in energies]
    excitation = all(_unique(lab[m]).size >= 2 for _, m in by_angle)
    angular = all(_unique(angle_cm[m]).size >= 2 for _, m in by_energy)
    choices = []
    if excitation:
        choices.append((len(by_angle), "energy"))
    if angular:
        choices.append((len(by_energy), "angle"))
    if not choices:
        raise _refuse(name, f"{len(lab)} points at {energies.size} energies and {fixed_angles.size} angles, not "
                      "one function of one variable", "points at fixed angles (each at two energies or more) or at "
                      "fixed energies (each at two angles or more)",
                      "an extrapolation record with the grid you want (Project.add('extrapolation', ...))")
    count, variable = min(choices)
    if count > MAX_BRANCHES:
        raise _refuse(name, f"{count} branches (one per fixed {'angle' if variable == 'energy' else 'energy'})",
                      f"at most {MAX_BRANCHES}", "an extrapolation record with the grid you want")
    out = []
    if variable == "energy":
        for a, m in by_angle:
            grid, at = _grid(lab[m], points)
            out.append((frame, float(a), grid, at))
    else:
        for e, m in by_energy:
            grid, at = _grid(angle_cm[m], points)
            out.append(("energy_lab", float(e), grid, at))
    return variable, out


def _columns_refusal(name, fields, roles):
    columns = str(fields.get("columns", "")).split()
    if columns == ["from_file"]:
        return _refuse(name, "its columns are read from the data file's header (from_file)", "the roles written in "
                       "the record's columns field", "Project.set(name, columns='E_lab[MeV] _ value error')")
    for token in columns:
        role = token.split("[")[0]
        if role not in CURVE_COLUMNS and role in roles:
            return _refuse(name, f"its column {role} gives each point its own value, which a grid point has none of",
                           "the columns " + ", ".join(CURVE_COLUMNS), "compare at the points: "
                           "Session.evaluate(x).dataset(name).model")
    return None


def plan_curves(session, names=None, x=None, points=DEFAULT_POINTS):
    """The extrapolation records of the curves of the data sets ``names`` (default: every active data set) of the
    data session ``session`` at its fit vector ``x`` (default ``parameters.x0``), ``points`` grid values per branch:
    a :class:`CurvePlan`.  A data set that cannot be drawn as a curve is in ``plan.refused``, with the three parts."""
    from . import kinds, schema
    from .project import Project
    if session.mode != "data":
        raise ParameterError("a curve is made from a session of the data", what="a session of the extrapolations",
                             allowed="pyazr.Session(file) (extrapolations=False)",
                             fix="the extrapolations are curves already: Session(file, extrapolations=True)"
                                 ".evaluate(x)")
    points = int(points)
    if points < 2:
        raise ParameterError(f"points={points}", what=f"points={points}", allowed="2 or more grid values per branch",
                             fix=f"model_curve(name, points={DEFAULT_POINTS})")
    ev = session.evaluate(x)
    by_name = {d.name: d for d in ev.data}
    names = list(by_name) if names is None else list(names)
    for n in names:
        if n not in by_name:
            raise ParameterKeyError(f"no active data set {n!r} in the session{close_names(n, by_name)}",
                                    what=f"the data set {n!r}", allowed="an active data set: " + listed(by_name),
                                    fix="Session.evaluate(x).data lists them")
    takes = set(kinds()["records"]["extrapolation"]["fields"]) - {"energies", "angles", "active", "orders"}
    roles = set(schema()["column_roles"])
    records_of = Project.open(session.file).records()
    records, layouts, refused, quantities = {}, {}, {}, {}
    for n in names:
        d = by_name[n]
        fields = records_of[n].fields
        quantities[n] = (d.quantity, d.unit)
        try:
            refusal = _columns_refusal(n, fields, roles)
            if refusal is not None:
                raise refusal
            own = [t.split("[")[0] for t in str(fields.get("columns", "")).split()]
            frame = "angle_lab" if "theta_lab" in own else "angle_cm" if "theta_cm" in own else None
            variable, branches = _branches(n, d, points, frame)
        except ParameterError as e:
            refused[n] = e
            continue
        if sum(len(g) for _, _, g, _ in branches) > MAX_POINTS:
            refused[n] = _refuse(n, f"{sum(len(g) for _, _, g, _ in branches)} grid points", f"at most {MAX_POINTS}",
                                 "fewer points per branch")
            continue
        base = {k: v for k, v in fields.items() if k in takes}
        recs = []
        for fixed, value, grid, _ in branches:
            f = dict(base)
            if variable == "energy":
                f["energies"] = f"list {_numbers(grid)} MeV lab"
                if fixed is not None:
                    f["angles"] = f"list {_numbers([value])} deg {'lab' if fixed == 'angle_lab' else 'cm'}"
            else:
                f["energies"] = f"list {_numbers([value])} MeV lab"
                f["angles"] = f"list {_numbers(grid)} deg cm"
            recs.append(f)
        records[n] = recs
        layouts[n] = (variable, [(fixed, value, grid, at) for fixed, value, grid, at in branches])
    text = "\n".join(f"{n} {r}" for n in names if n in records for r in records[n])
    return CurvePlan(tuple(names), records, layouts, refused, quantities, text)


class CurveSession:
    """The session of the extrapolations that evaluates a :class:`CurvePlan`: one ``.azr`` with the plan's records
    as its only extrapolations (``pyazr.bands`` trims the project), written with its output directory in a temporary
    directory of its own (the external-capture integrals of its grid apart from any other's; nothing beside the
    project), opened with the data session's options, the data session's directory for relative paths; the
    directory is removed when the session is closed.
    :meth:`evaluate` gives the curves at a fit vector of the data session; :meth:`close` (or a ``with`` block)
    closes the session."""

    def __init__(self, session, plan):
        self.data_session, self.plan = session, plan
        self.refused = dict(plan.refused)
        self.session = None
        # a record the engine refuses (a quantity it does not compute on a grid of that frame, ...) refuses its data
        # set with the engine's words; the others are made without it
        while True:
            self.order = [(n, i) for n in plan.names if n in plan.records and n not in self.refused
                          for i in range(len(plan.records[n]))]
            if not self.order:
                return
            try:
                self._open()
                return
            except FormatError as e:
                names = self._refused_by(e)
                if not names:
                    raise
                for n in names:
                    self.refused[n] = _refuse(n, f"the engine refuses its curve: {e.what or e}", e.allowed or
                                              "a record the engine computes on a grid", e.fix or
                                              "compare at the points: Session.evaluate(x).dataset(name).model")

    def _refused_by(self, error):
        """The data sets whose records ``error`` (the engine's diagnostics of the curves' project) names."""
        text = " ".join([str(error)] + [str(d.get("record") or "") for d in error.diagnostics])
        found = {int(k) for k in re.findall(r"\bcurve_(\d+)\b", text)}
        return sorted({self.order[k - 1][0] for k in found if 1 <= k <= len(self.order)})

    def _open(self):
        from .bands import _trimmed_project
        from .session import Session
        grids = [dict(self.plan.records[n][i], name=f"curve_{k + 1}") for k, (n, i) in enumerate(self.order)]
        project, _ = _trimmed_project(self.data_session.file, None, grids)
        # a directory of its own for the project and its output (the external-capture integrals of this grid apart
        # from any other's); nothing is written beside the project
        self.workdir = tempfile.mkdtemp(prefix="pyazr_curves_")
        project.set_config(output_dir=self.workdir + "/")
        self.path = os.path.join(self.workdir, "curves.azr")
        try:
            project.save(self.path)
            self.session = Session(self.path, cwd=self.data_session.cwd,
                                   **dict(self.data_session.options, extrapolations=True))
        except BaseException:
            shutil.rmtree(self.workdir, ignore_errors=True)
            raise

    def evaluate(self, x=None):
        """``{name: ModelCurve or ParameterError}`` at the data session's fit vector ``x`` (default its
        ``parameters.x0``), every name of the plan."""
        from .parameters import carry_model_parameters
        s = self.data_session
        x = np.array(s.parameters.x0 if x is None else x, float).ravel()
        out = dict(self.refused)
        if self.session is None:
            return out
        grids = self.session.evaluate(carry_model_parameters(x, s.parameters, self.session.parameters)).data
        parts = {}
        for (n, i), g in zip(self.order, grids):
            parts.setdefault(n, []).append(g)
        for n, evaluated in parts.items():
            variable, layout = self.plan.layouts[n]
            branches, start = [], 0
            for (fixed, value, grid, _), g in zip(layout, evaluated):
                if g.n != len(grid):
                    raise RuntimeError(f"the curve of {n!r}: the engine gave {g.n} points for a grid of {len(grid)}")
                branches.append(CurveBranch(start, start + len(grid), fixed, value))
                start += len(grid)

            def cat(key):
                return np.concatenate([np.asarray(getattr(g, key), float) for g in evaluated])
            # the beam's lab energy of each grid point: the grid of an energy branch, the fixed energy of an angle one
            lab = np.concatenate([grid if variable == "energy" else np.full(len(grid), value)
                                  for _, value, grid, _ in layout])
            quantity, unit = self.plan.quantities[n]
            out[n] = ModelCurve(n, quantity, unit, variable, cat("energies"), lab, cat("angles"),
                                cat("excitation_energies"), cat("model"),
                                np.concatenate([np.asarray(g.extra["sfactor_conversion"], float) for g in evaluated]),
                                np.concatenate([a for _, _, _, a in layout]), tuple(branches), x.copy())
        return out

    def close(self):
        """Close the session and remove its file."""
        if self.session is not None:
            self.session.close()
            self.session = None
            shutil.rmtree(self.workdir, ignore_errors=True)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def model_curves(session, names=None, x=None, points=DEFAULT_POINTS):
    """The model of the data sets ``names`` (default: every active one) of ``session`` on fine grids at the fit
    vector ``x`` (default ``parameters.x0``): ``{name: ModelCurve or ParameterError}`` (a refusal is returned, not
    raised, so that one data set does not stop the others).  One session of the extrapolations for all of them,
    closed before the call returns."""
    plan = plan_curves(session, names, x, points)
    with CurveSession(session, plan) as c:
        return c.evaluate(x)


def model_curve(session, name, x=None, points=DEFAULT_POINTS):
    """The model of the data set ``name`` on a fine grid (:class:`ModelCurve`); a data set that cannot be drawn as a
    curve raises :class:`pyazr.ParameterError` with the three parts."""
    out = model_curves(session, [name], x, points)[name]
    if isinstance(out, Exception):
        raise out
    return out
