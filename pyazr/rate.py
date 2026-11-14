"""The thermonuclear reaction rate N_A<sigma v> of a project, computed by the engine, with its uncertainty band.

The engine's rate (``src/ReactionRate.cpp``, the program's ``AZURE2 rate`` command) integrates

    N_A <sigma v> = N_A (8 / (pi mu))^(1/2) (kT)^(-3/2) int sigma(E) E exp(-E / kT) dE

over the c.m. energy of the entrance pair, from 1e-7 MeV to 50 MeV: sigma is computed once on a grid (48 points per
decade, and around each level of the parameters being evaluated its energy and 41 points on each side from 1/32 to 64
widths), a Steffen spline of ln sigma against ln E is integrated at every temperature (GSL's QAGP, relative accuracy
1e-7, with breakpoints at the grid's nodes). sigma is zero below the grid's lowest positive value and above 50 MeV; a
grid with fewer than three points of sigma > 0 falls back to a direct quadrature per temperature. For identical entrance
particles sigma carries (1 + delta_12) and the rate excludes the 1/(1 + delta_12) of the rate equation (the
compilations' convention). Exit -1 is the total radiative capture (every capture pair of the entrance's compound
nucleus).

:func:`reaction_rate` runs that computation in-process, in a session set up as the program's rate command sets up its
run, at any parameter values, on a temperature grid: the program's ``min + i step`` grid (``t_min``, ``t_max``,
``t_step``), a named grid of :data:`TEMPERATURE_GRIDS`, or a list.  :func:`rate_band` gives the band from a
covariance (linear propagation of ln N_A<sigma v>) or from drawn vectors (the 16th, 50th and 84th percentiles).
:func:`gamow_window` gives the Gamow peak and width per temperature, :func:`write_rate_table` the table as CSV with
a header, :func:`read_rate_table` a table to compare with.  No integral is done in Python.
"""
import csv
import io
import os
import shutil
import tempfile
import weakref
from dataclasses import dataclass, field

import numpy as np

from .rate_formats import STARLIB_T9

__all__ = ["TEMPERATURE_GRIDS", "UNITS", "Rate", "RateBand", "RateSession", "reaction_rate", "rate_band",
           "gamow_window", "write_rate_table", "read_rate_table", "table_at", "temperatures_of"]

#: Units of every array of a result.
UNITS = {"T9": "GK", "rate": "cm^3 mol^-1 s^-1", "E0": "MeV", "width": "MeV", "kT": "MeV"}

#: Boltzmann's constant as the engine has it (include/Constants.h, CODATA 2018), MeV / GK.
BOLTZMANN_MEV_PER_GK = 8.617333262e-2


def _decades(lo, hi, per_decade):
    n = int(round(np.log10(hi / lo) * per_decade))
    return [float(f"{v:.6g}") for v in lo * 10.0 ** (np.arange(n + 1) / per_decade)]


#: Named temperature grids (GK): name -> (what it is, the temperatures).  A new grid is one row.
TEMPERATURE_GRIDS = {
    "table": ("the rate tables' grid from 0.01 to 10 GK (NACRE's temperatures in that range)",
              [0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.018, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.07,
               0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6,
               0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]),
    "log20": ("20 temperatures per decade from 0.01 to 10 GK", _decades(0.01, 10.0, 20)),
    "starlib": ("STARLIB's 60 temperatures from 0.001 to 10 GK (pyazr.rate_formats.STARLIB_T9)", STARLIB_T9),
}


def temperatures_of(grid=None, temperatures=None, t_min=None, t_max=None, t_step=None):
    """``(kind, values)``: ``("grid", (t_min, t_max, t_step))`` for the program's grid, ``("list", array)`` for a named
    grid or a list.  Exactly one of the three ways."""
    ways = [grid is not None, temperatures is not None, t_min is not None or t_max is not None or t_step is not None]
    if sum(ways) != 1:
        raise ValueError("give the temperatures one way: grid=<name> (TEMPERATURE_GRIDS), temperatures=[...] (GK), "
                         "or t_min, t_max and t_step (GK, the program's grid)")
    if grid is not None:
        if grid not in TEMPERATURE_GRIDS:
            raise ValueError(f"no temperature grid {grid!r}: one of {', '.join(TEMPERATURE_GRIDS)}")
        return "list", np.asarray(TEMPERATURE_GRIDS[grid][1], float)
    if temperatures is not None:
        t = np.atleast_1d(np.asarray(temperatures, float))
        if t.size == 0 or not np.all(np.isfinite(t)) or np.any(t <= 0):
            raise ValueError("temperatures: finite numbers > 0 (GK)")
        return "list", t
    if None in (t_min, t_max, t_step) or not (t_min > 0 and t_max >= t_min and t_step >= 0):
        raise ValueError("the program's grid needs 0 < t_min <= t_max and t_step >= 0 (GK; t_step 0: t_min only)")
    return "grid", (float(t_min), float(t_max), float(t_step))


@dataclass
class Rate:
    """A rate on a temperature grid: ``T9`` (GK), ``rate`` (N_A<sigma v>, cm^3 mol^-1 s^-1), the pairs, the reaction
    label, where the parameters came from (``source``), the project file and :data:`UNITS`.  The engine's lines of
    the run (its progress, the grid it built) go to standard output, as the program's."""
    T9: np.ndarray
    rate: np.ndarray
    entrance: int
    exit: int
    reaction: str
    source: str
    project: str
    units: dict = field(default_factory=lambda: dict(UNITS))


@dataclass
class RateBand:
    """A rate's band: ``low``, ``median``, ``high`` (cm^3 mol^-1 s^-1) on ``T9``; ``best`` the rate at the
    parameters given (the centre of a linear band); ``how`` in words; ``method`` "covariance" or "draws"; ``varied``
    the parameters that moved; ``draws`` the number of vectors (draws)."""
    T9: np.ndarray
    low: np.ndarray
    median: np.ndarray
    high: np.ndarray
    best: np.ndarray
    method: str
    how: str
    varied: list
    draws: int = 0
    units: dict = field(default_factory=lambda: dict(UNITS))


class RateSession:
    """An engine session set up for the rate of ``entrance`` to ``exit`` (pair keys of the project ``file``), as the
    program's ``rate`` command sets up its run: the levels, no data.  Its parameters are the level parameters (no
    normalization or shift: the rate does not depend on them); ``names`` and ``values`` (the fit basis, every
    parameter) are the project's, with a parameter file's values when one is given (matched by name, as the program
    matches them; a name the rate has not, a data set's norm or shift, is left out).

    ``rate(values=None, ...)``: the rate at ``values`` ({name: value} over the session's, or a full vector) on the
    temperatures (:func:`temperatures_of`).  Options of the calculation as :class:`pyazr.Session` takes them
    (``options``: {keyword: bool}); ``coulomb_table`` as there.

    The session writes nothing beside the project: the engine's rate run writes its files (the cross-section grid
    ``cross_section_grid.dat``, an external-capture cache) into its output directory, which for this session is a
    temporary directory removed by :meth:`close` (or when the session is collected); the project's relative paths
    are read from ``cwd`` (default: the project's directory).  The engine's text (its progress bars) goes to the
    logger ``pyazr.engine`` (:mod:`pyazr.engine_text`), not to standard output."""

    def __init__(self, file, entrance, exit, *, parameter_file=None, cwd=None, options=None, coulomb_table=None):
        from . import _azure2, engine_text, session_checks, session_options
        from .errors import AZURE2Error
        from .project import Project, require_format3
        require_format3(file)
        self.file = os.path.abspath(file)
        self.cwd = os.path.abspath(cwd or os.path.dirname(self.file) or ".")
        self.entrance, self.exit = int(entrance), int(exit)
        if self.entrance == 0 or self.exit == 0 or self.exit < -1 or self.entrance == self.exit:
            raise ValueError(f"rate pairs: entrance {self.entrance}, exit {self.exit}: pair keys of <levels> (not 0), "
                             f"exit -1 for the total radiative capture, different from the entrance")
        opts = session_options.runtime_options(_azure2, dict(options or {}), False, coulomb_table, self.cwd)
        # the engine's output directory: a temporary one, the project's copy pointing at it (the relative paths are
        # read from cwd, opts.base_directory)
        self._scratch = tempfile.mkdtemp(prefix="pyazr_rate_")
        self._remove_scratch = weakref.finalize(self, shutil.rmtree, self._scratch, True)
        copy = Project.open(self.file)
        copy.set_config(output_dir=self._scratch.rstrip("/") + "/", checks_dir=self._scratch.rstrip("/") + "/")
        engine_file = copy.save(os.path.join(self._scratch, os.path.basename(self.file)))
        said = []
        try:
            with engine_text.held() as said:
                self._binding = _azure2.Session(engine_file, opts, rate_entrance=self.entrance, rate_exit=self.exit)
        except AZURE2Error as err:
            self.close()
            refused = session_checks.setup_refused(err, [x.replace(engine_file, self.file) for x in said], self.file)
            raise type(refused)(str(refused).replace(engine_file, self.file), what=refused.what.replace(
                engine_file, self.file), allowed=refused.allowed, fix=refused.fix) from None
        self.basis = 2 if self._option("use_park") else (1 if self._option("use_brune") else 0)
        n = len(self._binding.params_fixed())
        self.names = [self._binding.params_names(i) for i in range(n)]
        self.fixed = list(self._binding.params_fixed())
        self.values = np.asarray(self._binding.params_all_rwa(), float).copy()
        self.source = "the project's values"
        if parameter_file is not None:
            from .parameters import read_parameter_file
            from .errors import ParameterError
            tag, given = read_parameter_file(os.fspath(parameter_file))
            if tag is not None and tag != self.basis:
                raise ParameterError(f"{parameter_file} holds basis {tag} amplitudes (0 standard, 1 Brune, 2 Park); "
                                     f"the rate session uses basis {self.basis}")
            hit = [i for i, nm in enumerate(self.names) if nm in given]
            if not hit:
                raise ParameterError(f"{parameter_file}: no name matches a level parameter of this project")
            for i in hit:
                self.values[i] = given[self.names[i]]
            self.source = f"the parameter file {os.path.abspath(parameter_file)}"

    def _option(self, name):
        return bool(dict(self._binding.option_sources()).get(name, (False,))[0])

    def close(self):
        """Release the engine and remove the session's temporary output directory."""
        self._binding = None
        self._remove_scratch()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def vector(self, values=None):
        """The full parameter vector (fit basis) of ``values``: None (the session's), {name: value} over the
        session's (names the rate has not are left out: norms, shifts), or a vector of every parameter."""
        if values is None:
            return self.values.copy()
        if isinstance(values, dict):
            out = self.values.copy()
            pos = {n: i for i, n in enumerate(self.names)}
            for name, v in values.items():
                if name in pos:
                    out[pos[name]] = float(v)
            return out
        v = np.asarray(values, float)
        if v.shape != self.values.shape:
            raise ValueError(f"a vector of every rate parameter ({self.values.size} values, RateSession.names), or a "
                             f"dict by name; got {v.shape}")
        return v

    def rate(self, values=None, *, grid=None, temperatures=None, t_min=None, t_max=None, t_step=None, source=None):
        """The rate (:class:`Rate`) at ``values`` (:meth:`vector`) on the temperatures (:func:`temperatures_of`)."""
        from . import engine_text
        from .project import Project
        kind, t = temperatures_of(grid, temperatures, t_min, t_max, t_step)
        p = self.vector(values)
        if kind == "grid":
            with engine_text.held():
                pairs = self._binding.reaction_rate(p, t[0], t[1], t[2], "")
        else:
            fd, path = tempfile.mkstemp(prefix="temperatures_", suffix=".dat", dir=self._scratch)
            try:
                with os.fdopen(fd, "w") as fh:
                    fh.write("".join(f"{float(v)!r}\n" for v in t))
                with engine_text.held():
                    pairs = self._binding.reaction_rate(p, 0.0, 0.0, 0.0, path)
            finally:
                os.remove(path)
        reaction = Project.open(self.file).reaction(self.entrance, self.exit)
        return Rate(T9=np.array([a for a, _ in pairs]), rate=np.array([b for _, b in pairs]), entrance=self.entrance,
                    exit=self.exit, reaction=reaction, source=source or self.source, project=self.file)


def reaction_rate(file, entrance, exit, values=None, *, grid=None, temperatures=None, t_min=None, t_max=None,
                  t_step=None, parameter_file=None, cwd=None, options=None, coulomb_table=None):
    """N_A<sigma v> of ``entrance`` to ``exit`` (pair keys of the project ``file``; exit -1 the total radiative
    capture) at ``values`` (None: the project's or the parameter file's; {name: value}; a full vector of the rate's
    parameters), on a temperature grid (one of ``grid``, ``temperatures``, ``t_min``/``t_max``/``t_step``, GK): a
    :class:`Rate`, computed by the engine as the program's ``rate`` command computes it."""
    with RateSession(file, entrance, exit, parameter_file=parameter_file, cwd=cwd, options=options,
                     coulomb_table=coulomb_table) as s:
        return s.rate(values, grid=grid, temperatures=temperatures, t_min=t_min, t_max=t_max, t_step=t_step)


def rate_band(session, values=None, *, covariance=None, names=None, draws=None, grid=None, temperatures=None,
              t_min=None, t_max=None, t_step=None, step=1e-4, max_draws=None, seed=1, progress=None, stop=None):
    """The band of a rate (:class:`RateBand`) on a :class:`RateSession`.

    - From a covariance (``covariance``, n x n, over the parameters ``names``; a fit's ``covariance.dat`` with its
      ``#columns``, :func:`pyazr.read_covariance_dat`): ln N_A<sigma v> propagated linearly, sigma^2 = g^T C g with
      g = d ln R / d p by central differences of the engine's rate (steps of :func:`pyazr.bands.step_sizes`); the
      band is best x exp(-+sigma), the median the rate at ``values``.  Names the rate has not (norms, shifts) are
      left out: the marginal block, since the rate does not depend on them.
    - From draws (``draws``, m x n over ``names``: a sampler's chain after its burn-in, or vectors drawn from a
      covariance): the rate at each (at most ``max_draws``, chosen with ``seed``), and its 16th, 50th and 84th
      percentiles per temperature.

    ``progress(done, total)`` is called after each rate; ``stop()`` returning True ends the loop early (a band from
    the draws done)."""
    from .bands import step_sizes
    if (covariance is None) == (draws is None):
        raise ValueError("rate_band: give covariance= (with names=) or draws= (with names=)")
    if names is None:
        raise ValueError("rate_band: names= the parameters of the covariance's rows or of the draws' columns")
    names = list(names)
    pos = {n: i for i, n in enumerate(session.names)}
    keep = [j for j, n in enumerate(names) if n in pos]
    if not keep:
        raise ValueError("rate_band: none of the names is a parameter of the rate (level energies and widths)")
    sel = [names[j] for j in keep]
    base = session.vector(values)
    temps = dict(grid=grid, temperatures=temperatures, t_min=t_min, t_max=t_max, t_step=t_step)
    best = session.rate(base, **temps)
    tick = progress or (lambda done, total: None)
    if covariance is not None:
        C = np.asarray(covariance, float)[np.ix_(keep, keep)]
        cols = [pos[n] for n in sel]
        h = step_sizes(base, cols, step)
        g = np.zeros((len(cols), best.rate.size))
        for k, (c, hk) in enumerate(zip(cols, h)):
            up, down = base.copy(), base.copy()
            up[c] += hk
            down[c] -= hk
            g[k] = (np.log(session.rate(up, **temps).rate) - np.log(session.rate(down, **temps).rate)) / (2 * hk)
            tick(k + 1, len(cols))
            if stop and stop():
                raise InterruptedError("rate_band: stopped before every derivative was made")
        sigma = np.sqrt(np.einsum("it,ij,jt->t", g, C, g))
        return RateBand(T9=best.T9, low=best.rate * np.exp(-sigma), median=best.rate.copy(),
                        high=best.rate * np.exp(sigma), best=best.rate, method="covariance",
                        how=(f"linear propagation of ln N_A<sigma v> over {len(sel)} parameters (central "
                             f"differences, step {step:g} of each value); the band is the rate times exp(-+sigma)"),
                        varied=sel)
    D = np.atleast_2d(np.asarray(draws, float))
    if D.shape[1] != len(names):
        raise ValueError(f"rate_band: draws have {D.shape[1]} columns, names {len(names)}")
    if max_draws is not None and D.shape[0] > max_draws:
        D = D[np.random.default_rng(seed).choice(D.shape[0], size=int(max_draws), replace=False)]
    rates = []
    for k, row in enumerate(D):
        v = base.copy()
        for j, n in zip(keep, sel):
            v[pos[n]] = row[j]
        rates.append(session.rate(v, **temps).rate)
        tick(k + 1, D.shape[0])
        if stop and stop():
            break
    R = np.array(rates)
    lo, mid, hi = np.percentile(R, [16.0, 50.0, 84.0], axis=0)
    return RateBand(T9=best.T9, low=lo, median=mid, high=hi, best=best.rate, method="draws",
                    how=(f"the 16th, 50th and 84th percentiles of the rate over {len(R)} drawn vectors of "
                         f"{len(sel)} parameters"), varied=sel, draws=len(R))


def gamow_window(file, entrance, T9):
    """The Gamow peak E0 and its 1/e width (MeV, c.m.) of the entrance pair of the project ``file`` at the
    temperatures ``T9`` (GK), with ``E0 = (b kT / 2)^(2/3)``, ``width = 4 (E0 kT / 3)^(1/2)``, ``b = 2 pi eta(E)
    E^(1/2)`` from the engine's Sommerfeld parameter (non-relativistic, the pair's masses in u); the window
    ``E0 -+ width / 2`` is where a non-resonant charged-particle rate is made.  A neutral entrance has no Coulomb
    barrier: E0 and width are NaN and ``kT`` is the hint.  A dict of arrays (:data:`UNITS`)."""
    from . import physics
    from .project import Project
    pair = next((q for q in Project.open(file).pairs() if q.key == int(entrance)), None)
    if pair is None or pair.kind != 0:
        raise ValueError(f"gamow_window: {entrance} is not a particle pair of {file}")
    t = np.atleast_1d(np.asarray(T9, float))
    kT = BOLTZMANN_MEV_PER_GK * t
    z1, z2 = int(pair.light[0]), int(pair.heavy[0])
    if z1 * z2 == 0:
        nan = np.full_like(t, np.nan)
        return {"T9": t, "kT": kT, "E0": nan, "width": nan, "low": nan, "high": nan}
    b = 2.0 * np.pi * float(physics.sommerfeld_eta(z1, z2, pair.light_mass, pair.heavy_mass, 1.0))
    e0 = (b * kT / 2.0) ** (2.0 / 3.0)
    width = 4.0 * np.sqrt(e0 * kT / 3.0)
    return {"T9": t, "kT": kT, "E0": e0, "width": width, "low": e0 - width / 2.0, "high": e0 + width / 2.0}


def write_rate_table(path, rate, band=None, comparison=None, parameters=None):
    """Write the rate (and its band, and a comparison table interpolated to its temperatures) as CSV: a header of
    ``#`` lines naming the project, the reaction, the parameters' source and the units, then the columns T9, rate
    (and low, median, high; comparison, ratio).  Returns ``path``."""
    out = io.StringIO()
    out.write(f"# reaction rate N_A<sigma v> of {rate.reaction} (entrance pair {rate.entrance}, exit pair "
              f"{rate.exit}), computed by the AZURE2 engine\n")
    out.write(f"# project: {rate.project}\n# parameters: {parameters or rate.source}\n")
    out.write(f"# units: T9 in {UNITS['T9']}, rates in {UNITS['rate']}\n")
    cols = {"T9": rate.T9, "rate": rate.rate}
    if band is not None:
        out.write(f"# band: {band.how}\n")
        cols.update(low=band.low, median=band.median, high=band.high)
    if comparison is not None:
        out.write(f"# comparison: {comparison['path']} (interpolated in ln rate against ln T9)\n")
        ref = table_at(comparison, rate.T9)
        cols.update(comparison=ref, ratio=rate.rate / ref)
    w = csv.writer(out, lineterminator="\n")
    w.writerow(list(cols))
    for i in range(rate.T9.size):
        w.writerow([f"{cols[k][i]:.10g}" for k in cols])
    with open(path, "w") as fh:
        fh.write(out.getvalue())
    return path


def table_at(table, T9):
    """A rate table (:func:`read_rate_table`) at the temperatures ``T9``: interpolated linearly in ln rate against
    ln T9; NaN outside the table's range."""
    return np.exp(np.interp(np.log(np.asarray(T9, float)), np.log(table["T9"]), np.log(table["rate"]),
                            left=np.nan, right=np.nan))


def read_rate_table(path):
    """A rate table to compare with: two or three columns (T9 in GK, N_A<sigma v> in cm^3 mol^-1 s^-1, optionally
    its uncertainty in the same unit), separated by blanks or commas; ``#`` comments and a header line of words are
    skipped.  A dict with T9, rate, uncertainty (None without) and path."""
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.split("#", 1)[0].replace(",", " ").split()
            if not s:
                continue
            try:
                rows.append([float(v) for v in s[:3]])
            except ValueError:
                if rows:
                    raise ValueError(f"{path}: a line of words among the numbers: {line.strip()!r}") from None
    if not rows or min(len(r) for r in rows) < 2:
        raise ValueError(f"{path}: two or three columns (T9 GK, rate, optional uncertainty) are needed")
    a = np.array([r[:3] if len(r) >= 3 else r[:2] + [np.nan] for r in rows])
    order = np.argsort(a[:, 0])
    a = a[order]
    if np.any(a[:, 0] <= 0) or np.any(a[:, 1] <= 0):
        raise ValueError(f"{path}: temperatures and rates must be > 0")
    unc = None if np.all(np.isnan(a[:, 2])) else a[:, 2]
    return {"T9": a[:, 0], "rate": a[:, 1], "uncertainty": unc, "path": os.path.abspath(path)}
