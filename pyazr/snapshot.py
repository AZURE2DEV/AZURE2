"""A fit's values into its project: :func:`write_values` (a :class:`pyazr.Project` edited in memory) and
:func:`save_fit` (:meth:`pyazr.Session.save_fit`: a copy of the project on disk, read back and checked).

The project is edited through :class:`pyazr.Project`, the engine's edits of the level rows, the <thm> lines and the
records.  ``session`` (``self`` in the helpers) is the :class:`pyazr.Session` the fit comes from.
"""
import os
import shutil
import tempfile

import numpy as np

from . import session_checks
from .described import _engine as _engine_module
from .errors import ParameterError
from .locking import session_lock
from .paths import is_absolute


def write_values(session, project, x=None, norms="fitted"):
    """Write the values of the fit vector ``x`` (default: the start, ``session.parameters.x0``) of a
    :class:`pyazr.Session` ``session`` into ``project``, a :class:`pyazr.Project` of the session's file
    (``Project.open(session.file)``, or an edit of it), in memory: nothing is saved, the caller saves it (or keeps it
    as one edit of its history).  The values are those :meth:`pyazr.Session.save_fit` writes: the level parameters as
    the ``<levels>`` rows state them (energies in MeV, partial widths in eV, ANCs), fixed ones included; the THM
    coherent backgrounds; the session's potentials; with ``norms="fitted"`` every free nuisance into the record it
    comes from (``"nominal"``: the abundances only).  Returns ``project``."""
    if norms not in ("fitted", "nominal"):
        raise ParameterError(f"write_values: norms must be 'fitted' or 'nominal', not {norms!r}.")
    x = np.asarray(session.parameters.x0 if x is None else x, float).ravel()
    n = len(session.parameters.free)
    if x.size != n:
        raise ParameterError(f"write_values: x has {x.size} entries; the session has {n} free parameters "
                             f"(Session.parameters.x0 is such a vector, Session.parameters.free names its entries).")
    full = session._all_rwa(x)
    want = session.transform(x)
    place_values(session, project, x, want, full, norms)
    return project


def save_fit(session, path, x=None, param_sav=True, verify=True, norms="fitted",
             close_session=False, relative_to=None):
    """Snapshot a fit as a ``.azr`` you can reopen, plot, or hand over.

    ``path`` is required and is written to explicitly -- nothing is ever
    saved in place over the model you loaded.

    ``x`` is a fit vector (``parameters.x0`` is one); it defaults to the
    session's start. The conversion to the physical values a
    ``<levels>`` line holds (partial widths in eV, ANCs in fm^-1/2) is done
    here, so a caller cannot forget it.

    **Fixed parameters are written too.** A width marked fixed holds its
    reduced-width amplitude during a fit -- in pyazr as in the CLI's MIGRAD
    and LM -- so its *physical* value moves whenever the level's other
    widths do (the Brune denominator ``1 + sum_c gamma_c^2 dS_c/dE``
    contains them all).  The snapshot carries the physical value the fit
    actually had, so the fixed amplitude reads back unchanged.  (Before,
    the file kept the fixed width as it was read: a fit that tripled a
    level's free width reloaded as a different model, 4 % off in chi2 on
    tests/15N_p_a.)

    Two things a ``<levels>`` block cannot carry, which is why this is not
    just a write of the level values into ``<levels>``:

    - **Normalizations, energy shifts and the other nuisances are not in
      it.** They are values of the data-set and effect records (their
      nuisances: a start value, ``fit`` or ``fixed``, and a prior that
      states its centre).  ``norms`` chooses:

      ``"fitted"`` (default)
          every free parameter that is not a level parameter (norms,
          shift terms, shift-function parameters, abundances, ...) is
          written into the record it comes from, in that record's unit;
          the priors keep the centres they state.  The snapshot run on
          its own gives the fitted data chi-squared *and* the fit's
          priors (both checked by ``verify``).
      ``"nominal"``
          the records keep the loaded file's values (but the abundances,
          written at the fit's values); the fitted norms and shifts are
          in the ``.sav`` only (needs ``param_sav``).  Run with that
          ``.sav`` as the external parameter file (CLI: "use previous
          parameters") to reproduce the fit.

      ``param_sav`` (the default) writes the companion ``<name>.sav``
      beside the file, carrying every parameter, free and fixed.
    - **A written file is not a fit until it reads back as one.** With
      ``verify`` (the default) the result is reopened and checked two
      ways: its physical values transformed back must match, and the
      reduced-width amplitudes it reads for the particle and beta-decay
      channels of every level must be ``x``'s, up to a common sign of
      the level (the standard parametrization's eigenvector convention,
      no change of model). If either disagrees the files are removed and
      this raises, rather than leaving a snapshot that is quietly a
      mixture. The second check exists because the first cannot see a
      value two parameters map onto: the engine returned ``|B|`` for a
      beta-decay amplitude until 2026-10-12, and a snapshot with the
      sign lost passed. See :meth:`_snapshot_mismatch` for why photon
      channels are judged on physical values only.  With
      ``norms="fitted"`` (a session of the data) the reopened file's data chi-squared
      and its objective (:meth:`evaluate`: ``chi2``, ``objective``) must also equal the fit's
      (rel. 1e-6; 1e-4 in the standard parametrization, whose transform is
      inverted iteratively on reading).

    A ``.azr`` names its data files *relative to itself*: a snapshot written to
    another directory gets its paths rewritten relative to that directory
    (:func:`rebase_paths`), so it reopens there and reads the original's data;
    its output and checks directories are its own (relative to it).
    ``relative_to`` (a directory) names the paths relative to it instead,
    for a text opened with that working directory (the web page's server
    keeps the project's text so); the check then runs there.

    ``close_session`` closes this session once the files are written and
    what the check compares with has been computed, *before* the
    snapshot is reopened: one engine in memory instead of two (a 19F THM
    session is about 220 MB, 390 MB with the verifying one beside it).
    The session is unusable afterwards (:meth:`is_alive` is False); use
    it as the last call of a fit, as ``scripts/thm_model_average.py``
    does.

    Returns ``(azr_path, sav_path_or_None)``.
    """
    if norms not in ("fitted", "nominal"):
        raise ParameterError(f"save_fit: norms must be 'fitted' or 'nominal', not {norms!r} (see "
                             f"help(pyazr.Session.save_fit)).", what=f"save_fit: norms={norms!r}",
                             allowed="norms must be 'fitted' or 'nominal'", fix="help(pyazr.Session.save_fit)")
    if norms == "nominal" and not param_sav:
        raise ParameterError("save_fit: norms='nominal' keeps the fitted norms in the .sav "
                             "only; it needs param_sav=True.", what="save_fit: norms='nominal' with param_sav=False",
                             allowed="norms='nominal' keeps the fitted norms in the .sav only",
                             fix="it needs param_sav=True")
    x = np.asarray(session.parameters.x0 if x is None else x, float).ravel()
    free = [i for i, p in enumerate(session.parameters) if not p.fixed]
    if x.size != len(free):
        raise ParameterError(f"save_fit: x has {x.size} entries; the session has {len(free)} free parameters "
                             f"(Session.parameters.x0 is such a vector, Session.parameters.free names its entries).",
                             what=f"save_fit: x has {x.size} entries",
                             allowed=f"the session has {len(free)} free parameters",
                             fix="Session.parameters.x0 is such a vector, Session.parameters.free names its entries")
    # Every parameter, with the free ones from x, and every parameter's
    # physical value, fixed ones included (see above).
    full = session._all_rwa(x)
    want = session.transform(x)

    path = os.path.abspath(path)
    home = os.path.abspath(relative_to) if relative_to is not None else os.path.dirname(path)
    write_format3(session, path, x, want, full, norms, home)

    sav = None
    if param_sav:
        sav = session.parameters.write(os.path.splitext(path)[0] + ".sav", x)

    # What the check compares with, from this session before it may close.
    chi_want = pen_want = None
    if verify and norms == "fitted" and session.mode == "data":
        ev = session.evaluate(x)
        chi_want, pen_want = ev.chi2, ev.objective
    check_options = {k: v for k, v in session.options.items() if k not in ("extrapolations", "parameter_file")}
    check_options["extrapolations"] = session.mode == "extrap"
    if close_session:
        session.close()

    if verify:
        # Verify a *copy*, in the snapshot's own directory (its paths are
        # relative to it) and with its output redirected, so that opening it
        # does not overwrite a run's output/param.par there.
        probe = tmpdir = None
        try:
            tmpdir = tempfile.mkdtemp(prefix=".azr_verify_", dir=home)
            fd, probe = tempfile.mkstemp(suffix=".azr", dir=home)
            os.close(fd)
            from .project import Project
            from .session import Session
            copy = Project.open(path)
            copy.set_config(output_dir=tmpdir.rstrip("/") + "/")
            copy.save(probe)
            with Session(probe, cwd=home, **check_options) as check, session_lock(check):
                back_full = np.asarray(check._engine.params_all_rwa(), float)
                got = check.transform(check.parameters.x0)
                back = {p.name: float(back_full[i])
                        for i, p in enumerate(check.parameters)}
                chi_got = None
                if chi_want is not None:
                    ev = check.evaluate(check.parameters.x0)
                    chi_got, pen_got = ev.chi2, ev.objective
        except Exception as err:
            _discard(path, sav)
            raise RuntimeError(
                f"the snapshot was written but would not load back "
                f"({err}); it has been removed.") from err
        finally:
            if probe and os.path.exists(probe):
                os.remove(probe)
            if tmpdir and os.path.isdir(tmpdir):
                shutil.rmtree(tmpdir, ignore_errors=True)
        bad = _snapshot_mismatch(session, full, want, got, back)
        if bad:
            _discard(path, sav)
            raise RuntimeError(
                "the snapshot does not read back as the fit it was made "
                f"from ({'; '.join(bad)}); it has been removed. This means "
                "the .azr and the parameter set do not describe the same "
                "model.")
        if chi_got is not None:
            # The standard parametrization (basis 0) inverts its transform
            # iteratively on reading (energies to 1e-6 MeV), which leaves
            # ~1e-5 in the chi-squared of a steep spectrum
            # (tests/pyazr/beta_sign_test.py); Brune and Park read back
            # exactly.
            rtol = 1e-4 if session.basis == 0 else 1e-6
            if not np.isclose(chi_got, chi_want, rtol=rtol, atol=1e-9):
                _discard(path, sav)
                raise session_checks.snapshot_mismatch("chi2", chi_got, chi_want, rtol)
            if not np.isclose(pen_got, pen_want, rtol=rtol, atol=1e-9):
                _discard(path, sav)
                raise session_checks.snapshot_mismatch("objective", pen_got, pen_want, rtol)
    return path, sav

def _discard(*paths):
    """Remove a snapshot that failed verification, so nothing bad survives."""
    for f in paths:
        if f and os.path.exists(f):
            os.remove(f)

def _snapshot_mismatch(session, x, want, got, back):
    """Why a re-read snapshot is not the fit it was made from; ``[]`` if it is.

    ``x``, ``want``, ``got`` and ``back`` cover every parameter, fixed ones
    included (``params_all_rwa`` order; ``back`` by name).  Two
    comparisons, because neither alone is enough:

    - the **physical values** (:meth:`transform` of both sides), which
      is what the file holds.  A transform that maps two parameters onto
      one value passes this alone: the engine returned ``|B|`` for a
      beta-decay amplitude until 2026-10-12, and a snapshot with the
      sign lost read back as a different model (12C+alpha).
    - the **re-read reduced-width amplitudes** of the particle and
      beta-decay channels, and the level energies, against ``x`` --
      level by level, allowing a common sign of the whole level, which
      the standard parametrization's eigenvector convention may flip and
      which is no change of model.

    Photon channels are left to the first comparison on purpose: their
    amplitude is derived from the written width together with the
    external width of the level, a sum over the level's particle
    amplitudes, so it is only as exact as the written widths are, and the
    physical comparison is the direct one.
    """
    if got.shape != want.shape:
        return [f"{got.size} physical parameters instead of {want.size}"]
    bad = []
    levels = {}
    for i, p in enumerate(session.parameters):
        if p.kind not in ("energy", "width") or i >= want.size:
            continue
        if p.name not in back:
            bad.append(f"{p.name} is not in the re-read file")
            continue
        levels.setdefault((p.jgroup, p.level), []).append(
            (p, float(want[i]), float(got[i]), float(x[i]), back[p.name]))

    def same(pairs, rtol):
        """Equal, or equal with the level's common sign flipped."""
        if not pairs:
            return True
        a = np.array([v[0] for v in pairs])
        b = np.array([v[1] for v in pairs])
        return bool(np.allclose(a, b, rtol=rtol, atol=1e-12)
                    or np.allclose(a, -b, rtol=rtol, atol=1e-12))

    for (jg, lv), items in sorted(levels.items()):
        ok = all(abs(w - g) <= 1e-5 * max(1.0, abs(w))
                 for p, w, g, a, b in items if p.kind == "energy")
        # the written physical widths, every channel (sign carried)
        ok = ok and same([(w, g) for p, w, g, a, b in items
                          if p.kind == "width"], 1e-4)
        # the re-read amplitudes: particle and beta channels only, since
        # a photon channel's amplitude depends on the external-capture
        # integrals the re-read computed afresh
        ok = ok and same([(a, b) for p, w, g, a, b in items
                          if p.kind == "width"
                          and p.radiation_type in ("P", "F", "G")], 1e-3)
        ok = ok and all(abs(a - b) <= 1e-5 * max(1.0, abs(a))
                        for p, w, g, a, b in items if p.kind == "energy")
        if not ok:
            p0 = items[0][0]
            bad.append(f"level {p0.jpi} #{lv} (J-group {jg}, E = {p0.level_energy}) "
                       "re-reads to different parameters")
    return bad


def place_values(self, project, x, want, full, norms):
    """The fitted level values into the <levels> rows they come from (the engine's (group, level) and (pair key, l,
    2s) of each parameter, ``_azure2.level_value_edits``), the THM coherent backgrounds into their experiments, and
    the free nuisances (with ``norms="fitted"``; else the abundances) into the records they come from: edits of
    ``project`` in memory."""
    place_level_values(project, self.parameters, self.pairs, want)
    write_cbkg(self, project, np.asarray(full, float))
    project._write_parameter_values(
        {p.name: float(want[p.index]) for p in self.parameters.of_kind("abundance")} if norms != "fitted"
        else {p.name: float(x[p.free_index]) for p in self.parameters if p.free_index is not None})


def write_format3(self, path, x, want, full, norms, home):
    """save_fit of a format-3 project: :func:`place_values` on the session's project, written to ``path``."""
    from .project import Project
    project = Project.open(self.file)
    place_values(self, project, x, want, full, norms)
    rebase_paths(project, path, self.cwd, home)
    project.save(path)


def rebase_paths(project, path, base, home):
    """Make ``project`` read the same files from the directory ``home``: every relative path it names -- the record
    fields of type ``path`` in the engine's schema (a data set's file, covariance, THM weight table), the ``<config>``
    keys of kind path (``priorCovariance``) and the ``table:<file>`` values of ``<thm>`` experiments -- rewritten
    relative to the new directory; absolute paths and the output and checks directories (the snapshot's own) are
    kept.  ``base`` is the directory the paths are read from now (the session's working directory), ``home`` the one
    they are read from in the snapshot (its own directory, or save_fit's ``relative_to``).  The edits are validated
    as from ``home``; :meth:`Project.save` then writes the text to ``path``."""
    import os
    from .described import _engine
    old, new = os.path.abspath(base), os.path.abspath(home)
    project.path = os.path.join(new, os.path.basename(path))
    if os.path.normpath(old) == os.path.normpath(new):
        return
    refuse_moved_covariance(project, old, new)

    def moved(value):
        return value if is_absolute(value) else os.path.relpath(os.path.join(old, value), new)

    from .schema import schema
    fields = {k: {f for f, d in r["fields"].items() if d.get("type") == "path"} for k, r in schema()["records"].items()}
    for record in project.records().values():
        for f in fields.get(record.record_kind, ()):
            value = record.get(f)
            if isinstance(value, str) and value:
                project.set(record.record_name, **{f: moved(value)})
    _, keys = _engine().config_keys()
    config = project.config()
    for key, kind, _doc in keys:
        if kind == "path" and config.get(key):
            project.set_config(**{key: moved(config[key])})
    for name, entries in ((project.thm() or {}).get("experiments") or {}).items():
        for key, value in entries.items():
            if "table:" in value:
                head, file = value.split("table:", 1)
                project.set_thm_experiment(name, **{key: f"{head}table:{moved(file)}"})

def refuse_moved_covariance(project, old, new):
    """Refuse (before anything is rewritten) a snapshot read from ``new`` when a data set's covariance file names its
    data files by relative paths: the engine reads those names from the project's working directory, so from ``new``
    they name other files, and the snapshot would need a second copy of the covariance file.  The names come from the
    engine's reader of the file (``_azure2.covariance_file_names``).  A covariance file is the only file of a project
    that names other files (a prior covariance names parameters; THM, reference and weight tables hold numbers)."""
    import os
    from .described import _engine
    from .errors import ParameterError
    for record in project.records().values():
        file = record.get("covariance") if record.record_kind == "dataset" else None
        if not isinstance(file, str) or not file:
            continue
        names = _engine().covariance_file_names(file, old)
        relative = [n for n in names if not is_absolute(n)]
        if not relative:
            continue
        what = (f"save_fit to {new}: data set {record.record_name!r} has the covariance file {file}, which names its "
                f"data files relative to the project's directory ({', '.join(relative[:3])}"
                f"{', ...' if len(relative) > 3 else ''}); read from {new} they name other files")
        allowed = (f"a snapshot whose paths are read from the project's own directory {old}: a path there, or "
                   f"relative_to={old!r} for a text the caller reads with that directory (as the web page's server)")
        fix = f"call save_fit with a path in {old}"
        raise ParameterError(f"{what}.  Allowed: {allowed}.  To go on: {fix}.", what=what, allowed=allowed, fix=fix)


def place_level_values(project, parameters, pairs, want):
    """Write the level parameters of ``parameters`` (:class:`pyazr.Parameter`) at the values ``want`` (indexed by
    ``Parameter.index``; partial widths in eV, ANCs, energies in MeV) into the rows of ``project`` they come from:
    a level by the engine's (group, level), a channel by its pair KEY (the file's; ``Parameter.pair`` is the engine's
    pair number, mapped through ``pairs``), l and 2s.  A fixed width of a channel the level has no row for is the
    engine's zero input (a level of a J-group gets a fixed, zero width for every group channel it does not list):
    there is nothing to write.  Any other value without a row raises ValueError naming it, so that a snapshot is never
    a mixture of two fits."""
    pairkey = {p.number: p.key for p in pairs}
    level_params = [p for p in parameters if p.kind in ("energy", "width")]
    values = [(p.kind, p.jgroup, p.level, pairkey.get(p.pair, p.pair) if p.kind == "width" else 0,
               p.L if p.kind == "width" else 0, int(round(2 * (p.S or 0.0))) if p.kind == "width" else 0,
               float(want[p.index])) for p in level_params]
    unplaced = [f"{level_params[i].name} ({why})" for i, why in project._write_level_values(values)
                if not level_params[i].fixed]
    if unplaced:
        raise ValueError("save_fit could not place " + str(len(unplaced)) + " level parameter(s), so the result "
                         "would be a mixture of two fits:\n  " + "\n  ".join(unplaced[:8]))


def write_cbkg(self, model, values):
    """Write the THM coherent-background parameters (``cbkg_*``) back into
    the ``cbackground=`` values of their experiments, one explicit term
    per (entrance, exit) combination with its values (``f`` for fixed),
    so the snapshot reads back as the fit.  ``values`` is the full
    parameter vector in the session's order (the cbkg entries are not
    transformed).  Each parameter's experiment, combination and coefficient
    come from its record (:attr:`Parameter.thm_experiment` ...)."""
    cb = [(i, p) for i, p in enumerate(self.parameters) if p.kind == "cbkg"]
    if not cb:
        return
    names = list((model.thm() or {"experiments": {}})["experiments"])
    for name in names:
        combos = {}
        for i, p in cb:
            if p.thm_experiment != name:
                continue
            key = (p.thm_jpi, str(p.thm_exit_key), p.thm_channels)
            combos.setdefault(key, {})[p.thm_term] = (float(values[i]),
                                                      bool(self._fixed[i]))
        if not combos:
            continue
        terms = []
        for (jpi, exit_key, chans), v in combos.items():
            n = 4 if 2 in v or 3 in v else 2
            vals = ",".join(repr(v.get(k, (0.0, False))[0])
                            + ("f" if v.get(k, (0.0, False))[1] else "")
                            for k in range(n))
            terms.append(f"{jpi}:{exit_key}:{chans}"
                         + (":linear" if n == 4 else "") + "=" + vals)
        # the engine's canonical text of the terms, written by the engine into the experiment's line
        model.set_thm_experiment(name, cbackground=_engine_module().parse_thm_cbackground(";".join(terms))[1])
