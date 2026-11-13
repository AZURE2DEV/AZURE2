"""A fit's values into its project: :func:`write_values` (a :class:`pyazr.Project` edited in memory)
and the file writer of :meth:`pyazr.azure2.save_fit` (:func:`write_format3`, a copy of the project on disk).

The project is edited through :class:`pyazr.Project`, the engine's edits of the level rows, the <thm> lines and the
records.  ``self`` is the :class:`pyazr.azure2` session the fit comes from.
"""
import numpy as np

from .described import _engine as _engine_module
from .errors import ParameterError


def write_values(session, project, x=None, norms="fitted"):
    """Write the values of the fit vector ``x`` (default: the start, ``session.params_rwa``) of a
    :class:`pyazr.azure2` ``session`` into ``project``, a :class:`pyazr.Project` of the session's file
    (``Project.open(session.file)``, or an edit of it), in memory: nothing is saved, the caller saves it (or keeps it
    as one edit of its history).  The values are those :meth:`pyazr.azure2.save_fit` writes: the level parameters as
    the ``<levels>`` rows state them (energies in MeV, partial widths in eV, ANCs), fixed ones included; the THM
    coherent backgrounds; the session's potentials; with ``norms="fitted"`` every free nuisance into the record it
    comes from (``"nominal"``: the abundances only).  Returns ``project``."""
    if norms not in ("fitted", "nominal"):
        raise ParameterError(f"write_values: norms must be 'fitted' or 'nominal', not {norms!r}.")
    x = np.asarray(session.params_rwa if x is None else x, float).ravel()
    n = len(session.parameters.free)
    if x.size != n:
        raise ParameterError(f"write_values: x has {x.size} entries; the session has {n} free parameters "
                             f"(azure2.params_rwa is such a vector, azure2.parameters.free names its entries).")
    full = session._all_rwa(x)
    want = np.asarray(session.transform_all_rwa(full, include_fixed=True), float)
    place_values(session, project, x, want, full, norms)
    return project


def place_values(self, project, x, want, full, norms):
    """The fitted level values into the <levels> rows they come from (the engine's (group, level) and (pair key, l,
    2s) of each parameter, ``_azure2.level_value_edits``), the THM coherent backgrounds into their experiments, the
    session's potentials as records, and the free nuisances (with ``norms="fitted"``; else the abundances) into the
    records they come from: edits of ``project`` in memory."""
    place_level_values(project, self.parameters, self.pairs, want)
    write_cbkg(self, project, np.asarray(full, float))
    self._write_potential_records(project)
    project._write_parameter_values(
        {p.name: float(want[p.index]) for p in self.parameters.abundances} if norms != "fitted"
        else {p.name: float(x[p.free_index]) for p in self.parameters if p.free_index is not None})


def write_format3(self, path, x, want, full, norms):
    """save_fit of a format-3 project: :func:`place_values` on the session's project, written to ``path``."""
    from .project import Project
    project = Project.open(self.file)
    place_values(self, project, x, want, full, norms)
    project.save(path)

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
                                                      bool(self.fixed_params[i]))
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
