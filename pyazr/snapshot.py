"""The file writers of :meth:`pyazr.azure2.save_fit`: a fit's values into a copy of its project.

A format-3 project (:func:`write_format3`) is edited through :class:`pyazr.Project`, the engine's edits of the level
rows, the <thm> lines and the records; an older project (:func:`write_old_format`) through ``AzrModel``, deleted with
the old formats.  ``self`` is the :class:`pyazr.azure2` session the fit comes from.
"""
import numpy as np

from .described import _engine as _engine_module


def write_format3(self, path, x, want, full, norms):
    """save_fit of a format-3 project: the fitted level values into the <levels> rows they come from (the
    engine's (group, level) and (pair key, l, 2s) of each parameter, ``_azure2.level_value_edits``), the THM
    coherent backgrounds into their experiments, the session's potentials as records, and the free nuisances
    (with ``norms="fitted"``; else the abundances) into the records they come from; then written to ``path``."""
    from .project import Project
    project = Project.open(self.file)
    pairkey = {p.number: p.key for p in self.pairs}
    level_params = [p for p in self.parameters if p.kind in ("energy", "width")]
    values = [(p.kind, p.jgroup, p.level, pairkey.get(p.pair, p.pair) if p.kind == "width" else 0,
               p.L if p.kind == "width" else 0, int(round(2 * (p.S or 0.0))) if p.kind == "width" else 0,
               float(want[p.index])) for p in level_params]
    # A fixed width of a channel the level has no row for is the engine's zero input (a level of a J-group gets a
    # fixed, zero width for every group channel it does not list): there is nothing to write.
    unplaced = [f"{level_params[i].name} ({why})" for i, why in project._write_level_values(values)
                if not level_params[i].fixed]
    if unplaced:
        raise ValueError("save_fit could not place " + str(len(unplaced)) + " level parameter(s), so the result "
                         "would be a mixture of two fits:\n  " + "\n  ".join(unplaced[:8]))
    write_cbkg(self, project, np.asarray(full, float))
    self._write_potential_records(project)
    project._write_parameter_values(
        {p.name: float(want[p.index]) for p in self.parameters.abundances} if norms != "fitted"
        else {p.name: float(x[p.free_index]) for p in self.parameters if p.free_index is not None})
    project.save(path)

def write_old_format(self, path, x, want, full, norms):
    """save_fit of a format-1 or format-2 project (AzrModel; deleted with the old formats)."""
    from .azrfile import AzrModel
    model = AzrModel.from_file(self.file)
    model.apply_fit(self.parameters, want, physical=True, pairs=self.pairs,
                    include_fixed=True)
    for p in self.parameters.abundances:  # an <abundances> block: the fitted values
        model.set_abundance((p.nucleus_Z, p.nucleus_A), value=float(want[p.index]), sample=p.sample or "")
    write_cbkg(self, model, np.asarray(full, float))
    if norms == "fitted":
        values, centres = {}, {}
        datasets = self.datasets
        for p in self.parameters:
            if p.kind not in ("norm", "shift") or p.fixed or p.free_index is None:
                continue
            n, sh = values.get(p.segment_key, (None, None))
            v = float(x[p.free_index])
            values[p.segment_key] = (v, sh) if p.kind == "norm" else (n, v)
            # Keep the prior where the fit had it: an explicit centre for
            # every prior whose field now holds another value.
            d = datasets.by_key(p.segment_key)
            cn, cs = centres.get(p.segment_key, (None, None))
            if p.kind == "norm" and d.norm_error and v != d.nominal_norm:
                cn = d.nominal_norm
            elif (p.kind == "shift" and d.vary_shift and d.energy_shift_error
                  and d.nominal_shift != 0.0):
                # an explicit shift centre (the default, 0, needs no row)
                cs = d.nominal_shift
            if cn is not None or cs is not None:
                centres[p.segment_key] = (cn, cs)
        if values:
            model.set_segment_values(values)
        if centres:
            model.set_prior_centres(centres)
    # Electron-screening potentials are not written into the <targetInt> lines (their value is the prior's
    # centre; the fitted U_e travels in the .sav).
    self._write_potential_records(model).write(path)

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
    from .project import Project
    format3 = isinstance(model, Project)
    names = list((model.thm() or {"experiments": {}})["experiments"]) if format3 else model.thm_experiments()
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
        if format3:  # the engine's canonical text of the terms, written by the engine into the experiment's line
            model.set_thm_experiment(name, cbackground=_engine_module().parse_thm_cbackground(";".join(terms))[1])
        else:
            model.set_thm_cbackground(name, terms)
