"""The THM part of the web UI's server (docs/upgrade/web_ui.md section 6.7): the keys of the ``<thm>`` section as
record kinds of a schema (``GET /api/thm/keys``), and the page's own part of the session's THM diagnostics.

The keys are the engine's table (``src/ThmKeys.cpp``, served by ``pyazr.kinds()["thm"]``): the one the parser reads,
each key with its value type of the schema, choices, unit and doc, so that the page draws the forms with its form
generator.  A new key of ``<thm>`` is one row there and nothing here.

The diagnostics are rows of pyazr (``Session.describe()["functions"]``); :data:`PAGE` holds what is the page's own
(the THM tab draws them, the words of their status codes, the bound of one request, the fields left off the first
plot) and the two answers that are tables, not functions of an energy (:func:`_experiments`, :func:`_background`).
The answers on an energy axis are walked by :func:`pyazr.ui.physics.answer_parts` from the rows' ``path`` and
``lines``.
"""
import copy

TAB = "thm"
GROUP = "THM experiments"


def keys():
    """The ``<thm>`` keys as record kinds of a schema (``GET /api/thm/keys``): ``pyazr.kinds()["thm"]``."""
    import pyazr
    return copy.deepcopy(pyazr.kinds()["thm"])


def _experiments(row, out):
    """One table row per experiment: its data sets, background, points, chi-squared, status, the profiled norm and
    background coefficients, each with its sigma (``sigma`` is over (norm, b0, b1, b2)); the coherent background's
    parameters in a second table."""
    rows, cbkg = [], []
    for name, r in out.items():
        line = {"experiment": name, "segments": " ".join(map(str, r["segments"])), "background": r["background"],
                "points": r["points"], "chi2": r["chi2"], "status": r["status"], "norm": r["norm"],
                "norm_sigma": r["sigma"][0]}
        for k in range(3):
            line[f"b{k}"] = r["b"][k] if k < len(r["b"]) else None
            line[f"b{k}_sigma"] = r["sigma"][k + 1] if k + 1 < len(r["sigma"]) else None
        rows.append(line)
        cbkg += [{"experiment": name, "parameter": p, "value": v} for p, v in (r.get("cbkg") or {}).items()]
    columns = [{"key": k, "title": k, "unit": ""} for k in (
        "experiment", "segments", "background", "points", "chi2", "status", "norm", "norm_sigma", "b0", "b0_sigma",
        "b1", "b1_sigma", "b2", "b2_sigma")]
    tables = [{"title": "experiments", "columns": columns, "rows": rows}]
    if cbkg:
        tables.append({"title": "coherent background (fit parameters)", "rows": cbkg, "columns": [
            {"key": k, "title": k, "unit": ""} for k in ("experiment", "parameter", "value")]})
    return {"tables": tables}


def _background(row, out):
    """The profile (norm, b0, ...: value and sigma) and its covariance as tables; the data sets as a fact."""
    names = ["norm"] + [f"b{k}" for k in range(len(out["b"]))]
    values = [out["norm"]] + list(out["b"])
    cov = out["cov"]
    rows = [{"parameter": n, "value": v, "sigma": out["sigma"][k]} for k, (n, v) in enumerate(zip(names, values))]
    covariance = [dict({"parameter": n}, **{m: cov[k][j] for j, m in enumerate(names)}) for k, n in enumerate(names)]
    head = [{"key": "parameter", "title": "parameter", "unit": ""}]
    segments = next(r for r in row["returns"] if r["name"] == "segments")
    return {"facts": [{"name": "segments", "unit": "", "doc": segments["doc"],
                       "value": " ".join(map(str, out["segments"]))}], "tables": [
        {"title": "profile", "columns": head + [{"key": c, "title": c, "unit": ""} for c in ("value", "sigma")],
         "rows": rows},
        {"title": "covariance (from the profile)", "columns": head + [{"key": n, "title": n, "unit": ""}
                                                                      for n in names], "rows": covariance}]}


# The distortion factor costs about 17 ms per energy with an optical potential (0.3 s for 20 energies, 3.4 s for 200 on
# examples/o18_lacognata2010 without its DW vertex, one thread): 1000 energies keep a request near 20 s.
DISTORTION_ELEMENTS = 1000
NOT_REACHED = ("reached", -3, "beyond the reach of the spectator window or the DW vertex grid: the engine uses the "
                              "nearest data point's vertex there, not the vertex at E")

# The page's part of the THM diagnostics (pyazr.ui.physics.PAGE explains the fields).
_THM = {"tab": TAB, "group": GROUP, "axis_field": "E"}
PAGE = {
    "thm_experiments": dict(_THM, tables=_experiments),
    "thm_background": dict(_THM, tables=_background),
    "thm_lineshape": dict(_THM, hidden=("zeta", "eta_sb", "E_sF", "eta_0")),
    "thm_vertex": dict(_THM, status=NOT_REACHED, hidden=("M2_qf", "dw_q_delta", "dw_p_delta")),
    "thm_distortion": dict(_THM, max_elements=DISTORTION_ELEMENTS,
                           hidden=("M2", "M2_PW", "E_sF", "k_sF", "eta_sF", "theta_cm", "x", "q", "lmax")),
}
