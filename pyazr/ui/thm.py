"""The THM part of the web UI's server (docs/upgrade/web_ui.md section 6.7): the keys of the ``<thm>`` section as
record kinds of a schema (``GET /api/thm/keys``), and the session's THM diagnostics as rows of the shape of the
physics table (``GET /api/session/functions``, ``POST /api/session/functions/{fn}``).

The engine reads the ``<thm>`` keys with a grammar of its own (``ParseThmOptionLine``, ``ParseThmExperimentLine`` in
``src/ThmExperiment.cpp``) and describes none of them as data: :data:`KEYS` is written here from that grammar, in the
shape of a record kind of ``pyazr.schema()`` (a field per key: value type of the schema, choices, unit in the title,
doc), so that the page draws the forms with its form generator.  A key the table lacks is still shown and edited (as
text) by the page; the engine validates every edit.  ``tests/pyazr/ui/thm_routes_test.py`` compares the experiment
keys with the list the engine's own refusal names, so a key added to the engine fails the test until it has a row.

The session's THM calls return nested dictionaries (per channel, level or exit pair); :data:`THM_ROWS` describes
each as a row (arguments with unit, kind, domain, example and sweep; returned fields with units), and
:data:`ANSWERS` turns the adapter's plain answer (:class:`pyazr.ui.engine.Engine`) into the answer of a physics
function: fields on the energies (a field over several channels or levels has ``lines``, one label per line),
``facts`` (the call's scalars with units), ``tables`` (the answers that are not functions of an energy), and the
status per energy.  The page draws every row with the Physics tab's code (``physicsform.js``, ``functionview.js``).
"""
import numpy as np

from .physics import ARG, MAX_ELEMENTS

TAB = "thm"
GROUP = "THM experiments"

# -- the keys of <thm> ------------------------------------------------------------------------------------------------


def _key(type_, doc, *, title=None, values=None, default=None, grammar=None, applies_when=None, unit=None):
    field = {"type": type_, "doc": doc, "required": False}
    if title or unit:
        field["title"] = title or ""
    if unit:
        field["unit"] = unit
    if values is not None:
        field["values"] = list(values)
    if default is not None:
        field["default"] = default
    if grammar:
        field["grammar"] = grammar
    if applies_when:
        field["applies_when"] = applies_when
    return field


def _present(key):
    return {"if": "field_present", "field": key, "text": f"{key} is given"}


# The global keys (every line of <thm> but experiment[...]) and the keys of an experiment[<name>] line, each a record
# kind of a schema: {title, doc, fields}; a field's title carries its unit.  The key spectatorEnergy[<pair key>] (one
# pair's spectator energy) is written with its pair key: the page adds it as a key the table lacks.
KEYS = {
    "<thm>": {
        "title": "<thm> options", "doc": "The model options of every THM experiment of the project.",
        "fields": {
            "vertex": _key("enum", "boundary in the entrance vertex (B - 1) j_l - rho j_l': perlevel S_c(E_lambda) of "
                           "each level under Brune (alias real); constant B_c = S_c(E_1) of the lowest level of the J "
                           "group; onshell L_c(E) = S_c(E) + i P_c(E)",
                           values=("perlevel", "constant", "onshell"), default="constant"),
            "kinematics": _key("enum", "the kinematic factors the model carries, by how the HOES data were extracted: "
                               "lacognata exit k_f/mu_f (La Cognata's working formula); triple exit Gamma_f alone; "
                               "kf3body exit Gamma_f/(mu_f k_f); lambda32 exit Gamma_f, entrance 1/k_i on shell",
                               values=("lacognata", "triple", "kf3body", "lambda32"), default="lacognata"),
            "entranceL": _key("enum", "entrance partial waves of one channel spin summed coherently or incoherently "
                              "(for an observable integrated over the exit direction the l cross terms vanish)",
                              values=("incoherent", "coherent"), default="incoherent"),
            "coulombIntegral": _key("enum", "add the external Coulomb term 2 eta k int_a^inf O_l(kr)/O_l(ka) "
                                    "j_l(pr) dr to the vertex (Tribble eq. 2.79)", values=("0", "1"), default="0"),
            "spectatorEnergy": _key("number", "mean spectator kinetic energy <p_sx^2>/2mu_sx added to E + B in the "
                                    "half-off-shell momentum, for every THM pair (one pair: the key "
                                    "spectatorEnergy[<pair key>])", title="spectatorEnergy (MeV)", unit="MeV",
                                    default="0"),
        },
    },
    "experiment": {
        "title": "experiment", "doc": "One THM experiment: the data sets that name it share one profiled norm and "
                                      "its background; its kinematics, momentum window and distortion.",
        "fields": {
            "background": _key("enum", "the incoherent background b(E) added to the folded model, profiled with the "
                               "norm", values=("none", "const", "linear", "quadratic"), default="none"),
            "beam": _key("text", "the beam nucleus (a built-in name or Z,A,mass)", grammar="<nucleus>"),
            "target": _key("text", "the target nucleus, the Trojan horse a = x + s (a built-in name or Z,A,mass)",
                           grammar="<nucleus>"),
            "spectator": _key("text", "the spectator s (a built-in name or Z,A,mass)", grammar="<nucleus>"),
            "Ebeam": _key("number", "the lab beam energy, > 0", title="Ebeam (MeV)", unit="MeV"),
            "lineshape": _key("enum", "the Coulomb line shape N_C of the spectator (needs beam, target, spectator, "
                              "Ebeam)", values=("off", "on"), default="off"),
            "ps": _key("text", "the spectator-momentum window", grammar="delta | hulthen:pmin-pmax | "
                       "hulthen:a,b:pmin-pmax | gauss:FWHM:pmin-pmax | table:<file>"),
            "psNodes": _key("integer", "Gauss-Legendre nodes of the ps window, 1 to 64", applies_when=_present("ps")),
            "distortion": _key("text", "the distortion factor R(E) multiplying the model",
                               grammar="none | coulomb | optical | table:<file>"),
            "opticalAA": _key("text", "the a + A optical potential (distortion=optical or vertexModel=dw)",
                              grammar="plane | <global set>[:extrapolate] | <parameters>"),
            "opticalSF": _key("text", "the s + F optical potential (distortion=optical or vertexModel=dw)",
                              grammar="plane | <global set>[:extrapolate] | <parameters>"),
            "spectatorAngle": _key("text", "the spectator direction of the distortion factor and the DW vertex",
                                   grammar="qf | <lab degrees 0-180> | cm:<degrees 0-180>"),
            "distortionRef": _key("number", "the reference energy E_ref of R(E) = rho(E)/rho(E_ref) (c.m. of x + A)",
                                  title="distortionRef (MeV)", unit="MeV"),
            "distortionRatio": _key("enum", "rho = |M|^2/|M_PW|^2 (dwpw) or |M|^2 (dw, the published convention)",
                                    values=("dwpw", "dw"), default="dwpw"),
            "boundState": _key("text", "the s-x bound state of the distortion factor and the DW vertex",
                               grammar="whittaker | yukawa [:rmin in fm, 0-50]"),
            "theta": _key("text", "the c.m. angles of the exit pair relative to p_xA (degrees)",
                          grammar="all | thmin-thmax"),
            "vertexModel": _key("enum", "the entrance vertex: plane wave, or the surface term of the prior-form DWBA",
                                values=("pw", "dw"), default="pw"),
            "spectatorAngles": _key("text", "a window of spectator directions (lab polar angles to the beam, or c.m.)",
                                    grammar="[cm:]thmin-thmax | [cm:]table:<file>"),
            "spectatorAngleNodes": _key("integer", "Gauss-Legendre nodes of the spectatorAngles window, 1 to 64",
                                        applies_when=_present("spectatorAngles")),
            "cbackground": _key("text", "the coherent background: a THM-only amplitude c(E) = c0 + c1 E times the "
                                "entrance vertex, per J^pi and exit pair; terms separated by ';', a value followed "
                                "by f is fixed", grammar="<J><+|->:<exit pair>[:<s>,<l>,<s'>,<l'>][:const|:linear]"
                                                    "[=<v>,<v>[,<v>,<v>]]"),
        },
    },
}


def keys():
    """The ``<thm>`` keys as record kinds of a schema (``GET /api/thm/keys``)."""
    import copy
    return copy.deepcopy(KEYS)


# -- the session's THM diagnostics as rows ---------------------------------------------------------------------------

EXPERIMENT = dict(ARG, name="experiment", kind="experiment", doc="a THM experiment of the project",
                  domain="the name of an experiment[<name>] of <thm>")
ENERGIES = dict(ARG, name="energies", kind="real", unit="MeV", doc="c.m. energy E of x + A", domain="finite",
                example=0.5, sweep={"lower": 0.05, "upper": 2.0, "spacing": "linear"})


def _ret(name, unit, doc, lines=None, shown=True):
    """A returned field; ``shown`` False: left off the first plot (the page offers it among the fields)."""
    out = {"name": name, "unit": unit, "doc": doc}
    if lines:
        out["lines"] = lines
    if not shown:
        out["shown"] = False
    return out


def _row(name, summary, doc, arguments, returns, *, axis=None, facts=(), status=None, at_x=True,
         max_elements=MAX_ELEMENTS):
    """A row; ``at_x``: the call depends on the parameters (the page evaluates it at its working x);
    ``max_elements``: the energies one request may ask for (lower than the physics rows' where one costs more)."""
    return {"name": name, "group": GROUP, "tab": TAB, "op": name, "broadcast": False, "at_x": at_x,
            "max_elements": max_elements, "summary": summary,
            "doc": doc, "library": f"pyazr: Session.{name} (the session's)", "axis": axis,
            "arguments": [dict(a) for a in arguments], "returns": list(returns),
            "facts": [{"name": n, "unit": u, "doc": d} for n, u, d in facts],
            "status_words": dict({"0": "computed"}, **(status or {}))}


# The distortion factor costs about 17 ms per energy with an optical potential (0.3 s for 20 energies, 3.4 s for 200 on
# examples/o18_lacognata2010 without its DW vertex, one thread): 1000 energies keep a request near 20 s.
DISTORTION_ELEMENTS = 1000
NOT_FINITE = {"-1": "not a finite number at this energy"}
NOT_REACHED = {"-3": "beyond the reach of the spectator window or the DW vertex grid: the engine uses the nearest "
                     "data point's vertex there, not the vertex at E"}

THM_ROWS = {
    "thm_experiments": _row(
        "thm_experiments", "Every THM experiment as a chi-squared evaluation at the working x profiles it.",
        "One line per experiment: its data sets, the profiled norm n* and background b(E) = b0 + b1 E + b2 E^2 with "
        "their standard deviations from the profile (at fixed R-matrix parameters), its chi-squared and points; and "
        "the coherent background's fit parameters by name.",
        [], [_ret("segments", "", "the keys of its data sets"),
             _ret("background", "", "none, const, linear, quadratic"),
             _ret("points", "", "data points"), _ret("chi2", "", "the experiment's chi-squared"),
             _ret("status", "", "profiled, or what a degenerate profile fell back to"),
             _ret("norm", "", "the profiled norm n* (multiplies the data)"), _ret("norm_sigma", "", "its sigma"),
             _ret("b0", "", "background constant (the model's unit)"), _ret("b0_sigma", "", "its sigma"),
             _ret("b1", "", "background slope (the model's unit per MeV)"), _ret("b1_sigma", "", "its sigma"),
             _ret("b2", "", "background curvature (the model's unit per MeV^2)"), _ret("b2_sigma", "", "its sigma")]),
    "thm_background": _row(
        "thm_background", "The shared norm and background of one THM experiment at the working x.",
        "The profiled norm n* and background coefficients with their sigma and covariance from the profile (at fixed "
        "R-matrix parameters, not scaled by chi2/nu); E in b(E) is the c.m. energy of the THM entrance pair.",
        [EXPERIMENT], [_ret("parameter", "", "n* (norm), b0, b1, b2"), _ret("value", "", "the profiled value"),
                       _ret("sigma", "", "sqrt of the covariance's diagonal")],
        facts=[("chi2", "", "the experiment's chi-squared"), ("points", "", "data points"),
               ("background", "", "none, const, linear, quadratic"), ("status", "", "of the profile"),
               ("segments", "", "the keys of its data sets")]),
    "thm_lineshape": _row(
        "thm_lineshape", "The Coulomb line shape |N_C|^2 of an experiment with lineshape=on.",
        "|N_C|^2 = exp[2 zeta arctan(2 (E_l - E)/G_l)] per level of a J group coupling entrance and exit "
        "(Mukhamedzhanov, Kadyrov & Pang, EPJA 56 (2020) 233), zeta = eta_sB - eta_0 per exit pair, and eta_sb, the "
        "neglected s-b parameter (a validity check: it should be << 1).  The level poles are those of the working x.",
        [EXPERIMENT, ENERGIES],
        [_ret("NC2", "", "|N_C|^2", lines="exit pair and level"),
         _ret("zeta", "", "eta_sB - eta_0", lines="exit pair", shown=False),
         _ret("eta_sb", "", "the neglected s-b Sommerfeld parameter", lines="exit pair", shown=False),
         _ret("E_sF", "MeV", "E_aA - B - E", shown=False), _ret("eta_0", "", "Z_s Z_F alpha mu_sF/k_sF", shown=False)],
        axis="energies", facts=[("E_aA", "MeV", "c.m. energy of a + A"), ("B", "MeV", "binding energy of a = x + s"),
                                ("Zs", "", "spectator charge"), ("ZF", "", "charge of F"),
                                ("spectator", "", "the spectator")], status=NOT_FINITE),
    "thm_vertex": _row(
        "thm_vertex", "The THM entrance vertex |M_l|^2 per channel and level, averaged over the spectator window.",
        "<|M_l|^2> over the experiment's spectator-momentum window and |M_l|^2 at p_s = 0 per entrance channel and "
        "level (the vertex boundary of each level from the working x); with vertexModel=dw the surface term of the "
        "prior-form DWBA and its plane-wave limit.  An energy beyond the window's (or the DW grid's) reach has status "
        "-3: the engine uses the nearest data point's vertex there.",
        [EXPERIMENT, ENERGIES],
        [_ret("M2", "", "<|M_l|^2> over the window", lines="channel and level"),
         _ret("M2_qf", "", "|M_l|^2 at p_s = 0 (dw: at the spectatorAngle direction)", lines="channel and level",
              shown=False),
         _ret("M2_pw", "", "dw only: the plane-wave |M_l(p)|^2 at p = |k_aA - alpha k_sF|", lines="channel and level"),
         _ret("dw_q_delta", "MeV/c", "dw only: q of the spectatorAngle node", shown=False),
         _ret("dw_p_delta", "fm^-1", "dw only: p of the spectatorAngle node", shown=False)],
        axis="energies", facts=[("model", "", "pw or dw"), ("window", "", "the ps window (delta without one)"),
                                ("pair", "", "the entrance pair's key"), ("B", "MeV", "binding energy used"),
                                ("mu_sx", "MeV", "reduced mass of s + x"), ("radius", "fm", "channel radius")],
        status=NOT_REACHED),
    "thm_distortion": _row(
        "thm_distortion", "The distortion factor R(E) of an experiment with distortion=.",
        "Zero-range prior-form DWBA: R(E) = rho(E)/rho(E_ref), rho = |M|^2/|M_PW|^2 (dwpw) or |M|^2 (dw); R_model is "
        "the R the model uses (interpolated; the table's for distortion=table, then nothing else).  R does not depend "
        "on the R-matrix parameters.",
        [EXPERIMENT, ENERGIES],
        [_ret("R_model", "", "the R the model uses"), _ret("R", "", "R evaluated at E"),
         *(_ret(n, u, d, shown=False) for n, u, d in (
             ("M2", "fm^6", "|M|^2 (up to a constant)"), ("M2_PW", "fm^6", "|M_PW|^2 (up to a constant)"),
             ("E_sF", "MeV", "c.m. energy of s + F"), ("k_sF", "fm^-1", "its wave number"),
             ("eta_sF", "", "its Sommerfeld parameter"), ("theta_cm", "deg", "spectator c.m. angle to the beam"),
             ("x", "", "cos(k_sF, k_aA)"), ("q", "fm^-1", "|k_sF - beta k_aA|"),
             ("lmax", "", "the highest partial wave summed")))],
        axis="energies", facts=[("kind", "", "coulomb, optical or table"), ("description", "", "the settings"),
                                ("E_ref", "MeV", "reference energy"), ("E_aA", "MeV", "c.m. energy of a + A"),
                                ("B", "MeV", "binding energy of a = x + s"), ("k_aA", "fm^-1", "a + A wave number"),
                                ("eta_aA", "", "a + A Sommerfeld parameter"),
                                ("kappa", "fm^-1", "bound-state wave number"),
                                ("eta_b", "", "bound-state Sommerfeld parameter"), ("beta", "", "m_s/m_a")],
        status=NOT_FINITE, at_x=False, max_elements=DISTORTION_ELEMENTS),
}


def _facts(row, out):
    return [dict(f, value=out[f["name"]]) for f in row["facts"] if f["name"] in out]


def _finite_status(fields, n):
    """0 where every field (and every line of it) is finite at an energy, -1 elsewhere."""
    ok = np.ones(n, bool)
    for v in fields.values():
        a = np.asarray(v, float).reshape(-1, n)
        ok &= np.isfinite(a).all(axis=0)
    return [0 if o else -1 for o in ok]


def _sigma(cov, k):
    return float(np.sqrt(max(cov[k][k], 0.0)))


def _experiments(row, out):
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
    columns = [{"key": "experiment", "title": "experiment", "unit": ""}] + [
        {"key": c["name"], "title": c["name"], "unit": c["unit"]} for c in row["returns"]]
    tables = [{"title": "experiments", "columns": columns, "rows": rows}]
    if cbkg:
        tables.append({"title": "coherent background (fit parameters)", "rows": cbkg, "columns": [
            {"key": "experiment", "title": "experiment", "unit": ""}, {"key": "parameter", "title": "parameter",
                                                                        "unit": ""},
            {"key": "value", "title": "value", "unit": ""}]})
    return {"tables": tables}


def _background(row, out):
    names = ["norm"] + [f"b{k}" for k in range(len(out["b"]))]
    values = [out["norm"]] + list(out["b"])
    cov = out["cov"]
    rows = [{"parameter": n, "value": v, "sigma": out["sigma"][k]} for k, (n, v) in enumerate(zip(names, values))]
    covariance = [dict({"parameter": n}, **{m: cov[k][j] for j, m in enumerate(names)}) for k, n in enumerate(names)]
    head = [{"key": "parameter", "title": "parameter", "unit": ""}]
    facts = _facts(row, dict(out, segments=" ".join(map(str, out["segments"]))))
    return {"facts": facts, "tables": [
        {"title": "profile", "columns": head + [{"key": c, "title": c, "unit": ""} for c in ("value", "sigma")],
         "rows": rows},
        {"title": "covariance (from the profile)", "columns": head + [{"key": n, "title": n, "unit": ""}
                                                                      for n in names], "rows": covariance}]}


def _half(j):
    """A spin as AZURE2 writes it: 1, 3/2."""
    j = float(j)
    return str(int(j)) if j == int(j) else f"{int(round(2 * j))}/2"


def _spin(j, parity):
    return f"{_half(j)}{'+' if int(parity) > 0 else '-'}"


def _lineshape(row, out):
    n = len(out["E"])
    returns, lines = {"NC2": [], "zeta": [], "eta_sb": []}, {"NC2": [], "zeta": [], "eta_sb": []}
    for x in out["exits"]:
        for k in ("zeta", "eta_sb"):
            returns[k].append(x[k])
            lines[k].append(f"exit {x['pair']}")
        for lev in x["levels"]:
            returns["NC2"].append(lev["NC2"])
            lines["NC2"].append(f"exit {x['pair']}, level {lev['level']} ({_spin(lev['J'], lev['pi'])}, "
                                f"E {lev['E_level']:.6g} MeV, Gamma {lev['Gamma']:.4g} MeV)")
    returns.update(E_sF=out["E_sF"], eta_0=out["eta_0"])
    return {"returns": returns, "lines": lines, "status": _finite_status(returns, n), "facts": _facts(row, out)}


def _vertex(row, out):
    returns, labels = {"M2": [], "M2_qf": [], "M2_pw": []}, []
    for c in out["channels"]:
        for lev in c["levels"]:
            labels.append(f"{_spin(c['J'], c['pi'])} channel {c['channel']} (l {c['l']}, s {_half(c['s'])}), level "
                          f"{lev['level']}")
            for k in returns:
                if k in lev:
                    returns[k].append(lev[k])
    returns = {k: v for k, v in returns.items() if v}
    lines = {k: labels for k in returns}
    for k in ("dw_q_delta", "dw_p_delta"):
        if k in out:
            returns[k] = out[k]
    status = [0 if r else -3 for r in out["reached"]]
    return {"returns": returns, "lines": lines, "status": status, "facts": _facts(row, out)}


def _distortion(row, out):
    returns = {r["name"]: out[r["name"]] for r in row["returns"] if r["name"] in out}
    return {"returns": returns, "status": _finite_status(returns, len(out["E"])), "facts": _facts(row, out)}


# How each row's answer is made from the adapter's plain answer: returns (fields on the energies), lines, status,
# facts, tables (the parts of an answer of a physics function a THM call fills).
ANSWERS = {"thm_experiments": _experiments, "thm_background": _background, "thm_lineshape": _lineshape,
           "thm_vertex": _vertex, "thm_distortion": _distortion}


def answer(fn, row, out):
    """The parts of the answer of ``fn`` (a row of :data:`THM_ROWS`) from the adapter's answer ``out``."""
    parts = ANSWERS[fn](row, out)
    parts.setdefault("returns", {})
    parts.setdefault("status", [])
    return parts
