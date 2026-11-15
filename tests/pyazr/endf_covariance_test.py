#!/usr/bin/env python3
"""MF=32 of pyazr.endf (ENDF-102, 2023 edition, pp. 253-259; docs/upgrade/endf.md section 13.13): the covariance of
MF=2's formal parameters from HESSE in a formal-basis session at the transformed best fit, written and read back.

On T1 fitted by the engine's MIGRAD:
(a) the formal-basis session's objective at its start equals the fit's at x (the same model; to 1e-9);
(b) the covariance is symmetric and positive semi-definite; the fixed parameters of the listed resonances have zero
    rows; taken to the Brune basis through the Jacobian of the engine's door 2 (formal levels with B_c -> Brune
    levels, central differences), J C J^T is printed beside the Brune fit's own covariance (not compared: each
    session fixes the project's fixed parameters in its own basis, endf.md 13.13);
(c) T1 written with covariance="hesse": the own format check has no error; MF=32 read back (Material.covariance)
    equals the written matrix to the number form's digits (1e-6 of the largest element of a row) and its parameters
    are exactly MF=2's, each name a spin group and resonance of MF=2; MF=1 names MF=32;
(d) refused with what, allowed, fix: the fit's own covariance (Brune basis) as covariance=, a session of the
    extrapolations;
(e) with AZURE2_ENDF_UTILITIES set, CHECKR, FIZCON and PSYCHE on the written file: every finding printed.
Exits 77 (skipped) without the compiled engine module or iminuit.  Run from anywhere.
"""
import os
import shutil
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
try:
    import numpy as np
    import iminuit  # noqa: F401
    import pyazr
    from pyazr import _azure2  # noqa: F401
    import pyazr.endf as E
    from pyazr.endf.resonances import lab_factor
except ImportError as e:
    print(f"skip: {e}")
    sys.exit(77)

warnings.simplefilter("ignore")
TMP = tempfile.mkdtemp(prefix="endf_cov_")
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {str(detail)[:600]}"))
    if not ok:
        failures.append(name)


d = os.path.join(TMP, "T1")
shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T1"), d, ignore=shutil.ignore_patterns("output*", "checks*"))
for sub in ("output", "checks"):
    os.makedirs(os.path.join(d, sub), exist_ok=True)
os.chdir(d)
s = pyazr.Session("T1.azr")
fit = s.fit_with_engine("migrad")
print(f"T1 fitted: {fit.status}, objective {fit.objective:.10g}")
cov = E.formal_covariance(s, fit.x)
# (a)
rel = abs(cov.objective_formal - cov.objective_fit) / cov.objective_fit
check("the formal-basis session's objective = the fit's to 1e-9", rel < 1e-9, (cov.objective_fit, cov.objective_formal))
# (b)
M = cov.matrix
free = [k for k, f in enumerate(cov.fixed) if not f]
sub = M[np.ix_(free, free)] / np.outer(cov.sigma[free], cov.sigma[free])
eig = np.linalg.eigvalsh(sub)
zero_rows = [k for k, f in enumerate(cov.fixed) if f]
check("symmetric, positive semi-definite, the fixed parameters' rows zero",
      np.allclose(M, M.T, rtol=0, atol=0) and all(np.all(M[k] == 0) for k in zero_rows) and
      eig.min() > -1e-9 * eig.max(), (eig.min(), zero_rows))
record, _ = E.from_session(s, fit.x)
# the formal covariance taken to the Brune basis through the engine's door 2 (formal levels with B_c -> Brune levels):
# J C J^T against the Brune fit's own HESSE covariance of the level parameters (two numerical Hessians of one
# objective in two coordinates: equal to the step of the differences, not to rounding)
from pyazr.endf.resonances import amplitude_scale
e_scale, g_scale = lab_factor(record.pairs[0]) * 1e6, amplitude_scale(record.pairs[0])
ent = next(p.to_dict() for p in s.pairs if p.is_entrance)
s_ent = ent["sep_energy"] + ent["excitation"]
rows = [p.to_dict() for p in s.parameters]
pairs = [p.to_dict() for p in s.pairs]
groups = sorted({r["jgroup"] for r in rows if r["kind"] == "width" and r.get("radiation_type") == "P"})
brune_index = {p.name: k for k, p in enumerate(p for p in s.parameters if not p.fixed)}


def brune_of(gi, E, gam):
    """Door 2 on spin group gi: (E_x in MeV, gamma in MeV^1/2) per level, and the Brune session's names."""
    g = record.groups[gi]
    chans = [(record.pairs[c.pair - 1].mt and pairs_by_ppi[c.pair], c.l, int(round(2 * c.s))) for c in g.channels]
    B = [c.boundary for c in g.channels]
    out = pyazr._azure2.observed_levels_of_block(pairs, int(round(2 * g.J)), g.parity, chans, list(E), list(gam), B,
                                                 1e-12)
    return np.asarray(out["E"], float), np.asarray(out["gamma"], float)


entrance_key = ent["key"]
pairs_by_ppi = {1: entrance_key}
for k, p in enumerate([q for q in pairs if q["pair_kind"] == "particle" and q["key"] != entrance_key], start=2):
    pairs_by_ppi[k] = p["key"]
par_index = {p: k for k, p in enumerate(cov.parameters)}
names_of = {}  # (group, level, column) -> the Brune session's parameter name
for gi, gk in enumerate(groups):
    levels = sorted({r["level"] for r in rows if r["kind"] == "energy" and r.get("jgroup") == gk})
    chans = sorted({r["channel"] for r in rows if r["kind"] == "width" and r.get("jgroup") == gk and
                    r.get("radiation_type") == "P"})
    for p, r in zip(s.parameters, rows):
        if r.get("jgroup") != gk or r["kind"] not in ("energy", "width"):
            continue
        if r["kind"] == "energy":
            names_of[(gi, levels.index(r["level"]), 0)] = p.name
        elif r.get("radiation_type") == "P":
            names_of[(gi, levels.index(r["level"]), 1 + chans.index(r["channel"]))] = p.name
J = np.zeros((len(cov.parameters), len(cov.parameters)))
for gi in cov.groups:
    g = record.groups[gi]
    E0 = np.array([lv.energy / e_scale + s_ent for lv in g.levels])
    G0 = np.array([v / g_scale for lv in g.levels for v in lv.gam])
    nch = len(g.channels)
    for li in cov.levels[gi]:
        for col in range(nch + 1):
            a = par_index[(gi, li, col)]
            h = 1e-6 if col == 0 else 1e-6 * max(abs(G0[li * nch + col - 1]), 1e-3)
            outs = []
            for sign in (1, -1):
                E1, G1 = E0.copy(), G0.copy()
                if col == 0:
                    E1[li] += sign * h
                else:
                    G1[li * nch + col - 1] += sign * h
                outs.append(brune_of(gi, E1, G1))
            unit = e_scale if col == 0 else g_scale  # the covariance is in the file's units
            for lj in cov.levels[gi]:
                for cj in range(nch + 1):
                    b = par_index[(gi, lj, cj)]
                    hi = outs[0][0][lj] if cj == 0 else outs[0][1][lj * nch + cj - 1]
                    lo = outs[1][0][lj] if cj == 0 else outs[1][1][lj * nch + cj - 1]
                    J[b, a] = (hi - lo) / (2 * h) / unit
predicted = J @ cov.matrix @ J.T
worst = 0.0
for p_ in cov.parameters:
    for q_ in cov.parameters:
        na, nb = names_of.get(p_), names_of.get(q_)
        if na in brune_index and nb in brune_index:
            ref = fit.covariance[brune_index[na], brune_index[nb]]
            scale_ = np.sqrt(fit.covariance[brune_index[na], brune_index[na]] * fit.covariance[brune_index[nb],
                                                                                              brune_index[nb]])
            worst = max(worst, abs(predicted[par_index[p_], par_index[q_]] - ref) / scale_)
print(f"      the formal covariance through door 2's Jacobian against the Brune fit's HESSE: largest difference "
      f"{worst:.2e} of sigma_i sigma_j")
# not a check: the two covariances hold different parameters fixed (the project's fixed level energies and widths are
# fixed in the basis of each session), so they differ (docs/upgrade/endf.md 13.13); the number is printed for the report
# (c)
out = os.path.join(d, "T1.endf")
w = E.write(s, out, fit.x, covariance="hesse", header={"AUTH": "pyazr tests"})
report = E.check_format(out)
check("written with MF=32: the own format check has no error", not report.errors and (32, 151) in
      [(mf, mt) for mf, mt, _ in w.sections], [str(f) for f in report.findings])
m = E.read(out)
back = m.covariance()
res = m.resonances()
n = len(cov.parameters)
scale = np.maximum(np.max(np.abs(cov.matrix), axis=1), 1e-300)
dm = np.max(np.abs(back.matrix - cov.matrix) / scale[:, None]) if n else 0.0
dv = np.max(np.abs(back.values - cov.values) / np.maximum(np.abs(cov.values), 1.0))
mf2 = [(lv.energy, *lv.gam) for gi in cov.groups for li in cov.levels[gi] for lv in [res.groups[gi].levels[li]]]
same_mf2 = np.allclose(np.concatenate([np.atleast_1d(r) for r in mf2]), back.values, rtol=0, atol=0)
print(f"      MF=32 read back: {n} parameters ({len(back.names)} names, e.g. {back.names[:3]}), matrix to {dm:.1e}, "
      f"values to {dv:.1e}")
check("MF=32 read back = the written covariance to 1e-6 of each row's largest element (the number form keeps 7 to "
      "10 digits) and the values exactly MF=2's", len(back.names) == n and dm < 1e-6 and same_mf2, (dm, dv, same_mf2))
check("each name a spin group and a resonance of MF=2", all("#" in x and (x.endswith(".ER") or ".GAM" in x)
                                                            for x in back.names) and
      list(back.groups) == list(cov.groups), (back.names, back.groups))
text = " ".join(" ".join(m.header.text).split())
check("MF=1 names MF=32 and the fixed parameters", "MF=32 (LCOMP=1, one block)" in text and
      (not zero_rows or "zero rows" in text), text[-400:])
# (d)
for what, call in (("the Brune fit's covariance", lambda: E.write(s, os.path.join(d, "r.endf"), fit.x,
                                                                   covariance=fit.covariance)),
                   ("a session of the extrapolations", lambda: E.formal_covariance(
                       pyazr.Session("T1.azr", extrapolations=True)))):
    try:
        call()
        check(f"{what}: refused", False)
    except pyazr.FormatRefused as err:
        check(f"{what}: refused with what, allowed, fix", bool(err.what and err.allowed and err.fix), err)
# (e)
if os.environ.get("AZURE2_ENDF_UTILITIES") and os.environ.get("AZURE2_PREPRO"):
    print("\n".join("      " + line for line in str(E.check(out)).splitlines()))
s.close()
os.chdir(HERE)
shutil.rmtree(TMP, ignore_errors=True)
print(f"{len(failures)} failure(s)" + (": " + ", ".join(failures) if failures else ""))
sys.exit(1 if failures else 0)
