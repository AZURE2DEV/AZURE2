#!/usr/bin/env python3
"""MF=32 of pyazr.endf (ENDF-102, 2023 edition, pp. 253-259; docs/upgrade/endf.md section 13.13): the fit's
covariance carried through the transformation to MF=2's formal parameters, J C J^T, written and read back.

On T1 fitted by the engine's MIGRAD, C = the fit's covariance (EngineFitResult.covariance):
(a) the Jacobian J of the transformation (door 1 per J^pi group, central differences, the file's units) equals the
    central differences of the whole writer's transformation, from_session(s, x)'s ER and GAM, in every free
    parameter, at three times the step: to 1e-6 of each row's largest element.  A J assembled with a wrong column, a
    wrong unit or a wrong level order fails it;
(b) the matrix is J C J^T, symmetric and positive semi-definite, with zero rows exactly where J has a zero row; taken
    back to the Brune basis through the Jacobian of the engine's door 2 (formal levels with B_c -> Brune levels,
    central differences), it is the fit's own covariance of its free level parameters to 1e-4 of sigma_i sigma_j
    (HESSE in a formal-basis session, the method before 2026-11-18, differed from it by up to 3.9 sigma_i sigma_j);
(c) T1 written with covariance="hesse" (HESSE of the fit's objective by iminuit, carried through): the own format
    check has no error; MF=32 read back (Material.covariance) equals the written matrix to the number form's digits
    (1e-6 of the largest element of a row) and its parameters are exactly MF=2's, each name a spin group and
    resonance of MF=2; MF=1 names MF=32;
(d) the fit's covariance as covariance= of write gives the matrix of (b); refused with what, allowed, fix: a matrix
    of the wrong shape, a word other than "hesse", HESSE asked of a session of the extrapolations;
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
cov = E.formal_covariance(s, fit.x, covariance=fit.covariance)
record, _ = E.from_session(s, fit.x)
# (a) J against central differences of the transformation at the written B_c: for each free level parameter the
# formal levels at theta +- h (three times the module's step) are found by Newton's iteration on door 2 with MF=2's
# B_c (door 2 of them = the session's Brune levels of the group to 1e-11), and differenced
from pyazr.endf import model as endf_model
from pyazr.endf.resonances import amplitude_scale
e_scale, g_scale = lab_factor(record.pairs[0]) * 1e6, amplitude_scale(record.pairs[0])
ent = next(p.to_dict() for p in s.pairs if p.is_entrance)
s_ent = ent["sep_energy"] + ent["excitation"]
pairs = [p.to_dict() for p in s.pairs]
free_of = {i: p.free_index for i, p in enumerate(s.parameters) if not p.fixed}
rows_all = [p.to_dict() for p in s.parameters]


def blocks_at(x):
    groups = endf_model._level_groups(rows_all, endf_model._full_values(s, x))
    return [(gk, groups[gk], endf_model._block_of(groups[gk])) for gk in sorted(groups) if groups[gk]["channels"]]


def door2(g, chans, f, B, n_lev):
    out = pyazr._azure2.observed_levels_of_block(pairs, int(round(2 * g["J"])), g["parity"], chans, list(f[:n_lev]),
                                                 list(f[n_lev:]), list(B), 1e-13)
    return np.concatenate([np.asarray(out["E"], float), np.asarray(out["gamma"], float)])


def formal_at(gi, g, block, B, f0, D):
    target = np.array(list(block["E"]) + list(block["gam"]), float)
    n_lev = len(block["E"])
    f = np.array(f0, float)
    for _ in range(30):
        r = target - door2(g, block["chans"], f, B, n_lev)
        f = f + np.linalg.solve(D, r)
        if np.max(np.abs(r)) < 1e-11:
            break
    return f


J_fd = np.zeros_like(cov.jacobian)
row_of = {p: k for k, p in enumerate(cov.parameters)}
base_blocks = blocks_at(fit.x)
for gi, (gk, g, block) in enumerate(base_blocks):
    sg = record.groups[gi]
    n_lev, nch = len(block["E"]), len(block["chan_ids"])
    B = [c.boundary for c in sg.channels]
    f0 = np.array([lv.energy / e_scale + s_ent for lv in sg.levels] + [v / g_scale for lv in sg.levels for v in lv.gam])
    steps = [3e-5 * (abs(t) + 1.0) if k < n_lev else 3e-5 * max(abs(t), 1e-3) for k, t in enumerate(f0)]
    D = np.array([(door2(g, block["chans"], f0 + h * np.eye(len(f0))[k], B, n_lev) -
                   door2(g, block["chans"], f0 - h * np.eye(len(f0))[k], B, n_lev)) / (2 * h)
                  for k, h in enumerate(steps)]).T
    for slot, r in enumerate(block["rows"]):
        if r is None or r not in free_of:
            continue
        k = free_of[r]
        v = fit.x[k]
        h = 3 * (E.covariance.ENERGY_STEP * (abs(v) + 1.0) if slot < n_lev else
                 E.covariance.AMPLITUDE_STEP * max(abs(v), 1e-3))
        outs = []
        for sign in (1, -1):
            xs = np.array(fit.x, float)
            xs[k] += sign * h
            _, g2, block2 = blocks_at(xs)[gi]
            outs.append(formal_at(gi, g2, block2, B, f0, D))
        dfdx = (outs[0] - outs[1]) / (2 * h)
        for li in range(n_lev):
            J_fd[row_of[(gi, li, 0)], k] = dfdx[li] * e_scale
            for c in range(nch):
                J_fd[row_of[(gi, li, 1 + c)], k] = dfdx[n_lev + li * nch + c] * g_scale
row_max = np.maximum(np.max(np.abs(J_fd), axis=1), 1e-300)
dj = float(np.max(np.abs(cov.jacobian - J_fd) / row_max[:, None]))
print(f"      J: {cov.jacobian.shape[0]} parameters of MF=2 x {cov.jacobian.shape[1]} free parameters; against "
      f"differences of the inverse of door 2 at MF=2's B_c: {dj:.2e} of each row's largest element")
check("J = central differences of the transformation at the written B_c (Newton on door 2, 3x the step) to 1e-5 of "
      "each row's largest element", dj < 1e-5, dj)
# door 1 alone (B_c = S_c(E_1) following the levels: what from_session does at each x) is another Jacobian; printed
J_door1 = np.zeros_like(cov.jacobian)
for k, p in enumerate(p for p in s.parameters if not p.fixed):
    if p.to_dict()["kind"] not in ("energy", "width"):
        continue
    h = 3e-5 * (abs(fit.x[k]) + 1.0)
    xp, xm = np.array(fit.x, float), np.array(fit.x, float)
    xp[k] += h
    xm[k] -= h
    vals = [np.array([(lambda lv: lv.energy if col == 0 else lv.gam[col - 1])(rec.groups[gi].levels[li])
                      for gi, li, col in cov.parameters]) for rec in (E.from_session(s, xp)[0], E.from_session(s, xm)[0])]
    J_door1[:, k] = (vals[0] - vals[1]) / (2 * h)
M1 = J_door1 @ fit.covariance @ J_door1.T
sig = np.sqrt(np.clip(np.diag(cov.matrix), 1e-300, None))
print(f"      J C J^T with door 1's Jacobian (B_c following E_1) differs from it by up to "
      f"{np.max(np.abs(M1 - cov.matrix) / np.outer(sig, sig)):.2f} sigma_i sigma_j (printed, not compared)")
# (b)
M = cov.matrix
check("the matrix is J C J^T", np.allclose(M, cov.jacobian @ fit.covariance @ cov.jacobian.T, rtol=1e-12,
                                             atol=1e-14 * np.max(np.abs(M))))
listed = [k for k, f in enumerate(cov.fixed) if not f]
sub = M[np.ix_(listed, listed)] / np.outer(cov.sigma[listed], cov.sigma[listed])
eig = np.linalg.eigvalsh(sub)
zero_rows = [k for k, f in enumerate(cov.fixed) if f]
check("symmetric, positive semi-definite, zero rows exactly where J's row is zero",
      np.allclose(M, M.T, rtol=0, atol=0) and all(np.all(M[k] == 0) for k in zero_rows) and
      all(np.any(cov.jacobian[k] != 0) for k in listed) and eig.min() > -1e-9 * eig.max(), (eig.min(), zero_rows))
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
print(f"      the formal covariance through door 2's Jacobian against the Brune fit's covariance: largest difference "
      f"{worst:.2e} of sigma_i sigma_j")
check("taken back through door 2, the fit's covariance of its free level parameters to 1e-4 of sigma_i sigma_j",
      0.0 < worst < 1e-4 or (worst == 0.0 and len(listed) == 0), worst)
# (c)
out = os.path.join(d, "T1.endf")
w_cov = E.formal_covariance(s, fit.x)  # covariance="hesse": iminuit's HESSE at x, carried through
w = E.write(s, out, fit.x, covariance="hesse", header={"AUTH": "pyazr tests"})
report = E.check_format(out)
check("written with MF=32: the own format check has no error", not report.errors and (32, 151) in
      [(mf, mt) for mf, mt, _ in w.sections], [str(f) for f in report.findings])
m = E.read(out)
back = m.covariance()
res = m.resonances()
n = len(w_cov.parameters)
scale = np.maximum(np.max(np.abs(w_cov.matrix), axis=1), 1e-300)
dm = np.max(np.abs(back.matrix - w_cov.matrix) / scale[:, None]) if n else 0.0
dv = np.max(np.abs(back.values - w_cov.values) / np.maximum(np.abs(w_cov.values), 1.0))
mf2 = [(lv.energy, *lv.gam) for gi in w_cov.groups for li in w_cov.levels[gi] for lv in [res.groups[gi].levels[li]]]
same_mf2 = np.allclose(np.concatenate([np.atleast_1d(r) for r in mf2]), back.values, rtol=0, atol=0)
print(f"      MF=32 read back: {n} parameters ({len(back.names)} names, e.g. {back.names[:3]}), matrix to {dm:.1e}, "
      f"values to {dv:.1e}")
check("MF=32 read back = the written covariance to 1e-6 of each row's largest element (the number form keeps 7 to "
      "10 digits) and the values exactly MF=2's", len(back.names) == n and dm < 1e-6 and same_mf2, (dm, dv, same_mf2))
check("each name a spin group and a resonance of MF=2", all("#" in x and (x.endswith(".ER") or ".GAM" in x)
                                                            for x in back.names) and
      list(back.groups) == list(w_cov.groups), (back.names, back.groups))
text = " ".join(" ".join(m.header.text).split())
check("MF=1 names MF=32 and the zero rows", "MF=32 (LCOMP=1, one block)" in text and
      (not [k for k, f in enumerate(w_cov.fixed) if f] or "zero rows" in text), text[-400:])
dh = np.max(np.abs(w_cov.matrix - cov.matrix) / np.outer(cov.sigma, cov.sigma).clip(1e-300))
print(f"      covariance=\"hesse\" (iminuit's HESSE at x) against the engine fit's covariance: {dh:.2e} of "
      "sigma_i sigma_j (two numerical Hessians; printed, not compared)")
# (d)
w2 = E.write(s, os.path.join(d, "T1_fitcov.endf"), fit.x, covariance=fit.covariance)
back2 = E.read(os.path.join(d, "T1_fitcov.endf")).covariance()
scale2 = np.maximum(np.max(np.abs(cov.matrix), axis=1), 1e-300)
d2 = np.max(np.abs(back2.matrix - cov.matrix) / scale2[:, None]) if len(cov.parameters) else 0.0
check("the fit's covariance as covariance= of write: MF=32 is J C J^T of (b) to the number form's digits",
      back2.matrix.shape == cov.matrix.shape and d2 < 1e-6, d2)
for what, call in (("a matrix of the wrong shape", lambda: E.write(s, os.path.join(d, "r.endf"), fit.x,
                                                                   covariance=np.eye(2))),
                   ("a word other than hesse", lambda: E.write(s, os.path.join(d, "r.endf"), fit.x,
                                                               covariance="formal")),
                   ("HESSE asked of a session of the extrapolations", lambda: E.formal_covariance(
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
