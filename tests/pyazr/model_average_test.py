#!/usr/bin/env python3
"""Model averaging (pyazr/modelavg.py): weights, spreads and the .azr writer.

  1. Weights against closed forms: Akaike exp(-dAIC/2), BIC, chi2-only, flat,
     prior weights, chi2 rescaled by chi2/nu of the best variant.
  2. Spreads against closed forms: mean, stat (weighted mean of the
     variances), model (weighted variance of the means), total; partial
     coverage, a missing error, derived (value, sigma) quantities, the
     averaged covariance matrices.
  3. Refusals, the different-N warning, JSON/CSV output, Variant round trip.
  4. Names: parameter_label of Parameter-like records (pair numbers
     translated to file keys).
  5. The writer on tests/18O_p_a_thm: values land in the right <levels>
     fields, unknown names are skipped with a warning, nothing else changes,
     no values gives the template byte for byte, the template object is not
     modified.
  6. With the engine (skips without numpy or the compiled module): every free
     parameter of a session gets a name the writer finds in the file, and
     writing the session's own physical values gives back its model.

Pure Python except part 6: modelavg.py is loaded directly (no package
import), so 1-5 run where numpy is missing.

Run from anywhere:  python3 tests/pyazr/model_average_test.py
"""
import csv
import importlib.util
import json
import math
import os
import shutil
import sys
import tempfile
import warnings
from types import SimpleNamespace

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "18O_p_a_thm")
AZR = "18O_p_a_thm.azr"

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def close(a, b, tol=1e-12):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def load(name):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(ROOT, "pyazr", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ma = load("modelavg")
AzrModel = load("azrfile").AzrModel
Variant = ma.Variant


def refuses(what, fn):
    try:
        fn()
    except ValueError:
        check(f"refused: {what}", True)
        return
    check(f"refused: {what}", False, "no ValueError")


# ---------------------------------------------------------------------------
print("1. weights")
A = Variant("a", {"r": 4.1}, chi2=10.0, npoints=50, nfree=2,
            values={"x": 1.0, "y": 5.0}, errors={"x": 0.1, "y": 1.0})
B = Variant("b", {"r": 5.1}, chi2=12.0, npoints=50, nfree=2,
            values={"x": 2.0, "y": 6.0}, errors={"x": 0.2, "y": 1.0})
C = Variant("c", {"r": 6.1}, chi2=20.0, npoints=50, nfree=3,
            values={"x": 4.0}, errors={"x": 0.3})
V = [A, B, C]


def normalised(raw):
    s = sum(raw)
    return [r / s for r in raw]


def got(avg):
    return [w["weight"] for w in avg.weights]


aic = ma.model_average(V)                       # AIC 14, 16, 26
want = normalised([1.0, math.exp(-1.0), math.exp(-6.0)])
check("AIC = chi2 + 2k, w ~ exp(-dAIC/2)",
      all(close(g, w) for g, w in zip(got(aic), want)), f"{got(aic)} vs {want}")
check("AIC deltas", [w["delta"] for w in aic.weights] == [0.0, 2.0, 12.0],
      [w["delta"] for w in aic.weights])
bic = ma.model_average(V, weights="bic")
ln50 = math.log(50)
want = normalised([math.exp(-0.5 * (c + k * ln50 - (10 + 2 * ln50)))
                   for c, k in ((10, 2), (12, 2), (20, 3))])
check("BIC = chi2 + k ln N", all(close(g, w) for g, w in zip(got(bic), want)))
chi = ma.model_average(V, weights="chi2")
want = normalised([1.0, math.exp(-1.0), math.exp(-5.0)])
check("chi2 weights exp(-dchi2/2)", all(close(g, w) for g, w in zip(got(chi), want)))
flat = ma.model_average(V, weights="FLAT")
check("flat weights 1/3", all(close(g, 1 / 3) for g in got(flat)))
pr = ma.model_average(V, priors={"a": 0.5})
want = normalised([0.5, math.exp(-1.0), math.exp(-6.0)])
check("prior multiplies the weight", all(close(g, w) for g, w in zip(got(pr), want)))
A0 = Variant("a", chi2=10.0, npoints=50, nfree=2, values={"x": 1.0}, prior=0.0)
pz = ma.model_average([A0, B, C])
want = normalised([0.0, 1.0, math.exp(-5.0)])
check("prior 0 drops the variant, deltas from the live best",
      all(close(g, w) for g, w in zip(got(pz), want)) and pz.weights[1]["delta"] == 0.0)
check("weight(label)", close(aic.weight("b"), got(aic)[1]))
# A prior-0 variant takes no part at all: not as the 'best' of rescale, not
# through a missing error (stat), not in the min/max range.
Z = Variant("z", chi2=1.0, npoints=50, nfree=2, values={"x": 9.0}, prior=0.0)
pz2 = ma.model_average([A, B, Z], rescale="best")
xz = pz2.parameters["x"]
check("prior 0: rescale from the best live variant",
      pz2.rescale is None or close(pz2.rescale, max(1.0, 10.0 / 48)), str(pz2.rescale))
check("prior 0: its missing error leaves stat defined", xz.stat is not None, repr(xz))
check("prior 0: not in the min/max range", xz.min == 1.0 and xz.max == 2.0, f"{xz.min} {xz.max}")

H = [Variant("h1", chi2=96.0, npoints=50, nfree=2, values={"x": 1.0}),
     Variant("h2", chi2=100.0, npoints=50, nfree=2, values={"x": 2.0}),
     Variant("h3", chi2=110.0, npoints=50, nfree=2, values={"x": 3.0})]
rb = ma.model_average(H, rescale="best")        # best chi2/nu = 96/48 = 2
want = normalised([1.0, math.exp(-1.0), math.exp(-3.5)])
check("rescale='best' divides chi2 by chi2/nu of the best (2)",
      rb.rescale == 2.0 and all(close(g, w) for g, w in zip(got(rb), want)),
      f"{rb.rescale} {got(rb)}")
r3 = ma.model_average(H, rescale=4.0)
want = normalised([1.0, math.exp(-0.5), math.exp(-1.75)])
check("numeric rescale", all(close(g, w) for g, w in zip(got(r3), want)))
rl = ma.model_average(V, rescale="best")         # 10/48 < 1: no rescale
check("rescale='best' never sharpens (chi2/nu < 1 -> 1)", rl.rescale is None
      and all(close(g, w) for g, w in zip(got(rl), got(aic))))

# ---------------------------------------------------------------------------
print("\n2. spreads")
w = got(aic)
x = aic.parameters["x"]
m = sum(wi * v for wi, v in zip(w, (1.0, 2.0, 4.0)))
s2 = sum(wi * s * s for wi, s in zip(w, (0.1, 0.2, 0.3)))
M2 = sum(wi * (v - m) ** 2 for wi, v in zip(w, (1.0, 2.0, 4.0)))
check("mean = sum w m", close(x.mean, m))
check("stat = sqrt(sum w s^2)", close(x.stat, math.sqrt(s2)))
check("model = sqrt(sum w (m - mean)^2)", close(x.model, math.sqrt(M2)))
check("total = sqrt(stat^2 + model^2)", close(x.total, math.sqrt(s2 + M2)))
check("coverage 1, min/max", close(x.coverage, 1.0) and x.min == 1.0 and x.max == 4.0)
f = flat.parameters["x"]
check("flat: 1, 2, 4 -> mean 7/3, model^2 = 14/9",
      close(f.mean, 7 / 3) and close(f.model ** 2, 14 / 9)
      and close(f.stat ** 2, (0.01 + 0.04 + 0.09) / 3))
y = aic.parameters["y"]
wy = [w[0] / (w[0] + w[1]), w[1] / (w[0] + w[1])]
check("partial coverage: weights renormalised over the carriers",
      close(y.coverage, w[0] + w[1]) and close(y.mean, 5 * wy[0] + 6 * wy[1])
      and close(y.stat, 1.0), y)
D = Variant("d", chi2=10.0, npoints=50, nfree=2, values={"x": 3.0})
nd = ma.model_average([A, D], weights="flat")
check("a carrier without an error -> stat and total None, model still there",
      nd.parameters["x"].stat is None and nd.parameters["x"].total is None
      and close(nd.parameters["x"].model, 1.0))
one = ma.model_average([A])
check("one variant: mean = value, model 0, stat = error",
      one.parameters["x"].mean == 1.0 and one.parameters["x"].model == 0.0
      and close(one.parameters["x"].stat, 0.1))
P = Variant("p", chi2=1.0, npoints=10, nfree=1, values={"x": 1.0},
            derived={"wg": (2.0, 0.5), "S0": 7.0})
Q = Variant("q", chi2=1.0, npoints=10, nfree=1, values={"x": 1.0},
            derived={"wg": (4.0, 1.5), "S0": 9.0})
dq = ma.model_average([P, Q])
check("derived (value, sigma): mean 3, stat sqrt(1.25), model 1",
      close(dq.derived["wg"].mean, 3.0) and close(dq.derived["wg"].stat, math.sqrt(1.25))
      and close(dq.derived["wg"].model, 1.0))
check("derived without sigma: stat None", dq.derived["S0"].stat is None
      and close(dq.derived["S0"].model, 1.0))

CA = Variant("ca", chi2=10.0, npoints=20, nfree=2, values={"u": 1.0, "v": 2.0},
             covariance=(["u", "v"], [[0.04, 0.01], [0.01, 0.09]]))
CB = Variant("cb", chi2=10.0, npoints=20, nfree=2, values={"u": 3.0, "v": 1.0},
             errors={"u": 0.1, "v": 0.2})
cav = ma.model_average([CA, CB])
check("covariance fills errors from its diagonal",
      close(CA.errors["u"], 0.2) and close(CA.errors["v"], 0.3))
names, st, mo, to = cav.covariance()
check("covariance(): names", names == ["u", "v"], names)
check("stat matrix = sum w C_i (diagonal C for errors only)",
      close(st[0][0], 0.5 * 0.04 + 0.5 * 0.01) and close(st[0][1], 0.5 * 0.01)
      and close(st[1][1], 0.5 * 0.09 + 0.5 * 0.04), st)
check("model matrix = sum w d d^T", close(mo[0][0], 1.0) and close(mo[0][1], -0.5)
      and close(mo[1][1], 0.25), mo)
check("total diagonal = Averaged.total^2",
      close(to[0][0], cav.parameters["u"].total ** 2)
      and close(to[1][1], cav.parameters["v"].total ** 2))
check("values: the means", cav.values == {"u": 2.0, "v": 1.5}, cav.values)

# ---------------------------------------------------------------------------
print("\n3. refusals, warnings, output")
refuses("no variants", lambda: ma.model_average([]))
refuses("duplicate labels", lambda: ma.model_average([A, A]))
refuses("unknown method", lambda: ma.model_average(V, weights="dic"))
refuses("bad rescale text", lambda: ma.model_average(V, rescale="worst"))
refuses("rescale <= 0", lambda: ma.model_average(V, rescale=0))
refuses("prior for an unknown label", lambda: ma.model_average(V, priors={"z": 1}))
refuses("every prior 0", lambda: ma.model_average(V, priors={"a": 0, "b": 0, "c": 0}))
refuses("negative prior", lambda: Variant("n", chi2=1, prior=-1))
refuses("non-finite chi2", lambda: Variant("n", chi2=float("nan")))
refuses("covariance shape", lambda: Variant("n", chi2=1, covariance=(["a"], [[1, 2]])))
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    ma.model_average([A, Variant("e", chi2=5, npoints=40, nfree=2, values={"x": 1})])
check("different N warns", any("same number of points" in str(c.message) for c in caught))
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    ma.model_average([A, Variant("e", chi2=5, npoints=40, nfree=2)], weights="flat")
check("... but not for flat weights", not caught)

tmp = tempfile.mkdtemp(prefix="modelavg_test_")
try:
    js = aic.write_json(os.path.join(tmp, "avg.json"))
    d = json.load(open(js))
    check("JSON: weights and parameters",
          [r["label"] for r in d["weights"]] == ["a", "b", "c"]
          and close(d["parameters"]["x"]["model"], x.model))
    back = [Variant.from_dict(v) for v in d["variants"]]
    again = ma.model_average(back)
    check("Variant.to_dict / from_dict round trip",
          all(close(a, b) for a, b in zip(got(again), got(aic)))
          and close(again.parameters["x"].stat, x.stat))
    dd = ma.model_average([Variant.from_dict(v.to_dict()) for v in (P, Q)])
    check("derived errors survive the round trip",
          close(dd.derived["wg"].stat, dq.derived["wg"].stat))
    cc = ma.model_average([Variant.from_dict(CA.to_dict()), CB])
    check("covariance survives the round trip", cc.covariance()[1] == st)
    rows = list(csv.DictReader(open(aic.write_csv(os.path.join(tmp, "avg.csv")))))
    check("CSV: one row per quantity", [r["name"] for r in rows] == ["x", "y"]
          and close(float(rows[0]["total"]), x.total))
    text = aic.summary()
    check("summary has the weights and the table",
          "weight" in text and "model" in text and "c " in text)

    # -----------------------------------------------------------------------
    print("\n4. names")
    pairs = [SimpleNamespace(number=1, key=2), SimpleNamespace(number=2, key=1)]
    e = SimpleNamespace(kind="energy", J=1.5, parity=-1, level=2, name="energy_3")
    g = SimpleNamespace(kind="width", J=2.0, parity=1, level=1, pair=1, L=1,
                        S=0.5, name="width_1_1")
    n = SimpleNamespace(kind="norm", segment_key=3, name="segment_3_norm")
    cb = SimpleNamespace(kind="cbkg", name="cbkg_1+_6_x_re0")
    check("energy", ma.parameter_label(e) == "E[3/2-#2]", ma.parameter_label(e))
    check("width without pairs: engine number", ma.parameter_label(g) == "G[2+#1;p1;L1;S0.5]",
          ma.parameter_label(g))
    check("width with pairs: the file key", ma.parameter_label(g, pairs) == "G[2+#1;p2;L1;S0.5]")
    check("norm, cbkg", ma.parameter_label(n) == "norm[3]"
          and ma.parameter_label(cb) == "cbkg_1+_6_x_re0")

    # -----------------------------------------------------------------------
    print("\n5. the writer")
    work = os.path.join(tmp, "o18")
    shutil.copytree(SOURCE, work)
    template = os.path.join(work, AZR)
    original = open(template).read()
    model = AzrModel.from_file(template)
    slots = ma.azr_parameter_slots(model)
    keys = model.engine_level_keys()
    nE = len(keys)
    nG = sum(len(lv.channels) for lv in keys.values())
    check("one slot per level energy and per channel", len(slots) == nE + nG,
          f"{len(slots)} vs {nE} + {nG}")
    (jg, nlev), lv = sorted(keys.items())[0]
    ename = f"E[{ma.level_label(lv.J, lv.parity, nlev)}]"
    ch = lv.channels[-1]
    gname = ma.width_label(lv.J, lv.parity, nlev, ch.pair, ch.L, ch.S)
    check("names of the first level are slots", ename in slots and gname in slots,
          f"{ename} {gname} {sorted(slots)[:4]}")
    path = os.path.join(work, "avg.azr")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out, written, skipped = ma.write_averaged_azr(
            template, path, {ename: lv.energy + 0.0123, gname: 4321.5,
                             "norm[1]": 1.1, "G[9+#7;p1;L0;S0.5]": 1.0})
    check("written / skipped", sorted(written) == sorted([ename, gname])
          and sorted(skipped) == ["G[9+#7;p1;L0;S0.5]", "norm[1]"], (written, skipped))
    check("skipped names warn", any("skipped" in str(c.message) for c in caught))
    back = AzrModel.from_file(out)
    blv = back.engine_level_keys()[(jg, nlev)]
    check("energy in every line of the level",
          all(close(c.levelE, lv.energy + 0.0123) for c in blv.channels))
    check("width in its channel", close(blv.channels[-1].gamma, 4321.5))
    a_lines = original.splitlines()
    b_lines = open(out).read().splitlines()
    diff = [i for i, (p, q) in enumerate(zip(a_lines, b_lines)) if p != q]
    check("only that level's lines changed",
          len(a_lines) == len(b_lines) and len(diff) == len(lv.channels), diff)
    check("template file untouched", open(template).read() == original)
    ma.write_averaged_azr(template, path, {})
    check("no values: the template byte for byte", open(path).read() == original)
    mdl = AzrModel.from_file(template)
    before = mdl.to_text()
    elsewhere = os.path.join(tmp, "elsewhere") + "/"
    ma.write_averaged_azr(mdl, path, {gname: 7.0}, output_dir=elsewhere)
    check("an AzrModel template is not modified", mdl.to_text() == before)
    check("output_dir repoints <config>", elsewhere in open(path).read())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")              # u, v: not in the template
        ma.write_averaged_azr(template, path, cav)
    check("a ModelAverage is accepted (its means)", open(path).read() == original)

    # -----------------------------------------------------------------------
    print("\n6. the engine")
    try:
        import numpy as np
        sys.path.insert(0, ROOT)
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        from pyazr import azure2
    except Exception as err:                                   # engine not built
        print(f"skip the engine part: engine not available ({type(err).__name__}: {err})")
        sys.exit(1 if failures else 0)
    for junk in ("output", "checks"):
        os.makedirs(os.path.join(work, junk), exist_ok=True)
    free = AzrModel.from_file(template)
    for lv2 in free.levels:                     # free every energy and width
        lv2.set_fixed(False)
        for c in lv2.channels:
            c.channel_fixed = False
    run = free.write(os.path.join(work, "run.azr"))
    with azure2(run, cwd=work) as s:
        x0 = np.asarray(s.params_rwa, float)
        phys = np.asarray(s.transform_rwa(x0), float)
        rm = [p for p in s.parameters if not p.fixed and p.free_index is not None
              and p.kind in ("energy", "width")]
        labels = {ma.parameter_label(p, s.pairs): float(phys[p.free_index]) for p in rm}
        slots = ma.azr_parameter_slots(AzrModel.from_file(run))
        check(f"every free R-matrix parameter ({len(rm)}) has a slot",
              len(rm) > 0 and all(k in slots for k in labels),
              [k for k in labels if k not in slots])
        check("labels are unique", len(labels) == len(rm))
        x0chi = float(np.sum(s.calculate_chi2_rwa(x0)))
    out, written, skipped = ma.write_averaged_azr(run, os.path.join(work, "re.azr"), labels)
    check("all written", not skipped and len(written) == len(labels))
    with azure2(out, cwd=work) as s2:
        chi = float(np.sum(s2.calculate_chi2_rwa(np.asarray(s2.params_rwa, float))))
    check("the session's own values give back its chi2", close(chi, x0chi, 1e-6),
          f"{chi} vs {x0chi}")
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all model-averaging checks passed")
