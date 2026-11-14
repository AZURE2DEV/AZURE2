#!/usr/bin/env python3
"""The distorted-wave THM entrance vertex in pyazr (``vertexModel=dw``).

With ``vertexModel=dw`` on an experiment line the entrance vertex M_l of its
segments is the surface term of the prior-form DWBA built from the
experiment's distorted waves (ThmDwVertex.h; docs/source/theory/
thm_implementation.rst, "Distorted-wave entrance vertex"); R(E) is not
applied.  Checked on tests/18O_p_a_thm, 2H(18O,a15N)n at 54 MeV (neutron
spectator, one 1/2+ group with l = 0, so the HOES model is |M_0|^2 times a
level factor that does not depend on the vertex):

  1. Project: an experiment with vertexModel=dw round trip and canonical
     line; the combinations AZURE2 refuses are refused the same way (by the
     edit, or when a session is made), also set_thm(coulombIntegral /
     entranceL / spectatorEnergy) against a dw experiment.
  2. The engine against the CLI: evaluate().chi2 (rel 1e-9) and the output
     files of write_output_files (thm_experiments.out included).
  3. thm_vertex: model "dw"; plane waves in both channels: M2_qf == M2_pw
     (1e-8, the plane-wave vertex at p = |k_aA - alpha k_sF|) and M2 == M2_qf
     without a window; dw_p_delta from the kinematics.
  4. No folding: evaluate().data[i].model(dw)/evaluate().data[i].model(pw) == M2(dw)/M2(pw) of
     thm_vertex at the data energies, for point Coulomb and for a ps
     window (the nodes on its reachable part); 1e-8, the energies of
     s.energies being those of the points to about 1e-10.

Needs the compiled engine and an AZURE2 binary for 2-4; skips them cleanly.

Run from anywhere:  python3 tests/pyazr/thm_dw_vertex_test.py
"""
import glob
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "18O_p_a_thm")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thm_records import datasets, no_fold, refusal, set_experiment, thm_set  # noqa: E402
AZR = "18O_p_a_thm.azr"
KIN = "beam=18O target=d spectator=n Ebeam=54"

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
try:
    from pyazr import Project
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


def fresh_copy(dst):
    shutil.copytree(SOURCE, dst)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(dst, junk), ignore_errors=True)
        os.makedirs(os.path.join(dst, junk))
    for f in glob.glob(os.path.join(dst, "run.log")):
        os.remove(f)


with tempfile.TemporaryDirectory() as tmp:
    # -- 1. Project -----------------------------------------------------------
    print("1. Project: vertexModel")
    proj = os.path.join(tmp, "model")
    fresh_copy(proj)
    src = os.path.join(proj, AZR)
    base = dict(beam="18O", target="d", spectator="n", Ebeam=54)
    m = Project.open(src)
    set_experiment(m, "A", [1, 2], **base, ps="hulthen:0-30", psNodes=8, distortion="coulomb", vertexModel="dw")
    rec = {"beam": "18O", "target": "d", "spectator": "n", "Ebeam": "54", "ps": "hulthen:0-30", "psNodes": "8",
           "distortion": "coulomb", "vertexModel": "dw"}
    check("record carries vertexModel", m.thm()["experiments"] == {"A": rec}, m.thm())
    one = os.path.join(proj, "one.azr")
    m.save(one)
    line = f"experiment[A] {KIN} ps=hulthen:0-30 psNodes=8 distortion=coulomb vertexModel=dw\n"
    check("canonical line", line in open(one).read(), open(one).read()[-400:])
    check("reloaded: same record", Project.open(one).thm()["experiments"] == {"A": rec})
    bad = os.path.join(proj, "refused.azr")
    # Refused by the edit, or (the rules that need the pairs) by the engine when a session is made.
    for kwargs, frag in [
            (dict(vertexModel="dw"), "needs distortion=coulomb or distortion=optical"),
            (dict(distortion="coulomb", distortionRef=0.7, vertexModel="dw"),
             "belongs to the distortion factor R(E)"),
            (dict(distortion="coulomb", distortionRatio="dw", vertexModel="dw"),
             "belongs to the distortion factor R(E)"),
            (dict(distortion="coulomb", spectatorAngle=10, ps="hulthen:0-30", vertexModel="dw"),
             "spectatorAngle= (one direction) and a ps window"),
            (dict(distortion="coulomb", theta="30-60", vertexModel="dw"), "not available with vertexModel=dw"),
            (dict(distortion="coulomb", vertexModel="dwba"), "expected pw or dw")]:
        msg = refusal(lambda: set_experiment(Project.open(src), "B", [1], **base, **kwargs), bad)
        check(f"refused like AZURE2: {frag}", msg is not None and frag in msg, msg and msg[-400:])

    def with_option(project, key, value):
        project.set_thm({key: value})
        return project

    def dw_after_option(key, value):
        project = with_option(Project.open(src), key, value)
        return set_experiment(project, "A", [1, 2], **base, distortion="coulomb", vertexModel="dw")

    # the options against a dw experiment without a ps window (with one, spectatorEnergy is refused first for the
    # window: "a ps window and spectatorEnergy both set")
    dw_only = os.path.join(proj, "dw_only.azr")
    set_experiment(Project.open(src), "A", [1, 2], **base, distortion="coulomb", vertexModel="dw").save(dw_only)
    for key, value, frag in [("coulombIntegral", True, "coulombIntegral=1"),
                             ("entranceL", "coherent", "entranceL=coherent"),
                             ("spectatorEnergy", 0.3, "spectatorEnergy for entrance pair")]:
        msg = refusal(lambda: with_option(Project.open(dw_only), key, value), bad)
        check(f"set_thm({key}) refused like AZURE2 with a dw experiment", msg is not None and frag in msg,
              msg and msg[-400:])
        msg = refusal(lambda: dw_after_option(key, value), bad)
        check(f"dw experiment refused with {key}={value}", msg is not None and frag in msg, msg and msg[-400:])

    # -- the engine -----------------------------------------------------------
    import numpy as np
    from pyazr import Session
    binary = os.environ.get("AZURE2_BIN")
    if not binary:
        cands = [c for c in glob.glob(os.path.join(ROOT, "build*", "src", "AZURE2*"))
                 if os.path.isfile(c) and os.access(c, os.X_OK)]
        binary = max(cands, key=os.path.getmtime) if cands else None
    if binary is None:
        print("skip the engine part: no AZURE2 binary to compare against")
        sys.exit(1 if failures else 77)

    def project(name, block, fold=True):
        d = os.path.join(tmp, name)
        fresh_copy(d)
        text = open(os.path.join(d, AZR)).read()
        if not fold:
            text = no_fold(text)
        if block:  # the block's one experiment joined by both data sets
            text = thm_set(text, block, datasets(text))
        with open(os.path.join(d, "run.azr"), "w") as f:
            f.write(text)
        return d

    def output_files(work):
        out = os.path.join(work, "output")
        return {os.path.basename(f): open(f).read()
                for f in sorted(glob.glob(os.path.join(out, "*")))
                if os.path.isfile(f) and not os.path.basename(f).startswith("intEC")}

    block = f"experiment[A] {KIN} distortion=coulomb vertexModel=dw"
    cli_dir, py_dir = project("cli", block), project("py", block)
    proc = subprocess.run([binary, "--no-gui", "--no-readline", "run.azr"], cwd=cli_dir,
                          input="1\n\n\n7\n", text=True, capture_output=True, timeout=600)
    hit = re.findall(r"Total Chi-Squared:\s*([0-9.eE+-]+)", proc.stdout.replace("\r", "\n"))
    from_cli = output_files(cli_dir)
    if not hit or "thm_experiments.out" not in from_cli:
        print("FAIL: the CLI run produced no chi-squared / thm_experiments.out")
        print(proc.stdout[-2000:])
        sys.exit(1)
    cli_chi2 = float(hit[-1])

    print("2. the engine against the CLI")
    with Session(os.path.join(py_dir, "run.azr"), cwd=py_dir) as s:
        x = np.asarray(s.parameters.x0, float)
        chi2 = s.evaluate(x).chi2
        check(f"evaluate().chi2 == CLI (rel 1e-9): {chi2!r}", rel(chi2, cli_chi2) < 1e-9,
              f"{chi2} vs {cli_chi2}")
        s.write_output_files(x)
    from_py = output_files(py_dir)
    for name in sorted(n for n in from_cli if n.startswith("AZUREOut_") or n in
                       ("chiSquared.out", "normalizations.out", "thm_experiments.out")):
        check(f"write_output_files: {name} identical to the CLI's", from_py.get(name) == from_cli[name])
    check("thm_experiments.out has the DW vertex rows",
          from_cli["thm_experiments.out"].count("dw_vertex_point") == 3)

    print("3. thm_vertex")
    pl_dir = project("plane", f"experiment[A] {KIN} distortion=optical opticalAA=plane "
                              "opticalSF=plane vertexModel=dw")
    grid = np.linspace(0.45, 0.95, 26) + 0.003
    with Session(os.path.join(pl_dir, "run.azr"), cwd=pl_dir) as s:
        r = s.thm_vertex("A", grid)
        check(f"model {r['model']}", r["model"] == "dw")
        lv = r["channels"][0]["levels"][0]
        worst = float(np.max(np.abs(lv["M2_qf"] / lv["M2_pw"] - 1)))
        check(f"plane waves: M2_qf == M2_pw, the plane-wave vertex at p ({worst:.1e})", worst < 1e-8)
        worst = float(np.max(np.abs(lv["M2"] / lv["M2_qf"] - 1)))
        check(f"no window: M2 == M2_qf ({worst:.1e})", worst < 1e-12)
        # p = |k_aA - alpha k_sF| at qf (k_sF along k_aA), nuclear masses of the table.
        amu, hc = 931.49410242, 197.3269804
        md, mn, mp, m18 = 2.0135532134, 1.0086649159, 1.0072764675, 17.9947732059
        ueff = 931.494
        eaa = 54.0 * md / (m18 + md)
        bind = (mp + mn - md) * amu
        ka = math.sqrt(2 * md * m18 / (md + m18) * ueff * eaa) / hc
        musf = mn * (mp + m18) / (mn + mp + m18) * ueff
        alpha = m18 / (mp + m18)
        pk = [abs(ka - alpha * math.sqrt(2 * musf * (eaa - bind - e)) / hc) for e in grid]
        worst = max(rel(a, b) for a, b in zip(r["dw_p_delta"], pk))
        check(f"dw_p_delta == |k_aA - alpha k_sF| ({worst:.1e})", worst < 1e-6)

    print("4. no folding: the model ratio is the vertex ratio")
    for tag, extra in (("coulomb", "distortion=coulomb"),
                       ("window", "distortion=coulomb ps=hulthen:0-30 psNodes=8")):
        model, vert = {}, {}
        for name, blk in (("pw", f"experiment[A] {KIN}"),
                          ("dw", f"experiment[A] {KIN} {extra} vertexModel=dw")):
            d = project(f"nf_{tag}_{name}", blk, fold=False)
            with Session(os.path.join(d, "run.azr"), cwd=d) as s:
                x = np.asarray(s.parameters.x0, float)
                model[name] = np.concatenate([np.asarray(v) for v in [d.model for d in s.evaluate(x).data]])
                ecm = np.concatenate([np.asarray(v) for v in s.energies])
                vert[name] = s.thm_vertex("A", ecm, params=x)["channels"][0]["levels"][0]["M2"]
        worst = float(np.max(np.abs(model["dw"] / model["pw"] / (vert["dw"] / vert["pw"]) - 1)))
        ratio = vert["dw"] / vert["pw"]
        check(f"{tag}: evaluate().data[i].model(dw)/evaluate().data[i].model(pw) == M2(dw)/M2(pw) at {ecm.size} points "
              f"({ratio.min():.3f}-{ratio.max():.3f}, worst {worst:.1e})", worst < 1e-8)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all THM DW vertex checks passed")
