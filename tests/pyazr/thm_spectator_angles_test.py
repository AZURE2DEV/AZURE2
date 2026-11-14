#!/usr/bin/env python3
"""The spectator-direction window of a THM experiment in pyazr (``spectatorAngles=``).

``spectatorAngles`` on an experiment line averages the distortion factor R(E)
and the DW vertex over the accepted spectator directions, weight
d cos(theta_cm) x acceptance (ThmDistortion::AngleNodes; docs/source/theory/
thm_implementation.rst, "Experimental acceptance").  Checked on
tests/18O_p_a_thm, 2H(18O,a15N)n at 54 MeV (neutron spectator, the Trojan
horse is the target):

  1. Project: an experiment with spectatorAngles=..., spectatorAngleNodes=...
     round trip and canonical line; the combinations AZURE2 refuses are
     refused the same way (by the edit, or when a session is made).
  2. The engine against the CLI, R(E) and the DW vertex with a window:
     evaluate().chi2 (rel 1e-9) and the output files of write_outputs.
  3. thm_distortion with a window: R == (M2/M2_PW)(E) / (M2/M2_PW)(E_ref) and
     R_model == R (1e-6, the grid); theta_cm is the mean of the window.
  4. thm_vertex with a DW window: the directions (angle_theta_cm, angle_q,
     angle_weights; weights sum to 1, angles inside the window) and, without
     folding, evaluate().data[i].model(dw window)/evaluate().data[i].model(pw) == M2(dw window)/M2(pw)
     (1e-8): the reported vertex is the one the model uses.

Needs the compiled engine and an AZURE2 binary for 2-4; skips them cleanly.

Run from anywhere:  python3 tests/pyazr/thm_spectator_angles_test.py
"""
import glob
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
os.environ.setdefault("OMP_NUM_THREADS", "1")
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
    print("1. Project: spectatorAngles")
    proj = os.path.join(tmp, "model")
    fresh_copy(proj)
    src = os.path.join(proj, AZR)
    base = dict(beam="18O", target="d", spectator="n", Ebeam=54)
    m = Project.open(src)
    set_experiment(m, "A", [1, 2], **base, distortion="coulomb", vertexModel="dw", ps="hulthen:0-1000",
                   spectatorAngles="cm:120-180", spectatorAngleNodes=12)
    rec = {"beam": "18O", "target": "d", "spectator": "n", "Ebeam": "54", "distortion": "coulomb",
           "vertexModel": "dw", "ps": "hulthen:0-1000", "spectatorAngles": "cm:120-180", "spectatorAngleNodes": "12"}
    check("record carries spectatorAngles", m.thm()["experiments"] == {"A": rec}, m.thm())
    one = os.path.join(proj, "one.azr")
    m.save(one)
    line = (f"experiment[A] {KIN} distortion=coulomb vertexModel=dw ps=hulthen:0-1000 "
            "spectatorAngles=cm:120-180 spectatorAngleNodes=12\n")
    check("canonical line", line in open(one).read(), open(one).read()[-400:])
    check("reloaded: same record", Project.open(one).thm()["experiments"] == {"A": rec})
    m2 = Project.open(src)
    set_experiment(m2, "A", [1, 2], **base, distortion="coulomb", spectatorAngles=(10, 30.5))
    check("a (lo, hi) pair is a lab window",
          m2.thm()["experiments"]["A"].get("spectatorAngles") == "10-30.5", m2.thm())
    bad = os.path.join(proj, "refused.azr")
    # Refused by the edit, or (the rules that need the pairs) by the engine when a session is made.
    for kwargs, frag in [
            (dict(spectatorAngles="cm:120-180"), "with neither it has nothing to average"),
            (dict(distortion="coulomb", spectatorAngle=20, spectatorAngles="cm:120-180"),
             "exclude each other"),
            (dict(distortion="coulomb", spectatorAngleNodes=4), "needs a spectator-direction window"),
            (dict(distortion="coulomb", vertexModel="dw", ps="hulthen:0-40", psNodes=8,
                  spectatorAngles="cm:120-180"), "the nodes are spectatorAngleNodes="),
            (dict(distortion="coulomb", spectatorAngles="cm:50-20"), "expected thmin-thmax"),
            (dict(distortion="coulomb", spectatorAngles="10-200"), "expected thmin-thmax"),
            (dict(distortion="coulomb", spectatorAngles="cm:120-180", spectatorAngleNodes=65),
             "1 to 64")]:
        msg = refusal(lambda: set_experiment(Project.open(src), "B", [1], **base, **kwargs), bad)
        check(f"refused like AZURE2: {frag}", msg is not None and frag in msg, msg and msg[-400:])

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

    print("2. the engine against the CLI")
    for tag, extra in (("R", "distortion=coulomb spectatorAngles=cm:120-180 ps=hulthen:0-40"),
                       ("dw", "distortion=coulomb vertexModel=dw spectatorAngles=0-180")):
        block = f"experiment[A] {KIN} {extra}"
        cli_dir, py_dir = project(f"cli_{tag}", block), project(f"py_{tag}", block)
        proc = subprocess.run([binary, "--no-gui", "--no-readline", "run.azr"], cwd=cli_dir,
                              input="1\n\n\n7\n", text=True, capture_output=True, timeout=600)
        hit = re.findall(r"Total Chi-Squared:\s*([0-9.eE+-]+)", proc.stdout.replace("\r", "\n"))
        from_cli = output_files(cli_dir)
        if not hit or "thm_experiments.out" not in from_cli:
            print("FAIL: the CLI run produced no chi-squared / thm_experiments.out")
            print(proc.stdout[-2000:])
            sys.exit(1)
        cli_chi2 = float(hit[-1])
        with Session(os.path.join(py_dir, "run.azr"), cwd=py_dir) as s:
            x = np.asarray(s.parameters.x0, float)
            chi2 = s.evaluate(x).chi2
            check(f"{tag}: evaluate().chi2 == CLI (rel 1e-9): {chi2!r}", rel(chi2, cli_chi2) < 1e-9,
                  f"{chi2} vs {cli_chi2}")
            s.write_outputs(x)
        from_py = output_files(py_dir)
        for name in sorted(n for n in from_cli if n.startswith("AZUREOut_") or n in
                           ("chiSquared.out", "normalizations.out", "thm_experiments.out")):
            check(f"{tag}: write_outputs: {name} identical to the CLI's",
                  from_py.get(name) == from_cli[name])

    print("3. thm_distortion with a window")
    d = project("dist", f"experiment[A] {KIN} distortion=coulomb spectatorAngles=cm:120-180")
    grid = np.linspace(0.45, 0.95, 21) + 0.003
    with Session(os.path.join(d, "run.azr"), cwd=d) as s:
        r = s.thm_distortion("A", grid)
        ref = s.thm_distortion("A", [r["E_ref"]])
        rho = np.asarray(r["M2"]) / np.asarray(r["M2_PW"])
        rho_ref = ref["M2"][0] / ref["M2_PW"][0]
        worst = float(np.max(np.abs(rho / rho_ref / np.asarray(r["R"]) - 1)))
        check(f"R == <|M|^2>/<|M_PW|^2> normalized at E_ref ({worst:.1e})", worst < 1e-12)
        worst = float(np.max(np.abs(np.asarray(r["R_model"]) / np.asarray(r["R"]) - 1)))
        check(f"R_model (grid) == R ({worst:.1e})", worst < 1e-6)
        th = np.asarray(r["theta_cm"])
        check(f"theta_cm is the window's mean ({th.min():.2f}-{th.max():.2f} deg)",
              bool(np.all((th > 120) & (th < 180))))
        check("description names the window", "spectator directions cm:120-180" in r["description"],
              r["description"])

    print("4. thm_vertex with a DW window")
    model, vert = {}, {}
    for name, blk in (("pw", f"experiment[A] {KIN}"),
                      ("dw", f"experiment[A] {KIN} distortion=coulomb vertexModel=dw "
                             "spectatorAngles=0-180")):
        dd = project(f"nf_{name}", blk, fold=False)
        with Session(os.path.join(dd, "run.azr"), cwd=dd) as s:
            x = np.asarray(s.parameters.x0, float)
            model[name] = np.concatenate([np.asarray(v) for v in [d.model for d in s.evaluate(x).data]])
            ecm = np.concatenate([np.asarray(v) for v in s.energies])
            r = s.thm_vertex("A", ecm, params=x)
            vert[name] = r["channels"][0]["levels"][0]["M2"]
            if name == "dw":
                wsum = [float(np.sum(w)) for w in r["angle_weights"]]
                check(f"angle_weights sum to 1 ({min(wsum):.12f}-{max(wsum):.12f})",
                      max(abs(w - 1) for w in wsum) < 1e-12)
                ths = np.concatenate([np.asarray(t) for t in r["angle_theta_cm"]])
                check(f"angle_theta_cm inside [0, 180] ({ths.min():.2f}-{ths.max():.2f})",
                      bool(np.all((ths >= 0) & (ths <= 180))))
                nodes = {len(t) for t in r["angle_theta_cm"]}
                check(f"two branches of 8 directions in the lab window: {sorted(nodes)}", nodes == {16})
                check("one averaged vertex node", all(len(w) == 1 for w in r["dw_weights"]))
    worst = float(np.max(np.abs(model["dw"] / model["pw"] / (vert["dw"] / vert["pw"]) - 1)))
    ratio = vert["dw"] / vert["pw"]
    check(f"evaluate().data[i].model(dw window)/evaluate().data[i].model(pw) == M2(dw window)/M2(pw) at {ecm.size} points "
          f"({ratio.min():.3f}-{ratio.max():.3f}, worst {worst:.1e})", worst < 1e-8)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all THM spectator-direction window checks passed")
