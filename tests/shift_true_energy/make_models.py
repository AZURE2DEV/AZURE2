#!/usr/bin/env python3
"""Writes the projects of tests/shift_true_energy (run once; the outputs are committed).

    python3 tests/shift_true_energy/make_models.py BUILD     (a build directory: BUILD/stopping/azure2-stopping writes
                                                              the table, BUILD/src/AZURE2 the fit's data)

An energy shift is a calibration of the beam energy: the fold of a point measured at the nominal lab energy E_n with
the shift s is the fold of the same point without a shift at E_n + s(E_n) (include/ShiftedFold.h, the rule).  Every
data set of models/identity.azr comes in a pair: `<case>` with a shift and its nominal energies, `<case>_twin` without
a shift at the energies E_n + s(E_n) (17 significant digits), on the same target with the same beam and quadrature.

The physics is tests/stopping_yield's: 15N(p,a)12C across one isolated 1- level of 16O at E_r = 0.9 MeV c.m. (Gamma_p
60 eV, Gamma_a 40 eV; <config> and <levels> copied from tests/stopping_yield/models/stopping_yield.azr), E_r/f = 0.9605
MeV lab.  Four nominal energies: below the resonance, on its rising edge, on the plateau, at the back of the layer.

Cases, every combination of
  quantity and target  thin (quantity sigma on a layer: the target-averaged cross section), layer (quantity yield on
                       a layer, 1/particle), thick (quantity yield on an infinite target)
  stopping             formula  eps = 2.4e-20 x^-1.5 MeV cm^2 per active nucleus (d ln eps/d ln E = -1.5)
                       srim     `stopping evaluation` (SRIM-2013 for p in TiN, through the stopping module)
                       table    `stopping table`: SRIM-2013 for p in TiN from 1 keV to 1.2 MeV in 400 steps, written
                                by azure2-stopping (the thick target reaches 1 keV)
  straggling           none, bohr
  beam                 none, a Gaussian spread of 0.8 keV lab
and a shift, cycling through the three kinds: an offset of 3 keV with a prior (a varied shift, the planner's reach R > 0),
b sqrt(E) with b = 0.003 MeV^1/2 (fixed), and the function a0 + a1 x with a0 = 1 keV, a1 = 0.002 (fixed).  Layers: 1e18
active N atoms/cm^2 (about 24 keV for the formula and TiN), the thick target down to 1 keV.

models/fit.azr: the gradient of a fitted shift on a layer whose stopping depends on the energy.  The formula target
with 3e18 active N atoms/cm^2 (about 72 keV), a beam spread of 1.5 keV lab, a scan of 39 energies from 0.950 to 1.045
MeV lab across the rising edge, the plateau and the back of the layer (where E_out passes E_r and the window's width
enters the yield's slope);
the shift an offset fitted with the prior 0 +- 3 keV, starting at 0.  Its data are the engine's own model values at the
shift 2 keV (this script runs BUILD/src/AZURE2 on the scan with the shift fixed there), with 1 % errors.
"""
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(HERE, "models")
DATA = os.path.join(MODELS, "data")
SOURCE = os.path.join(HERE, "..", "stopping_yield", "models")

NOMINAL = [0.9585, 0.9620, 0.9750, 0.9840]  # MeV lab
OFFSET = 0.003                               # MeV
SQRT_B = 0.003                               # MeV^1/2
FN_A0, FN_A1 = 0.001, 0.002                  # MeV, 1


def shift(kind, e):
    """s(E) in MeV lab, as Segment::TotalEnergyShift forms it."""
    if kind == "offset":
        return OFFSET
    if kind == "sqrt":
        return SQRT_B * math.sqrt(e)
    return FN_A0 + FN_A1 * e


def section(text, name):
    a = text.index("<%s>" % name)
    b = text.index("</%s>" % name) + len("</%s>" % name)
    return text[a:b]


def main():
    os.makedirs(DATA, exist_ok=True)
    source = open(os.path.join(SOURCE, "stopping_yield.azr")).read()
    # the table: SRIM-2013 for p in TiN (azure2-stopping), the total stopping per atom of TiN, eV/(1e15 atoms/cm^2)
    build = os.path.abspath(sys.argv[1])
    tool = os.path.join(build, "stopping", "azure2-stopping")
    cmd = [tool, "stopping", "--ion-z", "1", "--ion-mass", "1.007825", "--material", "TiN", "--from", "0.001", "--to",
           "1.2", "--points", "400"]
    rows = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.splitlines()
    with open(os.path.join(DATA, "tin_p.tab"), "w") as f:
        for r in rows:
            if r.startswith("#"):
                f.write(r + "\n")
            elif r and r[0].isdigit():
                c = r.split()
                f.write("%s %s\n" % (c[0], c[3]))
        f.write("# written by tests/shift_true_energy/make_models.py: %s; columns E_MeV eps_total\n" % " ".join(cmd[1:]))

    stoppings = {
        "formula": ['    composition    Ti 1 N 1', '    stopping       formula "a0*x^a1" a0=2.4e-20 a1=-1.5 MeV*cm2/atom'],
        "srim": ['    composition    Ti 1 N 1', '    stopping       evaluation'],
        "table": ['    composition    Ti 1 N 1', '    stopping       table "data/tin_p.tab" MeV eV/1e15atoms/cm2'],
    }
    targets, beams, scales, datasets = [], [], [], []
    kinds = ["offset", "sqrt", "function"]
    k = 0
    for quantity in ("thin", "layer", "thick"):
        for stop in ("formula", "srim", "table"):
            for strag in ("none", "bohr"):
                tname = "%s_%s_%s" % ("inf" if quantity == "thick" else "lay", stop, strag)
                if not any(t[0] == tname for t in targets):
                    lines = ["target " + tname,
                             "    kind           " + ("infinite" if quantity == "thick" else "layer"),
                             "    active         15N"] + stoppings[stop]
                    if quantity != "thick":
                        lines.insert(3, "    areal_density  1e18 atoms/cm2 fixed")
                    lines.append("    straggling     " + strag)
                    targets.append((tname, lines))
                for beam in ("none", "spread"):
                    kind = kinds[k % 3]
                    k += 1
                    case = "%s_%s_%s_%s_%s" % (quantity, stop, strag, beam, kind)
                    if beam == "spread" and not beams:
                        beams.append(["beam b08", "    spread         0.0008 MeV lab"])
                    if kind == "offset":
                        shift_line = "    shift          offset %.17g MeV fit prior %.17g 0.001" % (OFFSET, OFFSET)
                    elif kind == "sqrt":
                        shift_line = "    shift          sqrt %.17g MeV^1/2 fixed" % SQRT_B
                    else:
                        scales.append(["energy_scale fn_" + case,
                                       '    function       "a0+a1*x" MeV',
                                       "    param          a0 %.17g fixed" % FN_A0,
                                       "    param          a1 %.17g fixed" % FN_A1])
                        shift_line = "    shift          fn_" + case
                    for twin in (False, True):
                        name = case + ("_twin" if twin else "")
                        with open(os.path.join(DATA, name + ".dat"), "w") as f:
                            f.write("# E_lab(MeV) unused value error (make_models.py)%s\n" %
                                    (": the nominal energies plus the shift of " + case if twin else ""))
                            for e in NOMINAL:
                                f.write("%.17g 0 1 1\n" % (e + shift(kind, e) if twin else e))
                        lines = ["dataset " + name,
                                 "    reaction       15N(p,a)12C",
                                 "    quantity       " + ("sigma" if quantity == "thin" else "yield"),
                                 "    unit           " + ("b" if quantity == "thin" else "1/particle"),
                                 "    file           data/%s.dat" % name,
                                 "    columns        E_lab[MeV] _ value error",
                                 "    norm           1 fixed"]
                        if not twin:
                            lines.append(shift_line)
                        lines.append("    target         " + tname)
                        if beam == "spread":
                            lines.append("    beam           b08")
                        lines.append("    quadrature     points 40")
                        datasets.append(lines)
    out = [section(source, "config"), section(source, "levels"), "<data>",
           "# make_models.py: every data set with a shift, then its twin without a shift at the shifted energies"]
    for s in scales:
        out += s + [""]
    for d in datasets:
        out += d
    out += ["</data>", "<targets>"]
    for _, t in targets:
        out += t + [""]
    for b in beams:
        out += b + [""]
    out += ["</targets>"]
    with open(os.path.join(MODELS, "identity.azr"), "w") as f:
        f.write("\n".join(out) + "\n")
    # the same pairs in two projects, so that the two members of a pair do not share the nodes of one mesh
    for part, keep in (("shifted", lambda d: not d[0].endswith("_twin")), ("twins", lambda d: d[0].endswith("_twin"))):
        body = [section(source, "config"), section(source, "levels"), "<data>"]
        if part == "shifted":
            for s_ in scales:
                body += s_ + [""]
        for d in datasets:
            if keep(d):
                body += d
        body += ["</data>", "<targets>"]
        for _, t in targets:
            body += t + [""]
        for b in beams:
            body += b + [""]
        body += ["</targets>"]
        with open(os.path.join(MODELS, part + ".azr"), "w") as f:
            f.write("\n".join(body) + "\n")
    print("wrote %d data sets, %d targets" % (len(datasets), len(targets)))
    write_fit(source, os.path.join(build, "src", "AZURE2"))


FIT_ENERGIES = [0.950 + 0.0025 * i for i in range(39)]
FIT_TRUE_SHIFT = 0.002


def fit_project(source, shift_line, data_file):
    return "\n".join([section(source, "config"), section(source, "levels"), "<data>",
                      "dataset scan",
                      "    reaction       15N(p,a)12C",
                      "    quantity       yield",
                      "    unit           1/particle",
                      "    file           data/%s" % data_file,
                      "    columns        E_lab[MeV] _ value error",
                      "    norm           1 fixed",
                      shift_line,
                      "    target         thick_layer",
                      "    beam           b15",
                      "    quadrature     points 40",
                      "</data>", "<targets>",
                      "target thick_layer",
                      "    kind           layer",
                      "    active         15N",
                      "    areal_density  3e18 atoms/cm2 fixed",
                      "    composition    Ti 1 N 1",
                      '    stopping       formula "a0*x^a1" a0=2.4e-20 a1=-1.5 MeV*cm2/atom',
                      "    straggling     none",
                      "",
                      "beam b15",
                      "    spread         0.0015 MeV lab",
                      "</targets>"]) + "\n"


def write_fit(source, engine):
    """models/fit.azr and its data: the engine's model values at the shift FIT_TRUE_SHIFT, 1 % errors."""
    import shutil
    import tempfile
    work = tempfile.mkdtemp(prefix="shift_true_energy_fit.")
    try:
        os.makedirs(os.path.join(work, "data"))
        os.makedirs(os.path.join(work, "output"))
        os.makedirs(os.path.join(work, "checks"))
        with open(os.path.join(work, "data", "scan.dat"), "w") as f:
            for e in FIT_ENERGIES:
                f.write("%.17g 0 1 1\n" % e)
        with open(os.path.join(work, "truth.azr"), "w") as f:
            f.write(fit_project(source, "    shift          offset %.17g MeV fixed" % FIT_TRUE_SHIFT, "scan.dat"))
        subprocess.run([engine, "--no-gui", "--no-readline", "truth.azr"], input="1\n\n\n7\n", cwd=work, check=True,
                       capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="1"))
        out = [l.split() for l in open(os.path.join(work, "output", "AZUREOut_aa=1_R=2.out"))
               if len(l.split()) >= 5 and not l.startswith("#")]
        values = [float(r[3]) for r in out]
    finally:
        shutil.rmtree(work)
    if len(values) != len(FIT_ENERGIES):
        raise SystemExit("the engine gave %d values for %d energies" % (len(values), len(FIT_ENERGIES)))
    with open(os.path.join(DATA, "fit_scan.dat"), "w") as f:
        f.write("# E_lab(MeV) unused yield error: the engine's yield at the shift %g MeV, 1 %% errors (make_models.py)\n"
                % FIT_TRUE_SHIFT)
        for e, y in zip(FIT_ENERGIES, values):
            f.write("%.17g 0 %.10e %.10e\n" % (e, y, 0.01 * y))
    with open(os.path.join(MODELS, "fit.azr"), "w") as f:
        f.write(fit_project(source, "    shift          offset 0 MeV fit prior 0 0.003", "fit_scan.dat"))
    print("wrote fit.azr: %d energies, data at the shift %g MeV" % (len(values), FIT_TRUE_SHIFT))


if __name__ == "__main__":
    main()
