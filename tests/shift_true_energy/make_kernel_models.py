#!/usr/bin/env python3
"""Writes the kernel models of tests/shift_true_energy (run once; the outputs are committed).

    python3 tests/shift_true_energy/make_kernel_models.py BUILD     (BUILD/src/AZURE2 writes the fits' data)

The rule of include/ShiftedFold.h for the kernels whose shape depends on the point's energy beyond a translation (the
owner's answer of 2026-11-18, docs/upgrade/OWNER_ANSWERS.md section 8, "What an energy shift means for three
kernels"): the fold of a point measured at E_n with the shift s is the fold of the same point without a shift at
E_n + s.  As in make_models.py, every data set comes in a pair: `<case>` with a shift and its nominal energies,
`<case>_twin` without a shift at E_n + s(E_n) (17 significant digits); a beam profile's reconstructed-energy window
moves with the point (window + s(E_n) in the twin's file).  The shift kinds in turn: an offset with a prior (a varied
shift: the planner's reach R > 0), b sqrt(E) (fixed), the function a0 + a1 x (fixed).

models/kernel_gauss.azr   15N(p,a)12C across the 1- level at E_r/f = 0.9605 MeV lab (tests/stopping_yield's <config>
                          and <levels>), a Gaussian beam spread sigma(E_lab) = a0 x^8 lab (2 keV at 0.96 MeV, 8 % more
                          at 3 keV higher: its 5 sigma support moves with the energy), without a target (`gauss`) and
                          on a layer of 1e18 15N/cm^2 (`gausslayer`: the target kernel's support reaches 25 sigma, so this
                          case holds with the rigid support too; it checks that the layer is not broken by the move),
                          and a free-gas Doppler broadening at 1e6 K alone (`doppler`, a window in velocity whose width
                          grows as sqrt(E))
models/kernel_profile.azr 12C(p,g)13N across the 1/2+ level at E_lab = 0.4588 MeV (tests/13N's <config> and <levels>),
                          a beam profile (one skewed Gaussian, xi 0.47 MeV, omega 15 keV, alpha -1.5) with the
                          detector's Gaussian resolution of 3 keV inside each point's +-4 keV window and the
                          detailed-balance weight (`detailed_balance yes`, whose norm is taken at the point's energy)
models/kernel_udr.azr     17O(n,a)14C (tests/udr_low_energy's <config> and <levels>) with a UDR kernel: the 0.2 ns
                          delay of tests/udr_low_energy, a 125 ns burst and a channel of 7.55 bins per decade, whose
                          width in time is t0(E) ln 10 / 15.1: the kernel is the table read at the point's energy
models/fit_<kernel>.azr   one scan of each kernel with an offset fitted (prior 0 +- 3 keV, start 0), data the engine's
                          own model values at the shift 2 keV, 1 % errors: the gradient of the shift (check 3)
"""
import math
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(HERE, "models")
DATA = os.path.join(MODELS, "data")
TESTS = os.path.join(HERE, "..")

OFFSET = 0.003                  # MeV
SQRT_B = 0.003                  # MeV^1/2
FN_A0, FN_A1 = 0.001, 0.002     # MeV, 1
FIT_TRUE_SHIFT = 0.002          # MeV


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


SOURCES = {
    "15N": os.path.join(TESTS, "stopping_yield", "models", "stopping_yield.azr"),
    "13N": os.path.join(TESTS, "13N", "13N.azr"),
    "17O": os.path.join(TESTS, "udr_low_energy", "model.azr"),
}

GAUSS_A0 = 0.002 / 0.9605 ** 8
TARGETS = [
    "target lay15n",
    "    kind           layer",
    "    active         15N",
    "    areal_density  1e18 atoms/cm2 fixed",
    "    composition    Ti 1 N 1",
    '    stopping       formula "a0*x^a1" a0=2.4e-20 a1=-1.5 MeV*cm2/atom',
    "    straggling     none",
    "",
    "target doppler1e6",
    "    kind           none",
    "    doppler        free_gas 1e6 K fixed",
    "",
    "beam bx8",
    '    spread         formula "a0*x^a1" a0=%.17g a1=8 MeV lab' % GAUSS_A0,
    "",
    "beam bprofile",
    "    profile          skewed_gaussian xi 0.47 MeV omega 0.015 MeV alpha -1.5 weight 1",
    "    truncation       0",
    "    detailed_balance yes",
    "",
    "detector dwin",
    "    energy_resolution 0.003 MeV lab",
    "",
    "detector dudr",
    "    udr            file data/udr_delta.dat flight_path 7 m burst 125 ns channel -7.55 ns centred yes",
    "",
]

# kernel -> (source, reaction, quantity, unit, effect lines, nominal lab energies, window half-width or None)
KERNELS = {
    "gauss": ("15N", "15N(p,a)12C", "sigma", "b", ["    beam           bx8"], [0.9480, 0.9530, 0.9585, 0.9620, 0.9700], None),
    "gausslayer": ("15N", "15N(p,a)12C", "sigma", "b", ["    target         lay15n", "    beam           bx8"],
              [0.9530, 0.9585, 0.9620, 0.9750, 0.9840], None),
    "doppler": ("15N", "15N(p,a)12C", "sigma", "b", ["    target         doppler1e6"],
                [0.9480, 0.9530, 0.9585, 0.9620, 0.9700], None),
    "profile": ("13N", "12C(p,g)13N", "sigma", "b", ["    beam           bprofile", "    detector       dwin"],
                [0.4450, 0.4520, 0.4570, 0.4620, 0.4700], 0.004),
    "udr": ("17O", "17O(n,a)14C", "sigma", "b", ["    detector       dudr"], [0.1000, 0.1500, 0.2000, 0.3000], None),
}
MODEL_OF = {"gauss": "kernel_gauss", "gausslayer": "kernel_gauss", "doppler": "kernel_gauss", "profile": "kernel_profile",
            "udr": "kernel_udr"}


def shift_line(kind, case, scales, prior=OFFSET):
    if kind == "offset":
        return "    shift          offset %.17g MeV fit prior %.17g 0.001" % (OFFSET, OFFSET)
    if kind == "sqrt":
        return "    shift          sqrt %.17g MeV^1/2 fixed" % SQRT_B
    scales.append(["energy_scale fn_" + case, '    function       "a0+a1*x" MeV', "    param          a0 %.17g fixed" % FN_A0,
                   "    param          a1 %.17g fixed" % FN_A1])
    return "    shift          fn_" + case


def write_points(name, energies, window, comment, values=None):
    with open(os.path.join(DATA, name + ".dat"), "w") as f:
        f.write("# E_lab(MeV) unused value error%s (make_kernel_models.py)%s\n" %
                (" window_lo window_hi" if window else "", comment))
        for i, (e, s) in enumerate(energies):
            v = values[i] if values else 1.0
            line = "%.17g 0 %.10e %.10e" % (e + s, v, 0.01 * v if values else 1.0)
            if window:
                line += " %.17g %.17g" % (e - window + s, e + window + s)
            f.write(line + "\n")


def dataset(name, reaction, quantity, unit, window, extra, effects):
    lines = ["dataset " + name,
             "    reaction       " + reaction,
             "    quantity       " + quantity,
             "    unit           " + unit,
             "    file           data/%s.dat" % name,
             "    columns        E_lab[MeV] _ value error" + (" window_lo window_hi" if window else ""),
             "    norm           1 fixed"] + extra + effects
    lines.append("    quadrature     points 40")
    return lines


def project(source, scales, datasets):
    out = [section(source, "config"), section(source, "levels"), "<data>",
           "# make_kernel_models.py: every data set with a shift, then its twin without a shift at the shifted energies"]
    for s in scales:
        out += s + [""]
    for d in datasets:
        out += d
    used = {l.split()[1] for d in datasets for l in d if l.split()[0] in ("target", "beam", "detector")}
    out += ["</data>", "<targets>"]
    for record in "\n".join(TARGETS).split("\n\n"):
        lines = [l for l in record.split("\n") if l]
        if lines and lines[0].split()[1] in used:
            out += lines + [""]
    out += ["</targets>"]
    return "\n".join(out) + "\n"


def main():
    os.makedirs(DATA, exist_ok=True)
    shutil.copy(os.path.join(TESTS, "udr_low_energy", "udr_delta.dat"), os.path.join(DATA, "udr_delta.dat"))
    texts = {k: open(p).read() for k, p in SOURCES.items()}
    kinds = ["offset", "sqrt", "function"]
    per_model = {}
    for kernel, (src, reaction, quantity, unit, effects, nominal, window) in KERNELS.items():
        model = MODEL_OF[kernel]
        scales, datasets = per_model.setdefault(model, (src, [], []))[1:]
        for kind in kinds:
            case = "%s_%s" % (kernel, kind)
            for twin in (False, True):
                name = case + ("_twin" if twin else "")
                pts = [(e, shift(kind, e) if twin else 0.0) for e in nominal]
                write_points(name, pts, window, ": the nominal energies plus the shift of " + case if twin else "")
                extra = [] if twin else [shift_line(kind, case, scales)]
                datasets.append(dataset(name, reaction, quantity, unit, window, extra, effects))
    for model, (src, scales, datasets) in per_model.items():
        with open(os.path.join(MODELS, model + ".azr"), "w") as f:
            f.write(project(texts[src], scales, datasets))
        print("wrote %s.azr: %d data sets" % (model, len(datasets)))
    engine = os.path.join(os.path.abspath(sys.argv[1]), "src", "AZURE2")
    for kernel in ("gauss", "profile", "udr"):
        write_fit(kernel, texts, engine)


def write_fit(kernel, texts, engine):
    """models/fit_<kernel>.azr and its data: the engine's model values at the shift FIT_TRUE_SHIFT, 1 % errors."""
    src, reaction, quantity, unit, effects, nominal, window = KERNELS[kernel]
    lo, hi = min(nominal), max(nominal)
    energies = [lo + (hi - lo) * i / 15.0 for i in range(16)]
    name = "fit_%s_scan" % kernel
    work = tempfile.mkdtemp(prefix="shift_kernel_fit.")
    try:
        for d in ("data", "output", "checks"):
            os.makedirs(os.path.join(work, d))
        shutil.copy(os.path.join(DATA, "udr_delta.dat"), os.path.join(work, "data", "udr_delta.dat"))
        write_to = os.path.join(work, "data", name + ".dat")
        with open(write_to, "w") as f:
            for e in energies:
                line = "%.17g 0 1 1" % e
                if window:
                    line += " %.17g %.17g" % (e - window, e + window)
                f.write(line + "\n")
        truth = [dataset(name, reaction, quantity, unit, window,
                         ["    shift          offset %.17g MeV fixed" % FIT_TRUE_SHIFT], effects)]
        with open(os.path.join(work, "truth.azr"), "w") as f:
            f.write(project(texts[src], [], truth))
        subprocess.run([engine, "calculate", "--no-gui", "truth.azr"], stdin=subprocess.DEVNULL, cwd=work, check=True,
                       capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="1"))
        outs = [p for p in os.listdir(os.path.join(work, "output")) if p.startswith("AZUREOut_aa=")]
        if len(outs) != 1:
            raise SystemExit("fit %s: expected one output file, found %s" % (kernel, outs))
        rows = [l.split() for l in open(os.path.join(work, "output", outs[0]))
                if len(l.split()) >= 5 and not l.startswith("#")]
        values = [float(r[3]) for r in rows]
    finally:
        shutil.rmtree(work)
    if len(values) != len(energies):
        raise SystemExit("fit %s: the engine gave %d values for %d energies" % (kernel, len(values), len(energies)))
    write_points(name, [(e, 0.0) for e in energies], window,
                 ": the engine's value at the shift %g MeV, 1 %% errors" % FIT_TRUE_SHIFT, values)
    fit = [dataset(name, reaction, quantity, unit, window, ["    shift          offset 0 MeV fit prior 0 0.003"], effects)]
    with open(os.path.join(MODELS, "fit_%s.azr" % kernel), "w") as f:
        f.write(project(texts[src], [], fit))
    print("wrote fit_%s.azr: %d energies, data at the shift %g MeV" % (kernel, len(values), FIT_TRUE_SHIFT))


if __name__ == "__main__":
    main()
