#!/usr/bin/env python3
"""Models for tests/transmission_yield, built at run time (check.sh) into a work directory (project format 3).

A transmission or yield data set (quantity transmission or neutron_yield) names a target of kind sample: its
areal_density is the sample thickness n (atoms/b), a parameter <target>.areal_density when it is fitted.

From tests/reaction_rate/models/n12C.azr (n + 12C -> 13C, elastic and capture, external capture off; the s-wave 1/2+
resonance at E_cm = 0.425 MeV with Gamma_n = 34 keV), data sets d (sigma_t, 12C(n,n)12C) and g (sigma_gamma,
12C(n,g)13C) on the data grid:

  bare     sigma_total and sigma, no effect
  thin     transmission T = exp(-N sigma_t) and neutron_yield Y0 = (sg/st)(1 - exp(-N st)), each on its sample target
  conv     thin with a Gaussian beam spread of SIGMA_LAB on both (the transform before the kernel)
  bconv    bare with the same spread (broadening sigma first: what the transform must NOT do)
  truth    extrapolations sigma_t and sigma_gamma on a 0.1 keV lab grid, no spread (mode 3)
  xbare, xthin   extrapolations (mode 3): sigma_t, sigma_gamma and their transmission and yield
  grad     thin + conv in one project, every width free, for the gradient check
  refuse_* data sets the engine must refuse

From tests/natural_neon/natural_neon.azr (three compound nuclei, natural sums):

  nat      only data set 1 (Vaughn, natural sum, no target effect) active
  natT     the same data set as a transmission of a natural sample NNAT thick: T = exp(-NNAT sigma_nat)
  ybare, yelast, yreac   n + 20Ne with 1 keV reaction widths: sigma_total and sigma(n,alpha); the elastic; a
                         neutron_yield to alpha + 17O (its sigma_t must be the total, which differs from the elastic)
  refuse_natsub   data set 7 (Junghans, natural sum under an energy-dependent spread) as a transmission
  refuse_natyield a natural reaction (n+20Ne -> n'+20Ne) as a neutron_yield

The fitted sample thickness ("areal_density n atoms/b fit [prior n sigma]", a target shared by several data sets):

  thin18     thin at N_TRUE = 0.18 atoms/b: its T and Y0 become the data of the fits below
  grad_fit   thin (1, 2), conv (3, 4; a 0.01 atoms/b prior), Doppler 300 K alone (5, 6; one target for a
             transmission and a yield) and Doppler + the Gaussian (7, 8), every thickness fitted, every width free,
             for the gradient check
  refuse_t_* malformed fit forms and a cross section on a sample target

    python3 make_models.py <n12C.azr> <natural_neon.azr> <workdir>

and, once thin18 has been calculated (mode 1),

    python3 make_models.py fitdata <n12C.azr> <thin18 dir> <workdir>

writes data/t18.dat, data/g18.dat (the thin18 values, 1 % errors) and data/t18w.dat (T, 20 % errors) and, with every
level parameter fixed,

  fitT       t18 as a transmission, "areal_density 0.15 atoms/b fit": Sd.areal_density free, no prior
  fitY       g18 as a yield
  fitG       t18 and g18 on one target S, 0.15 fit: one parameter from both
  fitP       t18w, "areal_density 0.15 atoms/b fit prior 0.15 0.01": a 0.01 atoms/b prior about 0.15
"""
import os
import re
import sys

N = 0.15            # sample thickness, atoms/barn (n sigma_t ~ 1 at the 1/2+ peak)
N_TRUE = 0.18       # the thickness of the fit data
NNAT = 0.05
SIGMA_LAB = 0.015   # Gaussian beam spread, MeV lab
E_DATA = [round(0.300 + 0.005 * k, 6) for k in range(71)]   # 0.300 .. 0.650 MeV lab


def split(text):
    """(the file up to <data>, <data> .. </targets> as the three section bodies, the file after </targets>)."""
    m = re.search(r"(.*?)<data>\n(.*?)</data>\n<extrapolations>\n(.*?)</extrapolations>\n<targets>\n(.*?)</targets>\n(.*)",
                  text, re.S)
    if not m:
        sys.exit("make_models: unexpected .azr layout (format 3 with <data> ... </targets>)")
    return m.groups()


def records(body):
    """The records of a section body, each a list of lines (comment lines dropped)."""
    out = []
    for block in re.split(r"\n\s*\n", body):
        lines = [l for l in block.splitlines() if l.strip() and not l.lstrip().startswith("#")]
        if lines:
            out.append(lines)
    return out


def text(recs):
    return "".join("\n".join(r) + "\n\n" for r in recs)


def setfield(rec, key, value):
    """The record with KEY's line replaced (value None: removed) or, absent, appended."""
    out, done = [], False
    for l in rec:
        if l.split()[0] == key and l.startswith("    "):
            done = True
            if value is not None:
                out.append("    %-14s %s" % (key, value))
        else:
            out.append(l)
    if not done and value is not None:
        out.append("    %-14s %s" % (key, value))
    return out


def record(header, fields):
    return [header] + ["    %-14s %s" % kv for kv in fields if kv[1] is not None]


def project(head, tail, data, extrapolations=(), targets=()):
    return (head + "<data>\n" + text(data) + "</data>\n<extrapolations>\n" + text(extrapolations) +
            "</extrapolations>\n<targets>\n" + text(targets) + "</targets>\n" + tail)


def write(path, content):
    with open(path, "w") as f:
        f.write(content)


def levels(head, edit):
    return re.sub(r"(<levels>\n)(.*?)(</levels>)",
                  lambda m: m.group(1) + "".join(edit(l) for l in m.group(2).splitlines(True)) + m.group(3),
                  head, flags=re.S)


def dataset(name, kind, dat, target=None, beam=None, points=None, window="0.29 0.66"):
    """kind: t (sigma_total), g (sigma of the capture), T (transmission), Y (neutron_yield of the capture)."""
    reaction = "12C(n,n)12C" if kind in "tT" else "12C(n,g)13C"
    quantity, unit = {"t": ("sigma_total", "b"), "g": ("sigma", "b"), "T": ("transmission", "1"),
                      "Y": ("neutron_yield", "1")}[kind]
    return record("dataset " + name, [("reaction", reaction), ("quantity", quantity), ("unit", unit),
                                      ("file", "data/" + dat), ("columns", "E_lab[MeV] _ value error"),
                                      ("energy_window", window + " MeV"), ("norm", "1 fixed"),
                                      ("target", target), ("beam", beam),
                                      ("quadrature", None if points is None else "points %d" % points)])


def sample(name, n, fit="fixed", doppler=None, active="12C"):
    """A sample target n atoms/b thick (fit: "fixed", "fit", "fit prior <n> <sigma>")."""
    return record("target " + name, [("kind", "sample"), ("active", active), ("areal_density", "%g atoms/b %s" % (n, fit)),
                                     ("doppler", doppler)])


def beam(name="B"):
    return record("beam " + name, [("spread", "%.6f MeV lab" % SIGMA_LAB)])


def old_grid(lo, hi, step):
    """The energies of an old test line: energy += step from lo while energy <= hi (EData::MakePoints)."""
    grid, e = [], lo
    while e <= hi:
        grid.append(repr(e))
        e += step
    return "list " + " ".join(grid) + " MeV lab"


def free_width(line):
    f = line.split()
    # a level line: column 10 (1-based 11) is the width's fix flag; free the open channels with a nonzero width
    if len(f) > 30 and f[10] == "1" and float(f[11]) != 0.0 and float(f[2]) < 10.0:
        f[10] = "0"
        return "  ".join(f) + "\n"
    return line


def n12c_models(src, work):
    head, _, _, _, tail = split(open(src).read())
    os.makedirs(os.path.join(work, "data"), exist_ok=True)
    with open(os.path.join(work, "data", "d.dat"), "w") as f:
        for e in E_DATA:
            f.write("%.6f 0 0.5 0.05\n" % e)
    with open(os.path.join(work, "data", "g.dat"), "w") as f:
        for e in E_DATA:
            f.write("%.6f 0 1.0e-7 1.0e-8\n" % e)

    def model(name, data, extrapolations=(), targets=(), h=head):
        write(os.path.join(work, name + ".azr"), project(h, tail, data, extrapolations, targets))

    bare = [dataset("d", "t", "d.dat"), dataset("g", "g", "g.dat")]
    model("bare", bare)
    model("thin", [dataset("d", "T", "d.dat", "Sd"), dataset("g", "Y", "g.dat", "Sg")], targets=[sample("Sd", N), sample("Sg", N)])
    model("conv", [dataset("d", "T", "d.dat", "Sd", "B", 200), dataset("g", "Y", "g.dat", "Sg", "B", 200)],
          targets=[sample("Sd", N), sample("Sg", N), beam()])
    model("bconv", [dataset("d", "t", "d.dat", None, "B", 200), dataset("g", "g", "g.dat", None, "B", 200)], targets=[beam()])

    def extrapolations(grid, kinds, targets=(None, None)):
        out = []
        for k, (kind, tgt) in enumerate(zip(kinds, targets), 1):
            reaction = "12C(n,n)12C" if kind in "tT" else "12C(n,g)13C"
            quantity, unit = {"t": ("sigma_total", "b"), "g": ("sigma", "b"), "T": ("transmission", "1"),
                              "Y": ("neutron_yield", "1")}[kind]
            out.append(record("extrapolation x%d" % k, [("reaction", reaction), ("quantity", quantity), ("unit", unit),
                                                       ("energies", grid), ("target", tgt)]))
        return out
    model("truth", [], extrapolations(old_grid(0.20, 0.75, 0.0001), "tg"))
    # extrapolations (mode 3): bare and transformed curves
    grid = old_grid(0.30, 0.65, 0.001)
    model("xbare", [], extrapolations(grid, "tg"))
    model("xthin", [], extrapolations(grid, "TY", ("Sx1", "Sx2")), [sample("Sx1", N), sample("Sx2", N)])
    # gradient: thin (points without sub-points) and conv (with) in one project, data sets 3, 4 under the spread
    g = levels(head, free_width)
    model("grad", [dataset("d", "T", "d.dat", "S1"), dataset("g", "Y", "g.dat", "S2"),
                   dataset("d_2", "T", "d.dat", "S3", "B", 60), dataset("g_2", "Y", "g.dat", "S4", "B", 60)],
          targets=[sample("S1", N), sample("S2", N), sample("S3", N), sample("S4", N), beam()], h=g)
    # the fitted thickness: data at N_TRUE; the gradient project with every thickness a parameter
    model("thin18", [dataset("d", "T", "d.dat", "Sd"), dataset("g", "Y", "g.dat", "Sg")],
          targets=[sample("Sd", N_TRUE), sample("Sg", N_TRUE)])
    prior = "fit prior %g 0.01" % N
    model("grad_fit",
          [dataset("d", "T", "d.dat", "S1"), dataset("g", "Y", "g.dat", "S2"),
           dataset("d_2", "T", "d.dat", "S3", "B", 60), dataset("g_2", "Y", "g.dat", "S4", "B", 60),
           dataset("d_3", "T", "d.dat", "dop", None, 60), dataset("g_3", "Y", "g.dat", "dop", None, 60),
           dataset("d_4", "T", "d.dat", "S7", "B", 60), dataset("g_4", "Y", "g.dat", "S8", "B", 60)],
          targets=[sample("S1", N, "fit"), sample("S2", N, "fit"), sample("S3", N, prior), sample("S4", N, prior),
                   sample("dop", N, "fit", "free_gas 300 K fixed"), sample("S7", N, "fit", "free_gas 300 K fixed"),
                   sample("S8", N, "fit", "free_gas 300 K fixed"), beam()], h=g)
    refuse_t = {
        "refuse_t_sigma": ([dataset("d", "T", "d.dat", "S")], [sample("S", N, "fit prior %g -0.01" % N)]),
        "refuse_t_nosigma": ([dataset("d", "T", "d.dat", "S")], [sample("S", N, "fit prior %g" % N)]),
        "refuse_t_zero": ([dataset("g", "Y", "g.dat", "S")], [sample("S", 0, "fit")]),
        # a cross section on the sample of a transmission (the old group member without its own block)
        "refuse_t_member": ([dataset("d", "T", "d.dat", "S"), dataset("g", "g", "g.dat", "S")], [sample("S", N, "fit")]),
    }
    for name, (data, targets) in refuse_t.items():
        model(name, data, targets=targets)
    # refusals (one data set each)
    model("refuse_yield_total", [setfield(dataset("d", "Y", "d.dat", "S"), "reaction", "12C(n,n)12C")], targets=[sample("S", N)])
    # charged entrance: the same data set on the p + 12C model next to it
    p12c = os.path.join(os.path.dirname(src), "p12C.azr")
    ph, _, _, _, pt = split(open(p12c).read())
    charged = setfield(dataset("d", "T", "d.dat", "S"), "reaction", "12C(p,p)12C")
    write(os.path.join(work, "refuse_charged.azr"), project(ph, pt, [charged], (), [sample("S", N)]))


def neon_models(src, work):
    head, data, _, targets, tail = split(open(src).read())
    data, targets = records(data), records(targets)

    def only(k, rec=None):
        """Data set K (1-based; REC in its place) active, every other one inactive."""
        out = []
        for i, r in enumerate(data):
            if i + 1 == k:
                out.append(rec if rec is not None else r)
            else:
                out.append([r[0], "    %-14s %s" % ("active", "no")] + r[1:])
        return out

    def natural_sample(name, n):
        """The <abundances> target 'natural' as a sample n atoms/b thick."""
        nat = [r for r in targets if r[0] == "target natural"][0]
        rec = setfield(["target " + name] + nat[1:], "kind", "sample")
        return rec[:2] + ["    %-14s %s" % ("active", "20Ne"), "    %-14s %s" % ("areal_density", "%g atoms/b fixed" % n)] + rec[2:]

    def transmission(rec, tgt):
        return setfield(setfield(setfield(rec, "quantity", "transmission"), "unit", "1"), "target", tgt)

    def model(name, d, t, h=head):
        write(os.path.join(work, name + ".azr"), project(h, tail, d, (), t))

    model("nat", only(1), targets)
    model("natT", only(1, transmission(data[0], "natS")), targets + [natural_sample("natS", NNAT)])
    # the yield's sigma_t is the TOTAL, not the elastic: n + 20Ne with a free 1 keV width on every open reaction
    # channel (as tests/total_optical), a yield to alpha + 17O next to sigma_t, sigma(n,alpha) and the elastic
    k = [0]

    def react(l):
        f = l.split()
        if len(f) > 30 and f[5] in ("1", "3", "5", "7") and float(f[2]) > float(f[22]) + float(f[25]) + 0.05:
            f[10] = "0"
            f[11] = "-1000" if k[0] % 2 else "1000"
            k[0] += 1
            return " ".join(f) + "\n"
        return l
    rhead = levels(head, react)

    def vaughn(name, reaction, quantity, unit, tgt=None):
        return record("dataset " + name, [("reaction", reaction), ("quantity", quantity), ("unit", unit),
                                          ("file", "data/Vaughn_ntot.dat"), ("columns", "E_lab[MeV] _ value error"),
                                          ("energy_window", "0 2.33 MeV"), ("norm", "1 fixed"), ("target", tgt)])
    model("ybare", [vaughn("tot", "20Ne(n,n)20Ne", "sigma_total", "b"), vaughn("na", "20Ne(n,a)17O", "sigma", "b")], [], rhead)
    model("yelast", [vaughn("el", "20Ne(n,n)20Ne", "sigma", "b")], [], rhead)
    model("yreac", [vaughn("y", "20Ne(n,a)17O", "neutron_yield", "1", "S")], [sample("S", NNAT, active="20Ne")], rhead)
    model("refuse_natsub", only(7, transmission(data[6], "natS")), targets + [natural_sample("natS", NNAT)])
    natyield = vaughn("y", "20Ne(n,n)20Ne[1.63365]", "neutron_yield", "1", "natS")
    model("refuse_natyield", only(1, natyield), targets + [natural_sample("natS", NNAT)])


def fix_all(head):
    """Every level energy and width fixed (fields 4 and 11 of a <levels> line set to 1)."""
    def level(line):
        f = line.split()
        if len(f) > 12:
            f[3] = "1"
            f[10] = "1"
            return "  ".join(f) + "\n"
        return line
    return levels(head, level)


def fitdata(src, thin18, work):
    """The fit models on the thin18 values (T and Y0 at N_TRUE), every level fixed."""
    head, _, _, _, tail = split(open(src).read())
    head = fix_all(head)
    for name, dat in (("t18", "1_R=1"), ("g18", "1_R=2")):
        rows = []
        for line in open(os.path.join(thin18, "output", "AZUREOut_aa=%s.out" % dat)):
            f = line.split()
            if len(f) >= 4:
                rows.append(float(f[3]))
        if len(rows) != len(E_DATA):
            sys.exit("make_models fitdata: %d values in %s, want %d" % (len(rows), dat, len(E_DATA)))
        for suffix, err in (("", 0.01), ("w", 0.20)):
            with open(os.path.join(work, "data", name + suffix + ".dat"), "w") as f:
                for e, v in zip(E_DATA, rows):
                    f.write("%.6f 0 %.10e %.10e\n" % (e, v, err * v))

    models = {
        "fitT": ([dataset("d", "T", "t18.dat", "Sd")], [sample("Sd", N, "fit")]),
        "fitY": ([dataset("g", "Y", "g18.dat", "Sg")], [sample("Sg", N, "fit")]),
        "fitG": ([dataset("d", "T", "t18.dat", "S"), dataset("g", "Y", "g18.dat", "S")], [sample("S", N, "fit")]),
        "fitP": ([dataset("d", "T", "t18w.dat", "Sd")], [sample("Sd", N, "fit prior %g 0.01" % N)]),
    }
    for name, (data, targets) in models.items():
        write(os.path.join(work, name + ".azr"), project(head, tail, data, (), targets))


if __name__ == "__main__":
    if sys.argv[1] == "fitdata":
        fitdata(sys.argv[2], sys.argv[3], sys.argv[4])
        sys.exit(0)
    n12c, neon, work = sys.argv[1:4]
    os.makedirs(work, exist_ok=True)
    n12c_models(n12c, work)
    neon_models(neon, work)
