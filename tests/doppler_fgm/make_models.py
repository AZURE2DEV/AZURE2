#!/usr/bin/env python3
"""Models for tests/doppler_fgm, built at run time (check.sh) into a work directory (project format 3).

From tests/reaction_rate/models/n12C.azr (n + 12C -> 13C, elastic and capture, external capture off), with the
s-wave 1/2+ resonance at E_cm = 0.4247 MeV narrowed from Gamma_n = 34 keV to 200 eV, so that free-gas Doppler
broadening at 300 K (Doppler width about 60 eV at E_lab = 0.46 MeV) changes the cross section by tens of percent.
Data sets t (sigma_t: 12C(n,n)12C as sigma_total) and g (sigma_gamma: 12C(n,g)13C) on the data grid:

  bare     no target, no beam
  dop      a target with "doppler free_gas T K" alone (sigma_D at the point)
  dopg     the Doppler target and a beam with a Gaussian spread SIGMA_LAB (Doppler first, then the Gaussian)
  gauss    the beam alone (to show the Doppler step matters under the Gaussian too)
  dopT     transmission (of sigma_t) and neutron_yield (of sigma_gamma) of a sample N atoms/b thick, the Doppler
           step on the sample's target record
  dopgT    the same with the beam too: Doppler, then the transform, then the Gaussian
  truth    extrapolations sigma_t and sigma_gamma on a 1 eV lab grid, no effect (mode 3)
  grad     dop (data sets 1, 2), dopgT (3, 4) and dopT (5, 6) in one project, for the gradient check
  grad_fit the same with fitted temperatures: targets D0 (1, 2) and D1 (5, 6) "fit grid_max 900", target D2
           "fit prior 300 20" (3, 4); D0.doppler, D1.doppler, D2.doppler free with the widths
  dop450   dop at T_TRUE = 450 K: its sigma_t and sigma_gamma become the data of the fits below
  refuse_* Doppler fields the engine must refuse

A target record holds the Doppler temperature and, for transmission and yield, the sample's areal density (kind
sample); a cross section takes a target of kind none or layer, so the plain data sets of grad have a target of their
own (the old model shared one temperature, doppler_1, between them and the transmission lines 5, 6).

    python3 make_models.py <n12C.azr> <workdir>

and, once dop450 has been calculated (mode 1),

    python3 make_models.py fitdata <n12C.azr> <dop450 dir> <workdir>

writes data/t450.dat, data/g450.dat (the dop450 values, 1 % errors), data/t450w.dat, data/g450w.dat (the same
values, 20 % errors) and, with every level parameter fixed,

  fit       the 1 % data, "doppler free_gas 300 K fit grid_max 900": start at 300 K, no prior, grids for 900 K
  fit_prior the 20 % data (T alone known to about 60 K), "fit prior 300 20": a 20 K prior about 300 K (grids for
            max(300 + 60, 1.21 * 300) = 363 K), so the fit lands between 300 K and 363 K
"""
import os
import re
import sys

N = 0.15            # sample thickness, atoms/barn
T_K = 300.0         # effective temperature, K
SIGMA_LAB = 1.0e-4  # Gaussian resolution, MeV lab (100 eV)
E_DATA = [round(0.4589 + 0.00005 * k, 7) for k in range(63)]   # 0.45890 .. 0.46200 MeV lab
GAMMA_N = "200"     # eV, the narrowed 1/2+ neutron width
T_TRUE = 450.0      # K, the temperature of the fit data


def split(text):
    """(the file up to <data>, the file after </targets>): configuration, levels and parameter settings."""
    a, b = text.find("<data>\n"), text.find("</targets>\n")
    if a < 0 or b < 0:
        sys.exit("make_models: unexpected .azr layout (format 3 with <data> ... </targets>)")
    return text[:a], text[b + len("</targets>\n"):]


def levels(head, edit):
    return re.sub(r"(<levels>\n)(.*?)(</levels>)",
                  lambda m: m.group(1) + "".join(edit(l) for l in m.group(2).splitlines(True)) + m.group(3),
                  head, flags=re.S)


def narrow(head):
    def level(line):
        f = line.split()
        if len(f) > 12 and f[0] == "0.5" and f[1] == "1" and f[11] == "34000":
            f[11] = GAMMA_N
            return "  ".join(f) + "\n"
        return line
    out = levels(head, level)
    if out.count(" " + GAMMA_N + " ") < 1:
        sys.exit("make_models: the 1/2+ level was not found")
    return out


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


def record(header, fields):
    """A format-3 record: the header and its "key value" lines (a value None is left out), then a blank line."""
    return header + "\n" + "".join("    %-14s %s\n" % kv for kv in fields if kv[1] is not None) + "\n"


def dataset(name, kind, dat, target=None, beam=None, points=None):
    """kind: t (sigma_t), g (sigma_gamma), T (transmission of sigma_t), Y (neutron_yield of sigma_gamma)."""
    reaction = "12C(n,n)12C" if kind in "tT" else "12C(n,g)13C"
    quantity, unit = {"t": ("sigma_total", "b"), "g": ("sigma", "b"), "T": ("transmission", "1"),
                      "Y": ("neutron_yield", "1")}[kind]
    return record("dataset " + name, [("reaction", reaction), ("quantity", quantity), ("unit", unit),
                                      ("file", "data/" + dat), ("columns", "E_lab[MeV] _ value error"),
                                      ("energy_window", "0.4585 0.4625 MeV"), ("norm", "1 fixed"),
                                      ("target", target), ("beam", beam),
                                      ("quadrature", None if points is None else "points %d" % points)])


def target(name, doppler, sample=False):
    """A target with a Doppler field (the text after the key) and, for transmission and yield, the sample."""
    fields = [("kind", "sample"), ("active", "12C"), ("areal_density", "%g atoms/b fixed" % N)] if sample else [("kind", "none")]
    return record("target " + name, fields + [("doppler", doppler)])


BEAM = record("beam B", [("spread", "%g MeV lab" % SIGMA_LAB)])


def project(head, tail, data, targets="", extrapolations=""):
    return (head + "<data>\n" + data + "</data>\n<extrapolations>\n" + extrapolations + "</extrapolations>\n<targets>\n" +
            targets + "</targets>\n" + tail)


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def plain(T):
    return "free_gas %g K fixed" % T


def main(src, work):
    head, tail = split(open(src).read())
    head = narrow(head)
    os.makedirs(os.path.join(work, "data"), exist_ok=True)
    with open(os.path.join(work, "data", "t.dat"), "w") as f:
        for e in E_DATA:
            f.write("%.7f 0 0.5 0.05\n" % e)
    with open(os.path.join(work, "data", "g.dat"), "w") as f:
        for e in E_DATA:
            f.write("%.7f 0 1.0e-6 1.0e-7\n" % e)

    # 1000 sub-points under the Gaussian: its window is wider than the Doppler one, and the 200 eV resonance must
    # still be resolved by the grid (400 gives 4e-4 at the dip of sigma_t, 1600 gives 3e-5).
    def pair(kinds, tgt=None, beam=None, points=None, suffix=""):
        return (dataset("t" + suffix, kinds[0], "t.dat", tgt, beam, points) +
                dataset("g" + suffix, kinds[1], "g.dat", tgt, beam, points))

    def model(name, data, targets=""):
        write(os.path.join(work, name + ".azr"), project(head, tail, data, targets))

    model("bare", pair("tg"))
    model("dop", pair("tg", "D", None, 400), target("D", plain(T_K)))
    model("dopg", pair("tg", "D", "B", 1000), target("D", plain(T_K)) + BEAM)
    model("gauss", pair("tg", None, "B", 1000), BEAM)
    model("dopT", pair("TY", "D", None, 400), target("D", plain(T_K), sample=True))
    model("dopgT", pair("TY", "D", "B", 1000), target("D", plain(T_K), sample=True) + BEAM)
    # the 1 eV grid of the old test lines (energy += step from 0.4560 while energy <= 0.4650)
    grid, e = [], 0.4560
    while e <= 0.4650:
        grid.append(repr(e))
        e += 0.000001
    energies = "list " + " ".join(grid) + " MeV lab"
    truth = "".join(record("extrapolation " + n, [("reaction", r), ("quantity", q), ("unit", "b"), ("energies", energies)])
                    for n, r, q in (("x1", "12C(n,n)12C", "sigma_total"), ("x2", "12C(n,g)13C", "sigma")))
    write(os.path.join(work, "truth.azr"), project(head, tail, "", "", truth))
    # one project, 120 sub-points: Doppler alone (1, 2), Doppler + Gaussian on a sample (3, 4), Doppler on a sample (5, 6)
    grad = pair("tg", "D0", None, 120) + pair("TY", "D2", "B", 120, "_2") + pair("TY", "D1", None, 120, "_3")
    model("grad", grad, target("D0", plain(T_K)) + target("D1", plain(T_K), sample=True) +
          target("D2", plain(T_K), sample=True) + BEAM)
    model("grad_fit", grad, target("D0", "free_gas %g K fit grid_max 900" % T_K) +
          target("D1", "free_gas %g K fit grid_max 900" % T_K, sample=True) +
          target("D2", "free_gas %g K fit prior %g 20" % (T_K, T_K), sample=True) + BEAM)
    model("dop450", pair("tg", "D", None, 400), target("D", plain(T_TRUE)))
    refuse = {
        "refuse_zero": "free_gas 0 K fixed",
        "refuse_missing": "free_gas K fixed",
        # the fit form: a negative or missing prior width, grid_max not above T, T <= 0
        "refuse_sigma": "free_gas 300 K fit prior 300 -5",
        "refuse_nosigma": "free_gas 300 K fit prior 300",
        "refuse_tmax": "free_gas 300 K fit grid_max 200",
        "refuse_negT": "free_gas -300 K fit grid_max 900",
    }
    for name, doppler in refuse.items():
        model(name, pair("tg", "D", None, 400), target("D", doppler))


def fitdata(src, dop450, work):
    """The fit models on the dop450 values (sigma_t, sigma_gamma at T_TRUE), 1 % errors, every level fixed."""
    head, tail = split(open(src).read())
    head = fix_all(narrow(head))
    for name, dat in (("t450.dat", "1_R=1"), ("g450.dat", "1_R=2")):
        rows = []
        for line in open(os.path.join(dop450, "output", "AZUREOut_aa=%s.out" % dat)):
            f = line.split()
            if len(f) >= 4:
                rows.append(float(f[3]))
        if len(rows) != len(E_DATA):
            sys.exit("make_models fitdata: %d values in %s, want %d" % (len(rows), dat, len(E_DATA)))
        for suffix, err in (("", 0.01), ("w", 0.20)):
            with open(os.path.join(work, "data", name.replace(".dat", suffix + ".dat")), "w") as f:
                for e, v in zip(E_DATA, rows):
                    f.write("%.7f 0 %.10e %.10e\n" % (e, v, err * v))

    def fit(t, g, doppler):
        return project(head, tail, dataset("t", "t", t, "D", None, 400) + dataset("g", "g", g, "D", None, 400),
                       target("D", doppler))
    write(os.path.join(work, "fit.azr"), fit("t450.dat", "g450.dat", "free_gas %g K fit grid_max 900" % T_K))
    write(os.path.join(work, "fit_prior.azr"), fit("t450w.dat", "g450w.dat", "free_gas %g K fit prior %g 20" % (T_K, T_K)))


if __name__ == "__main__":
    if sys.argv[1] == "fitdata":
        fitdata(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        os.makedirs(sys.argv[2], exist_ok=True)
        main(sys.argv[1], sys.argv[2])
