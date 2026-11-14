"""Turn reduced-width amplitudes into physical partial widths with the engine's own transformation.

AZURE2 fits reduced-width amplitudes but reports partial widths, and the step between them is not
``Gamma = 2 gamma^2 P``: it carries the level-shift normalisation ``1 + sum_c gamma_c^2 dS_c/dE`` (Brune), and
without the Brune parametrization the levels of a J-group mix.  ``Session.transform(x)`` is that step, as the
program writes ``<levels>`` from a fit; this script prints it beside the naive ``2 gamma^2 P``.

Print the widths at the start of the fit vector, or at a parameter file's values (from the repository's root)::

    python3 pyazr/examples/transform_widths.py tests/upgrade/T1/T1.azr [--params output/param.sav]

Map the samples of a chain onto physical values: a CSV with a header row whose columns are named after free
parameters (``Session.parameters.free.names``); a free parameter without a column keeps its start value::

    python3 pyazr/examples/transform_widths.py tests/upgrade/T1/T1.azr --chain samples.csv --show width_3_1

A naive ``2 gamma^2 P`` far above the engine's width means a strongly coupled level, where a width saturates at
``2 P_c / (dS_c/dE)``: a chain pressed against that ceiling says the data do not fix the coupling strength.
"""

import argparse
import os

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np  # noqa: E402

import pyazr  # noqa: E402


def level_table(s, x):
    """One line per free particle width: the amplitude, the engine's width and the naive ``2 gamma^2 P`` (eV)."""
    P = s.parameters
    phys = s.transform(x)
    pairs = {p.number: p for p in s.pairs}
    energy = {(p.jgroup, p.level): phys[i] for i, p in enumerate(P) if p.kind == "energy"}
    print(f"{'parameter':<12} {'E_x (MeV)':>10} {'pair':>4} {'l':>2} {'gamma (MeV^1/2)':>16} {'engine (eV)':>14} "
          f"{'2 gamma^2 P (eV)':>17}")
    for i, p in enumerate(P):
        if p.kind != "width" or pairs[p.pair].is_photon or p.free_index is None:
            continue
        ex = energy[(p.jgroup, p.level)]
        gamma = float(x[p.free_index])
        e_cm = ex - pairs[p.pair].threshold
        head = f"{p.name:<12} {ex:10.4f} {p.pair:4d} {p.L:2d} {gamma:16.6g} {abs(phys[i]):14.6g}"
        if e_cm <= 0:
            print(head + "  closed: an ANC in fm^-1/2")
            continue
        naive = 2e6 * gamma ** 2 * float(s.channel_functions(p.pair, e_cm, p.L)["P"])
        print(head + f" {naive:17.6g}")


def chain_values(s, path, show):
    """The 16th, 50th and 84th percentiles of the physical values of the parameters ``show`` over a chain."""
    P = s.parameters
    data = np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True))
    cols = [name for name in data.dtype.names if name in P.free.names]
    wanted = [P.names.index(P.by_name(name).name) for name in show]  # by_name names the close ones
    values = np.array([s.transform(P.put(P.x0, cols, [sample[c] for c in cols]))[wanted] for sample in data])
    print(f"{data.size} samples; percentiles 16, 50, 84 of the physical values (MeV, eV, fm^-1/2: Session.transform)")
    for j, name in enumerate(show):
        q = np.percentile(values[:, j], [16, 50, 84])
        print(f"{name:<12} {q[0]:14.6g} {q[1]:14.6g} {q[2]:14.6g}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", help="the format-3 project the parameters belong to")
    ap.add_argument("--params", help="a parameter file to take (param.sav, param.par)")
    ap.add_argument("--chain", help="samples: a CSV with a header row of free-parameter names")
    ap.add_argument("--show", default="", help="comma-separated parameter names to summarise over the chain")
    args = ap.parse_args()
    with pyazr.Session(args.azr, parameter_file=args.params) as s:
        if args.chain:
            chain_values(s, args.chain, [n.strip() for n in args.show.split(",") if n.strip()])
        else:
            level_table(s, s.parameters.x0)


if __name__ == "__main__":
    main()
