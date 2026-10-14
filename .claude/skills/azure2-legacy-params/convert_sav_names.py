#!/usr/bin/env python3
"""Rename an old-format AZURE2 parameter file (param.sav/.par/.fit) to the
current parameter names, and optionally check the result against the model.

Builds before dcfdb34 (2025-08-01) named level parameters
    j=<J>_la=<L>_energy      and   j=<J>_la=<L>_ch=<C>_rwa
Current builds (CNuc::FillMnParams) name them
    energy_<N>               and   width_<N>_<C>
N counts levels in J-group order (j ascending, then la ascending -- the same
double loop that wrote the old names); C is unchanged.  The values (level
energy, reduced-width amplitude) mean the same thing.  segment_* lines are
named the same in both and are copied untouched.

usage:
  convert_sav_names.py OLD NEW                  convert
  convert_sav_names.py OLD NEW --check PAR      convert, then compare NEW with
                                                PAR, a param.par written by
                                                the CURRENT build from the
                                                same .azr (mode 1, blank
                                                parameter file)
Exit status: 0 ok, 1 nothing to convert, 2 check failed.
"""
import argparse
import re
import sys

E_RE = re.compile(r"^j=(\d+)_la=(\d+)_energy$")
W_RE = re.compile(r"^j=(\d+)_la=(\d+)_ch=(\d+)_rwa$")


def read_vals(path):
    vals = {}
    for line in open(path):
        t = line.split()
        if len(t) >= 2 and not t[0].startswith("#") and t[0] != "parametrization":
            try:
                vals[t[0]] = float(t[1])
            except ValueError:
                pass
    return vals


def convert(src, dst):
    lines = open(src).read().splitlines()
    pairs = set()
    for line in lines:
        t = line.split()
        m = t and (E_RE.match(t[0]) or W_RE.match(t[0]))
        if m:
            pairs.add((int(m.group(1)), int(m.group(2))))
    if not pairs:
        return 0
    index = {p: n for n, p in enumerate(sorted(pairs), start=1)}
    out = []
    for line in lines:
        t = line.split()
        if not t:
            out.append(line)
            continue
        m = E_RE.match(t[0])
        if m:
            name = "energy_%d" % index[(int(m.group(1)), int(m.group(2)))]
        else:
            m = W_RE.match(t[0])
            name = ("width_%d_%s" % (index[(int(m.group(1)), int(m.group(2)))], m.group(3))
                    if m else t[0])
        # keep the value columns byte-for-byte; right-align the name as AZURE2 does
        rest = line[len(line) - len(line.lstrip()) + len(t[0]):]
        out.append("%20s" % name + rest)
    open(dst, "w").write("\n".join(out) + "\n")
    return len(index)


def check(new, par):
    a, b = read_vals(new), read_vals(par)
    lvl = lambda d: {k for k in d if k.startswith(("energy_", "width_"))}
    only_new, only_par = sorted(lvl(a) - lvl(b)), sorted(lvl(b) - lvl(a))
    zero = sorted(k for k in lvl(a) & lvl(b) if k.startswith("width_")
                  and (a[k] == 0.0) != (b[k] == 0.0))
    print("level parameters: %d in converted file, %d in model, %d in common"
          % (len(lvl(a)), len(lvl(b)), len(lvl(a) & lvl(b))))
    if only_new:
        print("  only in converted file:", only_new[:10], "..." if len(only_new) > 10 else "")
    if only_par:
        print("  only in model:", only_par[:10], "..." if len(only_par) > 10 else "")
    # a width that is exactly zero in one and not the other means the levels
    # or channels were paired up differently (structure changed since the fit)
    if zero:
        print("  zero/nonzero width mismatch:", zero[:10], "..." if len(zero) > 10 else "")
    ok = not (only_new or only_par or zero)
    print("check:", "OK" if ok else "FAILED")
    return ok


ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
ap.add_argument("old")
ap.add_argument("new")
ap.add_argument("--check", metavar="PAR")
args = ap.parse_args()

n = convert(args.old, args.new)
if n == 0:
    print("%s has no j=..._la=... names: already current format, nothing written" % args.old)
    sys.exit(1)
print("%d levels renamed; wrote %s" % (n, args.new))
if args.check and not check(args.new, args.check):
    sys.exit(2)
