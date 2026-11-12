#!/usr/bin/env bash
#
# The thin limit of a target that loses less than 1 eV (README.md): a layer whose support is narrower than
# MeshZeroWidth() has one node, and its fold is the kernel's thin limit (Kernel::thinLimit, experiment/src/Kernels.cpp),
# not 0, with or without an energy shift on its data set (ShiftedFold.h).  models/thin_limit.azr is calculated (menu 1)
# in three series of runs:
#   A. no energy shift, 33 areal densities N from 1e13 to 1e17 atoms/cm^2 (8 per decade; the 1 eV boundary lies near
#      1.7e15 at E_lab = 0.6 MeV);
#   B. N = 1e14 (one node) with an energy shift s on the two folded data sets: s = 0, +-2 keV, +-20 keV fixed, and
#      +-20 keV varied without a prior (its reach 3 x the support's width: a planned support still below 1 eV), and
#      +20 keV and -2 keV varied with a prior of 10 eV (a planned support of many nodes: the depth integral, as above
#      1 eV; README.md: -20 keV is left out, its <sigma> is 6e-6 off; a prior of keV is out of reach: the planned
#      support keeps the step of the layer's own grid);
#   C. N from 1e14 to 1e16 (17 densities) with s = +20 keV fixed, and again with s = +20 keV varied without a prior.
# In every run the bare data set is calculated at E0 + s and 1 keV either side (its data file is written with the
# shifted energies, the data set itself unshifted): the unfolded engine at the shifted energy is the reference.  At every
# N and each of three energies E0:
#   1. the yield on the layer Y and the target-averaged cross section <sigma> are finite and not 0;
#   2. Y / (1e-24 N) and <sigma> follow sigma(E0 + s), the bare cross section, within T |dsigma/dE| / sigma (T =
#      epsilon(E0) N the lab energy loss; the depth average lies T/2 below E0 + s, the factor 2 is margin) plus 1e-9
#      (the output's 11 digits), the tolerance of 3. at a single node; the slope from the bare sigma 1 keV either side;
#   3. at the densities with a single node the two equal sigma(E0 + s) to 1e-9 for a fixed shift (the node is the
#      shifted point's own energy), within 7 T |dsigma/dE| / sigma + 1e-9 for a varied shift (the node is the upper edge
#      of the planned support, 7 times the support's width: ShiftedFold.h);
#   4. a series with densities on both sides of the boundary crosses it continuously: neighbouring densities on the
#      two sides differ by less than the tolerance of 2. (of 3. on the single-node side).
# Single node: the planned support (c.m.) narrower than 1e-6 MeV: F x 1.02 T (layer: the grid margin) or F T (thin
# target), times 7 for a varied shift.  Before the thin limit, Y and <sigma> were exactly 0 below the boundary in series
# A; before the thin limit of a shifted single node, every shifted single-node point of B and C but s = 0 fixed was NaN
# ("the lower edge of the moved support ... is outside the point's planned nodes").
#
#   ./tests/thin_layer_limit/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/thin_layer_limit.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

# run NAME N S SHIFT: the model at areal density N with the record "shift SHIFT" on the folded data sets (none when
# SHIFT is empty), the bare data at E0 + S (MeV)
run() {
  local d="$WORK/$1"
  mkdir -p "$d/output" "$d/checks" "$d/data"
  python3 -I - "$HERE/models" "$d" "$2" "$3" "$4" <<'PY'
import sys
models, d, n, s, shift = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
text, out, done = open(models + "/thin_limit.azr").read().replace("areal_density  1e+14 ", "areal_density  %s " % n), [], 0
for line in text.split("\n"):
    out.append(line)
    if shift and line.startswith("    norm ") and done < 2:  # data sets thin and thin_sigma, not bare
        out.append("    shift          " + shift)
        done += 1
open(d + "/m.azr", "w").write("\n".join(out))
open(d + "/data/e0.dat", "w").write(open(models + "/data/e0.dat").read())
with open(d + "/data/bare.dat", "w") as f:
    for line in open(models + "/data/bare.dat"):
        if line.startswith("#"):
            f.write(line)
        else:
            e, *rest = line.split()
            f.write(" ".join(["%.9g" % (float(e) + s)] + rest) + "\n")
PY
  (cd "$d" && printf '1\n\n\n7\n' | $(guard_command 120) "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
}

densities="$(python3 -I -c 'print(" ".join("%.6e" % (10 ** (13 + k / 8)) for k in range(33)))')"
shifted_densities="$(python3 -I -c 'print(" ".join("%.6e" % (10 ** (14 + k / 8)) for k in range(17)))')"
cases=()
for n in $densities; do run "A_$n" "$n" 0 ""; cases+=("A no shift:A_$n:$n:0:1"); done
for s in 0 0.002 -0.002 0.02 -0.02; do run "B_fixed_$s" 1e+14 "$s" "offset $s MeV fixed"; cases+=("B shifts at 1e14:B_fixed_$s:1e14:$s:1"); done
for s in 0.02 -0.02; do run "B_fit_$s" 1e+14 "$s" "offset $s MeV fit"; cases+=("B shifts at 1e14:B_fit_$s:1e14:$s:7"); done
# a varied shift with a 10 eV prior: a planned support of 3 x 2 x 10 eV, many nodes (the ladder of the layer's support)
for s in 0.02 -0.002; do run "B_prior_$s" 1e+14 "$s" "offset $s MeV fit prior $s 1e-05"; cases+=("B shifts at 1e14:B_prior_$s:1e14:$s:inf"); done
for n in $shifted_densities; do
  run "Cfixed_$n" "$n" 0.02 "offset 0.02 MeV fixed"; cases+=("C fixed +20 keV:Cfixed_$n:$n:0.02:1")
  run "Cfit_$n" "$n" 0.02 "offset 0.02 MeV fit"; cases+=("C varied +20 keV:Cfit_$n:$n:0.02:7")
done

python3 -I - "$WORK" "${cases[@]}" <<'PY'
import math, os, sys

work, cases = sys.argv[1], [c.split(":") for c in sys.argv[2:]]
A0, A1 = 4.5e-22, -0.6              # epsilon = a0 x^a1, lab MeV cm^2 per active nucleus (models/thin_limit.azr)
F = 15.000109 / (1.007825 + 15.000109)  # c.m./lab factor of p + 15N
ZERO = 1.0e-6                       # MeshZeroWidth(), MeV
MARGIN = 1.02                       # TargetEffect::LayerGridMargin()
E0S = [0.6, 0.7, 0.8]
failures = []

def blocks(path):
    out, cur = [], []
    for line in open(path):
        if line.strip():
            cur.append([float(x) for x in line.split()])
        elif cur:
            out.append(cur); cur = []
    if cur:
        out.append(cur)
    return out

def check(name, ok, detail):
    print("  %s  %s: %s" % ("ok  " if ok else "FAIL", name, detail))
    if not ok:
        failures.append(name)

rows = []  # one per (case, E0)
for series, name, n, s, reach in cases:
    N, s, reach = float(n), float(s), float(reach)  # reach: the planned support over the support (inf: many nodes)
    path = os.path.join(work, name, "output", "AZUREOut_aa=1_R=2.out")
    if not os.path.exists(path):
        check("%s runs" % name, False, open(os.path.join(work, name, "run.log")).read()[-300:])
        continue
    thin, avg, bare = [[r[3] for r in b] for b in blocks(path)]
    for i, e0 in enumerate(E0S):
        s_lo, s0, s_hi = bare[3 * i:3 * i + 3]
        slope = abs(s_hi - s_lo) / 0.002 / s0          # |dsigma/dE_lab| / sigma at E0 + s, 1/MeV
        T = A0 * e0 ** A1 * N                          # lab MeV
        y = thin[i] / (1e-24 * N)
        rel = lambda v: abs(v / s0 - 1) if v and math.isfinite(v) else math.inf
        rows.append(dict(series=series, name=name, N=N, i=i, tol=T * slope + 1e-9, y=y, avg=avg[i],
                         one_layer=reach * F * MARGIN * T < ZERO, one_thin=reach * F * T < ZERO,
                         single_tol=1e-9 if reach == 1 else reach * T * slope + 1e-9,
                         rel_layer=rel(y), rel_thin=rel(avg[i])))
        r = rows[-1]
        # a varied shift's single node lies up to the planned width (7 T) above the moved support: that is the tolerance
        r["tol_layer"] = max(r["tol"], r["single_tol"]) if r["one_layer"] else r["tol"]
        r["tol_thin"] = max(r["tol"], r["single_tol"]) if r["one_thin"] else r["tol"]

for series in sorted(set(c[0] for c in cases)):
    rs = [r for r in rows if r["series"] == series]
    want = 3 * sum(1 for c in cases if c[0] == series)
    print("== %s" % series)
    finite = all(math.isfinite(r["y"]) and r["y"] != 0 and math.isfinite(r["avg"]) and r["avg"] != 0 for r in rs)
    check("%s 1. Y and <sigma> finite and not 0" % series, finite and len(rs) == want, "%d of %d points read" % (len(rs), want))
    if not rs:
        continue
    worst = max(rs, key=lambda r: max(r["rel_layer"] / r["tol_layer"], r["rel_thin"] / r["tol_thin"]))
    check("%s 2. Y/(1e-24 N) and <sigma> follow sigma(E0 + s) within T |sigma'|/sigma + 1e-9 (7 T for a varied shift's "
          "single node)" % series,
          all(r["rel_layer"] <= r["tol_layer"] and r["rel_thin"] <= r["tol_thin"] for r in rs),
          "worst %s, E0 = %g: %.2e / %.2e of a tolerance %.2e / %.2e" %
          (worst["name"], E0S[worst["i"]], worst["rel_layer"], worst["rel_thin"], worst["tol_layer"], worst["tol_thin"]))
    single = [r for r in rs if r["one_layer"] or r["one_thin"]]
    bad = [r for r in single if (r["one_layer"] and r["rel_layer"] > r["single_tol"]) or
                                (r["one_thin"] and r["rel_thin"] > r["single_tol"])]
    check("%s 3. single node: sigma(E0 + s) to 1e-9 (fixed shift), within 7 T |sigma'|/sigma + 1e-9 (varied)" % series,
          bool(single) and not bad,
          "%d points with one node, largest deviation %.2e%s" %
          (len(single), max([max(r["rel_layer"] if r["one_layer"] else 0, r["rel_thin"] if r["one_thin"] else 0)
                             for r in single] or [math.inf]),
           "" if not bad else "; e.g. %s E0 %g" % (bad[0]["name"], E0S[bad[0]["i"]])))
    if not series.startswith("B"):  # B has one density
        jumps = []
        for i in range(3):
            seq = sorted([r for r in rs if r["i"] == i], key=lambda r: r["N"])
            for a, b in zip(seq, seq[1:]):
                for k, key in enumerate(("one_layer", "one_thin")):
                    if a[key] != b[key]:  # the boundary lies between a and b
                        va, vb = (a["y"], b["y"]) if k == 0 else (a["avg"], b["avg"])
                        t = ("tol_layer", "tol_thin")[k]
                        jumps.append((abs(vb / va - 1), max(a[t], b[t]), a["N"], b["N"], i, k))
        check("%s 4. the scan crosses the boundary, continuous there" % series,
              len(jumps) == 6 and all(j[0] <= j[1] for j in jumps),
              "; ".join("%s E0 %g: %.2e (tol %.2e) between N %.3e and %.3e" % (("Y", "<sigma>")[j[5]], E0S[j[4]], j[0],
                                                                              j[1], j[2], j[3]) for j in jumps)
              or "no crossing")
    print("\n  case                     E0    one node (layer, thin)   |Y/(1e-24 N sigma) - 1|   |<sigma>/sigma - 1|   tolerance")
    for r in rs:
        if r["i"] == 0:
            print("  %-24s %.1f   %-5s %-5s               %.2e                  %.2e              %.2e" %
                  (r["name"], E0S[r["i"]], r["one_layer"], r["one_thin"], r["rel_layer"], r["rel_thin"], r["tol"]))
    print()
if failures:
    print("FAIL thin_layer_limit (%d)" % len(failures))
    sys.exit(1)
print("PASS thin_layer_limit")
PY
