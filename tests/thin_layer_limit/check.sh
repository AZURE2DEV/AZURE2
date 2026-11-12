#!/usr/bin/env bash
#
# The thin limit of a target that loses less than 1 eV (README.md): a layer whose support is narrower than
# MeshZeroWidth() has one node, and its fold is the kernel's thin limit (Kernel::thinLimit, experiment/src/Kernels.cpp),
# not 0.  models/thin_limit.azr is calculated (menu 1) at 33 areal densities N from 1e13 to 1e17 atoms/cm^2 (8 per
# decade; the 1 eV boundary lies near 1.7e15 at E_lab = 0.6 MeV), and at every N and each of three energies E0:
#   1. the yield on the layer Y and the target-averaged cross section <sigma> are finite and not 0;
#   2. Y / (1e-24 N) and <sigma> follow sigma(E0), the bare cross section at the point, within T |dsigma/dE| / sigma
#      (T = epsilon(E0) N the lab energy loss; the depth average lies T/2 below E0, the factor 2 is margin) plus 1e-9
#      (the output's 11 digits); the slope from the bare sigma 1 keV either side;
#   3. at the densities with a single node (T_cm below 1 eV, the layer's grid 1.02 T) the two equal sigma(E0) to 1e-9:
#      the node is the point's own energy;
#   4. both kinds of support occur in the scan (it crosses the boundary), and neighbouring densities on the two sides of
#      it differ by less than the tolerance of 2. (continuity).
# Before the thin limit, Y and <sigma> were exactly 0 below the boundary (1. and 3. fail).
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

densities="$(python3 -I -c 'print(" ".join("%.6e" % (10 ** (13 + k / 8)) for k in range(33)))')"
for n in $densities; do
  d="$WORK/$n"
  mkdir -p "$d/output" "$d/checks"
  ln -s "$HERE/models/data" "$d/data"
  sed "s/areal_density  1e+14 /areal_density  $n /" "$HERE/models/thin_limit.azr" > "$d/m.azr"
  (cd "$d" && printf '1\n\n\n7\n' | $(guard_command 120) "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
done

python3 -I - "$WORK" $densities <<'PY'
import math, os, sys

work, densities = sys.argv[1], sys.argv[2:]
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

rows = []  # (N, i, single_layer, single_thin, rel_layer, rel_thin, tol)
for n in densities:
    N = float(n)
    path = os.path.join(work, n, "output", "AZUREOut_aa=1_R=2.out")
    if not os.path.exists(path):
        check("N = %s runs" % n, False, open(os.path.join(work, n, "run.log")).read()[-300:])
        continue
    thin, avg, bare = [[r[3] for r in b] for b in blocks(path)]
    for i, e0 in enumerate(E0S):
        s_lo, s0, s_hi = bare[3 * i:3 * i + 3]
        slope = abs(s_hi - s_lo) / 0.002 / s0          # |dsigma/dE_lab| / sigma, 1/MeV
        T = A0 * e0 ** A1 * N                          # lab MeV
        tol = T * slope + 1e-9
        y = thin[i] / (1e-24 * N)
        rows.append((N, i, F * MARGIN * T < ZERO, F * T < ZERO, abs(y / s0 - 1) if y else math.inf,
                     abs(avg[i] / s0 - 1) if avg[i] else math.inf, tol, y, avg[i]))

finite = all(math.isfinite(r[7]) and r[7] != 0 and math.isfinite(r[8]) and r[8] != 0 for r in rows)
check("1. Y and <sigma> finite and not 0 at all %d points" % len(rows), finite and len(rows) == 99,
      "%d of 99 points read" % len(rows))
worst = max(rows, key=lambda r: max(r[4], r[5]) / r[6])
check("2. Y/(1e-24 N) and <sigma> follow sigma(E0) within T |sigma'|/sigma + 1e-9",
      all(r[4] <= r[6] and r[5] <= r[6] for r in rows),
      "worst at N = %.3e, E0 = %g: %.2e / %.2e of a tolerance %.2e" % (worst[0], E0S[worst[1]], worst[4], worst[5], worst[6]))
single = [r for r in rows if r[2] or r[3]]
check("3. single node: sigma(E0) to 1e-9",
      single and all((r[4] <= 1e-9 if r[2] else True) and (r[5] <= 1e-9 if r[3] else True) for r in single),
      "%d points with one node (layer or thin target), largest %.2e" %
      (len(single), max([max(r[4] if r[2] else 0, r[5] if r[3] else 0) for r in single] or [math.inf])))
jumps = []
for i in range(3):
    seq = [r for r in rows if r[1] == i]
    for a, b in zip(seq, seq[1:]):
        for k, (sa, sb) in enumerate(((a[2], b[2]), (a[3], b[3]))):
            if sa != sb:  # the boundary lies between a and b
                va, vb = (a[7], b[7]) if k == 0 else (a[8], b[8])
                jumps.append((abs(vb / va - 1), max(a[6], b[6]), a[0], b[0], i, k))
check("4. the scan crosses the boundary, continuous there",
      len(jumps) == 6 and all(j[0] <= j[1] for j in jumps),
      "; ".join("%s E0 %g: %.2e (tol %.2e) between N %.3e and %.3e" % (("Y", "<sigma>")[j[5]], E0S[j[4]], j[0], j[1], j[2], j[3])
                for j in jumps) or "no crossing")
print("\n  N           E0    one node (layer, thin)   |Y/(1e-24 N sigma) - 1|   |<sigma>/sigma - 1|   tolerance")
for r in rows:
    if r[1] == 0:
        print("  %.3e  %.1f   %-5s %-5s               %.2e                  %.2e              %.2e" %
              (r[0], E0S[r[1]], r[2], r[3], r[4], r[5], r[6]))
if failures:
    print("FAIL thin_layer_limit (%d)" % len(failures))
    sys.exit(1)
print("PASS thin_layer_limit")
PY
