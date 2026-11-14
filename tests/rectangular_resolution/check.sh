#!/usr/bin/env bash
#
# Rectangular energy resolution (the table kernel "rectangular", experiment/src/KernelRectangular.cpp): a point folded
# with a window of full width w = 20 keV lab must be the average of the unfolded cross section over
# [E0 - f w/2, E0 + f w/2] (c.m., f = E_cm/E_lab of the pair).
#
#   models/bare.azr : the level scheme of tests/convolution_window (a narrow resonance at 0.8508 MeV lab on a smooth
#                     background), unfolded on a 0.05 keV lab grid;
#   models/rect.azr : the points E_res + k w, k = -1 ... 1, each with `detector det_rect` (energy_resolution
#                     rectangular width 20 keV).  Both written by make_models.py.
#   (a) every folded point within 1e-3 of the trapezoidal average of the bare output over its window (the truth needs
#       no engine kernel: the bare run has none);
#   (b) the window: a point a full width from the peak sees little of it, the point on the peak averages it down.
#   The width as a fit parameter (kind `apparatus`, docs/upgrade/apparatus_parameters.md), models/fit.azr: two data
#   sets d1, d2 (two effects: their quadratures differ) name one detector det_fit, `width 14 keV fit`; the pseudo-data
#   (models/data/, make_models.py --data) are window averages of the bare output at w = 20 keV with 1 % noise.
#   (c) one parameter det_fit.energy_resolution.width for both data sets; the fit recovers 20 keV within 3 sigma;
#   (d) without noise (data written here from the bare run): w within 0.05 sigma of 20 keV, every data set's chi^2 < 0.1
#       (each data set is folded with the fitted width);
#   (e) with a prior (width 21 keV fit prior 21 0.05): the run's objective minus the data's chi^2 is ((w - 21)/0.05)^2,
#       chiSquared.out's Total-Kernel-Parameter-Chi-Squared too, and w lies between the prior centre and the data's 20;
#   (f) refusals: grid_max below the value and a width 0
#       (E_OUT_OF_RANGE), no mode (E_NUISANCE_MODE), `profiled` (E_MODE_NOT_ALLOWED), the width beside a target.
#   (g) the triangular resolution (experiment/src/KernelTriangular.cpp), base width w fitted the same way: noise-free
#       data made here as triangle-weighted averages of the bare output at w = 20 keV, models/fit.azr with
#       `triangular width 14 keV fit`: w within 0.05 sigma of 20 keV.
#
#   ./tests/rectangular_resolution/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/rectangular_resolution.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
OUT=output/AZUREOut_aa=1_R=3.extrap

run() {  # MODEL: mode 3 in $WORK/MODEL
  local dir="$WORK/$1"; mkdir -p "$dir/output" "$dir/checks"; cp "$HERE/models/$1.azr" "$dir/t.azr"
  (cd "$dir" && run_limited 600 "$AZURE2_BIN" extrapolate --no-gui t.azr > run.log 2>&1)
  [ -s "$dir/$OUT" ] || { echo "  FAIL  $1: no output"; tail -5 "$dir/run.log"; echo "rectangular_resolution: FAIL"; exit 1; }
}
fit() {  # NAME DATA [SED]: a fit (menu 2, one thread) of models/fit.azr, edited by SED, with the data of DATA
  local dir="$WORK/$1"; mkdir -p "$dir/output" "$dir/checks"; cp -r "$2" "$dir/data"
  sed -e "${3:-s/^//}" "$HERE/models/fit.azr" > "$dir/t.azr"
  (cd "$dir" && OMP_NUM_THREADS=1 run_limited 900 "$AZURE2_BIN" fit --no-gui t.azr > run.log 2>&1)
}
FAILS=0
refused() {  # WHAT CODE SED: --validate of models/fit.azr edited by SED prints CODE
  local dir="$WORK/refused"; mkdir -p "$dir"; cp -r "$HERE/models/data" "$dir/"
  sed -e "$3" "$HERE/models/fit.azr" > "$dir/t.azr"
  local msg; msg="$(cd "$dir" && "$AZURE2_BIN" --validate t.azr 2>&1)"
  if grep -q "$2" <<< "$msg"; then echo "  ok    (f) $1 is refused ($2)"
  else echo "  FAIL  (f) $1: no $2"; FAILS=$((FAILS + 1)); fi
}

echo "rectangular_resolution: $AZURE2_BIN"
run bare; run rect
python3 - "$WORK/bare/$OUT" "$WORK/rect/$OUT" <<'PY' || FAILS=$((FAILS + 1))
import sys

def rows(path):
    return [list(map(float, l.split())) for l in open(path) if len(l.split()) >= 4]

bare = rows(sys.argv[1])  # E_cm, ..., sigma (column 4)
rect = rows(sys.argv[2])
e = [r[0] for r in bare]
s = [r[3] for r in bare]
w_lab = 0.020
labs = [0.8508 + k * w_lab for k in (-1.0, -0.6, -0.3, 0.0, 0.25, 0.5, 1.0)]

def sigma(x):  # linear between the bare nodes
    for i in range(len(e) - 1):
        if e[i] <= x <= e[i + 1]:
            t = (x - e[i]) / (e[i + 1] - e[i])
            return (1 - t) * s[i] + t * s[i + 1]
    raise ValueError(x)

def average(lo, hi):  # trapezoid on the bare nodes inside [lo, hi] plus the two edges
    xs = [lo] + [x for x in e if lo < x < hi] + [hi]
    ys = [sigma(x) for x in xs]
    return sum(0.5 * (ys[i] + ys[i + 1]) * (xs[i + 1] - xs[i]) for i in range(len(xs) - 1)) / (hi - lo)

failures, worst = 0, 0.0
for r, lab in zip(rect, labs):
    f = r[0] / lab
    truth = average(r[0] - 0.5 * f * w_lab, r[0] + 0.5 * f * w_lab)
    d = abs(r[3] - truth) / truth
    worst = max(worst, d)
    print(f"        E_lab {lab:.4f}: folded {r[3]:.6e}  truth {truth:.6e}  rel {d:.1e}")
ok = worst <= 1e-3
print(("  ok    " if ok else "  FAIL  ") + f"(a) folded points within 1e-3 of the window average of the bare output (worst {worst:.1e})")
failures += not ok
peak, edge = rect[3][3], rect[0][3]
ok = peak > 10 * edge and peak < max(s)
print(("  ok    " if ok else "  FAIL  ") + f"(b) on the peak {peak:.3e} b (bare peak {max(s):.3e}), one width away {edge:.3e} b")
failures += not ok
sys.exit(1 if failures else 0)
PY

# ---- the width as a fit parameter
W='s/width 14 keV fit$/'
fit noisy "$HERE/models/data"
python3 "$HERE/make_models.py" --data "$WORK/bare/$OUT" "$WORK/exact-data" 0 > /dev/null
fit exact "$WORK/exact-data"
fit prior "$HERE/models/data" "${W}width 21 keV fit prior 21 0.05/"
python3 - "$WORK" <<'PY' || FAILS=$((FAILS + 1))
import sys
work = sys.argv[1]
W_TRUE = 0.020  # MeV

def sav(name):
    return {l.split()[0]: (float(l.split()[1]), float(l.split()[2])) for l in open(f"{work}/{name}/output/param.sav")
            if l.strip() and not l.startswith("#")}

def chi(name, key):
    for l in open(f"{work}/{name}/output/chiSquared.out"):
        t = l.split()
        if key in t:
            return float(t[t.index(key) + 1])
    return None

failures = 0
def report(ok, text):
    global failures
    print(("  ok    " if ok else "  FAIL  ") + text)
    failures += not ok

p = sav("noisy")
widths = [k for k in p if k.endswith(".width")]
w, dw = p.get("det_fit.energy_resolution.width", (float("nan"), float("nan")))
pull = (w - W_TRUE) / dw if dw > 0 else float("inf")
report(widths == ["det_fit.energy_resolution.width"] and abs(pull) <= 3,
       f"(c) {len(widths)} width parameter(s) for d1 and d2; w = {w * 1e3:.4f} +- {dw * 1e3:.4f} keV, pull "
       f"{pull:+.2f} (truth 20 keV, within 3 sigma)")
q = sav("exact")
w0, dw0 = q.get("det_fit.energy_resolution.width", (float("nan"), float("nan")))
segs = [float(l.split(",")[1]) for l in open(f"{work}/exact/output/chiSquared.out") if l.split(",")[0].strip().isdigit()]
pull0 = (w0 - W_TRUE) / dw0 if dw0 > 0 else float("inf")
report(abs(pull0) <= 0.05 and len(segs) == 2 and max(segs) < 0.1,
       f"(d) noise-free data: w = {w0 * 1e3:.6f} keV, {pull0:+.3f} sigma from 20 keV; chi^2 of d1, d2 "
       + ", ".join(f"{c:.2e}" for c in segs))
r = sav("prior")
wp = r.get("det_fit.energy_resolution.width", (float("nan"),))[0]
data = chi("prior", "Total-Chi-Squared:")
line = chi("prior", "Total-Kernel-Parameter-Chi-Squared:")
total = float([l.split(":")[1] for l in open(f"{work}/prior/run.log") if l.startswith("Total Chi-Squared:")][-1])
penalty = ((wp - 0.021) / 0.05e-3) ** 2
report(line is not None and abs(total - data - penalty) <= 1e-3 * max(1.0, penalty) and abs(line - penalty) <= 1e-3 * max(1.0, penalty)
       and min(w, 0.021) <= wp <= 0.021,
       f"(e) prior 21 +- 0.05 keV: w = {wp * 1e3:.4f} keV; objective - chi^2 = {total - data:.4f}, "
       f"((w - 21 keV)/0.05 keV)^2 = {penalty:.4f}, chiSquared.out {line}")
sys.exit(1 if failures else 0)
PY
refused "grid_max below the value" E_OUT_OF_RANGE "${W}width 14 keV fit grid_max 10/"
refused "a width of 0" E_OUT_OF_RANGE "${W}width 0 keV fit/"
refused "a width without its mode" E_NUISANCE_MODE "${W}width 14 keV/"
refused "a profiled width" E_MODE_NOT_ALLOWED "${W}width profiled/"
# beside a target: refused by the effect (TargetEffect::SelectKernels), before any parameter is made
dir="$WORK/target"; mkdir -p "$dir/output" "$dir/checks"; cp -r "$HERE/models/data" "$dir/"
sed -e 's/^    detector       det_fit$/    detector       det_fit\n    target         T1/' "$HERE/models/fit.azr" |
  sed -e 's|^</targets>|target T1\n    kind           layer\n    active         16O\n    areal_density  1e17 atoms/cm2 fixed\n    stopping       formula "1e-15" eV*cm2/atom\n    straggling     none\n</targets>|' > "$dir/t.azr"
msg="$(cd "$dir" && run_limited 120 "$AZURE2_BIN" calculate --no-gui t.azr 2>&1)"
if grep -q "cannot be combined with a target" <<< "$msg"; then echo "  ok    (f) the fitted width beside a target is refused"
else echo "  FAIL  (f) the fitted width beside a target: no refusal"; FAILS=$((FAILS + 1)); fi
# ---- (g) the triangular resolution's base width, fitted
python3 - "$WORK/bare/$OUT" "$WORK/tri-data" "$HERE/models/data" <<'PY'
import os, sys
rows = [list(map(float, l.split())) for l in open(sys.argv[1]) if len(l.split()) >= 4]
e, s = [r[0] for r in rows], [r[3] for r in rows]
f = e[0] / 0.79  # c.m./lab factor: the bare grid starts at 0.79 MeV lab
def at(x):
    for i in range(len(e) - 1):
        if e[i] <= x <= e[i + 1]:
            return s[i] + (x - e[i]) / (e[i + 1] - e[i]) * (s[i + 1] - s[i])
    raise ValueError(x)
os.makedirs(sys.argv[2], exist_ok=True)
for name in ("d1", "d2"):
    with open(os.path.join(sys.argv[2], name + ".dat"), "w") as out:
        for line in open(os.path.join(sys.argv[3], name + ".dat")):
            lab = float(line.split()[0])
            c, h = f * lab, 0.5 * f * 0.020
            xs = [c - h] + [x for x in e if c - h < x < c + h] + [c + h]
            if c not in xs:
                xs = sorted(xs + [c])
            t = lambda x: 1 - abs(x - c) / h
            num = sum(0.5 * (t(xs[i]) * at(xs[i]) + t(xs[i + 1]) * at(xs[i + 1])) * (xs[i + 1] - xs[i]) for i in range(len(xs) - 1))
            y = num / h  # the integral of t is h
            out.write(f"{lab:.6f} 0 {y:.6e} {0.01 * y:.6e}\n")
PY
fit triangular "$WORK/tri-data" "s/rectangular width 14 keV fit$/triangular width 14 keV fit/"
python3 - "$WORK/triangular/output/param.sav" <<'PY' || FAILS=$((FAILS + 1))
import sys
p = {l.split()[0]: (float(l.split()[1]), float(l.split()[2])) for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")}
w, dw = p.get("det_fit.energy_resolution.width", (float("nan"), 0.0))
pull = (w - 0.020) / dw if dw > 0 else float("inf")
ok = abs(pull) <= 0.05
print(("  ok    " if ok else "  FAIL  ") + f"(g) triangular, noise-free data at 20 keV: w = {w * 1e3:.6f} +- {dw * 1e3:.4f} keV ({pull:+.3f} sigma)")
sys.exit(0 if ok else 1)
PY
echo "rectangular_resolution: $([ $FAILS = 0 ] && echo PASS || echo "FAIL ($FAILS)")"
[ $FAILS = 0 ]
