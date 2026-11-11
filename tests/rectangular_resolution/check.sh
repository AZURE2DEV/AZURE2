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
  (cd "$dir" && printf '3\nn\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  [ -s "$dir/$OUT" ] || { echo "  FAIL  $1: no output"; tail -5 "$dir/run.log"; echo "rectangular_resolution: FAIL"; exit 1; }
}

echo "rectangular_resolution: $AZURE2_BIN"
run bare; run rect
python3 - "$WORK/bare/$OUT" "$WORK/rect/$OUT" <<'PY'
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
print("rectangular_resolution: " + ("PASS" if failures == 0 else f"FAIL ({failures})"))
sys.exit(1 if failures else 0)
PY
