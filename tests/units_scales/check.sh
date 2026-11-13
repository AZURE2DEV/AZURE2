#!/usr/bin/env bash
#
# Units of format 3 and the quantities that are scales of the engine's observables (docs/upgrade/project_format3.md
# sections 3.3, 7.2, 9.1; data_kinds.md sections 4.1, 4.6).  The chi^2 equalities against independent conversions are
# in tests/pyazr/units_scales_test.py (they need more digits than chiSquared.out prints); here, through the CLI:
#   1. the worked example of section 9.1 on T2's levels (models/yield_story.azr: counts per uC, charge state 1, a fixed
#      efficiency, a thick 15N target, keV energies): --describe and --print-predictions as committed;
#   2. its calculation at 1 and 2 threads: every output identical; a short fit lowers the chi^2;
#   3. an adaptive extrapolation grid of quantity sfactor refines on S, one of quantity sigma on sigma: below the
#      narrow level of tests/extrapolation_grids (E_cm 0.05 .. 0.25 MeV, the steep Coulomb slope of 12C(p,g)) the S grid
#      has fewer than a fifth of the sigma grid's points and still interpolates S within the tolerance against a
#      dense grid (with the level inside the range both grids are governed by its peak and coincide);
#   4. the S-factor extrapolation in keV*b equals the sigma extrapolation times E exp(2 pi eta) x 1000;
#   5. a yield on a target layer is built (tests/layer_yield tests it); a yield per uC without its charge state is refused.
#
#   ./tests/units_scales/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/units_scales.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# setup NAME FILE: $WORK/NAME with the models' data and FILE as m.azr
setup() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  ln -s "$HERE/models/data" "$WORK/$1/data"
  cp "$2" "$WORK/$1/m.azr"
}
# run NAME MENU THREADS [flags]: menu 1 (calculation), 2 (fit) or 3 (extrapolation) of $WORK/NAME/m.azr
run() {
  local name="$1" menu="$2" threads="$3"; shift 3
  (cd "$WORK/$name" && rm -rf output/* checks/* && printf '%s\nn\n\n\n7\n' "$menu" |
     OMP_NUM_THREADS="$threads" run_limited 600 "$AZURE2_BIN" --no-gui --no-readline "$@" m.azr > run.log 2>&1)
}
total() { awk '/Total-Chi-Squared:/ { print $2 }' "$1" | tail -1; }

# ---- 1. the worked example: --describe and --print-predictions
setup story "$HERE/models/yield_story.azr"
(cd "$WORK/story" && "$AZURE2_BIN" --describe m.azr 2> describe.err | sed "s|$WORK/story|<project>|g" |
  sed '/^settings of this run and machine (value, source):$/,$d' > describe.txt)  # the description; the rest is this machine's
if diff -u "$HERE/expected/story_describe.txt" "$WORK/story/describe.txt" > "$WORK/describe.diff"; then
  ok "story: --describe as committed"
else
  bad "story: --describe differs"; head -30 "$WORK/describe.diff"
fi
(cd "$WORK/story" && "$AZURE2_BIN" --print-predictions m.azr > predictions.txt 2> predictions.err)
if diff -u "$HERE/expected/story_predictions.txt" "$WORK/story/predictions.txt" > "$WORK/predictions.diff"; then
  ok "story: --print-predictions as committed ($(head -1 "$WORK/story/predictions.txt"))"
else
  bad "story: --print-predictions differs"; head -20 "$WORK/predictions.diff"
fi

# ---- 2. 1 and 2 threads; a short fit
run story 1 1 && cp -r "$WORK/story/output" "$WORK/story_1"
run story 1 2
start="$(total "$WORK/story_1/chiSquared.out")"
if [ -s "$WORK/story_1/chiSquared.out" ] && diff -r "$WORK/story_1" "$WORK/story/output" > /dev/null; then
  ok "story: calculation at 1 and 2 threads, every output identical (chi^2 $start)"
else
  bad "story: outputs differ between 1 and 2 threads (or none: $(grep -m1 ERROR "$WORK/story/run.log"))"
fi
run story 2 2 --max-calls=150
fitted="$(total "$WORK/story/output/chiSquared.out")"
if [ -n "$fitted" ] && awk -v a="$fitted" -v b="$start" 'BEGIN { exit !(a < b) }'; then
  ok "story: a fit of 150 calls lowers the chi^2 ($start -> $fitted)"
else
  bad "story: the fit did not lower the chi^2 ($start -> ${fitted:-none}; $(grep -m1 ERROR "$WORK/story/run.log"))"
fi

# ---- 3, 4. adaptive grids on S and on sigma below the narrow level; S in keV*b against sigma
narrow() {  # narrow NAME QUANTITY UNIT GRID
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  sed -e "s|@GRID@|$4|" -e "s|quantity       sigma|quantity       $2|" -e "s|unit           b$|unit           $3|" \
    "$HERE/../extrapolation_grids/narrow.azr.in" > "$WORK/$1/m.azr"
}
OUT=output/AZUREOut_aa=1_R=2.extrap
narrow s_adapt sfactor 'keV*b' "adaptive 0.05 0.25 MeV cm tolerance 1e-3 max_points 4000 start_step 0.05"
narrow x_adapt sigma b "adaptive 0.05 0.25 MeV cm tolerance 1e-3 max_points 4000 start_step 0.05"
narrow s_dense sfactor 'keV*b' "range 0.05 0.25 MeV cm step 0.0001"
narrow x_dense sigma b "range 0.05 0.25 MeV cm step 0.0001"
for n in s_adapt x_adapt s_dense x_dense; do run "$n" 3 "$OMP_NUM_THREADS"; done
if [ -s "$WORK/s_adapt/$OUT" ] && [ -s "$WORK/x_adapt/$OUT" ] && [ -s "$WORK/s_dense/$OUT" ] && [ -s "$WORK/x_dense/$OUT" ]; then
  # E_cm (column 1), sigma (column 4, b) and S (column 5, MeV b) of every grid, sigma or sfactor alike; S x 1000 in keV b
  res="$(python3 - "$WORK" "$OUT" <<'PY'
import sys
import numpy as np
work, out = sys.argv[1:3]
def load(name, col, scale):
    a = np.loadtxt(f"{work}/{name}/{out}", comments="#")
    return a[:, 0], a[:, col] * scale
es, ss = load("s_adapt", 4, 1000.0)
ex, sx = load("x_adapt", 4, 1000.0)
ed, sd = load("s_dense", 4, 1000.0)
_, s_from_sigma = load("x_dense", 4, 1000.0)
_, sigma_s = load("s_dense", 3, 1.0)
_, sigma_x = load("x_dense", 3, 1.0)
# 4: the S extrapolation = the sigma extrapolation, in its S column and in its sigma column
agree = max(np.max(np.abs(sd / s_from_sigma - 1)), np.max(np.abs(sigma_s / sigma_x - 1)))
worst = lambda e, s: np.max(np.abs(np.interp(ed, e, s) / sd - 1))
print(f"{len(es)} {len(ex)} {worst(es, ss):.2e} {worst(ex, sx):.2e} {agree:.2e}")
PY
)"
  read -r ns nx errs errx agree <<< "$res"
  if awk -v s="$errs" -v ns="$ns" -v nx="$nx" 'BEGIN { exit !(s <= 1e-3 && 5 * ns < nx) }'; then
    ok "below the level: the adaptive grid of S has $ns points and interpolates S to $errs; the grid of sigma $nx points ($errx)"
  else
    bad "adaptive S against sigma below the level: S grid $ns points $errs, sigma grid $nx points $errx"
  fi
  if awk -v a="$agree" 'BEGIN { exit !(a <= 1e-9) }'; then
    ok "sfactor extrapolation (keV*b) = the sigma extrapolation, sigma and S columns, to $agree"
  else
    bad "sfactor extrapolation against the S column of the sigma extrapolation: $agree"
  fi
else
  bad "adaptive grids: an extrapolation was not written ($(grep -h -m1 ERROR "$WORK"/*_*/run.log))"
fi

# ---- 5. refusals
sed -e 's|    kind           infinite|    kind           layer\n    areal_density  1e18 atoms/cm2 fixed|' \
  "$HERE/models/yield_story.azr" > "$WORK/story/layer.azr"  # beside its data: --describe reads them from there
msg="$("$AZURE2_BIN" --describe "$WORK/story/layer.azr" 2>&1)"
if echo "$msg" | grep -q "yield per incident particle on the layer"; then
  ok "a yield on a target layer is built (the layer kernel; its limits: tests/layer_yield)"
else
  bad "a yield on a target layer: $msg"
fi
grep -v '^    charge_state' "$HERE/models/yield_story.azr" > "$WORK/noq.azr"
msg="$("$AZURE2_BIN" --validate "$WORK/noq.azr" 2>&1)"
if echo "$msg" | grep -q "E_FIELD_REQUIRED" && echo "$msg" | grep -q "charge_state"; then
  ok "a yield per uC without charge_state is refused (E_FIELD_REQUIRED)"
else
  bad "the charge-state refusal: $msg"
fi

echo
if [ "$failures" -ne 0 ]; then
  echo "units_scales: $failures failure(s)"
  exit 1
fi
echo "units_scales: all checks passed"
