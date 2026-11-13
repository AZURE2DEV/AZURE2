#!/usr/bin/env bash
#
# Sub-points move rigidly with their data point under an energy-dependent shift (audit item 11, 2026-10-26).
#
# A sqrt(E) shift or a shift function used to move each sub-point of a target-integration grid by the shift at the
# sub-point's OWN energy (ESegment::UpdatePointEnergyWithShift), which stretched the integration window by 1 + ds/dE
# (a yield factor, ~1.9 % for b = 0.02 at 0.28 MeV) and left the resonance-anchored fine grid -- placed with the data
# point's shift (EData::ReplaceSubPointGrids) -- off a narrow resonance deep in the target by s' x depth.  A shift is a
# calibration of the measured beam energy: the depth grid moves with the data point and keeps its width.
# Models from make_models.py (deleted; README.md: how the models were made) (tests/thick_target_resonance, 14N(p,g) 259.7 keV, Gamma 0.97 keV, four thick-target points):
#
#   (a) sqrt (b = 0.02 on one four-point segment) == const (four one-point segments, each with the constant shift
#       b sqrt(E_k) of its point): the same calibration point by point, identical yields;
#   (b) shiftfn "a0*sqrt(x)" == sqrt (the other code path);
#   (c) sqrt == sqrt with AZURE_REFINE_QUADRATURE=8 (8 times the smooth intervals and the points per width) to 9.3e-5,
#       three times the 3.1e-5 measured (2026-11-12): the folded yields under the shift are converged on the mesh;
#   (d) the shift matters: sqrt differs from the unshifted model;
#   (e) (a) again with the resonance 100x narrower (10 eV), where the old grid missed the shifted resonance;
#   (f) narrow resonances are covered by the mesh: narrow_sqrt (a 9.7 eV level in a 3.2e19 at/cm^2 layer) == itself
#       refined eightfold to 5.4e-4, three times the 1.78e-4 measured (2026-11-12; x8 against x64: 3.2e-6).  The uniform
#       sub-point grid that useAdaptiveGrid=0 selected was a factor 100 low here and not converged at 64 times the
#       points; with the level anchors taken out of the mesh plan (MeshAnchorsOf empty) this check fails;
#   (g) the analytic plateau: narrow_sqrt (its shift and its formula stopping) minus narrow_sqrt_nores (the level's
#       proton width 0: the non-resonant part) = (lambda_r^2/2) omega Gamma_gamma (m+M)/(M eps(E_r) N) times the arctan
#       edge factor, at each point's true energy (plateau.py) to 5e-4, 2.4 times the 2.06e-4 measured (8.3e-5 refined eightfold): the
#       mesh's quadrature (f) and the R-matrix line shape against the Lorentzian.  With the stopping taken at the
#       nominal energies (the engine before 2026-11-13) it fails by 0.9 %.
#
#   ./tests/shift_sqrt_subpoints/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/shift_sqrt_subpoints.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
le() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r + 0 <= t + 0) }'; }

run() {  # MODEL [REFINE]: the run lands in $WORK/MODEL, or $WORK/MODEL.xREFINE under AZURE_REFINE_QUADRATURE=REFINE
  local model="$1" refine="${2:-}" dir="$WORK/$1"
  [ -n "$refine" ] && dir="$WORK/$1.x$refine"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$model.azr" "$dir/t.azr" && cp -r "$HERE/models/data" "$dir/"
  if [ -n "$refine" ]; then
    (cd "$dir" && printf '1\n\n\n7\n' | AZURE_REFINE_QUADRATURE="$refine" run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  else
    (cd "$dir" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  fi
  ls "$dir"/output/AZUREOut_aa=1_* > /dev/null 2>&1 || { bad "$1: no output"; tail -5 "$dir/run.log"; return 1; }
}
# yields MODEL -> "E_lab yield" per data point, every segment of the output file in order
yields() { cat "$WORK/$1"/output/AZUREOut_aa=1_*.out | awk 'NF >= 5 { printf "%.10e %.12e\n", $1, $4 }'; }
maxrel() {  # FILE1 FILE2
  paste "$1" "$2" | awk '{ d = $2 - $4; if (d < 0) d = -d; m = $4 < 0 ? -$4 : $4; if (m < 1e-300) m = 1e-300; r = d / m; if (r > w) w = r; n++ }
                        END { if (n == 0) print "nan"; else printf "%.3e\n", w }'
}

for m in sqrt const shiftfn noshift narrow_sqrt narrow_const narrow_sqrt_nores; do run "$m" || continue; done
for m in sqrt narrow_sqrt; do run "$m" 8 || continue; done
for m in sqrt const shiftfn noshift narrow_sqrt narrow_const sqrt.x8 narrow_sqrt.x8; do yields "$m" > "$WORK/$m.y" 2>/dev/null; done
n="$(grep -c . "$WORK/sqrt.y")"
[ "$n" -eq 4 ] || bad "expected 4 yield points, got $n"
[ "$(grep -c . "$WORK/const.y")" -eq 4 ] || bad "const: expected 4 points over 4 segments, got $(grep -c . "$WORK/const.y")"

r="$(maxrel "$WORK/sqrt.y" "$WORK/const.y")"
if le "$r" 1e-9; then ok "(a) sqrt shift = the same calibration as constant shifts point by point (max rel $r)"; else bad "(a) sqrt shift vs point-by-point constant shifts (max rel $r)"; fi
r="$(maxrel "$WORK/shiftfn.y" "$WORK/sqrt.y")"
if le "$r" 1e-12; then ok "(b) shiftfn a0*sqrt(x) = sqrtshift (max rel $r)"; else bad "(b) shiftfn vs sqrtshift (max rel $r)"; fi
r="$(maxrel "$WORK/sqrt.y" "$WORK/sqrt.x8.y")"
if le "$r" 9.3e-5; then ok "(c) sqrt = sqrt refined eightfold under the shift (max rel $r, tol 9.3e-5)"; else bad "(c) sqrt vs sqrt refined eightfold under the shift (max rel $r, tol 9.3e-5)"; fi
r="$(maxrel "$WORK/narrow_sqrt.y" "$WORK/narrow_const.y")"
if le "$r" 1e-9; then ok "(e) 10 eV resonance: sqrt shift = point-by-point constant shifts (max rel $r; the fine grid sits on the shifted resonance)"; else bad "(e) 10 eV resonance: sqrt shift vs point-by-point constant shifts (max rel $r)"; fi
r="$(maxrel "$WORK/sqrt.y" "$WORK/noshift.y")"
if le 5e-2 "$r"; then ok "(d) the shift moves the yields (max rel $r against b = 0)"; else bad "(d) the test is vacuous: the shift does nothing (max rel $r)"; fi
n="$(grep -c . "$WORK/narrow_sqrt.x8.y" 2>/dev/null)"
r="$(maxrel "$WORK/narrow_sqrt.y" "$WORK/narrow_sqrt.x8.y")"
if [ "${n:-0}" -eq 4 ] && le "$r" 5.4e-4; then ok "(f) 9.7 eV resonance: narrow_sqrt = itself refined eightfold (max rel $r, tol 5.4e-4)"
else bad "(f) 9.7 eV resonance: narrow_sqrt vs itself refined eightfold (max rel $r over ${n:-0} points, tol 5.4e-4)"; fi

if python3 -I "$HERE/plateau.py" "$WORK/narrow_sqrt/output" "$WORK/narrow_sqrt_nores/output" 5e-4 > "$WORK/plateau.txt" 2>&1; then
  ok "(g) 9.7 eV resonance: narrow_sqrt = the analytic plateau x arctan + the non-resonant part ($(tail -1 "$WORK/plateau.txt"))"
else bad "(g) the analytic plateau: $(tail -1 "$WORK/plateau.txt")"; cat "$WORK/plateau.txt"; fi

echo
if [ "$failures" -eq 0 ]; then echo "shift_sqrt_subpoints: PASS"; else echo "shift_sqrt_subpoints: FAIL ($failures)"; exit 1; fi
