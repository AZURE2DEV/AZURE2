#!/usr/bin/env bash
#
# Sub-points move rigidly with their data point under an energy-dependent shift (audit item 11, 2026-10-26).
#
# A sqrt(E) shift or a shift function used to move each sub-point of a target-integration grid by the shift at the
# sub-point's OWN energy (ESegment::UpdatePointEnergyWithShift), which stretched the integration window by 1 + ds/dE
# (a yield factor, ~1.9 % for b = 0.02 at 0.28 MeV) and left the resonance-anchored fine grid -- placed with the data
# point's shift (EData::ReplaceSubPointGrids) -- off a narrow resonance deep in the target by s' x depth.  A shift is a
# calibration of the measured beam energy: the depth grid moves with the data point and keeps its width.
# Models from make_models.py (tests/thick_target_resonance, 14N(p,g) 259.7 keV, Gamma 0.97 keV, four thick-target points):
#
#   (a) sqrt (b = 0.02 on one four-point segment) == const (four one-point segments, each with the constant shift
#       b sqrt(E_k) of its point): the same calibration point by point, identical yields;
#   (b) shiftfn "a0*sqrt(x)" == sqrt (the other code path);
#   (c) sqrt == a uniform 2000-point sub-point grid to 1 % (the fine grid sits on the shifted resonance);
#   (d) the shift matters: sqrt differs from the unshifted model;
#   (e) (a) again with the resonance 100x narrower (10 eV), where the old grid missed the shifted resonance.
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

run() {  # MODEL
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$1.azr" "$dir/t.azr" && cp -r "$HERE/models/data" "$dir/"
  (cd "$dir" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  ls "$dir"/output/AZUREOut_aa=1_* > /dev/null 2>&1 || { bad "$1: no output"; tail -5 "$dir/run.log"; return 1; }
}
# yields MODEL -> "E_lab yield" per data point, every segment of the output file in order
yields() { cat "$WORK/$1"/output/AZUREOut_aa=1_*.out | awk 'NF >= 5 { printf "%.10e %.12e\n", $1, $4 }'; }
maxrel() {  # FILE1 FILE2
  paste "$1" "$2" | awk '{ d = $2 - $4; if (d < 0) d = -d; m = $4 < 0 ? -$4 : $4; if (m < 1e-300) m = 1e-300; r = d / m; if (r > w) w = r; n++ }
                        END { if (n == 0) print "nan"; else printf "%.3e\n", w }'
}

for m in sqrt const shiftfn dense noshift narrow_sqrt narrow_const; do run "$m" || continue; done
for m in sqrt const shiftfn dense noshift narrow_sqrt narrow_const; do yields "$m" > "$WORK/$m.y" 2>/dev/null; done
n="$(grep -c . "$WORK/sqrt.y")"
[ "$n" -eq 4 ] || bad "expected 4 yield points, got $n"
[ "$(grep -c . "$WORK/const.y")" -eq 4 ] || bad "const: expected 4 points over 4 segments, got $(grep -c . "$WORK/const.y")"

r="$(maxrel "$WORK/sqrt.y" "$WORK/const.y")"
if le "$r" 1e-9; then ok "(a) sqrt shift = the same calibration as constant shifts point by point (max rel $r)"; else bad "(a) sqrt shift vs point-by-point constant shifts (max rel $r)"; fi
r="$(maxrel "$WORK/shiftfn.y" "$WORK/sqrt.y")"
if le "$r" 1e-12; then ok "(b) shiftfn a0*sqrt(x) = sqrtshift (max rel $r)"; else bad "(b) shiftfn vs sqrtshift (max rel $r)"; fi
r="$(maxrel "$WORK/sqrt.y" "$WORK/dense.y")"
if le "$r" 1e-2; then ok "(c) adaptive grid = uniform 2000-point grid under the shift (max rel $r)"; else bad "(c) adaptive vs uniform 2000-point grid under the shift (max rel $r)"; fi
r="$(maxrel "$WORK/narrow_sqrt.y" "$WORK/narrow_const.y")"
if le "$r" 1e-9; then ok "(e) 10 eV resonance: sqrt shift = point-by-point constant shifts (max rel $r; the fine grid sits on the shifted resonance)"; else bad "(e) 10 eV resonance: sqrt shift vs point-by-point constant shifts (max rel $r)"; fi
r="$(maxrel "$WORK/sqrt.y" "$WORK/noshift.y")"
if le 5e-2 "$r"; then ok "(d) the shift moves the yields (max rel $r against b = 0)"; else bad "(d) the test is vacuous: the shift does nothing (max rel $r)"; fi

echo
if [ "$failures" -eq 0 ]; then echo "shift_sqrt_subpoints: PASS"; else echo "shift_sqrt_subpoints: FAIL ($failures)"; exit 1; fi
