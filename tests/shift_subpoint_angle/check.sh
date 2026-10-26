#!/usr/bin/env bash
#
# Sub-points of a target/resolution integration are built at their parent's c.m. angle.  The energy-shift path
# (ESegment::UpdatePointEnergyWithShift) used to re-run the parent's lab-to-c.m. angle conversion on them, taking
# that c.m. angle for a lab angle: every sub-point moved by the full lab-c.m. difference as soon as a segment had
# any active shift (nonzero, varied, shiftgroup, sqrtshift, shiftfn).  16O(p,p) with a 0.1 keV convolution went from
# chi2 13,484 to 23,097 at a 1e-7 MeV shift.  Models from make_models.py (template: the 2026-10-25 audit model):
#
#   (a) a1 (shift 1e-7 MeV) equals a0 (shift 0) -- 16O(p,p) lab-frame differential, 11 angles, convolution;
#   (b) av (shift 0, varied) equals a0 -- a zero shift, varied or not, must not enter the conversion at all;
#   (c) g1 equals g0 on the capture segments -- a lab-frame differential 16O(p,g0) segment carrying the shift (gamma branch);
#   (d) a0 itself against expected/a0_chi2 (the unshifted value, unchanged by the fix).
#
#   ./tests/shift_subpoint_angle/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/shift_subpoint_angle.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
total() { awk '/^Total-Chi-Squared:/ { print $2 }' "$WORK/$1/output/chiSquared.out" 2>/dev/null; }
# segsum MODEL FROM TO -> chi2 summed over segments FROM..TO
segsum() { awk -F, -v a="$2" -v b="$3" '$1 ~ /^[0-9]+$/ && $1 >= a && $1 <= b { s += $2 } END { printf "%.6g", s }' "$WORK/$1/output/chiSquared.out" 2>/dev/null; }
# reldev A B -> |A-B|/|B|
reldev() { awk -v a="$1" -v b="$2" 'BEGIN { d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; if (m == 0) m = 1e-300; printf "%.2e", d / m }'; }
within() { awk -v d="$1" -v t="$2" 'BEGIN { exit !(d <= t) }'; }

run() {  # MODEL: mode 1 in $WORK/MODEL
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$1.azr" "$dir/t.azr"; cp -r "$HERE/data" "$dir/"
  (cd "$dir" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  [ -s "$dir/output/chiSquared.out" ] || { bad "$1: no chiSquared.out"; tail -5 "$dir/run.log"; }
}
same() {  # MODEL REF LABEL TOL
  local a b d; a=$(total "$1"); b=$(total "$2"); [ -n "$a" ] && [ -n "$b" ] || { bad "$3: missing chi2"; return; }
  d=$(reldev "$a" "$b")
  if within "$d" "$4"; then ok "$3: chi2 $a vs $b (rel. $d)"; else bad "$3: chi2 $a vs $b (rel. $d > $4)"; fi
}

echo "shift_subpoint_angle: $AZURE2_BIN"
for m in a0 a1 av g0 g1; do run $m; done
same a1 a0 "(a) 1e-7 MeV shift on a convolved lab-frame 16O(p,p) segment leaves chi2 unchanged" 1e-5
same av a0 "(b) a varied shift of exactly 0 leaves chi2 unchanged" 1e-6
g1=$(segsum g1 12 14); g0=$(segsum g0 12 14); d=$(reldev "$g1" "$g0")
if within "$d" 1e-5; then ok "(c) the shift on convolved lab-frame 16O(p,g0) segments leaves their chi2 unchanged: $g1 vs $g0 (rel. $d)"
else bad "(c) the shift on convolved lab-frame 16O(p,g0) segments changes their chi2: $g1 vs $g0 (rel. $d > 1e-5)"; fi
a0=$(total a0); exp=$(cat "$HERE/expected/a0_chi2"); d=$(reldev "$a0" "$exp")
if within "$d" 1e-3; then ok "(d) unshifted chi2 $a0 vs expected $exp (rel. $d)"; else bad "(d) unshifted chi2 $a0 vs expected $exp (rel. $d > 1e-3)"; fi

if [ "$failures" -eq 0 ]; then echo "shift_subpoint_angle: PASS"; else echo "shift_subpoint_angle: FAIL ($failures)"; exit 1; fi
