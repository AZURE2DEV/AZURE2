#!/usr/bin/env bash
#
# Shift groups ("shiftgroup <label>" at the end of a segment line): the segments
# carrying one label share ONE constant energy shift, the parameter of the
# group's first segment.
#
#   1. a follower takes its leader's shift, whatever its own value says:
#      chi2 identical to the same segment given that shift directly;
#   2. a fit recovers a common shift from two segments mislabelled by the same
#      amount, and writes the follower's parameter as the shared value
#      (param.sav, shifts.out) -- with Minuit's numerical gradient and with
#      --use-gradient;
#   3. two segments mislabelled by DIFFERENT amounts get one compromise shift,
#      strictly between the two (independent shifts would reach both).
#
#   ./tests/shift_groups/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/shift_groups.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
seg() { awk -F, -v k="$2" '$1 == k { print $2 }' "$1/output/chiSquared.out" 2>/dev/null; }
total() { awk '/^Total-Chi-Squared:/ { print $2 }' "$1/output/chiSquared.out" 2>/dev/null; }
par() { awk -v n="$2" '$1 == n { print $2 }' "$1/output/param.sav" 2>/dev/null; }
# relclose A B TOL -> 0 when |A-B| <= TOL*max(|B|,1e-300)
relclose() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; if (m == 0) m = 1e-300; exit !(d <= t * m) }'; }
absclose() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(d < t) }'; }

# mk NAME LEVELS_FIXED SEGMENT_LINE... -> $WORK/NAME/t.azr with those <segmentsData> lines
mk() {
  local name="$1" fix="$2"; shift 2
  local dir="$WORK/$name"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$HERE/data" "$dir/data"
  printf '%s\n' "$@" > "$dir/segs.txt"
  awk -v fix="$fix" -v segs="$dir/segs.txt" '
    /<levels>/ { inlev = 1; print; next }
    /<\/levels>/ { inlev = 0 }
    inlev && fix && NF > 12 { $4 = 1; $11 = 1; print; next }   # fix every energy and width
    /<segmentsData>/ { print; while ((getline l < segs) > 0) print l; skip = 1; next }
    /<\/segmentsData>/ { skip = 0 }
    skip { next }
    { print }' "$HERE/shift_groups.azr" > "$dir/t.azr"
}
run() {  # NAME MODE EXTRA_FLAGS
  local input
  if [ "$2" = 1 ]; then input='1\n\n\n7\n'; else input='2\n\n\n\n7\n'; fi
  (cd "$WORK/$1" && printf "$input" | "$AZURE2_BIN" $3 --no-gui --no-readline t.azr > run.log 2>&1)
  if [ ! -s "$WORK/$1/output/chiSquared.out" ]; then
    bad "$1: no chiSquared.out"; tail -15 "$WORK/$1/run.log"; return 1
  fi
}

# ---- 1. a follower takes the leader's shift --------------------------------
# Segment 2's own shift (0) is ignored: it uses segment 1's 0.004.
mk group 0 \
  "1 1 1 0 100 40 40 4 1 0 0 0.004 0 0 data/xs_40.dat 0 0 shiftgroup G" \
  "1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90.dat 0 0 shiftgroup G"
mk direct 0 \
  "1 1 1 0 100 40 40 4 1 0 0 0.004 0 0 data/xs_40.dat 0 0" \
  "1 1 1 0 100 90 90 4 1 0 0 0.004 0 0 data/xs_90.dat 0 0"
mk plain 0 \
  "1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40.dat 0 0" \
  "1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90.dat 0 0"
if run group 1 "" && run direct 1 "" && run plain 1 ""; then
  g1="$(seg "$WORK/group" 1)"; g2="$(seg "$WORK/group" 2)"; d1="$(seg "$WORK/direct" 1)"; d2="$(seg "$WORK/direct" 2)"; p2="$(seg "$WORK/plain" 2)"
  if relclose "$g1" "$d1" 1e-10; then ok "leader keeps its own shift (chi2 $g1 = $d1)"; else bad "leader: grouped $g1 != direct $d1"; fi
  if relclose "$g2" "$d2" 1e-10; then ok "follower uses the leader's shift (chi2 $g2 = $d2)"; else bad "follower: grouped $g2 != direct $d2"; fi
  if relclose "$g2" "$p2" 1e-3; then bad "follower unshifted (chi2 $g2 = unshifted $p2)"; else ok "follower really shifted (chi2 $g2 vs unshifted $p2)"; fi
  if grep -q "shift group 'G': segments #1 + 2 share one constant energy shift" "$WORK/group/run.log"; then ok "group noted in the log"; else bad "no shift-group note in the log"; fi
fi

# ---- 2. fit one common shift -----------------------------------------------
# Both files are labelled 6 keV high (L = E + 0.006) with the model's exact
# values as yields, so a = -0.006 is the chi2 = 0 minimum.  The follower's own
# line says "free from 0.01": ignored.
fit_common() {  # NAME EXTRA_FLAGS
  mk "$1" 1 \
    "1 1 1 0 100 40 40 4 1 0 0 0 1 1 data/xs_40_off6.dat 0 0 shiftgroup F" \
    "1 1 1 0 100 90 90 4 1 0 0 0.01 1 1 data/xs_90_off6.dat 0 0 shiftgroup F"
  run "$1" 2 "$2" || return
  local a1 a2 t dir="$WORK/$1"
  a1="$(par "$dir" segment_1_energy_shift)"; a2="$(par "$dir" segment_2_energy_shift)"; t="$(total "$dir")"
  if absclose "$a1" -0.006 1e-5; then ok "$1: shared shift a = $a1 (truth -0.006)"; else bad "$1: shared shift a = $a1, expected -0.006"; fi
  if relclose "$a2" "$a1" 1e-12; then ok "$1: follower written with the shared value ($a2)"; else bad "$1: follower written as $a2, leader $a1"; fi
  if awk -v t="$t" 'BEGIN { exit !(t + 0 < 1e-3) }'; then ok "$1: chi2 = 0 at the truth (chi2 $t)"; else bad "$1: chi2 $t, expected ~0"; fi
  if [ "$(awk -F, 'NR > 1 { print $1 }' "$dir/output/shifts.out" | tr '\n' ' ')" = "1 2 " ]; then ok "$1: shifts.out lists both members"; else bad "$1: shifts.out rows: $(awk -F, 'NR > 1 { print $1 }' "$dir/output/shifts.out" | tr '\n' ' ')"; fi
}
fit_common fit_numerical ""
fit_common fit_analytic "--use-gradient"

# ---- 3. different mislabels -> one compromise --------------------------------
mk compromise 1 \
  "1 1 1 0 100 40 40 4 1 0 0 0 1 1 data/xs_40_off6.dat 0 0 shiftgroup F" \
  "1 1 1 0 100 90 90 4 1 0 0 0 1 1 data/xs_90_off2.dat 0 0 shiftgroup F"
if run compromise 2 "--use-gradient"; then
  a1="$(par "$WORK/compromise" segment_1_energy_shift)"; a2="$(par "$WORK/compromise" segment_2_energy_shift)"
  if awk -v a="$a1" 'BEGIN { exit !(a > -0.0059 && a < -0.0021) }' && relclose "$a2" "$a1" 1e-12; then
    ok "one shift for both, between the two truths (a = $a1; truths -0.006 / -0.002)"
  else
    bad "expected one shift strictly between -0.006 and -0.002, got $a1 / $a2"
  fi
fi

if [ "$failures" -eq 0 ]; then echo "PASS: shift groups"; else echo "FAIL: $failures check(s)"; exit 1; fi
