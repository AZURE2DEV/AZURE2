#!/usr/bin/env bash
#
# A thick target (target integration with density 0) is an explicit mode, not a mutation of the effect
# (audit item 17, 2026-10-27).
#
# A <targetInt> line with target integration and a density of 0 means an infinitely thick target: the yield per
# incident particle, int_0^E0 sigma / epsilon dE, not divided by a density.  Up to dev 5df7024 this was done by
# overwriting the SHARED effect's density with 1e24 at the first point whose thickness came out < 1e-10, which
# (i) left that first point with a target thickness of 0, so with a beam spread or straggling on the line its
# yield was exactly 0 (the Gaussian branch integrates from E0 - thickness to E0), and (ii) also fired for a THIN
# target whose stopping cross section evaluated to <= 0 at some point, silently turning the rest of the line into a
# thick target.  Now TargetEffect::IsThickTarget is set from the density, each thick-target point integrates from
# its energy down to 1 keV, the yield is not divided, and a non-positive stopping power on a thin target is an ERROR.
#
# Models are made at run time from tests/13N/13N.azr (p + 12C; its line "1,2" is a thick target with no spread):
#   (a) the file as it is: the 13N total (expected/chiSquared.out);
#   (b) the thick line with a 2 keV beam spread: every yield of segments 1 and 2 is positive and the first point of
#       segment 1 is within 1 % of its value without the spread (the spread moves the neighbours by ~0.3 %);
#   (c) the thick line with straggling (coefficient 0.04): the same;
#   (d) a THIN target (density 1e18) with the stopping cross section a0 x + a1 = x - 0.5, negative below 0.5 MeV:
#       the run stops with an ERROR naming the segment and energy, and writes no chiSquared.out.
# dev 5df7024: (b) and (c) first point 0 (b: 0 vs 2.85e-14), (d) runs and writes negative yields (chi2 436,536).
#
#   ./tests/thick_target/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thick_target.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
near() { awk -v a="$1" -v b="$2" -v tol="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(d <= tol) }'; }

LINE_PREFIX='1  "1,2"  10  0  0  1  0  '
grep -qF "$LINE_PREFIX" "$BASE/13N.azr" || { echo "  FAIL  tests/13N/13N.azr no longer holds the thick-target line '$LINE_PREFIX...'"; exit 1; }

# make_model NAME AWK-PROGRAM : the 13N project with its thick-target line rewritten (fields: 1 active, 2 segments,
# 3 sub-points, 4 convolution, 5 sigma, 6 target integration, 7 density, 8 equation, 9 n, 10.. parameters,
# then Q coefficients, convolution coefficients, 23 straggling, 24 coefficient, 25 widths, 26 points per width)
make_model() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -r "$BASE/data" "$WORK/$1/"
  awk -v pre="$LINE_PREFIX" "{ sub(/\r\$/, \"\"); if (index(\$0, pre) == 1) { $2 } print }" "$BASE/13N.azr" > "$WORK/$1/t.azr"
}
run() { (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1); }
total() { awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$WORK/$1/output/chiSquared.out" 2> /dev/null; }
# seg1 NAME -> "E fit" rows of segment 1 (block 1 of the p -> gamma output file)
seg1() { awk 'BEGIN { RS = ""; FS = "\n" } NR == 1 { print }' "$WORK/$1/output/AZUREOut_aa=1_R=2.out" | awk 'NF >= 5 { print $1, $4 }'; }
seg2() { awk 'BEGIN { RS = ""; FS = "\n" } NR == 2 { print }' "$WORK/$1/output/AZUREOut_aa=1_R=2.out" | awk 'NF >= 5 { print $1, $4 }'; }

want="$(awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$BASE/expected/chiSquared.out")"

# ---- (a) the file as it is ----------------------------------------------------------------------------------------
make_model plain ''
run plain
ta="$(total plain)"
if [ -n "$ta" ] && near "$ta" "$want" 1.5; then ok "(a) thick target without spread: total $ta = $want"; else bad "(a) thick target without spread: total '$ta', expected $want"; tail -3 "$WORK/plain/run.log"; fi
first="$(seg1 plain | head -1 | awk '{ print $2 }')"
[ -n "$first" ] || { bad "(a) no segment-1 output"; first=0; }

# ---- (b) a beam spread on the thick line, (c) straggling ------------------------------------------------------------
check_spread() {  # TAG NAME
  local tag="$1" name="$2"
  run "$name"
  local n npos f
  n="$( { seg1 "$name"; seg2 "$name"; } | wc -l | tr -d ' ')"
  npos="$( { seg1 "$name"; seg2 "$name"; } | awk '$2 > 0' | wc -l | tr -d ' ')"
  f="$(seg1 "$name" | head -1 | awk '{ print $2 }')"
  if [ "$n" -ge 30 ] && [ "$npos" -eq "$n" ]; then ok "$tag every yield of segments 1 and 2 is positive ($n points)"; else bad "$tag $((n - npos)) of $n yields are not positive"; seg1 "$name" | head -2; fi
  if [ -n "$f" ] && awk -v a="$f" -v b="$first" 'BEGIN { exit !(b > 0 && a / b > 0.99 && a / b < 1.01) }'; then
    ok "$tag first point of segment 1: $f, within 1 % of $first without the spread"
  else
    bad "$tag first point of segment 1: '$f' vs $first without the spread"; tail -2 "$WORK/$name/run.log"
  fi
}
make_model conv '$4 = 1; $5 = 0.002;'
check_spread "(b) 2 keV beam spread:" conv
make_model strag '$23 = 1; $24 = 0.04;'
check_spread "(c) straggling:" strag

# ---- (d) a thin target with a negative stopping cross section --------------------------------------------------------
make_model negstop '$7 = 1e18; $8 = "\"a0*x+a1\""; $9 = 2; $10 = 1; $11 = -0.5; for (i = 12; i <= 17; i++) $i = ""; $0 = $0; $1 = $1;'
run negstop
if grep -q "^ERROR: target integration of segment #1 at E_lab = 0.19 MeV" "$WORK/negstop/run.log" && [ ! -s "$WORK/negstop/output/chiSquared.out" ]; then
  ok "(d) a thin target with a negative stopping cross section stops with an ERROR"
else
  bad "(d) a thin target with a negative stopping cross section: expected an ERROR and no output (total '$(total negstop)')"; grep -m2 "ERROR\|Chi" "$WORK/negstop/run.log"
fi

if [ "$failures" -eq 0 ]; then echo "thick_target: PASS"; else echo "thick_target: $failures check(s) FAILED"; fi
[ "$failures" -eq 0 ]
