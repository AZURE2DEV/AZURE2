#!/usr/bin/env bash
#
# The yield of a sum of reactions on a thick target (quantity yield with "r1 + r2").
#
# A thick-target yield counted over several exit groups (neutrons to every state of the residual nucleus,
# 19F(a,n)22Na on BaF2 and CaF2) is one data set whose reaction is a sum.  Format 3 refused it ("takes one
# reaction, not a sum") until 2026-11-25, so such a project of the old format did not convert.  The fold over the
# target is linear, so the yield of the sum must be the sum of the terms' yields on the same target: checked
# angle-integrated on an infinite target, and at a lab angle, where each term has its own c.m. angle and solid-angle
# factor, on a layer (1e19 atoms/cm2, beam energies from 2 MeV).  The layer because a thick-target yield at a lab
# angle is not computed for a reaction with a threshold, sum or not: its fold reaches energies where the lab angle
# is beyond the term's kinematic maximum, and the point fails with a warning (PointFold.cpp; the same on the old
# engine, "the lab angle cannot be converted").
#
# Model: tests/channel_space (12C + p with the inelastic exit 12C[0.8] and capture), a carbon target.
#
#   ./tests/yield_sum/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../channel_space"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/yield_sum.XXXXXX")"
trap "true" EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

INELASTIC="12C(p,p)12C[0.8]"
CAPTURE="12C(p,g)13N"
STOPPING='formula "((a0*(x*1e3)^a1+a2*(x*1e3)^a3)*(a4/((x*1e3)^a5)*ln(a6/(x*1e3)+a7*(x*1e3))))/((a4/((x*1e3)^a5)*ln(a6/(x*1e3)+a7*(x*1e3)))+(a0*(x*1e3)^a1+a2*(x*1e3)^a3))/1e+21" a0=2.40289 a1=0.00491497 a2=2.49101 a3=0.414939 a4=1858.36 a5=1.01581 a6=2504.17 a7=0.0181984 MeV*cm2/atom'

mkdir -p "$WORK/data"
# ten beam energies above the inelastic threshold; the values are not compared
awk 'BEGIN { for (i = 0; i < 10; i++) printf "%.3f 0 1e-9 1e-10\n", 1.2 + 0.4 * i }' > "$WORK/data/integrated.dat"
awk 'BEGIN { for (i = 0; i < 10; i++) printf "%.3f 60 1e-9 1e-10\n", 2.0 + 0.3 * i }' > "$WORK/data/angle.dat"

# model NAME REACTION KIND: channel_space with one data set, the yield of REACTION on a carbon target; KIND
# integrated (an infinite target) or angle (60 deg lab, a layer)
model() {
  awk -v reaction="$2" -v kind="$3" -v stopping="$STOPPING" '
    { sub(/\r$/, "") }
    /^<data>$/ {
      print
      print "dataset y"
      print "    reaction       " reaction
      print "    quantity       yield"
      if (kind == "angle") {
        print "    unit           1/particle/sr"
        print "    file           data/angle.dat"
        print "    columns        E_lab[MeV] theta_lab[deg] value error"
        print "    angle_window   0 180 deg"
      } else {
        print "    unit           1/particle"
        print "    file           data/integrated.dat"
        print "    columns        E_lab[MeV] _ value error"
      }
      print "    energy_window  0.5 7 MeV"
      print "    norm           1 fixed"
      print "    target         carbon"
      print "    quadrature     points 200"
      skip = 1; next
    }
    /^<\/data>$/ { skip = 0 }
    /^<targets>$/ {
      print
      print "target carbon"
      if (kind == "angle") {
        print "    kind           layer"
        print "    active         12C"
        print "    areal_density  1e+19 atoms/cm2 fixed"
      } else {
        print "    kind           infinite"
        print "    active         12C"
      }
      print "    stopping       " stopping
      print "    straggling     none"
      next
    }
    !skip { print }' "$SRC/channel_space.azr" > "$WORK/$1.azr"
}
run() {
  local dir="$WORK/run_$1"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$WORK/data" "$dir/data"
  cp "$WORK/$1.azr" "$dir/t.azr"
  (cd "$dir" && run_limited 600 "$AZURE2_BIN" calculate t.azr > run.log 2>&1)
  local n
  n="$(find "$dir/output" -name 'AZUREOut_*.out' | wc -l | tr -d ' ')"
  if [ "$n" != 1 ]; then bad "$1: $n output files, 1 expected"; tail -12 "$dir/run.log"; return 1; fi
  awk 'NF >= 4 && $1 + 0 == $1 { print $1, $4 }' "$(find "$dir/output" -name 'AZUREOut_*.out')" > "$WORK/$1.values"
  if [ "$(wc -l < "$WORK/$1.values" | tr -d ' ')" != 10 ]; then bad "$1: not 10 points"; return 1; fi
}
within() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r ~ /^ *[-+]?\.?[0-9]/ && r != "nan" && r + 0 <= t) }'; }

for kind in integrated angle; do
  model "sum_$kind" "$INELASTIC + $CAPTURE" "$kind"
  model "reversed_$kind" "$CAPTURE + $INELASTIC" "$kind"
  model "inelastic_$kind" "$INELASTIC" "$kind"
  model "capture_$kind" "$CAPTURE" "$kind"
  run "sum_$kind" && run "reversed_$kind" && run "inelastic_$kind" && run "capture_$kind" || continue
  # At a lab angle the numbers are written in the c.m. frame of the data set's first term: the sum is
  # a + rho b with rho = (dOmega_cm/dOmega_lab of term 2) / (that of term 1), and the reversed sum b + a / rho, which
  # gives rho without the kinematics.  Angle-integrated, rho = 1 and the reversed sum is the sum.
  line="$(paste "$WORK/sum_$kind.values" "$WORK/inelastic_$kind.values" "$WORK/capture_$kind.values" "$WORK/reversed_$kind.values" | awk -v kind="$kind" '
    { s = $2; a = $4; b = $6; v = $8
      rho = 1; if (kind == "angle") rho = a / (v - b)
      d = s - (a + rho * b); if (d < 0) d = -d; m = a + rho * b; if (m < 0) m = -m; if (m < 1e-300) m = 1e-300
      r = d / m; if (r > worst) worst = r
      d = v - (b + a / rho); if (d < 0) d = -d; m = v < 0 ? -v : v; if (m < 1e-300) m = 1e-300; if (d / m > worstReversed) worstReversed = d / m
      fa = a / (a + rho * b); fb = rho * b / (a + rho * b); if (fa > shareA) shareA = fa; if (fb > shareB) shareB = fb
      if (rho < lo || NR == 1) lo = rho; if (rho > hi || NR == 1) hi = rho; if (s > 0) positive++ }
    END { printf "%.3e %.3e %.3e %.3e %d %.4f %.4f\n", worst, worstReversed, shareA, shareB, positive, lo, hi }')"
  read -r worst worstReversed shareA shareB positive lo hi <<< "$line"
  if within "$worst" 1e-9; then ok "$kind: yield of the sum = sum of the yields (max rel diff $worst; frame factor $lo to $hi)"
  else bad "$kind: yield of the sum differs from the sum of the yields by $worst (limit 1e-9)"; fi
  if [ "$kind" = integrated ]; then
    if within "$worstReversed" 1e-9; then ok "$kind: the order of the terms does not matter (max rel diff $worstReversed)"
    else bad "$kind: the reversed sum differs by $worstReversed (limit 1e-9)"; fi
  elif awk -v lo="$lo" -v hi="$hi" 'BEGIN { exit !(lo ~ /^ *[-+]?\.?[0-9]/ && hi ~ /^ *[-+]?\.?[0-9]/ && lo + 0 > 0.5 && hi + 0 < 2.0 && (lo + 0 < 0.99 || hi + 0 > 1.01)) }'; then
    ok "$kind: the two terms have different frames (so the factor is tested)"
  else bad "$kind: frame factor $lo to $hi: not between 0.5 and 2, or 1"; fi
  if [ "$positive" = 10 ]; then ok "$kind: the yields are positive at all 10 energies"; else bad "$kind: only $positive of 10 yields positive"; fi
  # both terms matter (the smaller one is far above the limit of the comparison)
  if within "$shareA" 1e-5 || within "$shareB" 1e-5; then bad "$kind: a term never reaches 1e-5 of the sum (shares $shareA, $shareB)"
  else ok "$kind: both terms contribute (largest shares $shareA and $shareB)"; fi
done

# A ratio is still not a yield.
awk '{ sub(/\r$/, "") } $1 == "reaction" { print "    reaction       '"$INELASTIC"' / '"$CAPTURE"'"; next } { print }' "$WORK/sum_integrated.azr" > "$WORK/ratio.azr"
mkdir -p "$WORK/run_ratio/output" "$WORK/run_ratio/checks"; ln -s "$WORK/data" "$WORK/run_ratio/data"; cp "$WORK/ratio.azr" "$WORK/run_ratio/t.azr"
(cd "$WORK/run_ratio" && run_limited 120 "$AZURE2_BIN" calculate t.azr > run.log 2>&1)
if awk 'index($0, "E_REACTION_NOT_ALLOWED") && index($0, "ratio") { found = 1 } END { exit !found }' "$WORK/run_ratio/run.log"; then ok "a ratio as a yield is refused (E_REACTION_NOT_ALLOWED)"
else bad "a ratio as a yield was not refused"; tail -5 "$WORK/run_ratio/run.log"; fi

if [ "$failures" -eq 0 ]; then echo "PASSED"; else echo "FAILED: $failures"; fi
exit "$failures"
