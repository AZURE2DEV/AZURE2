#!/usr/bin/env bash
#
# A yield at a lab angle across the threshold of its reaction.
#
# A thick-target yield of an endothermic reaction at a lab angle is a fold from the beam energy down through the
# reaction's threshold.  Below the threshold there is no particle; just above it the ejectile is slower in the c.m.
# than the c.m. itself, so it goes into a forward cone only and every lab angle inside the cone is reached by two
# c.m. angles.  Until 2026-11-25 a fold that read such an energy failed ("the lab angle is beyond the kinematic
# maximum at an energy its fold reads"), so no thick-target yield at an angle could be computed for a reaction with
# a threshold (the neutrons of 19F(a,n)22Na on BaF2 and CaF2).  Now the fold takes the lab cross section at each of
# its energies: zero without a particle, the sum over both c.m. angles inside the cone (experiment::LabBranchesAt).
#
# The check does not use the kinematics it tests: the yields at 16 lab angles (Gauss-Legendre in cos psi), integrated
# over the sphere, must give the angle-integrated yield of the same target, which has no angle in it.
#
# Model: tests/channel_space (12C + p; the inelastic exit 12C[0.8], threshold 0.8672 MeV, two-valued up to 0.8733 MeV;
# capture), an infinite carbon target.
#
#   ./tests/lab_angle_threshold/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../channel_space"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/lab_threshold.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

INELASTIC="12C(p,p)12C[0.8]"
CAPTURE="12C(p,g)13N"
STOPPING='formula "((a0*(x*1e3)^a1+a2*(x*1e3)^a3)*(a4/((x*1e3)^a5)*ln(a6/(x*1e3)+a7*(x*1e3))))/((a4/((x*1e3)^a5)*ln(a6/(x*1e3)+a7*(x*1e3)))+(a0*(x*1e3)^a1+a2*(x*1e3)^a3))/1e+21" a0=2.40289 a1=0.00491497 a2=2.49101 a3=0.414939 a4=1858.36 a5=1.01581 a6=2504.17 a7=0.0181984 MeV*cm2/atom'
# 16-point Gauss-Legendre on [-1, 1]: the positive nodes and their weights
NODES="0.0950125098376374 0.2816035507792589 0.4580167776572274 0.6178762444026438 0.7554044083550030 0.8656312023878318 0.9445750230732326 0.9894009349916499"
WEIGHTS="0.1894506104550686 0.1826034150449236 0.1691565193950026 0.1495959888165767 0.1246289712555341 0.0951585116824926 0.0622535239386475 0.0271524594117542"

# model NAME REACTION KIND: channel_space with one data set of REACTION; KIND "angle" (a yield per steradian at the
# lab angles of data/NAME.dat on the infinite target), "integrated" (the yield) or "unfolded" (dsigma/dOmega at a
# lab angle, no target)
model() {
  awk -v reaction="$2" -v kind="$3" -v name="$1" -v stopping="$STOPPING" '
    { sub(/\r$/, "") }
    /^<data>$/ {
      print
      print "dataset y"
      print "    reaction       " reaction
      if (kind == "unfolded") { print "    quantity       dsigma_domega"; print "    unit           b/sr" }
      else { print "    quantity       yield"; print "    unit           " (kind == "angle" ? "1/particle/sr" : "1/particle") }
      print "    file           data/" name ".dat"
      print "    columns        E_lab[MeV] " (kind == "integrated" ? "_" : "theta_lab[deg]") " value error"
      if (kind != "integrated") print "    angle_window   0 180 deg"
      print "    energy_window  0.5 7 MeV"
      print "    norm           1 fixed"
      if (kind != "unfolded") { print "    target         carbon"; print "    quadrature     points 400" }
      skip = 1; next
    }
    /^<\/data>$/ { skip = 0 }
    /^<targets>$/ {
      print
      print "target carbon"
      print "    kind           infinite"
      print "    active         12C"
      print "    stopping       " stopping
      print "    straggling     none"
      next
    }
    !skip { print }' "$SRC/channel_space.azr" > "$WORK/$1.azr"
}
# run NAME: the project of model NAME; its output rows as "theta_cm value" in NAME.values; 1 when nothing was written
run() {
  local dir="$WORK/run_$1"
  mkdir -p "$dir/output" "$dir/checks" "$dir/data"
  cp "$WORK/$1.dat" "$dir/data/$1.dat"
  cp "$WORK/$1.azr" "$dir/t.azr"
  (cd "$dir" && run_limited 600 "$AZURE2_BIN" calculate t.azr > run.log 2>&1)
  local n
  n="$(find "$dir/output" -name 'AZUREOut_*.out' | wc -l | tr -d ' ')"
  [ "$n" = 1 ] || return 1
  awk 'NF >= 4 && $1 + 0 == $1 { print $3, $4 }' "$(find "$dir/output" -name 'AZUREOut_*.out')" > "$WORK/$1.values"
}
angles() {  # E_lab: the rows of a data file at the 16 lab angles
  awk -v e="$1" -v nodes="$NODES" 'BEGIN {
    n = split(nodes, x, " ")
    for (i = n; i >= 1; i--) printf "%s %.10f 1e-9 1e-10\n", e, atan2(sqrt(1 - x[i] * x[i]), -x[i]) * 45 / atan2(1, 1)
    for (i = 1; i <= n; i++) printf "%s %.10f 1e-9 1e-10\n", e, atan2(sqrt(1 - x[i] * x[i]), x[i]) * 45 / atan2(1, 1)
  }'
}
within() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r ~ /^ *[-+]?\.?[0-9]/ && r != "nan" && r + 0 <= t) }'; }

# 1. The yields at 16 lab angles integrate to the angle-integrated yield: far above the threshold, close to it, and
#    with the beam 7 keV above the two-valued range, so that the fold crosses all of it.
for e in 2.0 0.95 0.88; do
  angles "$e" > "$WORK/angle_$e.dat"
  echo "$e 0 1e-9 1e-10" > "$WORK/integrated_$e.dat"
  model "angle_$e" "$INELASTIC" angle
  model "integrated_$e" "$INELASTIC" integrated
  if ! run "angle_$e"; then bad "E = $e MeV: the yields at the lab angles were not written"; tail -5 "$WORK/run_angle_$e/run.log"; continue; fi
  if ! run "integrated_$e"; then bad "E = $e MeV: the angle-integrated yield was not written"; continue; fi
  if awk 'index($0, "WARNING: the calculation failed") { found = 1 } END { exit !found }' "$WORK/run_angle_$e/run.log"; then
    bad "E = $e MeV: a point failed"; else ok "E = $e MeV: no point fails"; fi
  # a row is written in the c.m. frame of the point: value / J with J = dOmega_lab/dOmega_cm =
  # (sin psi / sin theta)^2 cos(theta - psi) is the yield per lab steradian
  line="$(paste "$WORK/angle_$e.dat" "$WORK/angle_$e.values" | awk -v weights="$WEIGHTS" -v ref="$(awk '{ print $2 }' "$WORK/integrated_$e.values")" '
    BEGIN { n = split(weights, w, " "); rad = atan2(1, 1) / 45 }
    { k = NR <= n ? n + 1 - NR : NR - n; psi = $2 * rad; theta = $5 * rad
      j = (sin(psi) / sin(theta)) ^ 2 * cos(theta - psi); sum += w[k] * $6 / j; if ($6 > 0) positive++ }
    END { total = 8 * atan2(1, 1) * sum; d = total / ref - 1; if (d < 0) d = -d; printf "%.3e %.6e %.6e %d %d\n", d, total, ref, positive, NR }')"
  read -r diff total ref positive rows <<< "$line"
  if [ "$rows" = 16 ] && [ "$positive" = 16 ]; then ok "E = $e MeV: 16 positive yields"; else bad "E = $e MeV: $positive positive yields of $rows rows"; fi
  if within "$diff" 1e-5; then ok "E = $e MeV: the yields at the lab angles integrate to the angle-integrated yield ($total against $ref, rel diff $diff)"
  else bad "E = $e MeV: the integral over the lab angles $total differs from the angle-integrated yield $ref by $diff (limit 1e-5)"; fi
done

# 2. A sum with a term that has no particle at the point: below the inelastic threshold the yield of
#    capture + inelastic is the capture's, and the point is not refused.
echo "0.8 60 1e-9 1e-10" > "$WORK/sum_below.dat"; cp "$WORK/sum_below.dat" "$WORK/capture_below.dat"; cp "$WORK/sum_below.dat" "$WORK/reversed_below.dat"
model sum_below "$CAPTURE + $INELASTIC" angle
model capture_below "$CAPTURE" angle
model reversed_below "$INELASTIC + $CAPTURE" angle
if run sum_below && run capture_below; then
  r="$(paste "$WORK/sum_below.values" "$WORK/capture_below.values" | awk '{ d = $2 / $4 - 1; if (d < 0) d = -d; printf "%.3e\n", d }')"
  if within "$r" 1e-9; then ok "below the inelastic threshold capture + inelastic at 60 deg is the capture's yield (rel diff $r)"
  else bad "below the inelastic threshold capture + inelastic differs from capture by $r (limit 1e-9)"; fi
else bad "capture + inelastic below the inelastic threshold was not computed"; tail -5 "$WORK/run_sum_below/run.log"; fi
if run reversed_below && awk '$2 + 0 > 0 { found = 1 } END { exit !found }' "$WORK/reversed_below.values"; then
  ok "with the closed term first the point is computed too (in the lab frame)"
else bad "inelastic + capture below the inelastic threshold was not computed"; tail -5 "$WORK/run_reversed_below/run.log"; fi

# 3. A point that is not folded has no model where no particle reaches its lab angle: refused, as before.
echo "0.8 60 1e-9 1e-10" > "$WORK/unfolded.dat"
model unfolded "$INELASTIC" unfolded
run unfolded
if awk 'index($0, "ERROR") && index($0, "below the threshold of the exit pair") { found = 1 } END { exit !found }' "$WORK/run_unfolded/run.log"; then
  ok "a cross section at a lab angle below the threshold, without a target, is refused"
else bad "an unfolded point below the threshold was not refused"; tail -5 "$WORK/run_unfolded/run.log"; fi

if [ "$failures" -eq 0 ]; then echo "PASSED"; else echo "FAILED: $failures"; fi
exit "$failures"
