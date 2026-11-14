#!/usr/bin/env bash
#
# The components of composite data sets (SUM, RATIO, a second term at another angle) are points of their own at the
# energies of their base points: they share the base points' level-matrix solves through the energy tasks
# (docs/upgrade/energy_grid.md section 8), and the composite is the sum of its terms.  (Until 2026-11-11 each
# component point was mapped onto its base point, which calculated it; that mapping is deleted.)
# The 13N reference model is edited at run time into composites of both exit pairs (p+12C elastic, 1; capture to
# 13N g.s. with external capture, 2):
#
#   skowronski_luna_hpge  capture: SUM + 0.5 x capture
#   meyer_84              84.3 deg elastic, fixed energy shift 2 keV: SUM + 1000 x capture at the same angle
#   meyer_114             114.5 deg elastic: RATIO elastic / capture (unit 1)
#   kettner_0             0 deg capture, Q attenuation: SUM + capture at 90 deg lab (same exit, another angle)
#   kettner_55            55 deg capture: SUM + 1e-7 x elastic at 55 deg
#   extrapolations        0.3-3 MeV, 90 deg c.m.: SUM elastic + 1000 x capture, and SUM capture + 1000 x capture
#                         (Calculate Segments Without Data)
# and into the same model without the components ("parts": each data set its base reaction, meyer_84 with its shift;
# the extrapolation as two, elastic and capture at 90 deg c.m.; in the c.m. frame a SUM is the plain sum of its terms).
# (Format 3 cannot write two composites of the old model: a sum on the thick-target yield of Seagrave, and an
# angle-integrated elastic term of a charged pair in a sum; the capture data set skowronski_luna_hpge takes the first
# one's place, and the elastic term of kettner_55 stays differential.)
#
#   (a) the extrapolated SUM capture + 1000 x capture is 1001 x the parts' capture extrapolation, point by point
#       (1e-10, the 11 digits of the file; the component is a point of its own, calculated in its own plan task);
#   (b) the components add no level-matrix solution (AZURE_PROFILE "matrix fill"): the composite model solves as many
#       as the parts model, for the data (mode 1) and the extrapolation (mode 3) -- one per distinct compound energy;
#   (c) without the ratio: a fit with the analytic gradient starts from the chi2 of the calculation, the gradient
#       (one level-matrix adjoint per compound energy for all terms of a composite point, src/GradEnergyTasks.cpp)
#       agrees with finite differences, and the Jacobian path with the adjoint (AZURE_GRAD_CHECK);
#   (d) a capture pair that is only a component gets its external capture: the SUM elastic + 1000 x capture as the
#       only extrapolation is elastic + 1000 x capture of the parts, point by point (1e-10 of the sum; the capture
#       term is 1e-5 to 2e-2 of it), and the run writes intEC.extrap.  Until 2026-11-11 the command line switched
#       external capture on from the leading terms alone (src/ProjectSetup.cpp, ProjectReactionPairs), so this run
#       had none and its capture term was 0.54 to 2.9 times the right one; with the second extrapolation of (a) in
#       the project the same SUM was right.  Found by the agent that rewrote this check, 2026-11-11.
#
#   ./tests/component_mapping/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
# macOS has no timeout(1); run without a limit there
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/component_mapping.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# model RATIO OUT: the composites above (meyer_114's ratio only when RATIO = 1); RATIO = parts: no composite, the base
# reactions only (meyer_84 keeps its shift), and the extrapolation as its two terms
model() {
  awk -v ratio="$1" '
    { sub(/\r$/, "") }
    /^dataset / { ds = $2 }
    /^[a-z_]+ [A-Za-z0-9_-]+$/ && !/^dataset / { ds = "" }
    ratio == "parts" && $1 == "reaction" && ds == "meyer_84" { print; print "    shift          offset 0.002 MeV fixed"; next }
    ratio == "parts" && /^<extrapolations>$/ { print
      for (k = 1; k <= 2; k++) {
        print "extrapolation part" k; print "    reaction       " (k == 1 ? "12C(p,p)12C" : "12C(p,g)13N")
        print "    quantity       dsigma_domega"; print "    unit           b/sr"
        print "    energies       range 0.3 3.0 MeV lab step 0.01"; print "    angles         list 90 deg cm" }
      next }
    ratio == "parts" { print; next }
    $1 == "reaction" && ds == "skowronski_luna_hpge" { $0 = "    reaction       12C(p,g)13N + 0.5*12C(p,g)13N" }
    $1 == "reaction" && ds == "meyer_84" { $0 = "    reaction       12C(p,p)12C + 1000*12C(p,g)13N"; print; print "    shift          offset 0.002 MeV fixed"; next }
    $1 == "reaction" && ds == "meyer_114" && ratio == 1 { $0 = "    reaction       12C(p,p)12C / 12C(p,g)13N" }
    $1 == "unit" && ds == "meyer_114" && ratio == 1 { $0 = "    unit           1" }
    $1 == "reaction" && ds == "kettner_0" { $0 = "    reaction       12C(p,g)13N + 12C(p,g)13N @ 90 deg lab" }
    $1 == "reaction" && ds == "kettner_55" { $0 = "    reaction       12C(p,g)13N + 1e-7*12C(p,p)12C" }
    /^<extrapolations>$/ { print; print "extrapolation sum90"; print "    reaction       12C(p,p)12C + 1000*12C(p,g)13N"
      print "    quantity       dsigma_domega"; print "    unit           b/sr"
      print "    energies       range 0.3 3.0 MeV lab step 0.01"; print "    angles         list 90 deg cm"
      print "extrapolation cap1001"; print "    reaction       12C(p,g)13N + 1000*12C(p,g)13N"
      print "    quantity       dsigma_domega"; print "    unit           b/sr"
      print "    energies       range 0.3 3.0 MeV lab step 0.01"; print "    angles         list 90 deg cm"; next }
    { print }' "$REF/13N.azr" > "$WORK/$2"
}
model 1 t.azr
[ "$(grep -c '^    reaction       .* [+/] ' "$WORK/t.azr")" = 7 ] || { bad "model not built"; exit 1; }
model 0 g.azr
model parts p.azr
[ "$(grep -c '^extrapolation part' "$WORK/p.azr")" = 2 ] && [ "$(grep -c '^    reaction       .* [+/] ' "$WORK/p.azr")" = 0 ] ||
  { bad "parts model not built"; exit 1; }

run() {
  local name="$1" model="$2" command="$3"; shift 3
  local dir="$WORK/$name"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$REF/data" "$dir/data"
  cp "$WORK/$model" "$dir/t.azr"
  (cd "$dir" && for kv in "$@" AZURE_PROFILE=1; do export "$kv"; done &&
     run_limited 600 "$AZURE2_BIN" "$command" --no-gui t.azr > run.log 2>&1)
  if [ -z "$(ls "$dir/output" | grep -i 'AZUREOut')" ]; then bad "$name: no output written"; tail -15 "$dir/run.log"; return 1; fi
}
fills() { awk '/^\[profile\]   matrix fill/ { print $5 }' "$1/run.log"; }

run comp  t.azr calculate X=0 || exit 1
run parts p.azr calculate X=0 || exit 1
fc="$(fills "$WORK/comp")"; fp="$(fills "$WORK/parts")"
if [ -n "$fc" ] && [ -n "$fp" ] && [ "$fc" -eq "$fp" ]; then ok "(b) data: level-matrix solutions $fc, as many as without the components"
else bad "(b) data: level-matrix solutions '$fc' with the components, '$fp' without (must be equal)"; fi

run xcomp  t.azr extrapolate X=0 || exit 1
run xparts p.azr extrapolate X=0 || exit 1
fc="$(fills "$WORK/xcomp")"; fp="$(fills "$WORK/xparts")"
if [ -n "$fc" ] && [ -n "$fp" ] && [ "$fc" -eq "$fp" ]; then ok "(b) extrapolation: level-matrix solutions $fc, as many as its two terms alone"
else bad "(b) extrapolation: level-matrix solutions '$fc' with the components, '$fp' without (must be equal)"; fi
# (a) the extrapolated SUM capture + 1000 x capture (pair 2's file) against 1001 x the parts' capture, column 4
sum="$WORK/xcomp/output/AZUREOut_aa=1_R=2.extrap"; cap="$WORK/xparts/output/AZUREOut_aa=1_R=2.extrap"
if [ -s "$sum" ] && [ -s "$cap" ]; then
  res="$(awk '{ sub(/\r$/, "") } NF < 4 { next }
              FILENAME == ARGV[1] { s[FNR] = $4; e[FNR] = $1; n++; next }
              { if ($1 != e[FNR]) bad++; m++; x = 1001 * $4; d = s[FNR] - x; if (d < 0) d = -d
                r = d / (x < 0 ? -x : x); if (r > w) w = r }
              END { printf "%d %d %d %.2e", n, m, bad + 0, w }' "$sum" "$cap")"
  read -r n m nbad worst <<< "$res"
  if [ "$n" -ge 100 ] && [ "$n" = "$m" ] && [ "$nbad" = 0 ] && awk -v w="$worst" 'BEGIN { exit !(w ~ /^ *[-+]?\.?[0-9]/ && w + 0 <= 1e-10) }'; then
    ok "(a) extrapolated SUM capture + 1000 x capture = 1001 x capture at $n energies (worst $worst relative)"
  else bad "(a) extrapolated SUM against its terms: $n and $m rows, $nbad energies differ, worst $worst"; fi
else
  bad "(a) extrapolation files missing"; ls "$WORK/xcomp/output" "$WORK/xparts/output"
fi

# (d) the SUM elastic + 1000 x capture alone: no other extrapolation has the capture pair as its leading term
awk '/^extrapolation cap1001$/ { skip = 1 } /^<\/extrapolations>$/ { skip = 0 } !skip' "$WORK/t.azr" > "$WORK/alone.azr"
if [ "$(grep -c '^extrapolation ' "$WORK/alone.azr")" != 1 ]; then bad "(d) model not built"
elif run xalone alone.azr extrapolate X=0; then
  sum="$WORK/xalone/output/AZUREOut_aa=1_R=1.extrap"; el="$WORK/xparts/output/AZUREOut_aa=1_R=1.extrap"
  res="$(awk '{ sub(/\r$/, "") } NF < 4 { next }
              FILENAME == ARGV[1] { s[++n] = $4; e[n] = $1; next }
              FILENAME == ARGV[2] { a[++k] = $4; if ($1 != e[k]) bad++; next }
              { m++; if ($1 != e[m]) bad++; x = a[m] + 1000 * $4; d = s[m] - x; if (d < 0) d = -d
                r = d / (x < 0 ? -x : x); if (r > w) w = r }
              END { printf "%d %d %d %.2e", n, m, bad + 0, w }' "$sum" "$el" "$cap")"
  read -r n m nbad worst <<< "$res"
  if [ "$n" -ge 100 ] && [ "$n" = "$m" ] && [ "$nbad" = 0 ] && awk -v w="$worst" 'BEGIN { exit !(w ~ /^ *[-+]?\.?[0-9]/ && w + 0 <= 1e-10) }'; then
    ok "(d) SUM elastic + 1000 x capture alone = its terms at $n energies (worst $worst of the sum)"
  else bad "(d) SUM elastic + 1000 x capture alone against its terms: $n and $m rows, $nbad energies differ, worst $worst of the sum"; fi
  if [ -s "$WORK/xalone/output/intEC.extrap" ]; then ok "(d) the run calculated external capture (intEC.extrap written)"
  else bad "(d) no intEC.extrap: the run had no external capture"; fi
fi

# (c) the gradient check of tests/13N on the composite model
mkdir -p "$WORK/fit/output" "$WORK/fit/checks"; ln -s "$REF/data" "$WORK/fit/data"; cp "$WORK/g.azr" "$WORK/fit/t.azr"
run gcalc g.azr calculate X=0 || exit 1
(cd "$WORK/fit" && AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
   run_limited 1200 "$AZURE2_BIN" fit --use-gradient --no-gui t.azr > run.log 2>&1)
chi="$(awk '{ sub(/\r$/, "") } /^\[grad-check\] chi2 / { sub(",", "", $3); print $3 }' "$WORK/fit/run.log")"
want="$(awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$WORK/gcalc/output/chiSquared.out")"
if [ -n "$chi" ] && awk -v a="$chi" -v b="$want" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(d <= 1e-5) }'; then
  ok "(c) chi2 at the start of the fit $chi = calculation $want"
else bad "(c) chi2 at the start of the fit '$chi', calculation $want"; fi
awk '{ sub(/\r$/, "") }
     /^\[grad-check\] idx / { a = f = j = ""
       for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
       print a, f, j }' "$WORK/fit/run.log" > "$WORK/fit/rows"
nrow="$(awk 'END { print NR }' "$WORK/fit/rows")"
if [ "$nrow" -lt 5 ]; then bad "(c) gradient check listed $nrow parameters"
else
  fd="$(awk '{ a = $1 < 0 ? -$1 : $1; if (a > m) m = a; d = $1 - $2; if (d < 0) d = -d; if (d > w) w = d }
             END { printf "%.2e", w / m }' "$WORK/fit/rows")"
  awk -v w="$fd" 'BEGIN { exit !(w ~ /^ *[-+]?\.?[0-9]/ && w + 0 <= 1e-5) }' && ok "(c) gradient = finite differences for $nrow parameters ($fd of the largest)" \
    || bad "(c) gradient vs finite differences: $fd of the largest"
  jac="$(awk '{ a = $1 < 0 ? -$1 : $1; if (a > m) m = a; e = $1 - $3; if (e < 0) e = -e; if (e > v) v = e }
              END { printf "%.2e", v / m }' "$WORK/fit/rows")"
  awk -v w="$jac" 'BEGIN { exit !(w ~ /^ *[-+]?\.?[0-9]/ && w + 0 <= 1e-10) }' && ok "(c) Jacobian-path gradient = adjoint gradient ($jac of the largest)" \
    || bad "(c) Jacobian-path vs adjoint gradient: $jac of the largest"
fi

if [ "$failures" -eq 0 ]; then echo "PASS: components of composite data sets"; else echo "FAIL: $failures check(s)"; exit 1; fi
