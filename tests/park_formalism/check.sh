#!/usr/bin/env bash
#
# --use-park: Park's parametrization must reproduce the default (Brune) one.
#
# Park's level matrix (Phys. Rev. C 104, 064612, Eqs. 9, 21, 22, 26, 28) is
# Brune's alternative level matrix (Phys. Rev. C 66, 044611, Eq. 33) with each
# basis state rescaled: gamma_Park = gamma_Brune * sqrt(J), J = 1 - sum_c
# gamma_Park^2 dS_c/dE.  The collision matrix is invariant under that
# rescaling, so a project read from the same .azr (observed energies and
# widths) has to give the same cross sections in both modes -- while the two
# modes build their level matrices from different equations and different
# reduced width amplitudes.  Checked here:
#
#   1. mode-1 calculate: same chi2 from the same .azr, different param.par;
#   2. mode-2 fit (p+p): same minimum and same parameters.out, from
#      independent MIGRAD runs in each mode's own amplitudes;
#   3. parameter files: a Brune-mode param.sav read by a Park-mode run (and
#      the other way round) is converted on read and reproduces the chi2;
#   4. a width above the J > 0 bound is reported, and the penalty makes the
#      Park objective exceed the data chi2.
#
# Models: tests/13N (12C+p elastic scattering and capture, external capture
# included) and tests/identical_pp_res (a two-channel 3P2-3F2 group).
# Runs 13N three times, about three minutes.  The analytic derivatives are
# checked separately by tests/pyazr/park_gradient_test.py.
#
#   ./tests/park_formalism/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
# A hang fails the check instead of stalling the suite (tests/lib/guard.sh).
. "$HERE/../lib/guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-900}")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TOL="${TOL:-1e-5}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/park.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
close() { awk -v a="$1" -v c="$2" -v t="$3" 'BEGIN { d = a - c; if (d < 0) d = -d; s = a; if (s < 0) s = -s; exit !(d <= t * s) }'; }

# stage PROJECT NAME: a private copy of the project under $WORK/NAME.
stage() {
  local dir="$WORK/$2"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$HERE/../$1/data" "$dir/data"
  cp "$HERE/../$1/$1.azr" "$dir/"
}
# run NAME PROJECT MODE PARAMFILE [FLAG]: mode 1 (blank PARAMFILE = from the
# .azr) or mode 2 ("fit"); prints the total chi2.
run() {
  local dir="$WORK/$1" proj="$2" mode="$3" pfile="$4" flag="${5:-}"
  if [ "$mode" = fit ]; then
    (cd "$dir" && printf '2\nn\n%s\n\n7\n' "$pfile" | $RUN "$AZURE2_BIN" --no-gui --no-readline $flag "$proj.azr" > run.log 2>&1)
  else
    (cd "$dir" && printf '1\n%s\n\n7\n' "$pfile" | $RUN "$AZURE2_BIN" --no-gui --no-readline $flag "$proj.azr" > run.log 2>&1)
  fi
  awk '/^Total-Chi-Squared:/ { print $2 }' "$dir/output/chiSquared.out" 2>/dev/null
}

# ---- 1. same calculation, different amplitudes -----------------------------
for project in 13N identical_pp_res; do
  stage "$project" "$project.brune"; stage "$project" "$project.park"
  b="$(run "$project.brune" "$project" calc "")"
  p="$(run "$project.park" "$project" calc "" --use-park)"
  if [ -z "$b" ] || [ -z "$p" ]; then bad "$project: no chi-squared (Brune '$b', Park '$p')"
  elif close "$b" "$p" "$TOL"; then ok "$project: Brune $b = Park $p"
  else bad "$project: Brune $b, Park $p"; fi
  if cmp -s "$WORK/$project.brune/output/param.par" "$WORK/$project.park/output/param.par"; then
    bad "$project: param.par identical in both modes (--use-park had no effect)"
  else ok "$project: reduced width amplitudes differ between the modes"; fi
  grep -q "^ *#parametrization *1\." "$WORK/$project.brune/output/param.par" && \
  grep -q "^ *#parametrization *2\." "$WORK/$project.park/output/param.par" && \
    ok "$project: param.par tagged with its parametrization" || bad "$project: param.par not tagged"
done

# ---- 3. parameter files convert on read -------------------------------------
# A Brune-mode param.par carries Brune amplitudes; read into a Park run it must
# be rescaled by sqrt(J) and give the same chi2 (and the other way round).
for project in identical_pp_res 13N; do
  stage "$project" "$project.x1"; stage "$project" "$project.x2"
  cp "$WORK/$project.brune/output/param.par" "$WORK/$project.x1/brune.sav"
  cp "$WORK/$project.park/output/param.par" "$WORK/$project.x2/park.sav"
  ref="$(run "$project.brune" "$project" calc "" )"   # Brune, from the .azr
  c1="$(run "$project.x1" "$project" calc brune.sav --use-park)"
  c2="$(run "$project.x2" "$project" calc park.sav)"
  grep -q "Converted the parameter file from Brune to Park" "$WORK/$project.x1/run.log" && \
    ok "$project: Brune param file converted on read in Park mode" || bad "$project: no conversion message (Brune -> Park)"
  grep -q "Converted the parameter file from Park (observed) to Brune" "$WORK/$project.x2/run.log" && \
    ok "$project: Park param file converted on read in Brune mode" || bad "$project: no conversion message (Park -> Brune)"
  if close "$ref" "$c1" "$TOL" && close "$ref" "$c2" "$TOL"; then ok "$project: converted files reproduce chi2 $ref ($c1, $c2)"
  else bad "$project: converted files give $c1 (Brune file in Park) and $c2 (Park file in Brune), expected $ref"; fi
done

# ---- 2. independent fits land on the same physical parameters --------------
stage identical_pp_res fit.brune; stage identical_pp_res fit.park
fb="$(run fit.brune identical_pp_res fit "")"
fp="$(run fit.park identical_pp_res fit "" --use-park)"
if [ -z "$fb" ] || [ -z "$fp" ]; then bad "p+p fit: no chi-squared (Brune '$fb', Park '$fp')"
elif close "$fb" "$fp" 1e-4; then ok "p+p fit: Brune minimum $fb = Park minimum $fp"
else bad "p+p fit: Brune minimum $fb, Park minimum $fp"; fi
# parameters.out: every width G (keV) within 1e-3 relative
width_column() { awk '{ for (i = 1; i < NF - 1; i++) if ($i == "G" && $(i+1) == "=") print $(i+2) }' "$1"; }
paste <(width_column "$WORK/fit.brune/output/parameters.out") \
      <(width_column "$WORK/fit.park/output/parameters.out") | \
  awk 'BEGIN { bad = 0 } { d = $1 - $2; if (d < 0) d = -d; if (d > 1e-3 * ($1 < 0 ? -$1 : $1)) { bad++; print "        width " $1 " vs " $2 } } END { exit bad }'
[ $? -eq 0 ] && ok "p+p fit: physical widths agree to 1e-3 in both modes" || bad "p+p fit: physical widths differ between the modes"
grep -q "^ *#parametrization *2\." "$WORK/fit.park/output/param.sav" && ok "p+p fit: param.sav tagged Park" || bad "p+p fit: param.sav not tagged"

# ---- 4. the J > 0 wall ------------------------------------------------------
# Multiply the p+p 0+ level's width by 170 (300 keV -> 51 MeV, past the bound
# 2P/(dS/dE) at this radius).  Both modes warn; Park's chi2 then carries the
# penalty on top of the data chi2.
stage identical_pp_res wall
sed -i '0,/ 300000\.0 / s/ 300000\.0 / 51000000 /' "$WORK/wall/identical_pp_res.azr"
w="$(run wall identical_pp_res calc "" --use-park)"
grep -q "Park norm J = .* is not positive" "$WORK/wall/run.log" && ok "wall: J <= 0 reported (--use-park)" || bad "wall: no J <= 0 warning"
penw="$(awk '/^Total-Chi-Squared:/ { for (i = 1; i < NF; i++) if ($i == "Total-Park-Chi-Squared:") print $(i+1) }' "$WORK/wall/output/chiSquared.out")"
if [ -n "$penw" ] && awk -v p="$penw" 'BEGIN { exit !(p > 1.0) }'; then ok "wall: chiSquared.out reports the penalty, Total-Park-Chi-Squared = $penw"
else bad "wall: no positive Total-Park-Chi-Squared in chiSquared.out ('$penw')"; fi

[ "$failures" -eq 0 ] && echo "  PASS" || echo "  FAIL ($failures)"
[ "$failures" -eq 0 ]
