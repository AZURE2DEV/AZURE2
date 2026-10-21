#!/usr/bin/env bash
#
# Electron screening (a "screening Ue vary sigma" block on a <targetInt> line).  Models from
# make_models.py: tests/13N with only the low-energy capture segments 6, 7 (LUNA) and 8 (FELS) active.
#
#   1. screened / bare = exp(pi eta U_e / E) point by point on segments 6 and 7, segment 8 untouched;
#   2. the factor acts inside a target effect: a Gaussian beam spread of 1e-6 MeV gives the plain result;
#   3. the analytic d chi2 / d U_e (and the Jacobian column) agree with finite differences, with a beam
#      spread of 2 keV (sub-points) and a prior;
#   4. MIGRAD and --use-lm recover U_e = 1000 eV from data made with it (levels fixed), and a prior pulls
#      the fit towards its centre, with Total-Screening-Chi-Squared in chiSquared.out;
#   5. a screening block without its three numbers is refused.
#
#   ./tests/electron_screening/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/electron_screening.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME AZR MODE FLAGS [env...] -> $WORK/NAME (mode 1 calculate, 2 fit)
run() {
  local name="$1" azr="$2" mode="$3" flags="$4"; shift 4
  local dir="$WORK/$name"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$SRC/data" "$dir/data"
  [ -d "$WORK/fitdata" ] && ln -s "$WORK/fitdata" "$dir/fitdata"
  cp "$azr" "$dir/t.azr"
  local input='1\n\n\n7\n'
  [ "$mode" = 2 ] && input='2\nn\n\n\n'
  (cd "$dir" && for kv in "$@"; do export "$kv"; done
   printf "$input" | run_limited 900 "$AZURE2_BIN" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}
out() { echo "$WORK/$1/output/AZUREOut_aa=1_R=2.out"; }
# The value of parameter NAME in a run's param.sav (name value ... per line).
par() { awk -v n="$2" '$1 == n { print $2; exit }' "$WORK/$1/output/param.sav"; }

# ---- 1. the factor, point by point -------------------------------------------------------
run bare "$HERE/bare.azr" 1 "" OMP_NUM_THREADS="$OMP_NUM_THREADS"
run scr "$HERE/scr.azr" 1 "" OMP_NUM_THREADS="$OMP_NUM_THREADS"
if [ -s "$(out bare)" ] && [ -s "$(out scr)" ]; then
  # Rows in data order: segment 6 (32 points), 7 (18), then 8.  exp(pi eta Ue/E), eta of p + 12C with
  # AZURE2's constants (Constants.h) and the pair's masses (1.00728, 12).
  n6="$(grep -c . "$SRC/data/skowronski_luna_hpge.dat")"; n7="$(grep -c . "$SRC/data/skowronski_luna_bgo.dat")"
  read -r worst unscr nrow <<< "$(paste "$(out bare)" "$(out scr)" | awk -v n=$((n6 + n7)) '
    BEGIN { pi = 3.141592650; uconv = 931.4940880; fs = 1 / 137.0359996790; mu = 1.00728 * 12 / 13.00728 }
    NF >= 8 { r++; m = NF / 2; E = $1; b = $4; s = $(m + 4)
      if (r <= n) { eta = sqrt(uconv / 2) * fs * 6 * sqrt(mu / E); f = exp(pi * eta * 1.0e-3 / E)
                    d = s / (b * f) - 1; if (d < 0) d = -d; if (d > w) w = d }
      else { d = s / b - 1; if (d < 0) d = -d; if (d > u) u = d } }
    END { printf "%.3e %.3e %d\n", w, u, r }')"
  if awk -v w="$worst" 'BEGIN { exit !(w < 1e-9) }'; then ok "segments 6, 7: screened/bare = exp(pi eta U_e/E) (worst $worst, $nrow rows)"
  else bad "segments 6, 7: screened/bare differs from exp(pi eta U_e/E) by $worst"; fi
  if awk -v u="$unscr" 'BEGIN { exit !(u == 0) }'; then ok "segment 8 (no screening) unchanged"
  else bad "segment 8 changed by $unscr"; fi
else
  bad "bare/scr: no output"; tail -5 "$WORK/scr/run.log"
fi

# ---- 2. inside a target effect -------------------------------------------------------------
run conv "$HERE/conv.azr" 1 "" OMP_NUM_THREADS="$OMP_NUM_THREADS"
if [ -s "$(out conv)" ]; then
  d="$(paste "$(out scr)" "$(out conv)" | awk 'NF >= 8 { m = NF / 2; x = $(m + 4) / $4 - 1; if (x < 0) x = -x; if (x > w) w = x } END { printf "%.3e", w }')"
  if awk -v d="$d" 'BEGIN { exit !(d < 1e-5) }'; then ok "screening inside a 1e-6 MeV beam spread = plain screening (max rel diff $d)"
  else bad "screening inside a 1e-6 MeV beam spread differs from plain screening by $d"; fi
else
  bad "conv: no output"; tail -5 "$WORK/conv/run.log"
fi

# ---- 3. gradient ------------------------------------------------------------------------
run grad "$HERE/grad.azr" 2 "--use-gradient" OMP_NUM_THREADS="$OMP_NUM_THREADS" AZURE_GRAD_CHECK=exit
line="$(grep "components with relative disagreement" "$WORK/grad/run.log" | head -1)"
if echo "$line" | grep -q "screening 0/1;.*screening 0/1"; then ok "d chi2/d U_e: analytic = finite differences, adjoint and Jacobian paths"
else bad "gradient check: ${line:-no [grad-check] line}"; tail -5 "$WORK/grad/run.log"; fi
rel="$(sed -n 's/.*screening J-1.*analytic \([^ ]*\) *fd \([^ ]*\) *rel \([^ ]*\).*/\3/p' "$WORK/grad/run.log" | head -1)"
if [ -n "$rel" ] && awk -v r="$rel" 'BEGIN { exit !(r < 1e-6) }'; then ok "U_e gradient: relative difference $rel"
else bad "U_e gradient: relative difference ${rel:-missing}"; fi

# ---- 4. fits ----------------------------------------------------------------------------
# Data = the screened model (U_e = 1000 eV) with 1 % errors, at the LUNA energies; every level fixed.
if [ -s "$(out scr)" ]; then
  mkdir -p "$WORK/fitdata"
  awk 'NF >= 8' "$(out scr)" > "$WORK/scr_rows"
  i=0
  for f in skowronski_luna_hpge skowronski_luna_bgo; do
    awk -v skip="$i" 'NR == FNR { if (FNR > skip) v[FNR - skip] = $4; next }
                      NF >= 4 { k++; printf "%s %s %.12e %.12e\n", $1, $2, v[k], 0.01 * v[k] }' \
        "$WORK/scr_rows" "$SRC/data/$f.dat" > "$WORK/fitdata/$f.dat"
    i=$((i + $(grep -c . "$SRC/data/$f.dat")))
  done
  mkfit() {  # start nominal sigma vary -> azr on stdout
    awk -v ue="$1" -v sig="$2" '
      /<levels>/ { lv = 1 } /<\/levels>/ { lv = 0 }
      lv && NF > 20 { $4 = 1; $11 = 1 }
      { sub(/data\/skowronski_luna_hpge.dat/, "fitdata/skowronski_luna_hpge.dat")
        sub(/data\/skowronski_luna_bgo.dat/, "fitdata/skowronski_luna_bgo.dat")
        sub(/screening  1000  0  0/, "screening  " ue "  1  " sig) }
      { print }' "$HERE/scr.azr"
  }
  mkfit 300 0 > "$WORK/fit.azr"
  mkfit 300 50 > "$WORK/fit_prior.azr"
  for v in "fit:--use-gradient" "fit_lm:--use-lm"; do
    name="${v%%:*}"; flags="${v#*:}"
    run "$name" "$WORK/fit.azr" 2 "$flags" OMP_NUM_THREADS="$OMP_NUM_THREADS"
    ue="$(par "$name" screening_1)"
    if [ -n "$ue" ] && awk -v u="$ue" 'BEGIN { d = u / 1000 - 1; if (d < 0) d = -d; exit !(d < 1e-3) }'; then
      ok "$name ($flags): U_e $ue eV from 300 (truth 1000)"
    else bad "$name ($flags): U_e ${ue:-missing}, expected 1000"; tail -5 "$WORK/$name/run.log"; fi
  done
  run fit_prior "$WORK/fit_prior.azr" 2 "--use-gradient" OMP_NUM_THREADS="$OMP_NUM_THREADS"
  ue="$(par fit_prior screening_1)"
  sc="$(sed -n 's/.*Total-Screening-Chi-Squared: \([^ ]*\).*/\1/p' "$WORK/fit_prior/output/chiSquared.out" | head -1)"
  if [ -n "$ue" ] && [ -n "$sc" ] && awk -v u="$ue" -v s="$sc" 'BEGIN { p = ((u - 300) / 50)^2; d = s - p; if (d < 0) d = -d
         exit !(u > 300 && u < 1000 && d <= 1e-4 * (p + 1)) }'; then
    ok "prior 300 +- 50 eV pulls U_e to $ue eV; Total-Screening-Chi-Squared $sc = ((U_e - 300)/50)^2"
  else bad "prior fit: U_e ${ue:-missing}, Total-Screening-Chi-Squared ${sc:-missing}"; fi
else
  bad "fits skipped: no screened output to make data from"
fi

# ---- 5. a malformed block ---------------------------------------------------------------
run bad "$HERE/bad.azr" 1 "" OMP_NUM_THREADS="$OMP_NUM_THREADS"
if grep -q "a screening block is 'screening Ue_eV vary sigma_eV'" "$WORK/bad/run.log" && [ ! -s "$(out bad)" ]; then
  ok "a screening block without its sigma is refused"
else bad "a malformed screening block was not refused"; tail -5 "$WORK/bad/run.log"; fi

if [ "$failures" -eq 0 ]; then echo "PASS electron_screening"; exit 0; fi
echo "FAIL electron_screening ($failures)"; exit 1
