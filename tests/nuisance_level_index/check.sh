#!/usr/bin/env bash
#
# A nuisance prior on a level parameter penalizes THAT parameter.
#
# <parameterSettings> can mark a level energy or width as a nuisance parameter: the fit adds
# ((p - nominal)/sigma)^2 to chi^2 (AZURECalc::CalculateNuisanceChiSquared, and the matching term in
# AddNuisanceGradient).  The ParameterLimitsManager indexes the settings by the parameter's position among the
# NON-FIXED parameters, but the vector those functions receive is Minuit's FULL vector, fixed parameters included.
# Up to dev 89fa22f they read p[nonFixedIndex], so a prior on any level parameter preceded by a fixed level
# parameter was evaluated on a different parameter (in the 13N model, energy_1 is fixed: the prior on energy_3
# landed on width_1_1, and the prior on width_1_1 on energy_1), and the fit then pulled the wrong parameter
# toward the nominal value.  The least-squares drivers' penalty rows were right (audit 2026-10-25, item 14).
#
# Model: tests/13N/13N.azr (p + 12C; energy_1 fixed, width_1_1 free, energy_3 free; the file's
# <parameterSettings> block names every free level parameter with nuisance = 0).  Checks, all in mode 1
# ("calculate with data"; the run log's "Total Chi-Squared" includes the penalty, chiSquared.out does not):
#   (a) Level 3 Energy as a nuisance parameter with nominal 2.3882 = E_3 + 0.02 and sigma 0.01: the log total
#       is the 13N total (expected/chiSquared.out) + 4;
#   (b) Level 1 Channel 1 Width as a nuisance parameter at its own value (1.62 eV, sigma 0.162 eV): the log
#       total is the 13N total (penalty 0; this is the FIRST free parameter, preceded by the fixed energy_1);
#   (c) the same as (b) for Level 3 Energy at its own value with --no-brune (the nominal is converted to the
#       formal energy and compared with the formal energy parameter): the log total equals the --no-brune
#       calculation without the prior;
#   (d) a fit of model (a) with --use-gradient and AZURE_GRAD_CHECK=exit: the analytic gradient (penalty term
#       included) agrees with finite differences (1e-5 of the largest component);
#   (e) as (b) for the gamma width of level 3 (channel 6, -0.48 eV, sigma 0.048 eV): penalty 0 (the physical
#       nominal is converted to the reduced-width amplitude it corresponds to).
# dev 89fa22f: (a) +12460 (the prior evaluated on width_1_1 = 1.2719 in reduced units), (b) +56 (on
# energy_1 = 0), (c) +12018, (e) +273466, all by the same mechanism.
#
#   ./tests/nuisance_level_index/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/nuisance_level_index.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

E3_LINE='Level 3 Energy (MeV) 2.3682 0 0 0.01 0 0 level 1'
W11_LINE='Level 1 Channel 1 Width (eV) 1.62 0 0 0.162 0 0 level 0'
G36_LINE='Level 3 Channel 6 Width (eV) -0.48 0 0 -0.048 0 0 level 3'
for l in "$E3_LINE" "$W11_LINE" "$G36_LINE"; do
  grep -qF "$l" "$BASE/13N.azr" || { echo "  FAIL  tests/13N/13N.azr no longer holds the settings line '$l'"; exit 1; }
done

# make_model <name> <old settings line> <new settings line>
make_model() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -r "$BASE/data" "$WORK/$1/"
  awk -v old="$2" -v new="$3" '{ sub(/\r$/, ""); if ($0 == old) print new; else print }' "$BASE/13N.azr" > "$WORK/$1/t.azr"
}
# log_total <name> [flags...] : mode 1, prints the run log's total (data chi^2 + penalties)
log_total() {
  local name="$1"; shift
  (cd "$WORK/$name" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" "$@" --no-gui --no-readline t.azr > run.log 2>&1)
  awk '{ sub(/\r$/, "") } /^Total Chi-Squared:/ { v = $3 } END { print v }' "$WORK/$name/run.log"
}
near() { awk -v a="$1" -v b="$2" -v tol="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(d <= tol) }'; }

want="$(awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$BASE/expected/chiSquared.out")"
[ -n "$want" ] || { echo "  FAIL  no total in tests/13N/expected/chiSquared.out"; exit 1; }

# The log prints the total with six significant digits (~336,48x): compare to +-1.5.
make_model a "$E3_LINE" 'Level 3 Energy (MeV) 2.3882 0 0 0.01 0 1 level 1'
ta="$(log_total a)"
if [ -n "$ta" ] && near "$ta" "$(awk -v w="$want" 'BEGIN { print w + 4 }')" 1.5; then
  ok "(a) prior on energy_3 (nominal +0.02, sigma 0.01): total $ta = $want + 4"
else
  bad "(a) prior on energy_3: total '$ta', expected $want + 4"; tail -3 "$WORK/a/run.log"
fi

make_model b "$W11_LINE" 'Level 1 Channel 1 Width (eV) 1.62 0 0 0.162 0 1 level 0'
tb="$(log_total b)"
if [ -n "$tb" ] && near "$tb" "$want" 1.5; then
  ok "(b) prior on width_1_1 at its own value: total $tb = $want (penalty 0)"
else
  bad "(b) prior on width_1_1 at its own value: total '$tb', expected $want"; tail -3 "$WORK/b/run.log"
fi

make_model c0 "$E3_LINE" "$E3_LINE"
make_model c "$E3_LINE" 'Level 3 Energy (MeV) 2.3682 0 0 0.01 0 1 level 1'
tc0="$(log_total c0 --no-brune)"
tc="$(log_total c --no-brune)"
if [ -n "$tc0" ] && [ -n "$tc" ] && near "$tc" "$tc0" 1.5; then
  ok "(c) --no-brune, prior on energy_3 at its own value: total $tc = $tc0 without the prior"
else
  bad "(c) --no-brune, prior on energy_3 at its own value: total '$tc', without the prior '$tc0'"; tail -3 "$WORK/c/run.log"
fi

make_model e "$G36_LINE" 'Level 3 Channel 6 Width (eV) -0.48 0 0 0.048 0 1 level 3'
te="$(log_total e)"
if [ -n "$te" ] && near "$te" "$want" 1.5; then
  ok "(e) prior on the gamma width width_3_2 at its own value: total $te = $want (penalty 0)"
else
  bad "(e) prior on the gamma width width_3_2 at its own value: total '$te', expected $want"; tail -3 "$WORK/e/run.log"
fi

(cd "$WORK/a" && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
   run_limited 1200 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > fit.log 2>&1)
status=$?
# With a nuisance prior the analytic Gauss-Newton seed is off and MIGRAD takes numerical steps before its first
# gradient call, so the check runs at some later point of the fit, not at the start: only the agreement is tested.
if [ "$status" -eq 0 ] && grep -q '^\[grad-check\] chi2 ' "$WORK/a/fit.log"; then
  ok "(d) fit of (a) starts and reaches the gradient check"
else
  bad "(d) fit of (a): exit status $status, no gradient check in the log"; tail -3 "$WORK/a/fit.log"
fi
awk '{ sub(/\r$/, "") }
     /^\[grad-check\] idx / { a = f = ""
       for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1) }
       print a, f }' "$WORK/a/fit.log" > "$WORK/rows"
n="$(wc -l < "$WORK/rows" | tr -d ' ')"
if [ "$n" -lt 5 ]; then
  bad "(d) gradient check listed $n parameters"
else
  fd="$(awk '{ a = $1 < 0 ? -$1 : $1; if (a > m) m = a; d = $1 - $2; if (d < 0) d = -d; if (d > w) w = d }
             END { printf "%.2e", w / m }' "$WORK/rows")"
  awk -v w="$fd" 'BEGIN { exit !(w + 0 <= 1e-5) }' && ok "(d) gradient = finite differences for $n parameters ($fd of the largest)" \
    || bad "(d) gradient vs finite differences: $fd of the largest"
fi

if [ "$failures" -eq 0 ]; then echo "nuisance_level_index: PASS"; else echo "nuisance_level_index: $failures check(s) FAILED"; fi
[ "$failures" -eq 0 ]
