#!/usr/bin/env bash
#
# MIGRAD call limit (--max-calls=N) and the end-of-fit convergence report.
#
# Before 2026-10-31 every MIGRAD (and NLopt) fit had a fixed limit of 50,000 function calls and said nothing when it
# reached it: a fit of ~100+ free parameters (the 19F refit of 15N+a, 145 free) stopped there and wrote
# chiSquared.out and param.sav like a converged fit.  Now there is no limit unless --max-calls=N sets one, and every
# fit says how MIGRAD ended, on the screen and in output/fitStatus.out.  Checked on tests/identical_pp_res (p+p, a few
# seconds):
#
#   1. --max-calls=20: "MIGRAD did NOT converge: stopped at the call limit", the outputs named the lowest point
#      evaluated, fitStatus.out status call_limit, valid 0, and param.sav still written;
#   2. no option: "MIGRAD converged", status converged, valid 1, call_limit none, EDM below its goal;
#   3. the same with --use-gradient (the other MIGRAD call);
#   4. a malformed --max-calls is refused with an error (exit 2) and no fit is run (before November 2026: a warning,
#      and the fit ran without a limit).
#   5. --use-lm --max-calls=20: LM ignores the limit and says so (a "Note:" line); fitStatus.out is written with the
#      keys MIGRAD's has (one writer for every fitter), minimizer and fitter lm, call_limit "20 (--max-calls) ignored";
#   6. MINOS (menu 4): its progress lines say they are MINOS's ("MINOS evaluation:"), no "Iteration:" line after
#      "Performing parameter error analysis".
#
#   ./tests/fit_call_limit/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/calllimit.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# fit NAME [FLAGS...]: a mode-2 fit of a private copy of identical_pp_res under $WORK/NAME.
fit() {
  local dir="$WORK/$1"; shift
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$HERE/../identical_pp_res/data" "$dir/data"
  cp "$HERE/../identical_pp_res/identical_pp_res.azr" "$dir/"
  (cd "$dir" && printf "${MENU:-2\nn\n\n\n7\n}" |
     "$AZURE2_BIN" --no-gui --no-readline "$@" identical_pp_res.azr > run.log 2>&1)
}
field() { awk -v k="$2" '$1 == k { $1 = ""; sub(/^ /, ""); print; exit }' "$WORK/$1/output/fitStatus.out" 2>/dev/null; }

# ---- 1. a call limit that MIGRAD reaches ------------------------------------
fit limit --max-calls=20
if grep -q "WARNING: MIGRAD did NOT converge: stopped at the call limit" "$WORK/limit/run.log"; then
  ok "--max-calls=20: the run says MIGRAD stopped at the call limit"
else bad "--max-calls=20: no call-limit warning in the log"; fi
# the warning names the point written as the notice after it does: the lowest point evaluated, not the last one
if grep -q "outputs written are the lowest point it evaluated, not a minimum" "$WORK/limit/run.log" && \
   ! grep -q "outputs written are its last point" "$WORK/limit/run.log"; then
  ok "--max-calls=20: the warning says the outputs are the lowest point evaluated"
else bad "--max-calls=20: the warning does not say the outputs are the lowest point evaluated"; fi
[ "$(field limit status)" = call_limit ] && [ "$(field limit valid)" = 0 ] && \
  [ "$(field limit call_limit)" = "20 (--max-calls)" ] \
  && ok "--max-calls=20: fitStatus.out status call_limit, valid 0, limit 20" \
  || bad "--max-calls=20: fitStatus.out status '$(field limit status)' valid '$(field limit valid)' limit '$(field limit call_limit)'"
[ -s "$WORK/limit/output/param.sav" ] && ok "--max-calls=20: param.sav still written" \
  || bad "--max-calls=20: no param.sav"

# ---- 2./3. no limit: MIGRAD converges and says so ---------------------------
for mode in numeric gradient; do
  if [ $mode = gradient ]; then fit $mode --use-gradient; else fit $mode; fi
  if grep -q "^MIGRAD converged: valid minimum" "$WORK/$mode/run.log" && \
     ! grep -q "did NOT converge" "$WORK/$mode/run.log"; then
    ok "$mode: the run says MIGRAD converged"
  else bad "$mode: no convergence line in the log"; fi
  edm="$(field $mode edm)"; goal="$(field $mode edm_goal)"
  if [ "$(field $mode status)" = converged ] && [ "$(field $mode valid)" = 1 ] && \
     [ "$(field $mode call_limit)" = none ] && [ -n "$edm" ] && \
     awk -v e="$edm" -v g="$goal" 'BEGIN { exit !(e ~ /^ *[-+]?\.?[0-9]/ && e >= 0 && e < g) }'; then
    ok "$mode: fitStatus.out status converged, valid 1, no call limit, EDM $edm < $goal"
  else
    bad "$mode: fitStatus.out status '$(field $mode status)' valid '$(field $mode valid)' limit '$(field $mode call_limit)' EDM '$edm' goal '$goal'"
  fi
done

# ---- 4. a malformed value is refused: an error, exit 2, no fit (tests/cli_option_values: every case) ---------
fit malformed --max-calls=lots
rc=$?
grep -q "^ERROR: --max-calls needs a whole number" "$WORK/malformed/run.log" && [ "$rc" = 2 ] \
  && ok "--max-calls=lots: refused with an error, exit 2" || bad "--max-calls=lots: exit $rc, no error"
[ ! -e "$WORK/malformed/output/fitStatus.out" ] && ok "--max-calls=lots: no fit was run" \
  || bad "--max-calls=lots: a fit ran (status '$(field malformed status)')"

# ---- 5. LM ignores the limit, says so, and writes fitStatus.out with MIGRAD's keys -----------------------------
fit lm --use-lm --max-calls=20
keys() { awk '{ print $1 }' "$WORK/$1/output/fitStatus.out" 2>/dev/null | tr '\n' ' '; }
grep -q "^Note: lm does not honour a call limit: the limit 20 (--max-calls) was ignored" "$WORK/lm/run.log" \
  && ok "--use-lm --max-calls=20: the run says LM ignored the limit" || bad "--use-lm --max-calls=20: no notice"
if [ -n "$(keys lm)" ] && [ "$(keys lm)" = "$(keys limit)" ] && [ "$(field lm minimizer)" = lm ] &&
   [ "$(field lm fitter)" = lm ] && [ "$(field lm call_limit)" = "20 (--max-calls) ignored" ] &&
   [ "$(field lm edm)" = none ]; then
  ok "--use-lm: fitStatus.out with MIGRAD's keys ($(keys lm)), status $(field lm status), limit ignored"
else
  bad "--use-lm: fitStatus.out keys '$(keys lm)' (MIGRAD's '$(keys limit)'), minimizer '$(field lm minimizer)',\
 call_limit '$(field lm call_limit)'"
fi

# ---- 6. MINOS's progress lines are MINOS's ---------------------------------------------------------------------
MENU='4\n1\nn\n\n\n7\n' fit minos
after="$(tr '\r' '\n' < "$WORK/minos/run.log" | sed -n '/Performing parameter error analysis/,$p')"
if [ -n "$after" ] && grep -q "MINOS evaluation:" <<< "$after" && ! grep -q "Iteration:" <<< "$after"; then
  ok "MINOS: its progress lines say MINOS, no 'Iteration:' line during it"
else
  n_it=$(grep -c 'Iteration:' <<< "$after"); n_minos=$(grep -c 'MINOS evaluation:' <<< "$after")
  bad "MINOS: $n_it 'Iteration:' and $n_minos MINOS lines in it"
fi

echo
if [ "$failures" -gt 0 ]; then echo "FAIL: $failures check(s)"; exit 1; fi
echo "all checks passed"
