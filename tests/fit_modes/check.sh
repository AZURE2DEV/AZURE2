#!/usr/bin/env bash
# Every way of fitting on every kind of free-parameter set (2026-11-11).
#
# The four minimizer routes (MIGRAD with the numerical gradient, --use-gradient, --use-lm, --use-gsl-lm), each with
# --covariance-band, on four projects that differ in what is free:
#
#   none     tests/thick_target_resonance      no free parameter at all
#   levels   tests/upgrade/T1 (T1_start.azr)   level parameters and data-set norms
#   mixed    tests/layer_yield/models/fit.azr  level parameters and a target's areal density
#   target   tests/stopping_yield/models/fit.azr   two areal densities and a stopping scale, no level parameter
#
# Why it exists.  The head of 2026-11-11 19:25 ended with a segmentation fault in "none" with --use-lm (the row pass
# sized a buffer by a column count of 0) and passed the whole verification, because no test fitted a project without
# free parameters.  The routes share little code below the objective (scalar gradient, residual rows, central
# differences for the parameters of a target), so an edge of the parameter set is an edge of each route separately.
# Written, the check found a second one at once: "none" with MIGRAD and --covariance-band ended in a segmentation
# fault inside Minuit2 (MnHesse on a state without a variable parameter), in the code as it was at the branch point
# with dev; src/AZUREMain.cpp now skips HESSE there.
#
# Checked: every run exits 0 and prints a total chi-squared; in "none" it is the calculation's, digit for digit; in
# the others the four routes end at the same minimum (1e-5 relative; they agree to 1e-7 today); the band is written
# where a level parameter is free and skipped with a message that says why where none is.  (The --use-lm and
# --use-gsl-lm runs of "target" take half a minute each: only --use-gradient and --use-lm are run there.)
#
# Fails on the build of 19:25 (b8c970b0) in all four routes of "none", and on the build after the first fix
# (6cc1e76e) in "none" with MIGRAD (plain and --use-gradient).  About 95 s at 2 threads.
#
# Usage:
#   ./tests/fit_modes/check.sh path/to/AZURE2
set -u
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/fit_modes.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME SRCDIR PROJECT MENU [FLAGS...]: PROJECT of SRCDIR under $WORK/NAME with the menu answers MENU; sets $rc
# and $chi2 (the printed total chi-squared, empty if none).
run() {
  local name="$1" src="$2" proj="$3" menu="$4"; shift 4
  local dir="$WORK/$name"
  mkdir -p "$dir"
  cp -R "$src/data" "$dir/data"
  cp "$src/$proj" "$dir/"
  mkdir -p "$dir/output" "$dir/checks"
  (cd "$dir" && printf "$menu" | "$AZURE2_BIN" --no-gui --no-readline "$@" "$proj" > run.log 2>&1)
  rc=$?
  chi2="$(awk '/^Total Chi-Squared:/ { v = $3 } END { print v }' "$dir/run.log")"
}

# same A B: the two chi-squared agree to 1e-5 relative.
same() { awk -v a="$1" -v b="$2" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(a > 0 && b > 0 && d <= 1e-5 * a) }'; }

FIT='2\n\n\n\n7\n'
CALC='1\n\n\n7\n'
MODES="plain --use-gradient --use-lm --use-gsl-lm"
flag() { [ "$1" = plain ] || echo "$1"; }

echo "none: no free parameter (tests/thick_target_resonance)"
SRC="$HERE/../thick_target_resonance"
run none.calc "$SRC" thick_target_resonance.azr "$CALC"
calc="$chi2"
if [ "$rc" = 0 ] && [ -n "$calc" ]; then ok "calculation: chi-squared $calc"; else bad "calculation: exit $rc, chi-squared '$calc'"; fi
for m in $MODES; do
  run "none.$m" "$SRC" thick_target_resonance.azr "$FIT" $(flag "$m") --covariance-band
  if [ "$rc" = 0 ] && [ -n "$chi2" ] && [ "$chi2" = "$calc" ]; then ok "$m: exit 0, the calculation's chi-squared"
  else bad "$m: exit $rc, chi-squared '$chi2' (calculation $calc)"; fi
  if grep -q "Uncertainty band skipped" "$WORK/none.$m/run.log" && ! ls "$WORK/none.$m/output"/*.band > /dev/null 2>&1
  then ok "$m: the band is skipped with a message"; else bad "$m: no message that the band is skipped, or a band file"; fi
done

# fitted LABEL SRCDIR PROJECT BAND(yes|no) MODES...: the routes agree on the minimum; the band is there or refused.
fitted() {
  local label="$1" src="$2" proj="$3" band="$4" first="" m; shift 4
  for m in "$@"; do
    run "$label.$m" "$src" "$proj" "$FIT" $(flag "$m") --covariance-band
    if [ "$rc" != 0 ] || [ -z "$chi2" ]; then bad "$m: exit $rc, chi-squared '$chi2'"; continue; fi
    [ -n "$first" ] || first="$chi2"
    if same "$first" "$chi2"; then ok "$m: chi-squared $chi2"; else bad "$m: chi-squared $chi2, the first route gave $first"; fi
    local dir="$WORK/$label.$m" n
    n="$(find "$dir/output" -name '*.band' -size +0 | wc -l)"
    if [ "$band" = yes ]; then
      if [ "$n" -gt 0 ] && [ -s "$dir/output/covariance.dat" ]; then ok "$m: $n band file(s) and covariance.dat"
      else bad "$m: $n band files, covariance.dat $([ -s "$dir/output/covariance.dat" ] && echo there || echo missing)"; fi
    else
      if [ "$n" = 0 ] && grep -q "no free level parameter" "$dir/run.log"; then ok "$m: the band is skipped, the message says why"
      else bad "$m: $n band files; expected none and 'no free level parameter' in the log"; fi
    fi
  done
}

echo "levels: level parameters and norms (tests/upgrade/T1, T1_start.azr)"
fitted levels "$HERE/../upgrade/T1" T1_start.azr yes $MODES

echo "mixed: level parameters and an areal density (tests/layer_yield/models/fit.azr)"
fitted mixed "$HERE/../layer_yield/models" fit.azr yes $MODES

echo "target: two areal densities and a stopping scale (tests/stopping_yield/models/fit.azr)"
fitted target "$HERE/../stopping_yield/models" fit.azr no --use-gradient --use-lm

echo
if [ "$failures" = 0 ]; then echo "fit_modes: all checks passed"; else echo "fit_modes: $failures check(s) FAILED"; fi
[ "$failures" = 0 ]
