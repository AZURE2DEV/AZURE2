#!/usr/bin/env bash
#
# "nan" or "inf" in a <parameterSettings> field does not end the run (audit 2026-10-25, item 30; 2026-10-27).
#
# ParameterLimitsManager::ReadParameterSettings parses every numeric field with its own custom_stod, which threw
# std::invalid_argument for "nan" and "inf" -- which a fit error written from a singular covariance, or a GUI session
# that saved one, can leave in the block -- and nothing in AZUREMain caught it, so the CLI terminated with no output.
# The parser accepts them now (a fit error is a report, not an input).
#
#   (a) tests/13N as it is reproduces expected/chiSquared.out;
#   (b) the same with the fit_error field of a <parameterSettings> line set to "nan": the run completes with the same
#       total; (c) with "inf": the same.
# dev 2363ada: (b) and (c) terminate before any output.
#
#   ./tests/parameter_settings_nan/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/parameter_settings_nan.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
near() { awk -v a="$1" -v b="$2" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(a != "" && d <= 1e-5) }'; }

LINE='Level 1 Channel 1 Width (eV) 1.62 0 0 0.162 0 0 level 0'
grep -qF "$LINE" "$BASE/13N.azr" || { echo "  FAIL  tests/13N/13N.azr no longer holds the settings line '$LINE'"; exit 1; }
# model NAME FIT-ERROR-TOKEN: the 13N project with that line's fit_error field (the fourth token from the end) replaced
model() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -r "$BASE/data" "$WORK/$1/"
  awk -v old="$LINE" -v tok="$2" '{ sub(/\r$/, "") } $0 == old { $(NF - 3) = tok } { print }' "$BASE/13N.azr" > "$WORK/$1/t.azr"
}
run() { (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 900 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1); }
total() { awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$WORK/$1/output/chiSquared.out" 2> /dev/null; }
want="$(awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$BASE/expected/chiSquared.out")"

model plain 0; model withnan nan; model withinf inf
run plain; run withnan; run withinf
t="$(total plain)"
if near "$t" "$want"; then ok "(a) the plain 13N reproduces its expected total ($t)"; else bad "(a) plain 13N: '$t' vs $want"; fi
t="$(total withnan)"
if near "$t" "$want"; then ok "(b) fit_error = nan in <parameterSettings>: the run completes with the same total ($t)"; else bad "(b) fit_error = nan: total '$t' (expected $want); log: $(grep -i 'invalid_argument\|terminate\|custom_stod' "$WORK/withnan/run.log" | head -1)"; fi
t="$(total withinf)"
if near "$t" "$want"; then ok "(c) fit_error = inf: the run completes with the same total ($t)"; else bad "(c) fit_error = inf: total '$t' (expected $want)"; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
