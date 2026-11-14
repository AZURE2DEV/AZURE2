#!/usr/bin/env bash
#
# AZURE_PROFILE=1 reports the set-up and the evaluations apart (include/Profile.h, docs reference/command_line.rst).
# The set-up computes every point's energy-dependent values (Coulomb functions, penetrabilities, shifts) once; they
# used to be counted in the first evaluation's report of a fit, and a calculation (menu 1) labelled its report as an
# evaluation.  On tests/upgrade/T1:
#   (a) menu 1: a "set-up" report holding the energy-dependent values, then a "calculation" report holding the
#       evaluation and its stages and no energy-dependent values;
#   (b) menu 2 (two calls): the "set-up" report first, and no evaluation report holds the energy-dependent values;
#   (c) without AZURE_PROFILE nothing is reported.
#
#   ./tests/profile_report/check.sh path/to/AZURE2
set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/profile_report.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
. "$HERE/../lib/guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-300}")"
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME INPUT [ENV...]: T1 in $WORK/NAME; the profile (standard error) in err.log
run() {
  local dir="$WORK/$1" input="$2"; shift 2
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$HERE/../upgrade/T1/data" "$HERE/../upgrade/T1/T1.azr" "$dir/"
  (cd "$dir" && for kv in "$@"; do export "$kv"; done &&
     printf "$input" | $RUN "$AZURE2_BIN" --no-gui --no-readline --max-calls=2 T1.azr > run.log 2> err.log)
}
# reports FILE: one line per report, "<label>: <stages, comma-separated>"
reports() {
  awk '/^\[profile\] [^ ]/ { if (l != "") print l ": " s; l = $0; sub(/^\[profile\] /, "", l); sub(/ \(ms.*$/, "", l); s = ""; next }
       /^\[profile\]   / { x = $0; sub(/^\[profile\]   /, "", x); sub(/  +[0-9.]+( +[0-9]+)?$/, "", x); s = s (s == "" ? "" : ",") x }
       END { if (l != "") print l ": " s }' "$1"
}

run calc '1\n\n\n7\n' AZURE_PROFILE=1
r="$(reports "$WORK/calc/err.log")"
first="$(echo "$r" | sed -n 1p)"; second="$(echo "$r" | sed -n 2p)"
case "$first" in "set-up: "*"E-dependent values"*) ok "(a) menu 1: the set-up's report holds the energy-dependent values" ;;
  *) bad "(a) menu 1: first report '$first'" ;; esac
case "$second" in "calculation: evaluation,points,"*) grep -q "E-dependent" <<< "$second" \
    && bad "(a) menu 1: the calculation's report holds energy-dependent values" \
    || ok "(a) menu 1: a 'calculation' report with the evaluation's stages only" ;;
  *) bad "(a) menu 1: second report '$second'" ;; esac

run fit '2\n\n\n7\n' AZURE_PROFILE=1
r="$(reports "$WORK/fit/err.log")"
n="$(echo "$r" | grep -c '^evaluation: ')"
case "$(echo "$r" | sed -n 1p)" in "set-up: "*"E-dependent values"*)
  [ "$n" -ge 1 ] && ! grep -q "E-dependent" <<< "$(grep '^evaluation: ' <<< "$r")" \
    && ok "(b) menu 2: the set-up first, $n evaluation report(s) without the set-up's stages" \
    || bad "(b) menu 2: $n evaluation report(s); $(echo "$r" | grep '^evaluation: ' | sed -n 1p)" ;;
  *) bad "(b) menu 2: first report '$(echo "$r" | sed -n 1p)'" ;; esac

run plain '1\n\n\n7\n'
grep -q '^\[profile\]' "$WORK/plain/err.log" && bad "(c) a report without AZURE_PROFILE" || ok "(c) no report without AZURE_PROFILE"

echo
if [ "$failures" -gt 0 ]; then echo "FAIL: $failures check(s)"; exit 1; fi
echo "all checks passed"
