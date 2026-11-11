#!/usr/bin/env bash
#
# A parameter file given as a run's start is applied by name; its names are checked against the run's parameters
# (experiment/ParameterFileNames.h, AZUREParams::ApplyUserParameterFile).  A name the project has no parameter for was
# skipped and a parameter the file did not name kept the project's value, both without a word: the natural_neon
# example's .sav, still in the old names (segment_<k>_norm) after the project was converted to format 3, moved its
# chi2 from 25647 to 36346 silently.
#
# Model: tests/upgrade/T1 (format 3: 7 levels, 10 data sets, 5 extrapolations).
#   (a) its own param.par as the start of a data run: no diagnostic, the same chi2;
#   (b) energy_1 renamed bogus_1: W_PARAMETER_FILE_NAMES names bogus_1 (skipped) and energy_1 (not named); the run
#       goes on (exit 0);
#   (c) every name changed (old_<name>): E_PARAMETER_FILE_NO_MATCH, the run refused (exit 1, no chiSquared.out);
#   (d) the data run's param.par (data-set parameters, no extrapolation shifts) as the start of an extrapolation run:
#       no diagnostic (the data sets' names belong to the project's other run; the extrapolations' shifts are fixed).
#
#   ./tests/param_file_names/check.sh path/to/AZURE2
set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/param_file_names.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
. "$HERE/../lib/guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-300}")"
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME MENU-INPUT: T1 in $WORK/NAME with the given answers; sets $rc
run() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$HERE/../upgrade/T1/data" "$HERE/../upgrade/T1/T1.azr" "$dir/"
  [ -f "$WORK/start.par" ] && cp "$WORK"/*.par "$dir/"
  (cd "$dir" && printf "$2" | $RUN "$AZURE2_BIN" --no-gui --no-readline T1.azr > run.log 2>&1)
  rc=$?
}
total() { awk '/^Total-Chi-Squared:/ { print $2 }' "$WORK/$1/output/chiSquared.out" 2>/dev/null; }
said() { grep -q -- "$2" "$WORK/$1/run.log"; }

run plain '1\n\n\n7\n'
[ -s "$WORK/plain/output/param.par" ] || { bad "T1 wrote no param.par"; tail -5 "$WORK/plain/run.log"; exit 1; }
cp "$WORK/plain/output/param.par" "$WORK/start.par"
awk '$1 == "energy_1" { $1 = "bogus_1" } { print }' "$WORK/start.par" > "$WORK/partial.par"
awk '$1 != "#parametrization" { $1 = "old_" $1 } { print }' "$WORK/start.par" > "$WORK/none.par"

# (a)
run same '1\nstart.par\n\n7\n'
! said same "PARAMETER_FILE" && [ "$rc" = 0 ] && [ -n "$(total same)" ] && [ "$(total same)" = "$(total plain)" ] \
  && ok "(a) its own param.par: no diagnostic, chi2 $(total same) as without it" \
  || bad "(a) its own param.par: exit $rc, chi2 '$(total same)' vs '$(total plain)', $(grep -c PARAMETER_FILE "$WORK/same/run.log") diagnostic(s)"

# (b)
run partial '1\npartial.par\n\n7\n'
W="WARNING: partial.par: warning W_PARAMETER_FILE_NAMES: 47 of its 48 names applied; names this project has no parameter for (skipped): bogus_1; parameters it does not name (the project's values kept): energy_1"
said partial "$W" && [ "$rc" = 0 ] && [ -n "$(total partial)" ] \
  && ok "(b) one name changed: the warning names both sides, the run goes on" \
  || bad "(b) one name changed: exit $rc; $(grep PARAMETER_FILE "$WORK/partial/run.log" || echo 'no diagnostic')"

# (c)
run none '1\nnone.par\n\n7\n'
said none "ERROR: none.par: error E_PARAMETER_FILE_NO_MATCH: no name of the parameter file is a parameter of this project (old_energy_1, old_width_1_1," \
  && [ "$rc" = 1 ] && [ ! -s "$WORK/none/output/chiSquared.out" ] \
  && ok "(c) no name matches: refused, exit 1" \
  || bad "(c) no name matches: exit $rc; $(grep PARAMETER_FILE "$WORK/none/run.log" || echo 'no diagnostic')"

# (d)
run extrap '3\nn\nstart.par\n\n7\n'
said extrap "Reading User Parameter File" && ! said extrap "PARAMETER_FILE" && [ "$rc" = 0 ] \
  && ls "$WORK/extrap/output/"*.extrap > /dev/null 2>&1 \
  && ok "(d) a data run's param.par in an extrapolation run: no diagnostic" \
  || bad "(d) extrapolation run: exit $rc; $(grep PARAMETER_FILE "$WORK/extrap/run.log" || echo 'file not read or no output')"

echo
if [ "$failures" -gt 0 ]; then echo "FAIL: $failures check(s)"; exit 1; fi
echo "all checks passed"
