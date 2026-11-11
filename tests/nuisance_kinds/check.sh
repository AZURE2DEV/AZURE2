#!/usr/bin/env bash
#
# A <parameterSettings> nuisance row on a parameter that is neither a level width nor a segment parameter.
#
# The limits manager converts a width's settings (partial widths, eV) to reduced-width amplitudes through the
# level parameters (ParameterLimitsManager::ConvertPhysicalLimitToReduced).  It used to run that conversion for
# every nuisance row not on a segment, indexing the vector of level parameters and abundances
# (CNuc::GetTransformParams) with the index of, e.g., a screening potential or a THM coherent background: an
# out-of-bounds write (an _GLIBCXX_ASSERTIONS build aborts there).  Only a width is converted now (the parameter's
# record says what it is); any other kind keeps its settings as given.
#
#   screening   tests/screening_6li/fit.azr (screening_1 free), row "screening_1 300 0 1000 50 0 1 ..."
#               -> the run completes, chi2 equal to the run without the row (a screening potential carries its
#                  own prior; the nuisance flag adds nothing)
#   cbkg        tests/thm_coherent_background/toy with cbackground=1/2+:2=0.05,-0.08 (Re c0, Im c0 free), row
#               "cbkg_E_1/2+_2_1/2,0,1/2,0_re0 0 0 0 0.01 0 1 ..." -> chi2 of the run without the row
#               + ((0.05 - 0)/0.01)^2 = 25 (rel 1e-6): the prior is applied with its centre and width as given
# The "Total Chi-Squared:" line on stdout includes the priors.
#
#   ./tests/nuisance_kinds/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/thm_records.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-600}")"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/nuisance_kinds.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

# run NAME SOURCE_DIR AZR "EXTRA BLOCKS (appended) or empty" -> prints the total chi2; status in NAME/status
run() {
  local d="$WORK/$1"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$2/data" "$d/"
  { tr -d '\r' < "$2/$3"; [ -z "$4" ] || printf '%b' "$4"; } > "$d/run.azr"
  (cd "$d" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr > run.log 2>&1; echo $? > status)
  grep -oE 'Total Chi-Squared: [0-9.eE+-]+' "$d/run.log" | awk '{print $3}'
}

SCR="$HERE/../screening_6li"
# the toy with its data set in experiment E
mkdir -p "$WORK/toy"
cp -r "$HERE/../thm_coherent_background/toy/data" "$HERE/../thm_coherent_background/toy/toy.azr" "$WORK/toy/"
thm_join "$WORK/toy/toy.azr" E toy
TOY="$WORK/toy"
THM='<thm>\nexperiment[E] cbackground=1/2+:2=0.05,-0.08\n</thm>\n'

s0="$(run s0 "$SCR" fit.azr "")"
s1="$(run s1 "$SCR" fit.azr '<parameterSettings>\nscreening_1 300 0 1000 50 0 1 screening 0\n</parameterSettings>\n')"
t0="$(run t0 "$TOY" toy.azr "$THM")"
t1="$(run t1 "$TOY" toy.azr "$THM"'<parameterSettings>\ncbkg_E_1/2+_2_1/2,0,1/2,0_re0 0 0 0 0.01 0 1 level 0\n</parameterSettings>\n')"
echo "  chi2: screening $s0, with its nuisance row $s1;  cbkg $t0, with its nuisance row $t1"

bad=0
if [ "$(cat "$WORK/s1/status")" = 0 ] && [ -n "$s1" ] && [ "$s1" = "$s0" ]; then
  echo "  ok    screening: a nuisance row on screening_1 runs and adds nothing (its own prior)"
else
  echo "  FAIL  screening: status $(cat "$WORK/s1/status"), chi2 '$s1' (without the row '$s0')"; tail -3 "$WORK/s1/run.log"; bad=1
fi
if [ "$(cat "$WORK/t1/status")" = 0 ] && [ -n "$t0" ] && [ -n "$t1" ] &&
   awk -v a="$t1" -v z="$t0" 'BEGIN { d = (a - z) / 25 - 1; if (d < 0) d = -d; exit !(d < 1e-6) }'; then
  echo "  ok    cbkg: the nuisance prior on Re c0 adds ((0.05 - 0)/0.01)^2 = 25, centre and width as given"
else
  echo "  FAIL  cbkg: status $(cat "$WORK/t1/status"), chi2 '$t1' (without the row '$t0'; expected + 25)"; tail -3 "$WORK/t1/run.log"; bad=1
fi
exit "$bad"
