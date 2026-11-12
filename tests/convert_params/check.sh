#!/usr/bin/env bash
# A parameter file of an old project, renamed for its format-3 conversion (AZURE2 --convert-params, docs/upgrade/
# project_format3.md section 8.4; convert::ConvertParameterFile).
#
# The old project is the original of tests/upgrade/T2 (format 1, kept for the converter in tests/convert/fixtures);
# the parameter file is the param.sav of its fit by the engine that read old files (tests/upgrade/golden/T2/fit: 40
# names of data sets and target effects in the old grammar, segment_<k>_norm, density_4, ...).  The format-3 project
# is the committed conversion, tests/upgrade/T2/T2.azr.
#   (a) --convert-params writes the file with every old data-set name renamed (the converter's map; the level
#       parameters energy_<k>, width_<k>_<c> keep theirs) and refuses to overwrite;
#   (b) the format-3 project calculated (menu 1) from the renamed file gives the chi-squared recorded for the old
#       fit's parameters, expected/chi2_renamed (the same exact Coulomb functions, --no-coulomb-table), within 1e-5
#       relative (the tolerance of golden/T2's calculations; the file prints 6 digits).  The engine that read old files
#       gave the fit's own value there, 294.173 (golden/T2/fit/chiSquared.out); the shared mesh integrates the folded
#       segments 5-8 on other nodes (docs/upgrade/energy_grid.md section 7.1) and gives 294.164 at the same parameters
#       (today's engine refined eightfold: 293.649), so the recorded value is the calculation's, not the fit's;
#   (c) the old file itself, read by the format-3 project, applies only the level parameters: the run warns with
#       W_PARAMETER_FILE_NAMES and its chi-squared is not the fit's (the renaming is what makes (b) hold);
#   (d) a format-3 project as the map's source is refused, with the reason, exit 1.
# Shown failing: with one name of the map replaced by a wrong one (convert/src/Convert.cpp, RenameParameterFile), (b)
# fails (the commit that added the check).
#
# Usage:
#   ./tests/convert_params/check.sh path/to/AZURE2
set -u
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/convert_params.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

T2="$HERE/../upgrade/T2"
cp -R "$T2/data" "$T2/T2.azr" "$T2/cov_diff.txt" "$T2/prior_cov.txt" "$T2/prior_cov.f3.txt" "$WORK/"
cp "$HERE/../convert/fixtures/tests/upgrade/T2/T2.azr" "$WORK/old.azr"
cp "$HERE/../upgrade/golden/T2/fit/param.sav" "$WORK/old.sav"
mkdir -p "$WORK/output" "$WORK/checks"

# (a)
(cd "$WORK" && "$AZURE2_BIN" --convert-params old.sav new.sav --project old.azr > convert.log 2>&1); rc=$?
nold="$(grep -c '^ *segment_\|^ *density_\|^ *abundance_' "$WORK/old.sav")"
nren="$(grep -c ' -> ' "$WORK/convert.log")"
left="$(awk '$1 ~ /^(segment_|density_|abundance_)/' "$WORK/new.sav" 2>/dev/null | wc -l | tr -d ' ')"
[ "$rc" = 0 ] && [ "$nren" -gt 0 ] && [ "$nren" = "$nold" ] && [ "$left" = 0 ] \
  && [ "$(awk '{ print $2, $3 }' "$WORK/old.sav")" = "$(awk '{ print $2, $3 }' "$WORK/new.sav")" ] \
  && [ "$(grep -c '^ *energy_\|^ *width_' "$WORK/old.sav")" = "$(grep -c '^ *energy_\|^ *width_' "$WORK/new.sav")" ] \
  && ok "(a) $nren old data-set names renamed, the values and the level parameters kept" \
  || bad "(a) exit $rc, $nren renamed of $nold, $left old names left: $(tail -n 1 "$WORK/convert.log")"
(cd "$WORK" && "$AZURE2_BIN" --convert-params old.sav new.sav --project old.azr > again.log 2>&1) \
  && bad "(a) an existing new.sav was overwritten" || ok "(a) an existing new.sav is not overwritten"

# (b), (c): menu 1 with a parameter file
calc() {  # calc NAME FILE: menu 1 of T2.azr from FILE, outputs in $WORK/NAME
  rm -rf "$WORK/output" "$WORK/checks"; mkdir -p "$WORK/output" "$WORK/checks"
  (cd "$WORK" && printf '1\n%s\n\n7\n' "$2" | "$AZURE2_BIN" --no-gui --no-readline --no-coulomb-table T2.azr > "$1.log" 2>&1)
  rm -rf "$WORK/$1"; cp -R "$WORK/output" "$WORK/$1"
}
total() { awk '{ for (i = 1; i < NF; i++) if ($i == "Total-Chi-Squared:") print $(i + 1) }' "$1"; }
calc renamed new.sav
calc unrenamed old.sav
fit="$(total "$HERE/../upgrade/golden/T2/fit/chiSquared.out")"
want="$(cat "$HERE/expected/chi2_renamed")"
got="$(total "$WORK/renamed/chiSquared.out")"
old="$(total "$WORK/unrenamed/chiSquared.out")"
python3 -c "import sys; a, b = float('$got'), float('$want'); sys.exit(0 if abs(a - b) <= 1e-5 * b else 1)" 2>/dev/null \
  && ! grep -q "W_PARAMETER_FILE_NAMES\|E_PARAMETER_FILE_NO_MATCH" "$WORK/renamed.log" \
  && ok "(b) the format-3 project from the renamed file: chi2 $got, recorded at the old fit's parameters $want, no warning about names" \
  || bad "(b) chi2 '$got' against the recorded '$want' at the old fit's parameters ($(grep -m1 'PARAMETER_FILE' "$WORK/renamed.log"))"
grep -q "W_PARAMETER_FILE_NAMES" "$WORK/unrenamed.log" \
  && python3 -c "import sys; a, b = float('$old'), float('$fit'); sys.exit(0 if abs(a - b) > 1e-3 * b else 1)" 2>/dev/null \
  && ok "(c) the old file unrenamed: W_PARAMETER_FILE_NAMES, chi2 $old (not the fit's)" \
  || bad "(c) the old file unrenamed: chi2 '$old', warning: $(grep -c W_PARAMETER_FILE_NAMES "$WORK/unrenamed.log")"

# (d)
(cd "$WORK" && "$AZURE2_BIN" --convert-params old.sav new3.sav --project T2.azr > f3.log 2>&1); rc=$?
[ "$rc" = 1 ] && grep -q "ERROR: --convert-params: T2.azr: the project is format 3; give the old project" "$WORK/f3.log" \
  && [ ! -e "$WORK/new3.sav" ] \
  && ok "(d) a format-3 project as the map's source: refused, exit 1, no file" \
  || bad "(d) exit $rc, $(tail -n 1 "$WORK/f3.log")"

echo
if [ "$failures" -gt 0 ]; then echo "FAIL: $failures check(s)"; exit 1; fi
echo "all checks passed"
