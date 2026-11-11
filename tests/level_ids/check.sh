#!/usr/bin/env bash
#
# Lines of one J-group are read as one level when they carry the same level ID (column 9) within 1 keV, or the
# identical energy.  Until 2026-10-24 any lines within 1 keV were one level whatever their IDs (JGroup::IsLevel), so
# a narrow level put next to a broad one of the same J^pi was merged into it: the later lines' widths overwrote the
# earlier level's and the fit lost that level's parameters, with only a "lists channel N of its J-group twice"
# warning (17O+a/10-22-26_an3_jpiscan, 3/2- at 8.2900 MeV, ID 91, swallowed by 3/2- 8.2893 MeV, ID 22).
# Edits of the 13N reference model, which has the 3/2- level ID 5 at 3.5032 MeV (CLI defaults):
#
#   zero     a 3/2- level ID 13 at 3.5039 MeV (0.7 keV above ID 5) with zero widths, lines after all others:
#            chiSquared.out identical to the reference (before the fix its zeros replaced ID 5's widths);
#   narrow   the same level with 50 eV / -0.01 eV widths after all others, and
#   first    those lines before all others: the same chiSquared.out, different from the reference;
#            both print the "kept as two levels" warning;
#   same     ID 13 at exactly 3.5032 MeV with ID 5's widths: still one level (warning), chi2 = reference.
#
#   ./tests/level_ids/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
# macOS has no timeout(1); run without a limit there
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/level_ids.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# mk NAME WHERE ENERGY W1 W2 W3: copy of the reference with ID 5's three lines repeated as ID 13 at ENERGY with
# widths W1 W2 W3, placed at the END of <levels> or FIRST
mk() {
  awk -v e="$3" -v w1="$4" -v w2="$5" -v w3="$6" '
    function out() { for (i = 1; i <= n; i++) print held[i] }
    /<levels>/ { f = 1; print; next }
    /<\/levels>/ { f = 0; out() }
    f && $9 == 5 { n++; t = $0; split(t, a, " "); a[3] = e; a[9] = 13; a[12] = (n == 1 ? w1 : n == 2 ? w2 : w3)
                   s = a[1]; for (i = 2; i <= length(a); i++) s = s " " a[i]; held[n] = s }
    { print }' "$REF/13N.azr" > "$WORK/$1.tmp"
  if [ "$2" = first ]; then
    awk -v src="$WORK/$1.tmp" 'BEGIN { while ((getline l < src) > 0) if (split(l, a, " ") > 20 && a[9] == 13) h[++n] = l }
      { if (split($0, a, " ") > 20 && a[9] == 13) next; print; if ($0 == "<levels>") for (i = 1; i <= n; i++) print h[i] }' \
      "$WORK/$1.tmp" > "$WORK/$1.azr"
  else
    mv "$WORK/$1.tmp" "$WORK/$1.azr"
  fi
}
mk zero   end   3.5039 0 0 0
mk narrow end   3.5039 50 -0.01 0
mk first  first 3.5039 50 -0.01 0
mk same   end   3.5032 55200 -0.49 0.00072
n13() { awk 'split($0, a, " ") > 20 && a[9] == 13' "$1" | awk 'END { print NR }'; }
for m in zero narrow first same; do
  [ "$(n13 "$WORK/$m.azr")" = 3 ] || { bad "$m: model not built"; exit 1; }
done
[ "$(awk 'split($0, a, " ") > 20 { print a[9]; exit }' "$WORK/first.azr")" = 13 ] || { bad "first: ID 13 not first"; exit 1; }

# run NAME MODEL -> $WORK/NAME/output/chiSquared.out, $WORK/NAME/run.log
run() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$REF/data" "$dir/data"
  cp "$2" "$dir/t.azr"
  (cd "$dir" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  if [ ! -s "$dir/output/chiSquared.out" ]; then
    bad "$1: no chiSquared.out written"; tail -15 "$dir/run.log"; return 1
  fi
}
total() { grep -o 'Total-Chi-Squared: [^ ]*' "$1/output/chiSquared.out" | awk '{ print $2 }'; }
# awk only: the Windows runner has no diff
ndiff() { awk 'NR == FNR { a[FNR] = $0; n = FNR; next } { m = FNR; if (a[FNR] != $0) d++ } END { print d + (n > m ? n - m : m - n) }' "$1" "$2"; }
relclose() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; exit !(m > 0 && d <= t * m) }'; }

run ref "$REF/13N.azr" || exit 1
for m in zero narrow first same; do run "$m" "$WORK/$m.azr" || exit 1; done
ref="$(total "$WORK/ref")"

if [ "$(ndiff "$WORK/ref/output/chiSquared.out" "$WORK/zero/output/chiSquared.out")" = 0 ]; then
  ok "zero: a zero-width level 0.7 keV from ID 5 leaves chiSquared.out unchanged (total $ref)"
else bad "zero: total chi2 $(total "$WORK/zero"), reference $ref (ID 5's widths overwritten?)"; fi
tn="$(total "$WORK/narrow")"; tf="$(total "$WORK/first")"
if [ "$(ndiff "$WORK/narrow/output/chiSquared.out" "$WORK/first/output/chiSquared.out")" = 0 ]; then
  ok "narrow = first: line order does not matter (chiSquared.out identical, total $tn)"
else bad "narrow total $tn, first total $tf: the line order changes the model"; fi
if relclose "$tn" "$ref" 1e-6; then bad "narrow: total $tn equals the reference; the narrow level has no effect"
else ok "narrow: the narrow level changes the total ($ref -> $tn)"; fi
for m in zero narrow first; do
  if grep -q 'within 1 keV of level ID .*they are kept as two levels' "$WORK/$m/run.log"; then ok "$m: warning printed"
  else bad "$m: no 'kept as two levels' warning"; fi
  if grep -q 'of its J-group twice' "$WORK/$m/run.log"; then bad "$m: lines of two levels read as one"; fi
done
if [ "$(ndiff "$WORK/ref/output/chiSquared.out" "$WORK/same/output/chiSquared.out")" = 0 ]; then
  ok "same: ID 13 at ID 5's energy and widths is one level (chiSquared.out unchanged)"
else bad "same: total chi2 $(total "$WORK/same"), reference $ref"; fi
if grep -q 'have the same J^pi and energy; their lines are read as one level' "$WORK/same/run.log"; then
  ok "same: warning printed"; else bad "same: no same-energy warning"; fi
if grep -q 'WARNING: level' "$WORK/ref/run.log"; then bad "reference: unexpected level warning"; fi

if [ "$failures" -eq 0 ]; then echo "PASS: level IDs"; else echo "FAIL: $failures check(s)"; exit 1; fi
