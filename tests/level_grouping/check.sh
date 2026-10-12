#!/usr/bin/env bash
#
# How <levels> rows are grouped into levels (CNuc::Fill / JGroup::IsLevel).
#
# Each <levels> row is one channel of a level.  Rows used to be assigned to a level of their J^pi when the
# energies agreed within 1 keV, so two separate levels of one J^pi closer than that were silently merged into
# one level (the second one's channels appended to the first): a fit that drove two such levels together could
# not be written out and read back.  Rows are now grouped on the identical energy; the 1 keV tolerance is kept
# where a state typed separately in the particle pairs is matched to a level (the external-capture final state,
# ParseExternalCapture).  On copies of the 13N project:
#
#   1. a second 1/2+ level 0.6 keV above the 2.3682 MeV one is read as its own level, with a note;
#   2. the same level with all widths zero changes nothing: chi2 identical to the unmodified project;
#   3. the 1/2- ground state moved 0.4 keV away from the capture final state (E = 0) is still the level the
#      external capture attaches to: no extra 1/2- level, chi2 close to the base; moved 2 keV away it is not, and
#      the capture final state gets a level of its own.
#
#   ./tests/level_grouping/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
# macOS has no timeout(1); run without a limit there
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

SRC="$HERE/../13N"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/level_grouping.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# The 13N project without its target-integrated segments (1, 2, 8, 9) and their <targetInt> entries: the gas-
# target integration takes minutes and has nothing to do with level grouping.
awk '/^<targetInt>/ { inti = 1; print; next } /^<\/targetInt>/ { inti = 0 }
     /^<segmentsData>/ { inseg = 1; k = 0; print; next } /^<\/segmentsData>/ { inseg = 0 }
     inti && $1 == "1" { sub(/^1/, "0") }
     inseg && NF > 5 { k++; if (k == 1 || k == 2 || k == 8 || k == 9) sub(/^1/, "0") }
     { print }' "$SRC/13N.azr" > "$WORK/13N_noint.azr"

# mk NAME AWK_PROGRAM -> $WORK/NAME/13N.azr, that project with its <levels> rows passed through the awk
# program (which sees inlev=1 inside <levels>)
mk() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -R "$SRC/data" "$WORK/$1/"
  awk "/^<levels>/ { inlev = 1; print; next } /^<\/levels>/ { inlev = 0 } $2" "$WORK/13N_noint.azr" > "$WORK/$1/13N.azr"
}
# run NAME -> mode 1 (calculate with data) in $WORK/NAME
run() {
  (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline 13N.azr > run.log 2>&1)
  if [ ! -s "$WORK/$1/output/chiSquared.out" ]; then
    bad "$1: no chiSquared.out"; tail -15 "$WORK/$1/run.log"; return 1
  fi
}
chi() { awk '/^Total-Chi-Squared:/ { print $2 }' "$WORK/$1/output/chiSquared.out"; }
nlev() { grep -c "^J = $2  *E_level" "$WORK/$1/output/parameters.out"; }
relerr() { awk -v a="$1" -v b="$2" 'BEGIN { d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; printf "%.3e\n", d / m }'; }
below() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r + 0 <= t + 0) }'; }

# the unmodified project
mk base '{ print }'
run base || { echo "FAIL: base project did not run"; exit 1; }

# 1. a twin 1/2+ level at 2.3688 MeV (0.6 keV above 2.3682), level id 99, proton width 100 eV
TWIN='{ print } inlev && $1 == "0.5" && $2 == "1" && $3 == "2.3682" { n = split($0, f, " "); f[3] = "2.3688"; f[9] = 99; if (f[6] == "1") f[12] = 100; out = f[1]; for (i = 2; i <= n; i++) out = out " " f[i]; print out }'
mk twin "$TWIN"
if run twin; then
  n="$(nlev twin '0.5+')"
  if [ "$n" = 3 ]; then ok "twin 1/2+ level 0.6 keV away read as its own level (3 1/2+ levels)"; else bad "expected 3 1/2+ levels with the twin, got $n (merged?)"; fi
  if grep -q "less than 1 keV apart" "$WORK/twin/run.log"; then ok "close pair noted in the run log"; else bad "no note about the close pair"; fi
  if grep -q "E_level =  *2.3688" "$WORK/twin/output/parameters.out"; then ok "level at 2.3688 MeV in parameters.out"; else bad "no level at 2.3688 MeV in parameters.out"; fi
fi

# 2. the same twin with every width zero contributes nothing
ZERO='{ print } inlev && $1 == "0.5" && $2 == "1" && $3 == "2.3682" { n = split($0, f, " "); f[3] = "2.3688"; f[9] = 99; f[12] = 0; out = f[1]; for (i = 2; i <= n; i++) out = out " " f[i]; print out }'
mk zero "$ZERO"
if run zero; then
  n="$(nlev zero '0.5+')"; r="$(relerr "$(chi zero)" "$(chi base)")"
  if [ "$n" = 3 ]; then ok "zero-width twin read as its own level"; else bad "zero-width twin: expected 3 1/2+ levels, got $n"; fi
  if below "$r" 1e-10; then ok "zero-width twin leaves chi2 unchanged ($(chi zero) vs $(chi base), rel $r)"; else bad "zero-width twin changes chi2: $(chi zero) vs $(chi base) (rel $r)"; fi
fi

# 3. the 1/2- ground state (the capture final state, E = 0) written 0.4 keV off: still matched within 1 keV
OFF='inlev && $1 == "0.5" && $2 == "-1" && $3 == "0" { n = split($0, f, " "); f[3] = "0.0004"; out = f[1]; for (i = 2; i <= n; i++) out = out " " f[i]; print out; next } { print }'
mk off "$OFF"
if run off; then
  nb="$(nlev base '0.5-')"; n="$(nlev off '0.5-')"; r="$(relerr "$(chi off)" "$(chi base)")"
  if [ "$n" = "$nb" ]; then ok "final state 0.4 keV off still matched to the ground-state level ($n 1/2- levels, as in the base)"; else bad "final state 0.4 keV off: $n 1/2- levels, base has $nb (external capture made its own level)"; fi
  # moving the bound state by 0.4 keV is itself a small physics change (1.4e-3 here); a missed match would add a
  # level and drop the ANC to a placeholder
  if below "$r" 1e-2; then ok "chi2 close to the base ($(chi off) vs $(chi base), rel $r)"; else bad "chi2 moved: $(chi off) vs $(chi base) (rel $r)"; fi
fi

# 3b. ... and 2 keV off is outside the tolerance: the capture final state gets a level of its own
FAR='inlev && $1 == "0.5" && $2 == "-1" && $3 == "0" { n = split($0, f, " "); f[3] = "0.002"; out = f[1]; for (i = 2; i <= n; i++) out = out " " f[i]; print out; next } { print }'
mk far "$FAR"
if run far; then
  nb="$(nlev base '0.5-')"; n="$(nlev far '0.5-')"
  if [ "$n" = $((nb + 1)) ]; then ok "final state 2 keV off is not matched: the external capture adds its own 1/2- level ($n)"; else bad "final state 2 keV off: expected $((nb + 1)) 1/2- levels, got $n"; fi
fi

if [ "$failures" -eq 0 ]; then echo "PASS: level grouping"; else echo "FAIL: $failures check(s)"; exit 1; fi
