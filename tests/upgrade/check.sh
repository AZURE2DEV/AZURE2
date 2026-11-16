#!/usr/bin/env bash
#
# The upgrade test projects T1 and T3 (tests/upgrade/T1, T3) against today's engine:
#   (a) golden/golden.py: the calculate, extrapolate and rate outputs equal the golden files to 1e-6 relative;
#   (b) recover.py: a MIGRAD fit from the moved start model returns the true parameters (|pull| <= 4 on the
#       committed noisy data, <= 0.05 on noise-free data) with chi2(fit) <= chi2(truth);
#   (c) golden/T2/regen_compare.py: the T2 calculations at the start values and at the truth equal the golden
#       files (about 3 s). The T2 fit (about 10 min with exact Coulomb functions) is not run here; see
#       tests/upgrade/T2/README.md.
#   (w) a width entered as 1 keV prints in keV in parameters.out at ten level energies around T3's 17.64 MeV
#       (the unit was chosen by a comparison the last bit of the round trip decided; about 10 s).
# Needs python3 with numpy; without them it says so and passes (as the pyazr tests skip).
#
#   ./tests/upgrade/check.sh path/to/AZURE2
set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
export AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
PY="$(command -v python3 || true)"
if [ -z "$PY" ] || ! "$PY" -c 'import numpy' 2>/dev/null; then
  echo "skip: python3 with numpy not found"; exit 0
fi
failures=0
"$PY" "$HERE/golden/golden.py" || failures=$((failures + 1))
"$PY" "$HERE/recover.py" || failures=$((failures + 1))
"$PY" "$HERE/golden/T2/regen_compare.py" --binary "$AZURE2_BIN" --threads "$OMP_NUM_THREADS" \
  --only calc_start --only calc_truth || failures=$((failures + 1))

# (w) A width entered as 1 keV prints as 1.000000 keV in parameters.out at any level energy.  The round trip
# Gamma -> gamma -> Gamma lands one rounding step above or below 1e-3 MeV by the last bits of the Coulomb functions;
# the unit was chosen by the exact value, so T3's 1 keV width (17.64 MeV level, second channel) printed as
# "1000.000000 eV" on macOS (CI run of 0c6167c1) and, here with GCC, at 17.66 and 17.67 MeV.  Ten level energies
# around T3's, exact Coulomb functions as golden.py runs them.
. "$HERE/../lib/edit_guard.sh"
UNITS="$(mktemp -d "${TMPDIR:-/tmp}/upgrade_units.XXXXXX")"
trap 'rm -rf "$UNITS"' EXIT
cp -R "$HERE/T3/." "$UNITS/"
unit_bad=0
for e in 17.61 17.62 17.63 17.64 17.65 17.66 17.67 17.68 17.69 17.70; do
  if [ "$e" = 17.64 ]; then cp "$UNITS/T3.azr" "$UNITS/v.azr"
  elif ! edited "$UNITS/T3.azr" "$UNITS/v.azr" awk -v e="$e" '/^1\.0 1 17\.64 / { $3 = e } { print }'; then
    unit_bad=$((unit_bad + 1)); continue
  fi
  rm -rf "$UNITS/output" "$UNITS/checks"; mkdir -p "$UNITS/output" "$UNITS/checks"
  (cd "$UNITS" && "$AZURE2_BIN" calculate --no-coulomb-table v.azr < /dev/null > run.log 2>&1)
  line="$(awk -v e="$e" '/^J = / { lv = ($0 ~ "E_level = +" e "00 ") } lv && /R =  2  l =   1  s =  2.0/ { print; exit }' \
          "$UNITS/output/parameters.out" 2>/dev/null)"
  case "$line" in
    *"G  =     1.000000 keV"*) ;;
    *) echo "  FAIL  (w) level at $e MeV: the 1 keV width prints as '${line:-<no line>}'"; unit_bad=$((unit_bad + 1)) ;;
  esac
done
if [ "$unit_bad" -eq 0 ]; then echo "  ok    (w) a 1 keV width prints as 1.000000 keV at ten level energies"
else failures=$((failures + 1)); fi
if [ "$failures" -eq 0 ]; then echo "upgrade T1-T3: ok"; else echo "upgrade T1-T3: FAIL"; exit 1; fi
