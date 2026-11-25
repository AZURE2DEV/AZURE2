#!/bin/sh
# AZURE2 --convert does not depend on the user's locale: the program never calls setlocale, so a decimal-comma LC_ALL
# (de_DE, it_IT, pl_PL, ...) leaves strtod and the streams in the "C" locale.  Two history fixtures with decimals are
# converted with LC_ALL=C and with a decimal-comma locale of this machine; the outputs must be equal byte for byte.
#   locale_check.sh <AZURE2> <repository root>      exit 77 (skipped) where the machine has no decimal-comma locale
AZURE2_BIN=$1
ROOT=$2
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
COMMA=""
for l in de_DE.utf8 de_DE.UTF-8 de_DE it_IT.utf8 it_IT.UTF-8 pl_PL.utf8 pl_PL.UTF-8 fr_FR.utf8 fr_FR.UTF-8; do
  if locale -a 2>/dev/null | grep -qx "$l" && [ "$(LC_ALL=$l locale decimal_point 2>/dev/null)" = "," ]; then COMMA=$l; break; fi
done
if [ -z "$COMMA" ]; then echo "skip: no decimal-comma locale on this machine"; exit 77; fi
fail=0
for f in target_2025_01_coefficients energy_shift_2025_08 components_scaling_2026_05; do
  src="$ROOT/tests/convert/fixtures/history/$f/$f.azr"
  LC_ALL=C "$AZURE2_BIN" --convert "$src" "$WORK/$f.c.azr" > "$WORK/$f.c.log" 2>&1
  LC_ALL=$COMMA "$AZURE2_BIN" --convert "$src" "$WORK/$f.comma.azr" > "$WORK/$f.comma.log" 2>&1
  if [ -s "$WORK/$f.c.azr" ] && cmp -s "$WORK/$f.c.azr" "$WORK/$f.comma.azr"; then
    echo "ok    $f: the same bytes with LC_ALL=C and LC_ALL=$COMMA"
  else
    echo "FAIL  $f: the outputs differ (or none was written)"; fail=1
  fi
done
exit $fail
