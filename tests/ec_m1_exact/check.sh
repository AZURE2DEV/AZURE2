#!/usr/bin/env bash
#
# External capture, absolute, against an independent golden-rule calculation (truth.py -> expected/sigma*.dat),
# hard-sphere external capture only (models from make_models.py):
#
#   p + 16O -> 17F(0.495 MeV, 1/2+) + gamma
#   (a) E1 in the long-wavelength approximation;
#   (b) E1 with --no-long-wavelength: the charge term of the extended Siegert theorem (Friar and Fallieros, PRC 29,
#       1645 (1984)) plus the convection-current term (the default since 2026-10-23);
#   (c) E1 with --no-long-wavelength --no-ec-convection-current: the charge term alone;
#   (d) M1 (before 2026-10-23 every M1 term was 6.25 % low: sqrt(16 / 3) in integer arithmetic);
#   (e) M1 with --no-long-wavelength equals M1 without it (it used to be multiplied by the E1 charge function: x0.31);
#   alpha + 12C -> 16O(g.s.) + gamma, E1 + E2 (N = Z: the long-wavelength E1 charge nearly vanishes)
#   (f) long wavelength, (g) charge term, (h) charge + convection current;
#   saved intEC.dat files are not reused across
#   (i) the long-wavelength switch, (j) the convection-current switch.
#
#   ./tests/ec_m1_exact/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_m1_exact.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
OUT=output/AZUREOut_aa=1_R=3.out
TOL=2e-5

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run DIR MODEL [FLAGS] [EC FILE]: calculate with data in $WORK/DIR; EC FILE (e.g. output/intEC.dat) answers the
# "External Capture Amplitude File" prompt, so a saved file is offered for reuse
run() {
  local dir="$WORK/$1" model="$2" flags="${3:-}" ec="${4:-}"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$model.azr" "$dir/t.azr"
  cp -r "$HERE/models/data" "$dir/"
  # shellcheck disable=SC2086
  (cd "$dir" && printf '1\n\n%s\n7\n' "$ec" | run_limited 300 "$AZURE2_BIN" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}

# compare FILE COLUMN LABEL [TRUTH]: cross sections (column 4 of the output) against column COLUMN of the truth
# (expected/sigma.dat unless TRUTH names another file in expected/)
compare() {
  local file="$1" col="$2" label="$3" truth="${4:-sigma.dat}" res
  if [ ! -s "$file" ]; then bad "$label: no output"; return; fi
  res=$(awk -v col="$col" -v tol="$TOL" '
    NR == FNR { if ($1 !~ /^#/) t[++n] = $col; next }
    NF >= 4 && $4 > 0 { m++; d = ($4 - t[m]) / t[m]; if (d < 0) d = -d; if (d > w) w = d }
    END { if (m != n) { print "points " m " vs " n; exit } ; printf "%.2e%s", w, (w <= tol ? "" : " > tol") }
  ' "$HERE/expected/$truth" "$file")
  case "$res" in
    *tol*|points*) bad "$label: max rel. deviation $res" ;;
    *) ok "$label: max rel. deviation $res" ;;
  esac
}

echo "ec_m1_exact: $AZURE2_BIN"
run e1_lwa e1;  run e1_nlw e1 --no-long-wavelength;  run e1_chg e1 "--no-long-wavelength --no-ec-convection-current"
run m1_lwa m1;  run m1_nlw m1 --no-long-wavelength
compare "$WORK/e1_lwa/$OUT" 2 "(a) E1, long wavelength"
compare "$WORK/e1_nlw/$OUT" 4 "(b) E1, --no-long-wavelength (charge + convection current)"
compare "$WORK/e1_chg/$OUT" 3 "(c) E1, --no-long-wavelength --no-ec-convection-current (charge term)"
compare "$WORK/m1_lwa/$OUT" 5 "(d) M1, long wavelength"
compare "$WORK/m1_nlw/$OUT" 5 "(e) M1, --no-long-wavelength"
OUTA=output/AZUREOut_aa=1_R=2.out
run a_lwa a12c;  run a_chg a12c "--no-long-wavelength --no-ec-convection-current";  run a_nlw a12c --no-long-wavelength
compare "$WORK/a_lwa/$OUTA" 2 "(f) alpha+12C E1+E2, long wavelength" sigma_a12c.dat
compare "$WORK/a_chg/$OUTA" 3 "(g) alpha+12C E1+E2, charge term" sigma_a12c.dat
compare "$WORK/a_nlw/$OUTA" 4 "(h) alpha+12C E1+E2, charge + convection current" sigma_a12c.dat

# same4 A B: column 4 of two output files identical, row by row
same4() { awk 'NR == FNR { a[FNR] = $4; n = FNR; next } { if ($4 != a[FNR]) bad = 1; m = FNR } END { exit (bad || m != n) }' "$1" "$2"; }
SIGWARN="different integration grid or different inputs"

# (i) reuse across the long-wavelength switch: long-wavelength run first, then --no-long-wavelength in the same directory
run e1_reuse e1
run e1_reuse e1 --no-long-wavelength output/intEC.dat
if grep -q "$SIGWARN" "$WORK/e1_reuse/run.log" && [ -s "$WORK/e1_reuse/$OUT" ] &&
   same4 "$WORK/e1_nlw/$OUT" "$WORK/e1_reuse/$OUT"; then   # no cmp on the Windows runners
  ok "(i) intEC.dat of a long-wavelength run is refused (signature warning) and recomputed with --no-long-wavelength"
else
  bad "(i) intEC.dat of a long-wavelength run was reused with --no-long-wavelength"
fi
# (j) the same across the convection-current switch
run a_reuse a12c --no-long-wavelength
run a_reuse a12c "--no-long-wavelength --no-ec-convection-current" output/intEC.dat
if grep -q "$SIGWARN" "$WORK/a_reuse/run.log" && [ -s "$WORK/a_reuse/$OUTA" ] && same4 "$WORK/a_chg/$OUTA" "$WORK/a_reuse/$OUTA"; then
  ok "(j) intEC.dat with the convection current is refused and recomputed with --no-ec-convection-current"
else
  bad "(j) intEC.dat with the convection current was reused with --no-ec-convection-current"
fi

if [ "$failures" -eq 0 ]; then echo "ec_m1_exact: PASS"; else echo "ec_m1_exact: FAIL ($failures)"; exit 1; fi
