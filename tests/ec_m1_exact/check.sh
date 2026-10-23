#!/usr/bin/env bash
#
# External capture, absolute: p + 16O -> 17F(0.495 MeV, 1/2+) + gamma, hard-sphere external capture only (models from
# make_models.py), against an independent golden-rule calculation (truth.py -> expected/sigma.dat):
#
#   (a) E1 in the long-wavelength approximation;
#   (b) E1 with --no-long-wavelength (the extended-Siegert charge term, Friar and Fallieros, PRC 29, 1645 (1984));
#   (c) M1 (before 2026-10-23 every M1 term was 6.25 % low: sqrt(16 / 3) in integer arithmetic);
#   (d) M1 with --no-long-wavelength equals M1 without it (it used to be multiplied by the E1 charge function: x0.31);
#   (e) a saved intEC.dat is not reused across the long-wavelength switch: a --no-long-wavelength run in the
#       output directory of a long-wavelength run recomputes and gives the fresh answer.
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

# compare FILE COLUMN(of expected/sigma.dat) LABEL: cross sections (column 4 of the output) against the truth
compare() {
  local file="$1" col="$2" label="$3" res
  if [ ! -s "$file" ]; then bad "$label: no output"; return; fi
  res=$(awk -v col="$col" -v tol="$TOL" '
    NR == FNR { if ($1 !~ /^#/) t[++n] = $col; next }
    NF >= 4 && $4 > 0 { m++; d = ($4 - t[m]) / t[m]; if (d < 0) d = -d; if (d > w) w = d }
    END { if (m != n) { print "points " m " vs " n; exit } ; printf "%.2e%s", w, (w <= tol ? "" : " > tol") }
  ' "$HERE/expected/sigma.dat" "$file")
  case "$res" in
    *tol*|points*) bad "$label: max rel. deviation $res" ;;
    *) ok "$label: max rel. deviation $res" ;;
  esac
}

echo "ec_m1_exact: $AZURE2_BIN"
run e1_lwa e1;  run e1_nlw e1 --no-long-wavelength
run m1_lwa m1;  run m1_nlw m1 --no-long-wavelength
compare "$WORK/e1_lwa/$OUT" 2 "(a) E1, long wavelength"
compare "$WORK/e1_nlw/$OUT" 3 "(b) E1, --no-long-wavelength"
compare "$WORK/m1_lwa/$OUT" 4 "(c) M1, long wavelength"
compare "$WORK/m1_nlw/$OUT" 4 "(d) M1, --no-long-wavelength"

# (e) reuse across the switch: long-wavelength run first, then --no-long-wavelength in the same directory
run e1_reuse e1
run e1_reuse e1 --no-long-wavelength output/intEC.dat
if grep -q "different integration grid or different inputs" "$WORK/e1_reuse/run.log" && [ -s "$WORK/e1_reuse/$OUT" ] && cmp -s <(awk '{print $4}' "$WORK/e1_reuse/$OUT") <(awk '{print $4}' "$WORK/e1_nlw/$OUT"); then
  ok "(e) intEC.dat of a long-wavelength run is refused (signature warning) and recomputed with --no-long-wavelength"
else
  bad "(e) intEC.dat of a long-wavelength run was reused with --no-long-wavelength"
fi

if [ "$failures" -eq 0 ]; then echo "ec_m1_exact: PASS"; else echo "ec_m1_exact: FAIL ($failures)"; exit 1; fi
