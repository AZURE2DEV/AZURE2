#!/usr/bin/env bash
#
# Lab -> c.m. conversion of capture (gamma-ray) angular distributions (EPoint::ConvertLabAngleGammas and
# ConvertCrossSectionGammas), models from make_models.py (12C(a,g)16O(g.s.), E1 + E2 external capture, strongly
# asymmetric):
#
#   1. cm.azr gives the c.m. differential cross section sigma*(theta*) at E_lab x theta*;
#   2. those are boosted to the lab here, independently of AZURE2: with beta the c.m. velocity
#      (sqrt(T (T + 2 m1 c^2)) / ((m1 + m2) c^2 + T)), a photon at theta* appears at
#      cos theta = (cos theta* + beta)/(1 + beta cos theta*) with
#      dsigma/dOmega_lab = sigma* (1 + beta cos theta*)^2 / (1 - beta^2)   (the inverse of the solid-angle Jacobian);
#   3. lab.azr reads them as lab-frame data (isDiff 1), errors 1e-3 relative:
#      (a) chi2 = 0 -- AZURE2 takes every point back to the same theta* and sigma*;
#      (b) the lab-frame model column of the output equals the lab data (the fit is taken back with the same factor).
# Before 2026-10-23 the cross-section factor was (1 - beta)^2/(1 - beta cos theta): (a) fails with chi2/N ~ 1100 (errors 0.1 %, conversion off by up to 3 %).
#
#   ./tests/gamma_lab_frame/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/gamma_lab_frame.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
OUT=output/AZUREOut_aa=1_R=2.out

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

run() {   # run MODEL: calculate with data (mode 1) in $WORK/run_MODEL
  local dir="$WORK/run_$1"
  mkdir -p "$dir/output" "$dir/checks" "$dir/data"
  cp "$HERE/models/$1.azr" "$dir/t.azr"
  cp "$HERE/models/data/"*.dat "$dir/data/" 2> /dev/null
  [ -f "$WORK/lab.dat" ] && cp "$WORK/lab.dat" "$dir/data/"
  [ "$1" = lab ] && echo 1 > "$dir/kinflag.dat"   # lab-frame output files (EData::WriteOutputFiles)
  (cd "$dir" && printf '1\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
}

echo "gamma_lab_frame: $AZURE2_BIN"
run cm
if [ ! -s "$WORK/run_cm/$OUT" ]; then bad "c.m. model gave no output"; echo "gamma_lab_frame: FAIL"; exit 1; fi

# boost: output columns E_cm Ex angle sigma ... ; the data file's lab energies are in models/data/cm.dat
paste "$HERE/models/data/cm.dat" "$WORK/run_cm/$OUT" | awk -v m1=4.002603 -v m2=12.0 -v u=931.494088 '
  NF >= 8 {
    T = $1; ts = $2 * atan2(0, -1) / 180; s = $8
    b = sqrt(T * (T + 2 * m1 * u)) / ((m1 + m2) * u + T)
    c = (cos(ts) + b) / (1 + b * cos(ts)); th = atan2(sqrt(1 - c * c), c) * 180 / atan2(0, -1)
    sl = s * (1 + b * cos(ts))^2 / (1 - b * b)
    printf "%.12g  %.12g  %.12e  %.12e\n", T, th, sl, 1e-3 * sl
  }' > "$WORK/lab.dat"
n=$(wc -l < "$WORK/lab.dat")
[ "$n" -gt 0 ] || { bad "no boosted points"; exit 1; }

run lab
chi2=$(awk -F, 'NR > 1 && $1 + 0 == 1 { print $2 }' "$WORK/run_lab/output/chiSquared.out")
if [ -z "$chi2" ]; then
  bad "(a) no chi2 for the lab segment"
elif awk -v c="$chi2" -v n="$n" 'BEGIN { exit !(c / n < 1e-6) }'; then
  ok "(a) lab-frame data boosted from the c.m. model: chi2 = $chi2 over $n points"
else
  bad "(a) lab-frame data boosted from the c.m. model: chi2 = $chi2 over $n points (expected ~0)"
fi
# (b) output columns: E_lab Ex angle model_lab S_model data_lab err ...
dev=$(awk 'NF >= 7 && $6 != 0 { d = ($4 - $6) / $6; if (d < 0) d = -d; if (d > w) w = d } END { printf "%.2e", w }' "$WORK/run_lab/$OUT")
if awk -v d="$dev" 'BEGIN { exit !(d < 1e-6) }'; then
  ok "(b) lab-frame model column equals the lab data (max rel. deviation $dev)"
else
  bad "(b) lab-frame model column differs from the lab data (max rel. deviation $dev)"
fi

if [ "$failures" -eq 0 ]; then echo "gamma_lab_frame: PASS"; else echo "gamma_lab_frame: FAIL ($failures)"; exit 1; fi
