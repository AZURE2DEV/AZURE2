#!/usr/bin/env bash
#
# Gaussian kernel window: a point whose beam-spread Gaussian reaches a narrow resonance 10^4 x the local cross
# section must get the full Gaussian average.  Before 2026-10-26 every kernel was cut at +-3 sigma and renormalised
# (TargetEffect::convolutionRange = 3): the 0.27 % of weight beyond 3 sigma, and the resonance with it, was dropped --
# a point 3.5 sigma from the peak lost the peak entirely, one at 2.5 sigma lost 2-4 % (17O+a ODeSA targets).  The
# window is 5 sigma now (lost weight 5.7e-7).
#
#   bare.azr : 16O(p,g) angle-integrated, no external capture, a 1 keV-wide 1/2- resonance at E_lab 0.8508 MeV on a
#              ~1e-7 b background (fine grid);  conv.azr : the points E_res + k sigma_lab, k = +-1.5 ... +-4, with a
#              3 keV (lab) Gaussian beam spread.  expected/truth.dat = full-Gaussian convolution of the bare output
#              (truth.py; the bare output does not depend on the kernel).
#   (a) every point within 2e-3 of the truth;  (b) the bare run itself is unchanged (peak and background).
#
#   ./tests/convolution_window/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/convolution_window.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
OUT=output/AZUREOut_aa=1_R=3.extrap

run() {  # MODEL: mode 3 in $WORK/MODEL
  local dir="$WORK/$1"; mkdir -p "$dir/output" "$dir/checks"; cp "$HERE/models/$1.azr" "$dir/t.azr"
  (cd "$dir" && printf '3\nn\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  [ -s "$dir/$OUT" ] || { bad "$1: no output"; tail -5 "$dir/run.log"; }
}

echo "convolution_window: $AZURE2_BIN"
run bare; run conv
# (b) the bare excitation function: peak height and position against expected/bare_peak
read pk pe bg < "$HERE/expected/bare_peak"
res=$(awk -v pk="$pk" -v pe="$pe" -v bg="$bg" 'NF >= 4 && $4 > m { m = $4; e = $1 } NF >= 4 && NR < 50 { s += $4; n++ }
  END { d1 = (m - pk) / pk; if (d1 < 0) d1 = -d1; d2 = e - pe; if (d2 < 0) d2 = -d2; d3 = (s / n - bg) / bg; if (d3 < 0) d3 = -d3;
        printf "%.1e %.1e %.1e", d1, d2, d3 }' "$WORK/bare/$OUT")
set -- $res
if awk -v a="$1" -v b="$2" -v c="$3" 'BEGIN { exit !(a <= 1e-4 && b <= 2e-5 && c <= 1e-4) }'; then ok "(b) bare resonance unchanged (peak $1, position $2 MeV, background $3)"
else bad "(b) bare resonance moved (peak $1, position $2 MeV, background $3)"; fi
# (a) convolved points against the full-Gaussian truth
worst=$(paste <(grep -v '^#' "$HERE/expected/truth.dat") <(awk 'NF >= 4 { print $4 }' "$WORK/conv/$OUT") | awk '
  NF >= 4 { d = ($4 - $3) / $3; if (d < 0) d = -d; printf "k=%s %.2e\n", $1, d; if (d > w) { w = d; wk = $1 } }
  END { printf "WORST %s %.2e\n", wk, w }')
echo "$worst" | grep -v WORST | sed 's/^/        /' | paste -sd' ' | fold -w 110 | sed 's/^/        /'
w=$(echo "$worst" | awk '/^WORST/ { print $3 }'); wk=$(echo "$worst" | awk '/^WORST/ { print $2 }')
if awk -v w="$w" 'BEGIN { exit !(w <= 2e-3) }'; then ok "(a) convolved points within 2e-3 of the full-Gaussian truth (worst $w at k=$wk)"
else bad "(a) convolved points off the full-Gaussian truth: worst $w at k=$wk (> 2e-3)"; fi
if [ "$failures" -eq 0 ]; then echo "convolution_window: PASS"; else echo "convolution_window: FAIL ($failures)"; exit 1; fi
