#!/usr/bin/env bash
#
# UDR kernel at low energy: identities that any time-of-flight resolution
# function must satisfy, on the Koehler & Graff 17O(n,a) points (25 meV - 1 MeV,
# a 1/v cross section below a few keV) in the Oliva 2024 17O+n level scheme.
#
#   * a 0.2 ns delta kernel returns the bare cross section at every point (no
#     zeros: a kernel narrower than the sub-point grid must still integrate);
#   * a symmetric rectangular channel in time -- 7.55 bins per decade, i.e.
#     +-7.6 % of the flight time, 0.24 ms at 0.1 eV -- returns the bare value on
#     the 1/v region, because sigma ~ t is linear in the delay and the hat-function
#     weights integrate a linear function exactly;
#   * a 64 ns channel under a 125 ns burst is indistinguishable from bare below 10 keV.
#
# Before 2026-11-07 the sub-point window's low edge was clamped to 1 keV, so every
# point below 1 keV (and every kernel wide enough to reach below it) got an
# INVERTED grid: a bias of half the channel width on 1/v (0.924 at 0.1 eV) and
# exactly zero whenever the kernel's support ended at the first sub-point.
#
#   ./tests/udr_low_energy/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/udr_low_energy.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME "targetInt line or empty"  -> $WORK/NAME/output/AZUREOut_aa=1_R=2.out
run() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -r "$HERE/data" "$dir/"; cp "$HERE/udr_delta.dat" "$dir/"
  if [ -n "$2" ]; then
    awk -v line="$2" '{print} /^<targetInt>/{print line}' "$HERE/model.azr" > "$dir/t.azr"
  else
    cp "$HERE/model.azr" "$dir/t.azr"
  fi
  (cd "$dir" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --use-brune --no-gui --no-readline t.azr > run.log 2>&1)
  if [ ! -s "$dir/output/AZUREOut_aa=1_R=2.out" ]; then bad "$1: no output written"; tail -15 "$dir/run.log"; return 1; fi
}
# ratio NAME EMAX TOL -> fails if any point with E_cm <= EMAX has |fit/bare - 1| > TOL, or a zero/non-finite value anywhere
ratio() {
  paste "$WORK/bare/output/AZUREOut_aa=1_R=2.out" "$WORK/$1/output/AZUREOut_aa=1_R=2.out" | awk -v emax="$2" -v tol="$3" -v name="$1" '
    NF >= 18 { e = $1; b = $4; f = $13; n++
      if (f == 0 || f != f || f == "inf" || f == "-inf") { zeros++; if (zeros <= 3) printf "    %s: zero/non-finite at E_cm %.3e (bare %.4e)\n", name, e, b }
      else if (e <= emax) { r = f / b - 1; if (r < 0) r = -r; if (r > worst) { worst = r; we = e } } }
    END { printf "    %s: %d points, %d zero/non-finite, worst |ratio-1| = %.2e at E_cm %.3e (E <= %g, tol %g)\n", name, n, zeros, worst, we, emax, tol
          exit !(zeros == 0 && worst <= tol && n == 67) }'
}
udr() { printf '1 "1" 200 0 0 0 0 "" 0 0 0 0 "" 0 0 0.04 20 50 udr "udr_delta.dat" 7 %s %s 1' "$1" "$2"; }

run bare "" || exit 1
run delta "$(udr 0 0)" && { ratio delta 10 1e-4 && ok "0.2 ns delta kernel reproduces the bare cross section at all 67 points" || bad "delta kernel"; }
run wide  "$(udr 0 -7.55)" && { ratio wide 1e-5 3e-3 && ok "7.55 bins/decade channel returns the bare value on the 1/v region" || bad "wide channel"; }
run chan  "$(udr 125 64)" && { ratio chan 1e-2 1e-3 && ok "64 ns channel under a 125 ns burst is bare below 10 keV" || bad "64 ns channel"; }

echo "udr_low_energy: $failures failure(s)"
[ "$failures" -eq 0 ]
