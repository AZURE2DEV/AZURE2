#!/usr/bin/env bash
#
# A time-of-flight resolution kernel (UDR) at low energy: identities that any such kernel must satisfy, on the
# Koehler & Graff 17O(n,a) points (25 meV - 1 MeV, a 1/v cross section below a few keV) in the Oliva 2024 17O+n
# level scheme (model.azr; README.md).
#
#   (a) a 0.2 ns delta kernel returns the bare cross section at every point (no zeros: a kernel narrower than the
#       mesh resolves must still integrate);
#   (b) a symmetric rectangular channel in time -- 7.55 bins per decade, i.e. +-7.6 % of the flight time, 0.24 ms at
#       0.1 eV -- returns the bare value on the 1/v region, because sigma ~ t is linear in the delay and the
#       hat-function weights integrate a linear function exactly;
#   (c) a 64 ns channel under a 125 ns burst is indistinguishable from bare below 10 keV.
#
# Before the repair (docs/upgrade/energy_grid.md, "Time-of-flight windows at low energy") the window's low edge was
# clamped to 1 keV and a window narrower than 1 eV was a single node at its high edge: zeros at 17, 7 and 3 of the 67
# points in (a), (b), (c).  From dev's tests/udr_low_energy (a68474f0), model converted with --convert.
#
#   ./tests/udr_low_energy/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/edit_guard.sh"
RUN="$(guard_command 600)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/udr_low_energy.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
OUT="output/AZUREOut_aa=1_R=2.out"

# with_udr BURST_NS CHANNEL FILE: the model with a detector that folds the data set with the kernel (a negative
# CHANNEL is |CHANNEL| bins per decade of energy, as on the old <targetInt> line)
with_udr() {
  awk -v burst="$1" -v channel="$2" '
    { print }
    /^    norm / { print "    detector       det_koehler1991"; print "    quadrature     points 200 width_multiplier 20 per_width 50" }
    /^<targets>/ { print "detector det_koehler1991"
                   print "    udr            file udr_delta.dat flight_path 7 m burst " burst " ns channel " channel " ns centred yes" }' "$3"
}
# run NAME [BURST CHANNEL] -> $WORK/NAME/$OUT
run() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -r "$HERE/data" "$dir/"; cp "$HERE/udr_delta.dat" "$dir/"
  if [ "$#" -eq 3 ]; then
    edited "$HERE/model.azr" "$dir/t.azr" with_udr "$2" "$3" || { bad "$1: the model was not edited"; return 1; }
  else
    cp "$HERE/model.azr" "$dir/t.azr"
  fi
  (cd "$dir" && $RUN "$AZURE2_BIN" calculate t.azr < /dev/null > run.log 2>&1)
  if [ ! -s "$dir/$OUT" ]; then bad "$1: no output written"; tail -15 "$dir/run.log"; return 1; fi
}
# ratio NAME EMAX TOL: fails if a point with E_cm <= EMAX has |folded/bare - 1| > TOL, or any point is zero or not
# finite, or the 67 points are not all there
ratio() {
  paste "$WORK/bare/$OUT" "$WORK/$1/$OUT" > "$WORK/$1.both"
  awk -v emax="$2" -v tol="$3" -v name="$1" '
    NF >= 18 { e = $1; b = $4; f = $13; n++
      if (!(f ~ /^ *[-+]?\.?[0-9]/) || f + 0 == 0 || !(b ~ /^ *[-+]?\.?[0-9]/) || b + 0 == 0) { zeros++; if (zeros <= 3) printf "    %s: zero or not a number at E_cm %.3e (bare %s)\n", name, e, b }
      else if (e <= emax) { r = f / b - 1; if (r < 0) r = -r; if (r > worst) { worst = r; we = e } } }
    END { printf "    %s: %d points, %d zero or not finite, worst |ratio-1| = %.2e at E_cm %.3e (E <= %g, tol %g)\n", name, n, zeros, worst, we, emax, tol
          exit !(zeros == 0 && worst <= tol && n == 67) }' "$WORK/$1.both"
}

run bare || exit 1
run delta 0 0 && { ratio delta 10 1e-4 && ok "(a) 0.2 ns delta kernel reproduces the bare cross section at all 67 points" || bad "(a) delta kernel"; }
run wide 0 -7.55 && { ratio wide 1e-5 3e-3 && ok "(b) 7.55 bins/decade channel returns the bare value on the 1/v region" || bad "(b) wide channel"; }
run chan 125 64 && { ratio chan 1e-2 1e-3 && ok "(c) 64 ns channel under a 125 ns burst is bare below 10 keV" || bad "(c) 64 ns channel"; }

echo "udr_low_energy: $failures failure(s)"
[ "$failures" -eq 0 ]
