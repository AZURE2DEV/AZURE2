#!/usr/bin/env bash
#
# A convolution across a resonance whose widths are given as reduced-width
# amplitudes must converge as its integration grid is refined.
#
# The optional 33rd field of a level line (gammaIsRWA, used by THM projects)
# marks the value as an amplitude in MeV^1/2.  The adaptive grid used to sum it
# as a width in eV, sized a 46 eV resonance as 0.7 eV, packed its points into a
# sliver around the peak and bridged the flanks with keV-wide linear
# interpolation -- so the convolved cross section depended on the number of
# sub-points (tests/17O: chi2 13.8 -> 53 -> 88 at 150/600/2000 points).
#
# Model: tests/17O (a THM segment with a 21 keV Gaussian, RWA-flagged levels).
# The check runs it at 150 and 600 sub-points and requires the total chi^2 to
# agree to 1%.
#
#   ./tests/thm_rwa_grid/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../17O"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_rwa_grid.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
. "$HERE/../lib/guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-900}")"

run_with_points() {   # run_with_points N -> prints total chi2
  local d="$WORK/n$1"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$SRC/data" "$d/"
  local azr; azr="$(ls "$SRC"/*.azr | sed -n 1p)"
  # the number of sub-points of the fold: "points" of the quadrature field (the data set's and the extrapolation's)
  awk -v n="$1" '$1 == "quadrature" { sub(/points [0-9]+/, "points " n) } { print }' "$azr" > "$d/run.azr"
  [ "$(grep -c "^ *quadrature  *points $1 " "$d/run.azr")" -ge 1 ] || { echo "  FAIL  no quadrature field set to $1 points" >&2; exit 1; }
  (cd "$d" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr >/dev/null 2>&1)
  grep -oE 'Total-Chi-Squared: [0-9.eE+-]+' "$d/output/chiSquared.out" | awk '{print $2}'
}

c150="$(run_with_points 150)"
c600="$(run_with_points 600)"
echo "  chi2 at 150 points: $c150   at 600 points: $c600"
if [ -z "$c150" ] || [ -z "$c600" ]; then
  echo "  FAIL  a run produced no chi-squared"; exit 1
fi
if awk -v a="$c150" -v b="$c600" 'BEGIN{d=(a-b)/b; if(d<0)d=-d; exit !(d<0.01)}'; then
  echo "  ok    convolution converged (agreement < 1%)"; exit 0
else
  echo "  FAIL  convolution depends on the grid density"; exit 1
fi
