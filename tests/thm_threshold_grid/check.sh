#!/usr/bin/env bash
#
# A THM excitation function evaluated on a fine grid through the entrance
# threshold is finite and smooth.
#
# THM (HOES) segments may carry points at and just above E_cm = 0, and the
# Gaussian folding of a THM segment spans +-5 sigma, so its sub-point grid
# crosses the entrance threshold (EData.cpp lets THM sub-points go below it).
# A sub-point landing at round-off distance from E = 0 (E ~ -1e-16) used to
# get a NaN shift function (the plain-double Whittaker function underflows;
# for |E| < ~1e-5 MeV it even returned S = 0), and the folded cross section
# of every point whose grid hit it came out at 1e16-1e20 against neighbours
# of ~0.005 -- on 6Li_d every multiple of 9 keV below 0.11 MeV.  See
# ShftFunc::operator() and IsCoulombThreshold in EPoint.cpp.
#
# Models: tests/6Li_d (6Li+d, 30 keV sigma) and tests/7Li_p_a (p+7Li, 30 keV
# sigma), with their data replaced by a 1 keV grid over E_cm = 0-0.3 MeV (the
# first point at E_cm = 0 exactly).  Passes when every calculated value is
# finite and positive and no two neighbours differ by more than MAXRATIO.
#
#   ./tests/thm_threshold_grid/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_threshold_grid.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
if command -v timeout >/dev/null 2>&1; then RUN="timeout ${TEST_TIMEOUT:-900}"; else RUN=""; fi

# The smooth curve moves by < 0.5 % per keV here; the spike was 1e21.
MAXRATIO=1.02
fail=0

# grid NAME PROJECT DATAFILE LABPERCM
grid() {
  local d="$WORK/$1"
  mkdir -p "$d/output" "$d/checks" "$d/data"
  # E_cm = 0, 0.001, ..., 0.3 MeV in the lab convention of the data files.
  awk -v f="$4" 'BEGIN { for (i = 0; i <= 300; i++) printf "%.10e 0.0 1.0 0.1\n", i * 0.001 * f }' \
    > "$d/data/grid.dat"
  sed "s#$3#data/grid.dat#" "$HERE/../$2/$2.azr" > "$d/run.azr"
  (cd "$d" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr > log 2>&1)
  local out
  out="$(ls "$d"/output/AZUREOut_*.out 2>/dev/null | head -1)"
  if [ -z "$out" ]; then
    echo "  FAIL  $1: no AZUREOut file"; tail -5 "$d/log" | sed 's/^/        /'; fail=1; return
  fi
  # column 1 E_cm, column 4 the calculated (folded) cross section
  if awk -v maxr="$MAXRATIO" -v name="$1" '
    NF >= 4 { n++; e[n] = $1; y[n] = $4 }
    END {
      bad = 0; worst = 1; at = 0
      if (n < 250) { printf "  FAIL  %s: only %d points\n", name, n; exit 1 }
      for (i = 1; i <= n; i++) {
        if (!(y[i] > 0) || y[i] != y[i] || y[i] > 1e300) {
          printf "  FAIL  %s: E_cm = %s gives %s\n", name, e[i], y[i]; bad = 1
        }
        if (i > 1 && y[i] > 0 && y[i-1] > 0) {
          r = y[i] / y[i-1]; if (r < 1) r = 1 / r
          if (r > worst) { worst = r; at = e[i] }
        }
      }
      if (worst > maxr) {
        printf "  FAIL  %s: neighbours differ by x%.4g at E_cm = %s\n", name, worst, at; bad = 1
      }
      if (!bad)
        printf "  ok    %s: %d points from E_cm = %s, finite, max neighbour ratio %.5f\n", name, n, e[1], worst
      exit bad
    }' "$out"; then :; else fail=1; fi
}

grid 6Li_d 6Li_d data/pizzone_thm.dat 1.3333333333
grid 7Li_p_a 7Li_p_a data/tumino_thm.dat 1.1428571429

if [ "$fail" -eq 0 ]; then
  echo "  ok    THM grids through the entrance threshold are finite and smooth"; exit 0
else
  echo "  FAIL  THM threshold grid"; exit 1
fi
