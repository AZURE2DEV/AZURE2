#!/usr/bin/env bash
#
# An energy shift of a THM segment moves every point and every folding
# sub-point, below the entrance threshold too, and the sub-point grid follows
# the narrow level in the shifted frame.
#
# Model: tests/17O (17O(n,a) THM, points from E_lab = -54 keV, the 33 eV 5-
# level at E_cm = 74 keV, 21 keV Gaussian).  For a shift d (lab, MeV) the
# project is run twice:
#   A  segment_1_energy_shift = d in an external parameter file;
#   B  the data file's energies moved by d, no shift.
# The two must give the same model at every point (column 4 of the output;
# column 1 is the shifted E_cm in both) to 1e-4.  Two shifts: +6 and -4 keV.
#
# Before, ESegment::UpdatePointEnergiesWithShift moved only points and
# sub-points with E > 0: the four points below threshold kept their energies
# (27 % off) and the sub-points at E <= 0 stayed behind, tearing the fold of
# every point within reach of threshold; and EPoint::RefreshSubPointGrid did
# not rebuild the grid of a shifted segment, so the lattice built around the
# narrow level at the unshifted energies sat 6 keV off it (30 % off).
#
#   ./tests/thm_energy_shift/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../17O"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-600}")"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_energy_shift.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
OUT="AZUREOut_aa=1_R=2.out"
TOL=1e-4
fail=0

# run NAME [PARFILE]: calculate a copy of the project in $WORK/NAME, whose
# data/ may have been edited; PARFILE is named relative to the run directory
# (read from stdin, a POSIX path would not reach a native Windows binary).
setup() {
  local d="$WORK/$1"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$SRC/data" "$d/"
  cp "$SRC/17O.azr" "$d/run.azr"
}
run() {
  (cd "$WORK/$1" && printf '1\n%s\n\n7\n' "${2:-}" | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr > log 2>&1)
  [ -s "$WORK/$1/output/$OUT" ] || { echo "  FAIL  $1: no $OUT"; tail -5 "$WORK/$1/log" | sed 's/^/        /'; exit 1; }
}

setup base; run base
for d in 0.006 -0.004; do
  setup "A$d"
  awk -v d="$d" '{ sub(/\r$/, "") } $1 == "segment_1_energy_shift" { $2 = d } { print }' \
    "$WORK/base/output/param.par" > "$WORK/A$d/shift.par"
  grep -q "^segment_1_energy_shift $d " "$WORK/A$d/shift.par" ||
    { echo "  FAIL  no segment_1_energy_shift in param.par"; exit 1; }
  run "A$d" shift.par
  setup "B$d"
  for f in "$WORK/B$d"/data/*.dat; do
    awk -v d="$d" '/^#/ || !NF { print; next } { $1 = sprintf("%.9f", $1 + d); print }' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
  done
  run "B$d"
  paste "$WORK/A$d/output/$OUT" "$WORK/B$d/output/$OUT" | awk -v tol="$TOL" -v d="$d" '
    { gsub(/\r/, "") }
    NF < 8 { next }
    { h = NF / 2; e = $1; a = $4; b = $(h + 4); n++
      if (e < 0) below++
      de = e - $(h + 1); if (de < 0) de = -de
      if (de > 1e-6) { bad++; if (bad <= 5) printf "  point %d: E_cm %s shifted, %s moved in the file\n", n, e, $(h + 1); next }
      r = a / b - 1; if (r < 0) r = -r
      if (r > worst) { worst = r; at = e }
      if (r > tol) { bad++; if (bad <= 5) printf "  point at E_cm %8.2f keV: shifted %.6e, moved data %.6e (ratio %.4g)\n", e * 1000, a, b, a / b } }
    END {
      printf "  shift %+g MeV: %d points (%d below threshold), max |rel diff| %.2e at E_cm %.1f keV\n", d, n, below, worst, at * 1000
      exit !(n >= 20 && below >= 3 && !bad) }' || { echo "  FAIL  shift $d: the shifted segment differs from the moved data"; fail=1; }
done

if [ "$fail" -eq 0 ]; then echo "  ok    a THM energy shift moves every point and sub-point, and the grid follows"; fi
exit "$fail"
