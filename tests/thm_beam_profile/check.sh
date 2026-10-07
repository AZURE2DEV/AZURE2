#!/usr/bin/env bash
#
# The beam-profile kernel (<targetInt> ... beamprofile ...) folds a THM (HOES)
# point below the entrance threshold the way the Gaussian fold does.
#
# Model: examples/o18_lacognata2008 (18O(p,a)15N THM, 20/90/144 keV levels,
# 31 points from E_cm = -45 keV).  The same data are folded twice:
#   G  the Gaussian fold, sigma = 18 keV (lab);
#   B  a beam-profile kernel with a flat beam (one component, omega = 100 MeV)
#      and a 2 eV energy window around each point (data columns 5-6), with the
#      same 18 keV resolution: K(E) is then the Gaussian to 1e-9.
# B must equal G at every point to 1e-4 (it does to ~1e-6).  The beam-profile window was clipped
# at E_cm = +1 keV, not at the THM floor: B was 66 times G at -45 keV, and it
# cut its resolution tails at 4 sigma, not the 5 sigma of the THM fold (-2 %
# at 225 keV, 81 keV above a 289 eV level).
#
#   ./tests/thm_beam_profile/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../../examples/o18_lacognata2008"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-600}")"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_beam_profile.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
OUT="AZUREOut_aa=1_R=2.out"
SIGMA=0.018
TOL=1e-4

# project NAME TARGETINT-LINE [window]: a copy of the example with the given
# <targetInt> line; "window" adds the +-1 eV window columns to the data.
project() {
  local d="$WORK/$1"
  mkdir -p "$d/data" "$d/output" "$d/checks"
  awk -v w="${3:-}" '/^#/ || !NF { next }
       { if (w != "") printf "%s %s %s %s %.9f %.9f\n", $1, $2, $3, $4, $1 - 1e-6, $1 + 1e-6
         else print }' "$SRC/data/thm_lc2008_hoes_bkgsub.dat" > "$d/data/thm_lc2008_hoes_bkgsub.dat"
  awk -v t="$2" '/^<targetInt>/ { print; print t; skip = 1; next }
       /^<\/targetInt>/ { skip = 0 } skip { next } { print }' \
    "$SRC/o18_lacognata2008.azr" > "$d/run.azr"
  (cd "$d" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr > log 2>&1)
  [ -s "$d/output/$OUT" ] || { echo "  FAIL  $1: no $OUT"; tail -5 "$d/log" | sed 's/^/        /'; exit 1; }
}

project gauss "1  \"1\"  50  1  $SIGMA  0  0  \"\" 0  0  0  0  \"\" 0  0 0.04 5 50"
project beam "1  \"1\"  50  0  0  0  0  \"\"  0  0  0  0  \"\"  0  0  0.04  5  50  beamprofile  1  0.1  100  0  1  $SIGMA  0  0" window

# Column 1 E_cm, column 4 the folded model.  (\r: Windows line ends.)
paste "$WORK/gauss/output/$OUT" "$WORK/beam/output/$OUT" | awk -v tol="$TOL" '
  { gsub(/\r/, "") }
  NF < 8 { next }
  { h = NF / 2; e = $1; g = $4; b = $(h + 4)
    if (e != $(h + 1) || !(g > 0)) { printf "  FAIL  point %d: energies %s / %s, model %s\n", n + 1, e, $(h + 1), g; bad++; next }
    n++; r = b / g - 1; if (r < 0) r = -r
    if (e < 0) below++
    if (r > worst) { worst = r; at = e }
    if (r > tol) { bad++; if (bad <= 5) printf "  point at E_cm %7.2f keV: beam profile %.6e, Gaussian %.6e (ratio %.4g)\n", e * 1000, b, g, b / g } }
  END {
    printf "  %d points (%d below threshold): max |beam profile / Gaussian - 1| = %.2e at E_cm %.1f keV\n", n, below, worst, at * 1000
    if (n < 30 || below < 4) { print "  FAIL  too few points compared"; exit 1 }
    if (bad) { print "  FAIL  the beam-profile kernel does not fold THM points as the Gaussian does"; exit 1 }
    print "  ok    beam-profile kernel == Gaussian fold at every THM point, below threshold too" }'
