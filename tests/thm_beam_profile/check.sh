#!/usr/bin/env bash
#
# The beam-profile kernel (a beam record's profile) folds a THM (HOES)
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
. "$HERE/../lib/edit_guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-600}")"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_beam_profile.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
OUT="AZUREOut_aa=1_R=2.out"
SIGMA=0.018
TOL=1e-4

# project NAME AWK [window]: a copy of the example with its records edited by the awk program AWK (the edit must change
# the file: tests/lib/edit_guard.sh); "window" adds the +-1 eV window columns to the data.
project() {
  local d="$WORK/$1"
  mkdir -p "$d/data" "$d/output" "$d/checks"
  awk -v w="${3:-}" '/^#/ || !NF { next }
       { if (w != "") printf "%s %s %s %s %.9f %.9f\n", $1, $2, $3, $4, $1 - 1e-6, $1 + 1e-6
         else print }' "$SRC/data/thm_lc2008_hoes_bkgsub.dat" > "$d/data/thm_lc2008_hoes_bkgsub.dat"
  edited "$SRC/o18_lacognata2008.azr" "$d/run.azr" awk "$2" || { echo "  FAIL  $1: the edit of the records hit nothing"; exit 1; }
  (cd "$d" && $RUN "$AZURE2_BIN" calculate --no-gui run.azr > log 2>&1)
  [ -s "$d/output/$OUT" ] || { echo "  FAIL  $1: no $OUT"; tail -5 "$d/log" | sed 's/^/        /'; exit 1; }
}

# G: the beam's spread set to SIGMA (lab).
GAUSS='$1 == "spread" { $0 = "    spread         '"$SIGMA"' MeV lab" } { print }'
# B: the beam a profile of one flat component (xi 0.1 MeV, omega 100 MeV), a detector of resolution SIGMA (lab), and
# the data's window columns (what AZURE2 --convert writes for the old line "... beamprofile 1 0.1 100 0 1 SIGMA 0 0").
BEAM='/^beam beam_thm_lc2008_hoes_bkgsub$/ { beamrec = 1; print
        print "    profile          skewed_gaussian xi 0.1 MeV omega 100 MeV alpha 0 weight 1"
        print "    truncation       0"; print "    detailed_balance no"; next }
      beamrec && !/^ / { beamrec = 0 }
      beamrec { next }
      $1 == "columns" { $0 = "    columns        E_lab[MeV] _ value error window_lo window_hi" }
      /^    beam +beam_thm_lc2008_hoes_bkgsub$/ { print; print "    detector       det_window"; next }
      /^<\/targets>/ { print ""; print "detector det_window"; print "    energy_resolution '"$SIGMA"' MeV lab" }
      { print }'
project gauss "$GAUSS"
project beam "$BEAM" window
grep -q '^    profile ' "$WORK/beam/run.azr" && grep -q '^    detector       det_window' "$WORK/beam/run.azr" &&
  ! grep -q '^    spread' "$WORK/beam/run.azr" || { echo "  FAIL  beam: the beam record is not the profile"; exit 1; }

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
