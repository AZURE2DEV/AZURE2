#!/usr/bin/env bash
#
# E1-M1 interference of external capture in a differential cross section, against an independent m-scheme
# calculation (audit 2026-10-25, 04_external_capture_photons S5 / checklist t3-26; 2026-10-27).
#
# The E external-capture amplitude carries -i i^(li + L - lf) and the M1 one +i (ECIntegral.cpp); angle-integrated
# cross sections never see the relative phase, so the E-only and M-only tests of tests/ec_m1_exact cannot tell whether
# E1-M1 interference in an angular distribution has the right sign.  truth_dsdo.py computes dsigma/dOmega of
# p + 16O capture from the dipole interaction -d.E - mu.B between the Coulomb + hard-sphere scattering state (beam
# along z, each partial wave with its Coulomb and hard-sphere phase) and the ANC-normalized bound state, summed over
# magnetic substates and photon helicities: no multipole phase convention enters.  Two final states
# (make_models.py):
#   s : 17F(0.4953, 1/2+), E1 from l = 1 and M1 (spin term) from l = 0.  With a spin-0 target and an unpolarized beam
#       these CANNOT interfere (E1 spin-diagonal, the spin operator traceless): dsigma/dOmega(E1 + M1) = E1 + M1 exactly.
#   p : a fictitious 1/2- state (l_f = 1), E1 from l = 0 and 2, M1 (orbital + spin) from l = 1: the orbital M1
#       amplitude interferes with E1 in the odd Legendre orders (a forward/backward asymmetry).
# expected/dsdo_{s,p}.dat hold the truth at E_cm = 0.3, 1.0, 2.0 MeV and seven angles, E1 only, M1 only, E1 + M1.
#
# For each case, with the photon pair's multipolarity mask 1 (E1), 2 (M1), 3 (both) and one c.m. differential capture
# segment (isDiff 4) at those energies and angles:
#   (a) E1 only and (b) M1 only equal the truth at every point (2e-5 of the largest value of the row's energy: the
#       differential capture machinery and units);
#   (c) E1 + M1 equals the truth at every point (2e-5): the relative phase of the E and M amplitudes.
#
#   ./tests/ec_e1m1_phase/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_e1m1_phase.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
TOL=2e-5

run() {  # run CASE MASK
  local d="$WORK/$1_$2"
  mkdir -p "$d/output" "$d/checks" "$d/data"
  cp "$HERE/models/data/dsdo.dat" "$d/data/"
  sed 's/\r$//' "$HERE/models/$1_$2.azr" > "$d/t.azr"
  (cd "$d" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  awk '{ sub(/\r$/, "") } NF >= 5 { printf "%.4f %.1f %.10e\n", $1, $3, $4 }' "$d/output/AZUREOut_aa=1_R=3.out" 2> /dev/null > "$d.fit"
}
# compare CASE MASK COLUMN LABEL: every row of the run against column COLUMN of the case's truth; the difference is
# relative to the largest truth value at that energy (a point where the truth vanishes, E1 at 0 and 180 deg, is held
# to the same absolute tolerance)
compare() {
  local case="$1" mask="$2" col="$3" label="$4" exp="$HERE/expected/dsdo_$1.dat" fit="$WORK/$1_$2.fit"
  local n; n="$(wc -l < "$fit" | tr -d ' ')"
  [ "$n" -eq 21 ] || { bad "$label: the run wrote $n points, expected 21"; tail -3 "$WORK/$1_$2/run.log"; return; }
  local worst; worst="$(awk -v col="$col" '
    NR == FNR { if ($0 !~ /^#/) { key = sprintf("%.4f %.1f", $1, $2); t[key] = $col; e = sprintf("%.4f", $1); v = $col < 0 ? -$col : $col; if (v > big[e]) big[e] = v } next }
    { key = sprintf("%.4f %.1f", $1, $2); if (!(key in t)) { print "missing"; exit 2 }
      e = sprintf("%.4f", $1); d = ($3 - t[key]) / big[e]; if (d < 0) d = -d; if (d > m) m = d }
    END { printf "%.3e", m }' "$exp" "$fit")"
  if [ -n "$worst" ] && awk -v m="$worst" -v tol="$TOL" 'BEGIN { exit !(m <= tol) }'; then ok "$label: AZURE2 = truth at all 21 points (largest difference $worst of the row maximum)"; else bad "$label: largest difference '$worst' of the row maximum (tolerance $TOL)"; paste "$fit" <(awk '!/^#/' "$exp") | awk -v col="$col" '{ printf "    E %s theta %s: AZURE2 %s truth %s\n", $1, $2, $3, $(3 + col) }' | head -7; fi
}

for case in s p; do
  [ -f "$HERE/expected/dsdo_$case.dat" ] || { bad "no expected/dsdo_$case.dat (run truth_dsdo.py)"; continue; }
  for mask in 1 2 3; do run "$case" "$mask"; done
  compare "$case" 1 3 "case $case (a) E1 only"
  compare "$case" 2 4 "case $case (b) M1 only"
  compare "$case" 3 5 "case $case (c) E1 + M1 (the relative phase)"
done

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
