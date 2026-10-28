#!/usr/bin/env bash
#
# The excitation-energy column includes the pair's excitation energy (audit 2026-10-25, 01 B3; 2026-10-28).
#
# EPoint::ConvertLabEnergy / ConvertExcitationEnergy / ConvertDecayEnergy set excitation_energy_ = E_cm + S, while the
# physics uses E_cm + S + E* of the pair (CalcEDependentValues).  The E_x column of AZUREOut_* (and the API's
# excitation energies) was therefore low by E* for a segment whose entrance pair is an excited-state pair -- the
# inverse-channel analyzing-power segment the docs recommend for a measured outgoing polarization, for one.
#
#   tests/7Li_p_ay (pair 3 = n + 7Be*, E* = 0.429 MeV, S = 18.898 MeV): a test segment 3 -> 3, angle-integrated,
#   E_lab 1.0-2.0 MeV.  (a) column 2 of AZUREOut_aa=3_R=3.extrap equals E_cm + 18.898 + 0.429 at every energy
#   (old build: + 18.898 only, 0.429 low); (b) the ground-state pair 5 (E* = 0) is unchanged: E_x = E_cm + 17.25.
#
#   ./tests/excitation_energy_column/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LI="$HERE/../7Li_p_ay"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/excitation_energy_column.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

mkdir -p "$WORK/output" "$WORK/checks"; cp -r "$LI/data" "$WORK/"
awk '/<segmentsTest>/{print; print "1 3 3 1.0 2.0 0.25 0 0 0 0"; print "1 5 5 1.0 2.0 0.25 0 0 0 0"; skip=1; next} /<\/segmentsTest>/{skip=0} !skip{print}' "$LI/7Li_p_ay.azr" > "$WORK/m.azr"
(cd "$WORK" && printf '3\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
# column 1 = E_cm, column 2 = E_x: check E_x - E_cm against S + E*
colcheck() {  # FILE EXPECTED-OFFSET LABEL
  local f="$1" want="$2" label="$3"
  if [ -s "$f" ] && awk -v w="$want" 'NF { d = $2 - $1 - w; if (d < 0) d = -d; if (d > 1e-6) bad = 1; n++ } END { exit (bad || n == 0) }' "$f"; then
    ok "$label: E_x - E_cm = $want MeV at every energy ($(awk 'NF' "$f" | wc -l | tr -d ' ') energies)"
  else bad "$label: E_x - E_cm = $(awk 'NF { printf "%.4f ", $2 - $1 }' "$f" 2> /dev/null | cut -c1-40) (expected $want)"; fi
}
colcheck "$WORK/output/AZUREOut_aa=3_R=3.extrap" 19.327 "(a) excited-state pair n + 7Be* (S 18.898 + E* 0.429)"
colcheck "$WORK/output/AZUREOut_aa=5_R=5.extrap" 17.25 "(b) ground-state pair p + 7Li (S 17.25)"

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
