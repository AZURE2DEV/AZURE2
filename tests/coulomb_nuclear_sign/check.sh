#!/usr/bin/env bash
#
# Sign of the Coulomb-nuclear interference in elastic scattering, against an independent partial-wave sum.
#
# p + 12C hard-sphere scattering (includeAllJPi, l <= 4, every level width zero) at E_lab = 1-3 MeV and
# theta = 60-160 deg.  expected/dsdo.dat is |f_C + f_N|^2 from truth.py (mpmath Coulomb functions, textbook
# amplitudes, no AZURE2 convention).  Reversing the interference term 2 Re[f_C^* f_N] changes the cross section
# by 7-74 % on this grid (tolerance 1e-6), so a reversed sign of the Coulomb amplitude or of the i in the
# interference term fails it; so does any error in the hard-sphere or Coulomb phases (sigma_l - sigma_0).
# Checked on the A-matrix and R-matrix routes, with and without --use-brune (a zero-width model does not depend
# on the parametrization, so all four must agree with the truth).
#
#   ./tests/coulomb_nuclear_sign/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/coulomb_nuclear_sign.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

TOL=1e-6  # relative; the two Coulomb-function routines agree to ~1e-9 here

for route in true false; do
  for flags in "--use-brune" "--no-brune"; do
    tag="$([ "$route" = true ] && echo A-matrix || echo R-matrix) $flags"
    dir="$WORK/run_${route}${flags}"
    mkdir -p "$dir/output" "$dir/checks"
    sed "1,2s/^\(true\|false\) #Perform A-Matrix/$route #Perform A-Matrix/" "$HERE/models/hard_sphere_pc.azr" > "$dir/t.azr"
    (cd "$dir" && printf '3\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" $flags --no-gui --no-readline t.azr > run.log 2>&1)
    out="$dir/output/AZUREOut_aa=1_R=1.extrap"
    if [ ! -s "$out" ]; then bad "$tag: no extrapolation written"; tail -15 "$dir/run.log"; continue; fi
    # rows with theta > 0, in file order, against expected/dsdo.dat (same order)
    res="$(awk 'NF >= 4 && $3 > 0 { print $1, $3, $4 }' "$out" | paste - <(grep -v '^#' "$HERE/expected/dsdo.dat") | awk -v tol="$TOL" '
      { n++; if (NF != 6) { bad++; next }
        de = $1 - $4; if (de < 0) de = -de; dt = $2 - $5; if (dt < 0) dt = -dt;
        if (de > 1e-6 * $4 || dt > 1e-6 * $5) { grid++ }
        r = ($3 - $6) / $6; if (r < 0) r = -r; if (r > worst) worst = r }
      END { printf "%d %d %d %.3e\n", n, bad + 0, grid + 0, worst }')"
    read -r n nbad ngrid worst <<< "$res"
    nexp="$(grep -vc '^#' "$HERE/expected/dsdo.dat")"
    if [ "$n" != "$nexp" ] || [ "$nbad" != 0 ] || [ "$ngrid" != 0 ]; then
      bad "$tag: grid differs from expected ($n rows vs $nexp, $nbad malformed, $ngrid energy/angle mismatches)"
    elif awk -v w="$worst" -v t="$TOL" 'BEGIN { exit !(w ~ /^ *[-+]?\.?[0-9]/ && w + 0 <= t + 0) }'; then
      ok "$tag: $n points, max rel diff $worst"
    else
      bad "$tag: max rel diff $worst > $TOL"
    fi
  done
done

if [ "$failures" -eq 0 ]; then echo "PASS: Coulomb-nuclear interference sign"; else echo "FAIL: $failures check(s)"; exit 1; fi
