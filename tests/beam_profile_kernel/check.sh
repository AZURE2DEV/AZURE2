#!/usr/bin/env bash
#
# External capture is the same in the formal and the Brune basis (cnuc_split.md section 8).
# The two bases are two parametrizations of one R-matrix, so every observable must agree.  An external-capture
# final state enters the capture amplitude through its asymptotic normalization, which follows from its Brune
# amplitude (the residue of its own pole); its formal amplitude is a mixture of the poles of its J^pi group.  Here
# the 6.049 MeV 0+ final state shares its group with the ground state and two higher levels, and the formal basis
# used to take its formal amplitude: the capture to it differed from the Brune basis by a factor of order 1.
#
#   12C(a,g) to the 6.049 MeV 0+ (exit key 8) and to the ground state (key 4) at 90 degrees, 10 lab energies
#   E_alpha = 1.0-2.8 MeV, calculated without data in the Brune basis and with --no-brune; the cross sections
#   must agree to 1e-9 relative (printed with 11 digits).
#
# Both with --no-ignore-externals and with the default --ignore-externals, which drops the channel capture through a
# level whose photon width to the final state is given as zero.  That mask is defined on the input (Brune) levels and
# carried to the formal levels by the projector T^-1 P T (docs/upgrade/external_capture.md section 3); before step 4
# it was a test |gamma| < 1e-8 on the amplitudes of the basis in use, and the default flags differed by 9.7e-2.
#
#   ./tests/beam_profile_kernel/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/beam_profile_kernel.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

for ext in no-ignore-externals ignore-externals; do
for basis in brune no-brune; do
  D="$WORK/$ext/$basis"
  mkdir -p "$D/output" "$D/checks"
  cp -r "$HERE/data" "$D/"
  # the records the converter writes for the old test lines "1 1 8 1.0 3.0 0.2 90 90 0 1" and "1 1 4 ...": the
  # energies as the old grid computed them, folded by the beams and detector of the data sets 1 and 2
  awk '/^<extrapolations>/ { print; while ((getline l < recs) > 0) print l; skip = 1; next }
       /^<\/extrapolations>/ { skip = 0 } !skip' recs="$HERE/extrapolations.txt" "$HERE/beam_profile_kernel.azr" > "$D/t.azr"
  flag=""
  [ "$basis" = no-brune ] && flag=--no-brune
  (cd "$D" && run_limited 600 "$AZURE2_BIN" extrapolate --no-gui --$ext $flag t.azr > run.log 2>&1)
done

for r in 8 4; do
  f="AZUREOut_aa=1_R=$r.extrap"
  B="$WORK/$ext/brune/output/$f"; F="$WORK/$ext/no-brune/output/$f"
  if [ ! -s "$B" ] || [ ! -s "$F" ]; then
    bad "--$ext R=$r: no output (see run.log)"
    continue
  fi
  rel=$(paste "$B" "$F" | awk -v n="$(awk 'NR==1{print NF}' "$B")" '
    NF == 0 { next }
    { a = $4; b = $(n + 4); d = a - b; if (d < 0) d = -d; m = (a < 0 ? -a : a); k++
      if (m > 0 && d / m > w) w = d / m; else if (m == 0 && d > 0) w = 1 }
    END { printf "%.2e %d", w, k }')
  read -r worst count <<< "$rel"
  if [ "${count:-0}" -ge 10 ] && awk -v w="$worst" 'BEGIN { exit !(w ~ /^ *[-+]?\.?[0-9]/ && w <= 1e-9) }'; then
    ok "--$ext R=$r: formal and Brune capture agree, largest relative difference $worst over $count energies"
  else
    bad "--$ext R=$r: formal and Brune capture differ by $worst (relative) over $count energies"
  fi
done
done

[ "$failures" -eq 0 ] && echo "PASS" || echo "FAIL ($failures)"
exit $((failures > 0))
