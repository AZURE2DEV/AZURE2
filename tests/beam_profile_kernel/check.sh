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
# --no-ignore-externals: the default drops the external resonant capture of a level whose internal photon amplitude
# is zero, which is defined on the levels of the basis in use and so differs between the bases (an open item of
# step 2, HANDOFF).
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

for basis in brune no-brune; do
  mkdir -p "$WORK/$basis/output" "$WORK/$basis/checks"
  cp -r "$HERE/data" "$WORK/$basis/"
  awk '/<segmentsTest>/ { print; print "1 1 8 1.0 3.0 0.2 90 90 0 1"; print "1 1 4 1.0 3.0 0.2 90 90 0 1"; skip = 1; next }
       /<\/segmentsTest>/ { skip = 0 } !skip' "$HERE/beam_profile_kernel.azr" > "$WORK/$basis/t.azr"
  flag=""
  [ "$basis" = no-brune ] && flag=--no-brune
  (cd "$WORK/$basis" && printf '3\nn\n\n\n7\n' |
    run_limited 600 "$AZURE2_BIN" --no-gui --no-readline --no-ignore-externals $flag t.azr > run.log 2>&1)
done

for r in 8 4; do
  f="AZUREOut_aa=1_R=$r.extrap"
  if [ ! -s "$WORK/brune/output/$f" ] || [ ! -s "$WORK/no-brune/output/$f" ]; then
    bad "R=$r: no output (see run.log)"
    continue
  fi
  rel=$(paste "$WORK/brune/output/$f" "$WORK/no-brune/output/$f" | awk -v n="$(awk 'NR==1{print NF}' "$WORK/brune/output/$f")" '
    NF == 0 { next }
    { a = $4; b = $(n + 4); d = a - b; if (d < 0) d = -d; m = (a < 0 ? -a : a); k++
      if (m > 0 && d / m > w) w = d / m; else if (m == 0 && d > 0) w = 1 }
    END { printf "%.2e %d", w, k }')
  read -r worst count <<< "$rel"
  if [ "${count:-0}" -ge 10 ] && awk -v w="$worst" 'BEGIN { exit !(w <= 1e-9) }'; then
    ok "R=$r: formal and Brune capture agree, largest relative difference $worst over $count energies"
  else
    bad "R=$r: formal and Brune capture differ by $worst (relative) over $count energies"
  fi
done

[ "$failures" -eq 0 ] && echo "PASS" || echo "FAIL ($failures)"
exit $((failures > 0))
