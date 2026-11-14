#!/usr/bin/env bash
#
# External capture finds the exit to the photon pair by the pair's number, not by taking the number as a position in
# the list of exits (docs/upgrade/cnuc_split.md section 11.9, defect 1).
#
# The exits of an entrance skip every pair it has no pathway to (a beta pair, a pair of another compound nucleus).
# Until 2026-11-11 the external-capture code read the exit list at the photon pair's number, so with such a pair
# numbered below the photon pair it read another exit or past the end of the list.  Two models of the same physics
# (models/, README.md): tests/13N's levels, four capture data sets and one elastic one, and a 13C level with an
# n + 12C channel (another compound nucleus, so no pathway from p + 12C):
#   neutron_pair_first.azr  the 13C level first: pairs n + 12C (1), p + 12C (2), photon (3); the photon exit is the
#                           second exit of p + 12C;
#   photon_pair_first.azr   the 13C level last: p + 12C (1), photon (2), n + 12C (3); position = number.
# Checks:
#   (a) the pathways check file of each model shows where the photon exit is: position 2 of pair 3, and 2 of 2;
#   (b) neutron_pair_first runs (the base read past the end of the exit list: an empty intEC.dat and SIGSEGV);
#   (c) its intEC.dat holds the amplitudes of every capture data set and component (1560 lines) and equals the other
#       model's, line for line, and so do its signature, every AZUREOut file and chiSquared.out (the outputs that list
#       the levels in file order -- param.par, parameters.out, scatteringLengths.out -- differ by that order only and
#       are not compared).
#
#   ./tests/ec_exit_position/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
RUN="$(guard_command 900)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_exit_position.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

for m in neutron_pair_first photon_pair_first; do
  mkdir -p "$WORK/$m/output" "$WORK/$m/checks"
  cp -r "$HERE/models/data" "$WORK/$m/"
  cp "$HERE/models/$m.azr" "$WORK/$m/"
  (cd "$WORK/$m" && $RUN "$AZURE2_BIN" calculate --no-gui "$m.azr" > run.log 2>&1)
  echo "$?" > "$WORK/$m/status"
done
A="$WORK/neutron_pair_first"
B="$WORK/photon_pair_first"

# (a) the external pathways of the check file: entrance pair, exit position ("Decay #"), exit pair ("Decay Pair #")
exits() { awk '/External Reaction Pathways/ { f = 1; next } f && $1 ~ /^[0-9]+$/ { print $1, $2, $3 }' "$1/checks/pathways.chk" | sort -u; }
[ "$(exits "$A")" = "2 2 3" ] && ok "neutron_pair_first: the photon pair (3) is the second exit of p + 12C (pair 2)" \
  || bad "neutron_pair_first: external pathways (entrance, position, pair): '$(exits "$A" | tr '\n' ';')', want '2 2 3'"
[ "$(exits "$B")" = "1 2 2" ] && ok "photon_pair_first: the photon pair (2) is the second exit of p + 12C (pair 1)" \
  || bad "photon_pair_first: external pathways (entrance, position, pair): '$(exits "$B" | tr '\n' ';')', want '1 2 2'"

# (b)
for d in "$A" "$B"; do
  if [ "$(cat "$d/status")" = 0 ] && [ -s "$d/output/chiSquared.out" ]; then
    ok "$(basename "$d") runs"
  else
    bad "$(basename "$d"): exit status $(cat "$d/status"), no chiSquared.out"; tail -3 "$d/run.log"
  fi
done

# (c)
lines="$(wc -l < "$A/output/intEC.dat" 2>/dev/null | tr -d ' ')"
[ "${lines:-0}" = 1560 ] && ok "neutron_pair_first: intEC.dat holds 1560 amplitudes" \
  || bad "neutron_pair_first: intEC.dat holds '${lines:-none}' amplitudes, want 1560 ((32 + 32 + 32 + 32 + 67) points x 8 pathways)"
compared=0
for f in "$B"/output/AZUREOut_* "$B/output/chiSquared.out" "$B/output/intEC.dat" "$B/output/intEC.dat.sig"; do
  name="$(basename "$f")"
  if cmp -s "$f" "$A/output/$name"; then
    compared=$((compared + 1))
  else
    bad "$name differs between the two orders of the pairs"
  fi
done
[ "$compared" -ge 6 ] && ok "$compared output files equal in the two orders of the pairs" || bad "only $compared files compared"

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
