#!/usr/bin/env bash
# External capture into an UNBOUND final state (E_f > 0 above the entrance pair's threshold), treated as a narrow
# resonance (ECIntegral::IntegrateUnbound): the final state's exterior is its outgoing wave O(k_f r)/O(k_f a),
# O = G + iF, its strength the level's observed width in that channel, and the N_f factor (hbar^2/2 mu a) dS/dE at E_f.
#
# Models (models/, written by AZURE2 --convert from the models dev's test builds): 12C + p (13N) with a capture pair
# 3 into the 1/2+ level at E_x = 2.3682 MeV, 424.69 keV above the proton threshold (S = 1.94351 MeV): the 12C(p,g1)
# transition of Rolfs and Azuma, NPA 227, 291 (1974).  E1 lines from the 1/2- and 3/2- levels; no internal gamma
# widths into it.
#   (a) full: a NOTE names the unbound final state; below the photon threshold the capture is 0, never NaN;
#   (b) pure direct capture (every entrance-side particle width 0, final width 1 eV so N_f = 1; pure_u): divided by
#       the same transition into a mirror BOUND state at E_x = S - 0.42469 MeV (pure_b; with E_gamma^3 and g^2 taken
#       out), it is |I_O|^2 / I_W^2, computed independently (a damping factor exp(-eps r), Richardson in eps, instead
#       of the integration by parts of the code; the numbers REF below): within 0.5 % (measured 0.18 %);
#   (c) continuity through the threshold: the same final state 10, 20 and 40 keV on either side (c10b ... c40u), same
#       reduced width: sigma / (E_gamma^3 g^2) unbound over bound, minus 1, extrapolates to 0 at the threshold;
#   (d) the formal parametrization (--no-brune) gives the Brune capture;
#   (e) fit --use-gradient (grad: pseudo-data of (a), the final state's width moved): the analytic gradient is used
#       and agrees with finite differences.
#
#   ./tests/ec_unbound_final/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_unbound_final.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
S=1.94351
EF=2.3682

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME MODEL COMMAND [VAR=VALUE ...]: models/MODEL.azr in $WORK/NAME -> run.log
run() {
  local name="$1" model="$2" command="$3"
  shift 3
  mkdir -p "$WORK/$name/output" "$WORK/$name/checks"
  cp "$HERE/models/$model.azr" "$WORK/$name/t.azr"
  cp -R "$HERE/models/data" "$WORK/$name/data"
  (cd "$WORK/$name" && env "$@" $(guard_command 900) "$AZURE2_BIN" $command t.azr < /dev/null > run.log 2>&1)
}
CAP3="output/AZUREOut_aa=1_R=3.extrap"
# the reported reduced width of the proton channel of the 1/2+ level at EX
gint() { awk -v ex="$2" '{ sub(/\r$/, "") } $1 == "J" && $3 == "0.5+" { on = (($6 - ex)^2 < 1e-8) ; next }
                          on && $1 == "R" && $3 == "1" { for (i = 1; i < NF; i++) if ($i == "g_int") { print $(i + 2); exit } }' "$WORK/$1/output/parameters.out"; }

# (a) ----------------------------------------------------------------------------------------------------------
run full full extrapolate
if grep -q "NOTE: external capture into pair 3: its final state E_x = 2.3682 MeV is UNBOUND, 424.69 keV above" "$WORK/full/run.log"; then
  ok "NOTE names the unbound final state (424.69 keV above the 12C + p threshold)"
else bad "no NOTE on the unbound final state ($(grep -a "NOTE\|ERROR" "$WORK/full/run.log" | sed -n 1p | cut -c1-100))"; fi
if [ ! -s "$WORK/full/$CAP3" ]; then bad "no capture output into pair 3 ($(tail -1 "$WORK/full/run.log" | cut -c1-100))"
else
  r="$(awk -v t=0.42469 'NF >= 4 { x = $(NF - 1); if (x != x || x ~ /nan/) nan++; else if ($1 < t && x + 0 != 0) below++; else if ($1 > t + 0.005 && x + 0 > 0) above++; n++ }
          END { printf "%d %d %d %d", n, nan + 0, below + 0, above + 0 }' "$WORK/full/$CAP3")"
  set -- $r
  if [ "$2" -eq 0 ] && [ "$3" -eq 0 ] && [ "$4" -gt 40 ]; then ok "capture: $1 points, none NaN, 0 below the photon threshold, positive above ($4)"
  else bad "capture: $1 points, $2 NaN, $3 nonzero below the photon threshold, $4 positive above"; fi
fi

# (b) ----------------------------------------------------------------------------------------------------------
# E_lab 0.57, 0.87, 1.30, 1.73 MeV = E_cm 0.5259, 0.8026, 1.1993, 1.5960 MeV
MIRROR=1.51882
run pure_u pure_u extrapolate; run pure_b pure_b extrapolate
gu="$(gint pure_u "$EF")"; gb="$(gint pure_b "$MIRROR")"
# |I_O|^2 / I_W^2 of the independent integral (five damping factors, Richardson orders 3 and 4 agree to 2e-5 or better)
REF="0.5259 1641.158 0.8026 25.0766 1.1993 6.95073 1.5960 3.92679"
res="$(paste <(awk 'NF >= 4 { print $1, $(NF - 1) }' "$WORK/pure_u/$CAP3" 2>/dev/null) <(awk 'NF >= 4 { print $(NF - 1) }' "$WORK/pure_b/$CAP3" 2>/dev/null) |
  awk -v S="$S" -v eu="$EF" -v eb="$MIRROR" -v gu="$gu" -v gb="$gb" -v ref="$REF" '
    BEGIN { n = split(ref, r, " "); for (i = 1; i < n; i += 2) want[r[i]] = r[i + 1] }
    NF == 3 { e = sprintf("%.4f", $1); if (!(e in want) || !($3 + 0 > 0) || !(gu ~ /^ *[-+]?\.?[0-9]/ && gu + 0 > 0) || !(gb ~ /^ *[-+]?\.?[0-9]/ && gb + 0 > 0)) next
              got = ($2 / $3) * ((($1 + S - eb) / ($1 + S - eu))^3) * (gb / gu)^2
              d = got / want[e] - 1; if (d < 0) d = -d; if (d > worst) worst = d; out = out sprintf(" %s:%.4g/%.4g", e, got, want[e]); m++ }
    END { printf "%d %.3e%s", m, worst, out }')"
set -- $res
if [ "${1:-0}" -eq 4 ] && awk -v w="$2" 'BEGIN { exit !(w ~ /^ *[-+]?\.?[0-9]/ && w + 0 < 0.005) }'; then ok "direct capture into the unbound state / mirror bound state = |I_O|^2/I_W^2 of the reference within 0.5 % (worst $2;${res#* * })"
else bad "direct capture ratio vs reference: ${1:-0} of 4 energies, worst ${2:-?} (g $gu, $gb;${res#* * })"; fi

# (c) ----------------------------------------------------------------------------------------------------------
# the same reduced width on either side, g ~ 0.01 (so N_f = 1 to 1e-4 and the exterior integrals alone are compared;
# N_f itself is continuous because dS/dE is).  The difference d(delta) = ratio - 1 at delta = 10, 20, 40 keV is fitted
# with c0 + a delta + b delta^2; continuity means c0 = 0.  Measured on dev 2026-11-13: d(10 keV) 0.03-0.12, c0 3e-4
# to 2.6e-3 (at most 2 % of d(10 keV), the size the six printed digits of g allow).  Required: |c0| < 0.005.
for m in c40b c40u c20b c20u c10b c10u; do run "$m" "$m" extrapolate; done
ratio() {  # BOUND UNBOUND EXB EXU -> per energy, sigma/(E_g^3 g^2) unbound over bound
  local g1 g2; g1="$(gint "$1" "$3")"; g2="$(gint "$2" "$4")"
  paste <(awk 'NF >= 4 { print $1, $(NF - 1) }' "$WORK/$1/$CAP3" 2>/dev/null) <(awk 'NF >= 4 { print $(NF - 1) }' "$WORK/$2/$CAP3" 2>/dev/null) |
    awk -v S="$S" -v xb="$3" -v xu="$4" -v gb="$g1" -v gu="$g2" 'NF == 3 && $2 + 0 > 0 && $3 + 0 > 0 && gb + 0 > 0 && gu + 0 > 0 {
      printf "%.4f %.6f\n", $1, ($3 / (($1 + S - xu)^3 * gu^2)) / ($2 / (($1 + S - xb)^3 * gb^2)) }'
}
cres="$(paste <(ratio c10b c10u 1.93351 1.95351) <(ratio c20b c20u 1.92351 1.96351) <(ratio c40b c40u 1.90351 1.98351) |
  awk 'NF == 6 { d1 = $2 - 1; d2 = $4 - 1; d4 = $6 - 1
                 b = ((d4 - d2) - 2 * (d2 - d1)) / 600; a = (d2 - d1 - 300 * b) / 10; c0 = d1 - 10 * a - 100 * b
                 c = c0 < 0 ? -c0 : c0; if (c > wc) wc = c; dd = d1 < 0 ? -d1 : d1; if (dd > w1) w1 = dd; n++ }
       END { printf "%d %.2e %.3f", n, wc, w1 }')"
set -- $cres
if [ "${1:-0}" -ge 6 ] && awk -v c="$2" 'BEGIN { exit !(c ~ /^ *[-+]?\.?[0-9]/ && c + 0 < 0.005) }'; then
  ok "continuous through the threshold: unbound/bound - 1 extrapolates to |c0| <= $2 at delta = 0 over $1 energies (it is up to $3 at 10 keV)"
else bad "threshold continuity: ${1:-0} energies, extrapolated difference at delta = 0 up to ${2:-?} (limit 0.005)"; fi

# (d) ----------------------------------------------------------------------------------------------------------
run full_nb full "extrapolate --no-brune"
d="$(paste <(awk 'NF >= 4 { print $(NF - 1) }' "$WORK/full/$CAP3" 2>/dev/null) <(awk 'NF >= 4 { print $(NF - 1) }' "$WORK/full_nb/$CAP3" 2>/dev/null) |
     awk '$1 + 0 > 0 { r = ($2 - $1) / $1; if (r < 0) r = -r; if (r > m) m = r; n++ } END { if (n < 40) print "nan"; else printf "%.3e", m }')"
if awk -v d="$d" 'BEGIN { exit !(d ~ /^ *[-+]?\.?[0-9]/ && d + 0 <= 1e-7) }'; then ok "--no-brune gives the Brune capture into the unbound state (max rel $d)"
else bad "--no-brune vs Brune capture into the unbound state: max rel '$d'"; fi

# (e) ----------------------------------------------------------------------------------------------------------
run grad grad "fit --use-gradient" AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000
line="$(awk '{ sub(/\r$/, "") } /^\[grad-check\] chi2 / { print; exit }' "$WORK/grad/run.log")"
case "$line" in
  *"analytic block used"*"energy 0/"*"gamma 0/"*) ok "--use-gradient: analytic gradient used, agrees with finite differences ($(echo "$line" | grep -o 'gamma 0/[0-9]*' | sed -n 1p))" ;;
  *) bad "--use-gradient: '$(echo "$line" | cut -c1-160)'"; grep -a "grad-check\|ERROR" "$WORK/grad/run.log" | sed -n '1,4p' | cut -c1-200 ;;
esac

if [ "$failures" -eq 0 ]; then echo "PASS: external capture into an unbound final state"; else echo "FAIL: $failures check(s)"; exit 1; fi
