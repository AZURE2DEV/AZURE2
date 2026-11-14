#!/usr/bin/env bash
# External capture into an UNBOUND final state (E_f > 0 above the entrance pair's threshold), treated as a narrow
# resonance since 2026-11-13 (ECIntegral::IntegrateUnbound): the final state's exterior is its outgoing wave
# O(k_f r)/O(k_f a), O = G + iF, its strength the level's observed width in that channel, and the N_f factor
# (hbar^2/2 mu a) dS/dE at E_f.  Before, such a state was taken as BOUND by |E_f| without a word, and every point
# below the photon threshold (E < E_f) came out NaN.
#
# Model: tests/formal_ec_final_state/models/A_orig.azr (12C + p, 13N) with a capture pair 3 into the 1/2+ level at
# E_x = 2.3682 MeV, 424.69 keV above the proton threshold (S = 1.94351 MeV): the 12C(p,g1) transition of Rolfs
# and Azuma, NPA 227, 291 (1974).  E1 lines from the 1/2- and 3/2- levels; no internal gamma widths into it.
#   (a) a NOTE names the unbound final state; below the photon threshold the capture is 0, never NaN;
#   (b) pure direct capture (every entrance-side particle width 0, final width 1 eV so N_f = 1): divided by the
#       same transition into a mirror BOUND state at E_x = S - 0.42469 MeV (with E_gamma^3 and g^2 taken out),
#       it is |I_O|^2 / I_W^2, which reference_integral.py computes independently (damping factor exp(-eps r),
#       Richardson in eps, instead of the integration by parts of the code): within 0.5 % (measured 0.18 %);
#   (c) continuity through the threshold: the same final state 10, 20 and 40 keV on either side, same reduced
#       width: sigma / (E_gamma^3 g^2) unbound over bound, minus 1, extrapolates to 0 at the threshold;
#   (d) the formal parametrization (--no-brune) gives the Brune capture;
#   (e) --use-gradient: the analytic gradient is used and agrees with finite differences (final width free).
#
#   ./tests/ec_unbound_final/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_unbound_final.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
BASE="$HERE/../formal_ec_final_state/models/A_orig.azr"
S=1.94351
EF=2.3682

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

run() {  # NAME FLAGS: mode 3 (extrapolation) of $WORK/NAME/t.azr
  (cd "$WORK/$1" && printf '3\nn\n\n\n7\n' | run_limited 900 "$AZURE2_BIN" $2 --no-gui --no-readline t.azr > run.log 2>&1)
}
CAP3="output/AZUREOut_aa=1_R=3.extrap"
# model NAME EX WIDTH ZERO SEGMENTS...: A_orig with capture pair 3 into the 1/2+ level moved to EX (MeV; the level's
# proton line gets WIDTH -- eV if unbound, ANC if bound), every other 1/2-, 3/2- proton width zeroed if ZERO = 1,
# and the given <segmentsTest> lines.  Capture into pair 3 is E1 only (mask 1): with E2 as well, d-wave capture
# runs through the 5/2+ resonance at E_cm = 1.60 MeV and the pure-direct comparison of (b) is no longer E1 alone.
model() {
  local name="$1" ex="$2" w="$3" zero="$4"; shift 4
  mkdir -p "$WORK/$name/output" "$WORK/$name/checks"
  awk -v ex="$ex" -v w="$w" -v zero="$zero" -v segs="$(printf '%s\n' "$@")" '
    /^<levels>/ { lev = 1; print; next }
    /^<\/levels>/ { lev = 0; print; next }
    /^<segmentsTest>/ { print; printf "%s\n", segs; skip = 1; next }
    skip && /^<\/segmentsTest>/ { skip = 0; print; next }
    skip { next }
    lev && NF >= 31 && $1 == "0.5" && $2 == "1" && $3 == "2.3682" {
      $3 = ex; if ($27 == 0) $12 = w; print " " $0; next }
    lev && NF >= 31 && $6 == 1 && $27 == 0 && zero == 1 && ($1 == "0.5" || $1 == "1.5") && $2 == "-1" { $12 = 0; print " " $0; next }
    lev && NF >= 31 && $6 == 2 && $8 == 2 && (($1 == "0.5" && $2 == "-1" && $3 == "20") || ($1 == "1.5" && $2 == "-1")) {
      print; $6 = 3; $12 = 0; $15 = "0.5"; $16 = "1"; $17 = ex; $31 = 1; print " " $0; next }
    { print }' "$BASE" > "$WORK/$name/t.azr"
}
# the reported reduced width of the moved level's proton channel
gint() { awk -v ex="$2" '{ sub(/\r$/, "") } $1 == "J" && $3 == "0.5+" { on = (($6 - ex)^2 < 1e-8) ; next }
                          on && $1 == "R" && $3 == "1" { for (i = 1; i < NF; i++) if ($i == "g_int") { print $(i + 2); exit } }' "$WORK/$1/output/parameters.out"; }

# (a) ----------------------------------------------------------------------------------------------------------
model full "$EF" 34000 0 "1 1 3 0.40 1.0 0.01 0 0 0 0"
run full "--use-brune"
if grep -q "NOTE: external capture into pair 3: its final state E_x = 2.3682 MeV is UNBOUND, 424.69 keV above" "$WORK/full/run.log"; then
  ok "NOTE names the unbound final state (424.69 keV above the 12C + p threshold)"
else bad "no NOTE on the unbound final state ($(grep -a "NOTE\|ERROR" "$WORK/full/run.log" | head -1 | cut -c1-100))"; fi
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
SEGS=("1 1 3 0.57 0.57 0.01 0 0 0 0" "1 1 3 0.87 0.87 0.01 0 0 0 0" "1 1 3 1.30 1.30 0.01 0 0 0 0" "1 1 3 1.73 1.73 0.01 0 0 0 0")
MIRROR=1.51882
model pure_u "$EF" 1 1 "${SEGS[@]}"
model pure_b "$MIRROR" 0.017 1 "${SEGS[@]}"
run pure_u "--use-brune"; run pure_b "--use-brune"
gu="$(gint pure_u "$EF")"; gb="$(gint pure_b "$MIRROR")"
# |I_O|^2 / I_W^2 from reference_integral.py (2026-11-13; R-matrix/AZURE2_code_audit/11-13-26_ec_unbound_final/readme)
# (five damping factors, Richardson orders 3 and 4 agree to 2e-5 or better)
REF="0.5259 1641.158 0.8026 25.0766 1.1993 6.95073 1.5960 3.92679"
res="$(paste <(awk 'NF >= 4 { print $1, $(NF - 1) }' "$WORK/pure_u/$CAP3" 2>/dev/null) <(awk 'NF >= 4 { print $(NF - 1) }' "$WORK/pure_b/$CAP3" 2>/dev/null) |
  awk -v S="$S" -v eu="$EF" -v eb="$MIRROR" -v gu="$gu" -v gb="$gb" -v ref="$REF" '
    BEGIN { n = split(ref, r, " "); for (i = 1; i < n; i += 2) want[r[i]] = r[i + 1] }
    NF == 3 { e = sprintf("%.4f", $1); if (!(e in want) || $3 + 0 <= 0 || gu + 0 <= 0 || gb + 0 <= 0) next
              got = ($2 / $3) * ((($1 + S - eb) / ($1 + S - eu))^3) * (gb / gu)^2
              d = got / want[e] - 1; if (d < 0) d = -d; if (d > worst) worst = d; out = out sprintf(" %s:%.4g/%.4g", e, got, want[e]); m++ }
    END { printf "%d %.3e%s", m, worst, out }')"
set -- $res
if [ "${1:-0}" -eq 4 ] && awk -v w="$2" 'BEGIN { exit !(w < 0.005) }'; then ok "direct capture into the unbound state / mirror bound state = |I_O|^2/I_W^2 of the reference within 0.5 % (worst $2;${res#* * })"
else bad "direct capture ratio vs reference: ${1:-0} of 4 energies, worst $2 (g $gu, $gb;${res#* * })"; fi

# (c) ----------------------------------------------------------------------------------------------------------
# the same reduced width on either side, g ~ 0.01 (so N_f = 1 to 1e-4 and the exterior integrals alone are
# compared; N_f itself is continuous because dS/dE is): ANC (bound) and width in eV (unbound) scaled from runs at
# ANC 1 or 1e4 and widths 1e-12 / 1e-24 eV (g ~ ANC, g ~ sqrt(width)); the ratio uses each run's reported g.
# The difference d(delta) = ratio - 1 at delta = 10, 20, 40 keV is fitted with c0 + a delta + b delta^2; continuity
# means c0 = 0.  Measured 2026-11-13: d(10 keV) 0.03-0.12, c0 3e-4 to 2.6e-3 (at most 2 % of d(10 keV), the size
# the six printed digits of g allow).  Required: |c0| < 0.005 at every energy.
CSEG="1 1 3 0.6 2.0 0.2 0 0 0 0"
model c40b 1.90351 1.105 1 "$CSEG"; model c40u 1.98351 7.3e-9 1 "$CSEG"
model c20b 1.92351 31.6 1 "$CSEG"; model c20u 1.96351 5.05e-14 1 "$CSEG"
model c10b 1.93351 8898 1 "$CSEG"; model c10u 1.95351 2.74e-21 1 "$CSEG"
for m in c40b c40u c20b c20u c10b c10u; do run "$m" "--use-brune" & done; wait
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
if [ "${1:-0}" -ge 6 ] && awk -v c="$2" 'BEGIN { exit !(c < 0.005) }'; then
  ok "continuous through the threshold: unbound/bound - 1 extrapolates to |c0| <= $2 at delta = 0 over $1 energies (it is up to $3 at 10 keV)"
else bad "threshold continuity: ${1:-0} energies, extrapolated difference at delta = 0 up to $2 (limit 0.005)"; fi

# (d) ----------------------------------------------------------------------------------------------------------
mkdir -p "$WORK/full_nb/output" "$WORK/full_nb/checks"; cp "$WORK/full/t.azr" "$WORK/full_nb/"
run full_nb "--no-brune"
d="$(paste <(awk 'NF >= 4 { print $(NF - 1) }' "$WORK/full/$CAP3" 2>/dev/null) <(awk 'NF >= 4 { print $(NF - 1) }' "$WORK/full_nb/$CAP3" 2>/dev/null) |
     awk '$1 + 0 > 0 { r = ($2 - $1) / $1; if (r < 0) r = -r; if (r > m) m = r; n++ } END { if (n < 40) print "nan"; else printf "%.3e", m }')"
if awk -v d="$d" 'BEGIN { exit !(d ~ /^[0-9.e+-]+$/ && d + 0 <= 1e-7) }'; then ok "--no-brune gives the Brune capture into the unbound state (max rel $d)"
else bad "--no-brune vs Brune capture into the unbound state: max rel '$d'"; fi

# (e) ----------------------------------------------------------------------------------------------------------
# pseudo-data: the Brune capture of (a), 5 %, lab frame; then the final state's width moved 34 -> 30 keV
if [ -s "$WORK/full/$CAP3" ]; then
  mkdir -p "$WORK/grad/output" "$WORK/grad/checks" "$WORK/grad/data"
  awk 'NF >= 5 && $1 > 0.45 && (++n % 3) == 1 { printf "%.10e 0 %.10e %.10e\n", $1 * (1.00728 + 12) / 12, $(NF - 1), 0.05 * $(NF - 1) }' \
    "$WORK/full/$CAP3" > "$WORK/grad/data/pseudo.dat"
  awk -v seg="1 1 3 0 5 0 180 0 1 0 0 0 0 0 data/pseudo.dat 0 0" '
    /^<segmentsData>/ { print; print seg; next }
    !moved && $3 == "2.3682" && sub(/ 34000 /, " 30000 ") { moved = 1 }
    { print }' "$WORK/full/t.azr" > "$WORK/grad/t.azr"
  (cd "$WORK/grad" && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
     run_limited 900 "$AZURE2_BIN" --use-brune --use-gradient --no-gui --no-readline t.azr > grad.log 2>&1)
  line="$(awk '{ sub(/\r$/, "") } /^\[grad-check\] chi2 / { print; exit }' "$WORK/grad/grad.log")"
  case "$line" in
    *"analytic block used"*"energy 0/"*"gamma 0/"*) ok "--use-gradient: analytic gradient used, agrees with finite differences ($(echo "$line" | grep -o 'gamma 0/[0-9]*'))" ;;
    *) bad "--use-gradient: '$(echo "$line" | cut -c1-160)'" ;;
  esac
else bad "--use-gradient: no pseudo-data (no capture output in (a))"; fi

if [ "$failures" -eq 0 ]; then echo "PASS: external capture into an unbound final state"; else echo "FAIL: $failures check(s)"; exit 1; fi
