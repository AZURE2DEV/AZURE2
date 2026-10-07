#!/usr/bin/env bash
#
# Coulomb line shape of the THM spectator: `lineshape=on` on an
# `experiment[<name>]` line of the <thm> block (ThmLineshape.h;
# docs/source/theory/thm_implementation.rst, "Coulomb line shape").  Each
# level's exit amplitude is multiplied by
#   N_C = exp(pi zeta/2) (E_l - E - i G_l/2)^(-i zeta),
#   |N_C|^2 = exp[2 zeta arctan(2 (E_l - E)/G_l)]    (Mukhamedzhanov et al. 2020 eq. 62)
# with zeta = eta_sB - eta_0 = Z_s alpha (Z_B mu_sB - Z_F mu_sF)/k_sF and
# E_sF = E_aA - B - E.
#
# Model: tests/18O_p_a_thm (19F: p+18O -> a+15N, two 1/2+ levels, two THM
# segments with free norms).  Its real reaction, 2H(18O,a15N)n, has a neutron
# spectator; a charged one is made up as 18O(3He,a15N)d at 115 MeV (x = p,
# s = d, B(3He) = 5.4934 MeV; E_sF ~ 10.4 MeV, zeta ~ -0.14).
#
#   (a) lineshape=off changes nothing: byte-identical to the same line without
#       the key, and one-segment experiments with lineshape=off to no block;
#   (b) a neutral spectator (zeta = 0): lineshape=on is byte-identical to off;
#   (c) one level, no folding: model(on)/model(off) = exp[2 zeta arctan(2 (E_l
#       - E)/G_l)] at every point to 1e-6, zeta from the masses here and
#       G_l from parameters.out; the same with a second level of the same
#       J^pi that interferes but does not decay to a+15N (the factor belongs
#       to the decaying level, inside the coherent sum);
#   (d) zeta < 0 (Z_B = 7 < Z_F = 9): the peak moves up in energy, the model
#       is lowered below E_l and raised above it (fine grid, one level);
#   (e) refused: lineshape=on without kinematics, a value other than on/off,
#       E_sF <= 0 at a data point, the formal (non-Brune) parameterization.
#
#   ./tests/thm_lineshape/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../18O_p_a_thm"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_lineshape.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/check_common.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-300}")"
OUT="AZUREOut_aa=1_R=2.out"
KIN="beam=18O target=3He spectator=d Ebeam=115"
fail=0

# run NAME BLOCK [EDITS [FLAGS]] -- a "calculate with data" run of a copy of
# the project with BLOCK as its <thm> block ("" = none), the .azr first passed
# through the awk programs named in EDITS (variable names, in order), extra
# AZURE2 FLAGS.  Status in NAME/status.
run() {
  local d="$WORK/$1"
  rm -rf "$d"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$SRC/data" "$d/"
  cp "$SRC/18O_p_a_thm.azr" "$d/run.azr"
  local edit
  for edit in ${3:-}; do awk "${!edit}" "$d/run.azr" > "$d/tmp.azr" && mv "$d/tmp.azr" "$d/run.azr"; done
  [ -z "$2" ] || printf '\n<thm>\n%s\n</thm>\n' "$2" >> "$d/run.azr"
  (cd "$d" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline ${4:-} run.azr 2>&1 | head -c 1000000 > log;
   echo "${PIPESTATUS[1]}" > status)
}

# awk programs editing the .azr
NOFOLD='/<targetInt>/ { print; T = 1; next } /<\/targetInt>/ { T = 0 } T { next } { print }'
ONELEVEL='/<levels>/ { L = 1 } /<\/levels>/ { L = 0 } L && $3 == "8.805800" { next } { print }'
# level 2 kept, its a+15N width set to zero
NOEXIT2='/<levels>/ { L = 1 } /<\/levels>/ { L = 0 } L && $3 == "8.805800" && $6 == 2 { $12 = 0 } { print }'

# (a) ----------------------------------------------------------------------
echo "(a) lineshape=off == no key"
run nokey "experiment[A] segments=1,2 $KIN"
run off "experiment[A] segments=1,2 $KIN lineshape=off"
ran nokey && ran off
for f in chiSquared.out normalizations.out parameters.out thm_experiments.out "$OUT"; do
  same "$WORK/nokey/output/$f" "$WORK/off/output/$f" && ok "$f byte-identical" || bad "$f differs"
done
run base ""
run offone "experiment[A] segments=1 lineshape=off
experiment[B] segments=2 lineshape=off"
ran base && ran offone
for f in chiSquared.out normalizations.out parameters.out "$OUT"; do
  same "$WORK/base/output/$f" "$WORK/offone/output/$f" && ok "one-segment, lineshape=off: $f == no <thm> block" \
    || bad "one-segment, lineshape=off: $f differs from no block"
done

# (b) ----------------------------------------------------------------------
echo "(b) neutral spectator: N_C = 1"
run noff "experiment[A] segments=1,2 beam=18O target=d spectator=n Ebeam=54"
run non "experiment[A] segments=1,2 beam=18O target=d spectator=n Ebeam=54 lineshape=on"
ran noff && ran non
for f in chiSquared.out normalizations.out parameters.out "$OUT"; do
  same "$WORK/noff/output/$f" "$WORK/non/output/$f" && ok "spectator n: $f byte-identical" || bad "spectator n: $f differs"
done
z="$(tr -d '\r' < "$WORK/non/output/thm_experiments.out" | awk '$1 == "zeta[2]" { print $2, $3 }')"
[ "$z" = "0.0000000000e+00 0.0000000000e+00" ] && ok "thm_experiments.out: zeta = 0" || bad "zeta row '$z'"
grep -q "the spectator is neutral, so N_C = 1" "$WORK/non/log" && ok "log says N_C = 1" || bad "no neutral-spectator note"

# (c) ----------------------------------------------------------------------
# zeta(E) from the masses: nuclear masses (u) of the built-in table (3He, d,
# p, 18O) and of the .azr pairs (p + 18O, a + 15N); B = m_p + m_d - m_3He.
ratio_check() {  # ratio_check TAG ON OFF
  local on="$WORK/$2/output" off="$WORK/$3/output"
  local G EL
  # Level 1 (8.6026 MeV): total width from parameters.out (keV -> MeV).
  G="$(tr -d '\r' < "$on/parameters.out" | awk '/E_level/ { L = ($6 == "8.6026") } L && / G  =/ { s += $12 } END { printf "%.12e", s / 1000 }')"
  EL="$(awk 'BEGIN { printf "%.12e", 8.6026 - 7.9936 }')"
  paste -d' ' <(tr -d '\r' < "$on/$OUT" | awk 'NF > 4 { print $1, $4 }') <(tr -d '\r' < "$off/$OUT" | awk 'NF > 4 { print $1, $4 }') |
    awk -v G="$G" -v EL="$EL" '
    BEGIN { u = 931.4940954; hbarc = 197.3269804; alpha = 1 / 137.035999084
            m3 = 3.0149322434; md = 2.0135532134; mp = 1.0072764675; m18 = 17.9947732059
            mpair = 1.00727647 + 17.99477097; m15 = 14.99626884
            eaa = 115 * m3 / (m18 + m3); bind = (mp + md - m3) * 931.49410242
            musf = md * mpair / (md + mpair) * u; musb = md * m15 / (md + m15) * u }
    { if ($1 != $3) { print "ENERGY", $1, $3; bad = 1; exit }
      esf = eaa - bind - $1; k = sqrt(2 * musf * esf) / hbarc
      zeta = alpha * (7 * musb - 9 * musf) / (hbarc * k)
      want = exp(2 * zeta * atan2(2 * (EL - $1), G))
      d = ($2 / $4) / want - 1; if (d < 0) d = -d; if (d > w) w = d
      if ($1 < lo || n == 0) lo = $1; if ($1 > hi) hi = $1; n++
      if (want < rmin || n == 1) rmin = want; if (want > rmax) rmax = want; zl = zeta }
    END { if (!bad) printf "%d %.3e %.4f %.4f %.4f %.4f %.5f\n", n, w, lo, hi, rmin, rmax, zl }' > "$WORK/$1.ratio"
  read -r n w lo hi rmin rmax zl < "$WORK/$1.ratio"
  if [ -n "${w:-}" ] && awk -v w="$w" -v n="$n" 'BEGIN { exit !(n > 30 && w < 1e-6) }'; then
    ok "$1: on/off = exp[2 zeta arctan(2(E_l-E)/G_l)] at $n points, E $lo-$hi MeV, G_l = $G MeV, zeta ~ $zl, ratio $rmin-$rmax: worst rel $w"
  else
    bad "$1: ratio check failed ($(cat "$WORK/$1.ratio"))"
  fi
}
echo "(c) isolated level: model ratio = |N_C|^2"
run c_off "experiment[A] segments=1,2 $KIN" "NOFOLD ONELEVEL"
run c_on "experiment[A] segments=1,2 $KIN lineshape=on" "NOFOLD ONELEVEL"
ran c_off && ran c_on && ratio_check one c_on c_off
run c2_off "experiment[A] segments=1,2 $KIN" "NOFOLD NOEXIT2"
run c2_on "experiment[A] segments=1,2 $KIN lineshape=on" "NOFOLD NOEXIT2"
ran c2_off && ran c2_on && ratio_check interfering c2_on c2_off
same "$WORK/c_off/output/$OUT" "$WORK/c2_off/output/$OUT" && bad "the second level changes nothing" \
  || ok "the second level does change the model (it interferes through the level matrix)"
# ... and with both levels decaying, the factor is not a common one.
run c3_off "experiment[A] segments=1,2 $KIN" "NOFOLD"
run c3_on "experiment[A] segments=1,2 $KIN lineshape=on" "NOFOLD"
if ran c3_off && ran c3_on; then
  keep=$fail
  ratio_check two c3_on c3_off > "$WORK/two.txt" 2>&1  # expected to fail
  fail=$keep
  grep -q FAIL "$WORK/two.txt" && ok "two decaying levels: the ratio is no longer the one-level factor" \
    || bad "two decaying levels: ratio equals the one-level factor"
fi

# (d) ----------------------------------------------------------------------
echo "(d) zeta < 0: the peak moves up"
GRID='/<segmentsData>/ { S = 1; print; next } /<\/segmentsData>/ { S = 0 } S && NF > 8 { S++; if (S == 2) { $4 = 0.3; $5 = 1.0; $15 = "data/grid.dat" } else $1 = 0 } { print }'
run d_off "experiment[A] segments=1 $KIN" "NOFOLD ONELEVEL GRID"
run d_on "experiment[A] segments=1 $KIN lineshape=on" "NOFOLD ONELEVEL GRID"
for r in d_off d_on; do
  awk 'BEGIN { for (e = 0.40; e <= 0.8200001; e += 0.0005) printf "%.4f 0.0 1.0 0.1\n", e }' > "$WORK/$r/data/grid.dat"
  (cd "$WORK/$r" && rm -rf output && mkdir output && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr 2>&1 | head -c 1000000 > log)
done
peak() { tr -d '\r' < "$1" | awk 'NF > 4 && $4 > m { m = $4; e = $1 } END { printf "%.4f", e }'; }
if [ -f "$WORK/d_off/output/$OUT" ] && [ -f "$WORK/d_on/output/$OUT" ]; then
  p0="$(peak "$WORK/d_off/output/$OUT")"; p1="$(peak "$WORK/d_on/output/$OUT")"
  awk -v a="$p0" -v b="$p1" 'BEGIN { exit !(b > a + 0.002) }' && ok "peak $p0 -> $p1 MeV (up)" || bad "peak $p0 -> $p1 MeV"
  read -r below above <<< "$(paste -d' ' <(tr -d '\r' < "$WORK/d_on/output/$OUT" | awk 'NF > 4 { print $1, $4 }') \
      <(tr -d '\r' < "$WORK/d_off/output/$OUT" | awk 'NF > 4 { print $4 }') |
    awk '{ r = $2 / $3 } $1 < 0.600 && r >= 1 { b++ } $1 > 0.618 && r <= 1 { a++ } END { print b + 0, a + 0 }')"
  [ "$below" = 0 ] && [ "$above" = 0 ] && ok "on/off < 1 below E_l = 0.609 MeV, > 1 above" \
    || bad "ratio on the wrong side: $below points below, $above above"
else
  bad "grid runs failed"
fi

# (e) ----------------------------------------------------------------------
echo "(e) refused"
refuse() {  # refuse NAME MESSAGE-FRAGMENT BLOCK [EDITS [FLAGS]]
  run "$1" "$3" "${4:-}" "${5:-}"
  if [ "$(cat "$WORK/$1/status")" != 0 ] && grep -q "ERROR: <thm> experiment\[" "$WORK/$1/log" &&
     grep -q -- "$2" "$WORK/$1/log"; then
    ok "$1: $(grep -m1 'ERROR' "$WORK/$1/log" | tr -d '\r')"
  else
    bad "$1: not refused with '$2' (status $(cat "$WORK/$1/status"))"; tail -3 "$WORK/$1/log" | sed 's/^/        /'
  fi
}
refuse no_kinematics "lineshape=on needs the kinematics" "experiment[A] segments=1,2 lineshape=on"
refuse bad_value "expected on or off" "experiment[A] segments=1,2 $KIN lineshape=1"
refuse no_energy "the spectator has no energy left" "experiment[A] segments=1,2 beam=18O target=3He spectator=d Ebeam=40 lineshape=on"
refuse formal "needs the Brune parameterization" "experiment[A] segments=1,2 $KIN lineshape=on" "" "--use-rmc"

echo
if [ "$fail" -eq 0 ]; then echo "PASS: THM Coulomb line shape"; else echo "FAIL"; exit 1; fi
