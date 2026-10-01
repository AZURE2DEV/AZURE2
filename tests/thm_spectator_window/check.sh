#!/usr/bin/env bash
#
# Spectator-momentum window of a THM experiment: `ps=` on an
# `experiment[<name>]` line of the <thm> block (ThmLineshape.h
# ThmSpectatorWindow; docs/source/theory/thm_implementation.rst,
# "Spectator-momentum window").  At fixed E the spectator direction fixes
# q = |p_s| = |k_sF - beta k_aA|, and the three-body phase space is
# d cos(theta_cm) = q dq/(beta k_sF k_aA), so the HOES model at E is
#   sigma(E) = Int |phi(q)|^2 sigma(E; q) q dq / Int |phi(q)|^2 q dq
# over the q in the window that the kinematics reach at E, where node q adds
# T_s = q^2/2mu_sx to E + B in the entrance vertex (as the scalar
# spectatorEnergy does), by Gauss-Legendre in cos(theta_cm).  (Until October
# 2026 the weight was |phi|^2 p^2 dp, the measure of events integrated over E
# as well; the pins of (c) and (d) changed with it.)
#
# Model: tests/18O_p_a_thm (19F: p+18O -> a+15N, l = 0 entrance, B = 2.2246
# MeV, two THM segments with free norms), whose reaction 2H(18O,a15N)n has a
# deuteron Trojan horse: x = p, s = n, mu_sx = m_p m_n/(m_p + m_n); at 54 MeV
# the reachable q runs from |k_sF - beta k_aA| (<= 8 MeV/c over the data) to
# beyond 130 MeV/c.
#
#   (a) ps=delta changes nothing: byte-identical to the same line without the
#       key (with and without kinematics, and with spectatorEnergy set);
#   (b) a window shrinking to a point reproduces spectatorEnergy = p^2/2mu_sx:
#       pmin = pmax (one node) and a 2e-4 MeV/c wide window, to 1e-8;
#   (c) a flat table (|phi|^2 constant on [20, 40] MeV/c) equals the average
#       of 41 single-point spectatorEnergy runs with the weight q (Simpson's
#       rule), to 1e-6; and the 16-node Hulthen window [0, 40] agrees with 32
#       nodes; the node tables (ps_table in thm_experiments.out) grow
#       linearly with the nodes and stay below 4 kB per point at 32 nodes;
#   (d) near a node of the vertex M_0 (a point p_s = 24 MeV/c puts it at
#       E ~ 0.70 MeV, no folding) the window fills it;
#   (e) refusals;
#   (f) one acceptance: ps= alone is spectatorAngles=cm:0-180 with the same
#       cut (byte-identical with as many nodes), and a c.m. angle window
#       narrows it.
#
#   ./tests/thm_spectator_window/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../18O_p_a_thm"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_spectator_window.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
if command -v timeout >/dev/null 2>&1; then RUN="timeout ${TEST_TIMEOUT:-300}"; else RUN=""; fi
OUT="AZUREOut_aa=1_R=2.out"
KIN="beam=18O target=d spectator=n Ebeam=54"
# mu_sx (MeV) from the built-in nuclear masses of p and n, as the engine.
MUSX="$(awk 'BEGIN { mp = 1.0072764675; mn = 1.0086649159; printf "%.15e", mp * mn / (mp + mn) * 931.49410242 }')"
fail=0

# run NAME BLOCK [EDITS] -- a "calculate with data" run of a copy of the
# project with BLOCK as its <thm> block ("" = none), the .azr first passed
# through the awk programs named in EDITS.  Status in NAME/status.
run() {
  local d="$WORK/$1"
  rm -rf "$d"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$SRC/data" "$d/"
  cp "$SRC/18O_p_a_thm.azr" "$d/run.azr"
  local edit
  for edit in ${3:-}; do awk "${!edit}" "$d/run.azr" > "$d/tmp.azr" && mv "$d/tmp.azr" "$d/run.azr"; done
  [ -z "$2" ] || printf '\n<thm>\n%s\n</thm>\n' "$2" >> "$d/run.azr"
  [ -z "${TABLE:-}" ] || printf '%b' "$TABLE" > "$d/ps.dat"
  (cd "$d" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr 2>&1 | head -c 1000000 > log;
   echo "${PIPESTATUS[1]}" > status)
}
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; fail=1; }
same() { [ -f "$1" ] && [ -f "$2" ] && [ "$(tr -d '\r' < "$1")" = "$(tr -d '\r' < "$2")" ]; }
ran() { [ "$(cat "$WORK/$1/status")" = 0 ] && [ -f "$WORK/$1/output/$OUT" ] || { bad "run $1 failed"; tail -5 "$WORK/$1/log" | sed 's/^/        /'; return 1; }; }
model() { tr -d '\r' < "$WORK/$1/output/$OUT" | awk 'NF > 4 { print $1, $4 }'; }
# worst relative difference of the model columns of two runs (energies must agree)
worst() {
  paste -d' ' <(model "$1") <(model "$2") |
    awk '{ if ($1 != $3) { print "ENERGY"; exit } d = $2 / $4 - 1; if (d < 0) d = -d; if (d > w) w = d; n++ }
         END { printf "%d %.3e\n", n, w }'
}
T_of() { awk -v p="$1" -v mu="$MUSX" 'BEGIN { printf "%.17g", p * p / (2 * mu) }'; }

NOFOLD='/<targetInt>/ { print; T = 1; next } /<\/targetInt>/ { T = 0 } T { next } { print }'

# (a) ----------------------------------------------------------------------
echo "(a) ps=delta == no key"
run nokey "experiment[A] segments=1,2 $KIN"
run delta "experiment[A] segments=1,2 $KIN ps=delta"
ran nokey && ran delta
for f in chiSquared.out normalizations.out parameters.out thm_experiments.out "$OUT"; do
  same "$WORK/nokey/output/$f" "$WORK/delta/output/$f" && ok "$f byte-identical" || bad "$f differs"
done
run nokin "experiment[A] segments=1,2"
run deltanokin "experiment[A] segments=1,2 ps=delta"
ran nokin && ran deltanokin
same "$WORK/nokin/output/$OUT" "$WORK/deltanokin/output/$OUT" && ok "no kinematics: ps=delta byte-identical" \
  || bad "no kinematics: ps=delta differs"
run se "spectatorEnergy=0.4
experiment[A] segments=1,2 $KIN"
run sedelta "spectatorEnergy=0.4
experiment[A] segments=1,2 $KIN ps=delta"
ran se && ran sedelta
same "$WORK/se/output/$OUT" "$WORK/sedelta/output/$OUT" && ok "spectatorEnergy=0.4: ps=delta byte-identical" \
  || bad "spectatorEnergy=0.4: ps=delta differs"
same "$WORK/nokey/output/$OUT" "$WORK/se/output/$OUT" && bad "spectatorEnergy=0.4 changes nothing" || true

# (b) ----------------------------------------------------------------------
echo "(b) a window shrinking to a point == spectatorEnergy = p^2/2mu_sx (mu_sx = $MUSX MeV)"
T30="$(T_of 30)"
run pt_se "spectatorEnergy=$T30
experiment[A] segments=1,2 $KIN"
run pt_one "experiment[A] segments=1,2 $KIN ps=hulthen:30-30"
run pt_narrow "experiment[A] segments=1,2 $KIN ps=hulthen:29.9999-30.0001"
run pt_gauss "experiment[A] segments=1,2 $KIN ps=gauss:60:29.9999-30.0001"
if ran pt_se && ran pt_one && ran pt_narrow && ran pt_gauss; then
  for r in pt_one pt_narrow pt_gauss; do
    read -r n w <<< "$(worst $r pt_se)"
    awk -v w="$w" -v n="$n" 'BEGIN { exit !(n > 30 && w < 1e-8) }' \
      && ok "$r vs spectatorEnergy=$T30: $n points, worst rel $w" || bad "$r vs spectatorEnergy: $n points, worst rel $w"
  done
  grep -q "one node" "$WORK/pt_one/log" && ok "pmin = pmax is one node" || bad "pmin = pmax: no 'one node' in the log"
fi

# (c) ----------------------------------------------------------------------
echo "(c) flat |phi|^2 on [20, 40] MeV/c == Simpson average, weight q, of 41 spectatorEnergy runs"
TABLE='# p_s (MeV/c)  w\n20 1\n40 1\n' run flat "experiment[A] segments=1,2 $KIN ps=table:ps.dat"
ran flat
: > "$WORK/simpson.txt"
for i in $(seq 0 40); do
  p="$(awk -v i=$i 'BEGIN { printf "%.1f", 20 + i / 2 }')"
  run "se$i" "spectatorEnergy=$(T_of $p)
experiment[A] segments=1,2 $KIN"
  ran "se$i" || continue
  c=$(( i == 0 || i == 40 ? 1 : (i % 2 == 1 ? 4 : 2) ))
  model "se$i" | awk -v c=$c -v p="$p" '{ printf "%d %s %.15e\n", NR, $1, c * p * $2 }' >> "$WORK/simpson.txt"
done
# Int q sigma dq / Int q dq = (h/3) sum c_i q_i sigma_i / 600, h = 0.5 MeV/c.
awk '{ s[$1] += $3; e[$1] = $2; n = ($1 > n ? $1 : n) } END { for (i = 1; i <= n; i++) printf "%s %.12e\n", e[i], s[i] / 3600 }' \
  "$WORK/simpson.txt" > "$WORK/simpson.avg"
read -r n w <<< "$(paste -d' ' <(model flat) "$WORK/simpson.avg" |
  awk '{ if ($1 != $3) { print "0 ENERGY"; exit } d = $2 / $4 - 1; if (d < 0) d = -d; if (d > w) w = d; n++ } END { printf "%d %.3e\n", n, w }')"
awk -v w="$w" -v n="$n" 'BEGIN { exit !(n > 30 && w < 1e-6) }' \
  && ok "table (16 Gauss-Legendre nodes) vs Simpson (41 runs, h = 0.5 MeV/c): $n points, worst rel $w" \
  || bad "table vs Simpson: $n points, worst rel $w"
same "$WORK/flat/output/$OUT" "$WORK/se20/output/$OUT" && bad "the flat window equals its midpoint" || true
read -r n w <<< "$(worst flat se20)"
ok "(size) flat window vs its midpoint p_s = 30 MeV/c: worst rel $w"
run h16 "experiment[A] segments=1,2 $KIN ps=hulthen:0-40"
run h8 "experiment[A] segments=1,2 $KIN ps=hulthen:0-40 psNodes=8"
run h32 "experiment[A] segments=1,2 $KIN ps=hulthen:0-40 psNodes=32"
if ran h16 && ran h8 && ran h32; then
  read -r n w16 <<< "$(worst h16 h32)"
  read -r n w8 <<< "$(worst h8 h32)"
  awk -v w="$w16" 'BEGIN { exit !(w < 1e-6) }' && ok "hulthen [0, 40]: 16 vs 32 nodes worst rel $w16 (8 vs 32: $w8)" \
    || bad "hulthen [0, 40]: 16 vs 32 nodes worst rel $w16"
  read -r n w <<< "$(worst h16 nokey)"
  ok "(size) hulthen [0, 40] vs the quasi-free point p_s = 0: worst rel $w"
  grep -q "^ps: hulthen a=0.2317 b=1.202 fm^-1, p_s in \[0, 40\] MeV/c, 16 Gauss-Legendre nodes in cos theta_cm" \
    "$WORK/h16/output/thm_experiments.out" && ok "thm_experiments.out lists the window" || bad "no ps: line in thm_experiments.out"
  # Memory of the node tables (EPoint::ThmPsTable): entrance channels only,
  # flat, so linear in the nodes and a few hundred bytes per point.  The
  # per-node [J group][channel] vectors of every channel they replace held
  # ~4.5 kB per point and node on the 19F THM model (2.2 GB with 16 nodes).
  pst() { awk '$1 == "ps_table" { printf "%.0f %.0f", $2, $3 }' "$WORK/$1/output/thm_experiments.out"; }
  read -r np8 b8 <<< "$(pst h8)"; read -r np16 b16 <<< "$(pst h16)"; read -r np32 b32 <<< "$(pst h32)"
  if [ -n "${b32:-}" ] && [ "$np8" = "$np16" ] && [ "$np16" = "$np32" ] && [ "$np32" -gt 0 ]; then
    awk -v n="$np32" -v a="$b8" -v b="$b16" -v c="$b32" 'BEGIN { exit !(c / n < 4096 && (c - b) * 8 == (b - a) * 16 && b > a) }' \
      && ok "ps_table: $np32 points (with sub-points), $b8 / $b16 / $b32 bytes at 8 / 16 / 32 nodes ($(awk -v n="$np32" -v c="$b32" 'BEGIN { printf "%.0f", c / n }') B per point at 32)" \
      || bad "ps_table not bounded/linear: $np32 points, $b8 / $b16 / $b32 bytes at 8 / 16 / 32 nodes"
  else
    bad "no ps_table line in thm_experiments.out (h8/h16/h32: '$np8 $b8' '$np16 $b16' '$np32 ${b32:-}')"
  fi
fi

# (d) ----------------------------------------------------------------------
echo "(d) a vertex node: p_s = 24 MeV/c vs the window [0, 40] (no folding)"
run node "experiment[A] segments=1,2 $KIN ps=hulthen:24-24" NOFOLD
run fill "experiment[A] segments=1,2 $KIN ps=hulthen:0-40" NOFOLD
if ran node && ran fill; then
  read -r en mn mx wn wx <<< "$(paste -d' ' <(model node) <(model fill) |
    awk 'NR == 1 || $2 < mn { mn = $2; en = $1; wn = $4 } $2 > mx { mx = $2 } $4 > wx { wx = $4 } END { print en, mn, mx, wn, wx }')"
  awk -v mn="$mn" -v mx="$mx" -v wn="$wn" -v wx="$wx" 'BEGIN { exit !(mn / mx < 1e-3 && wn / wx > 0.05) }' \
    && ok "at E = $en MeV: point p_s = 24: model/max = $(awk -v a="$mn" -v b="$mx" 'BEGIN { printf "%.2e", a / b }'), window: $(awk -v a="$wn" -v b="$wx" 'BEGIN { printf "%.3f", a / b }') (filled by a factor $(awk -v a="$wn" -v b="$mn" 'BEGIN { printf "%.2e", a / b }'))" \
    || bad "node not filled: E = $en, point $mn/$mx, window $wn/$wx"
fi

# (e) ----------------------------------------------------------------------
echo "(e) refused"
refuse() {  # refuse NAME MESSAGE-FRAGMENT BLOCK
  run "$1" "$3"
  if [ "$(cat "$WORK/$1/status")" != 0 ] && grep -q "ERROR: <thm> experiment\[" "$WORK/$1/log" &&
     grep -q -- "$2" "$WORK/$1/log"; then
    ok "$1: $(grep -m1 'ERROR' "$WORK/$1/log" | tr -d '\r')"
  else
    bad "$1: not refused with '$2' (status $(cat "$WORK/$1/status"))"; tail -3 "$WORK/$1/log" | sed 's/^/        /'
  fi
}
refuse no_kinematics "needs the kinematics of the reaction" "experiment[A] segments=1,2 ps=hulthen:0-40"
refuse reversed "expected delta, hulthen" "experiment[A] segments=1,2 $KIN ps=hulthen:40-20"
refuse negative "expected delta, hulthen" "experiment[A] segments=1,2 $KIN ps=gauss:50:-10-40"
refuse unknown_kind "expected delta, hulthen" "experiment[A] segments=1,2 $KIN ps=lorentz:0-40"
refuse bad_ab "Hulthen a,b" "experiment[A] segments=1,2 $KIN ps=hulthen:1.2,0.2:0-40"
refuse bad_fwhm "FWHM" "experiment[A] segments=1,2 $KIN ps=gauss:0:0-40"
refuse no_table "cannot read the ps table" "experiment[A] segments=1,2 $KIN ps=table:missing.dat"
TABLE='20 1\n10 1\n' refuse table_order "strictly increasing" "experiment[A] segments=1,2 $KIN ps=table:ps.dat"
TABLE='20 0\n40 0\n' refuse table_zero "zero in every row" "experiment[A] segments=1,2 $KIN ps=table:ps.dat"
TABLE='20 1\n' refuse table_short "at least two rows" "experiment[A] segments=1,2 $KIN ps=table:ps.dat"
refuse nodes_range "1 to 64" "experiment[A] segments=1,2 $KIN ps=hulthen:0-40 psNodes=65"
refuse nodes_alone "psNodes= needs a ps window" "experiment[A] segments=1,2 $KIN psNodes=8"
refuse twice "given twice" "experiment[A] segments=1,2 $KIN ps=hulthen:0-40
experiment[A] ps=delta"
refuse with_se "spectatorEnergy both set" "spectatorEnergy=0.4
experiment[A] segments=1,2 $KIN ps=hulthen:0-40"
refuse with_se_pair "spectatorEnergy both set" "spectatorEnergy[1]=0.4
experiment[A] segments=1,2 $KIN ps=gauss:60:0-40"
refuse theta "expected all or thmin-thmax" "experiment[A] segments=1,2 theta=10-0"
refuse out_of_reach "do not overlap the ps window" "experiment[A] segments=1,2 $KIN ps=hulthen:200-300"
refuse with_angle "spectatorAngle= (one direction) and a ps window" \
  "experiment[A] segments=1,2 $KIN ps=hulthen:0-40 distortion=coulomb spectatorAngle=10"
refuse angle_nodes "the nodes are spectatorAngleNodes" \
  "experiment[A] segments=1,2 $KIN ps=hulthen:0-40 psNodes=8 spectatorAngles=cm:120-180"
refuse angles_alone "it has nothing to average" "experiment[A] segments=1,2 $KIN spectatorAngles=cm:120-180"

# (f) ----------------------------------------------------------------------
echo "(f) one acceptance: ps= alone == spectatorAngles=cm:0-180 with the cut"
run all16 "experiment[A] segments=1,2 $KIN ps=hulthen:0-40 spectatorAngles=cm:0-180 spectatorAngleNodes=16"
run cm160 "experiment[A] segments=1,2 $KIN ps=hulthen:0-40 spectatorAngles=cm:160-180 spectatorAngleNodes=16"
if ran h16 && ran all16 && ran cm160; then
  same "$WORK/h16/output/$OUT" "$WORK/all16/output/$OUT" && ok "ps=hulthen:0-40 == + spectatorAngles=cm:0-180 (16 nodes): byte-identical" \
    || { read -r n w <<< "$(worst h16 all16)"; bad "ps alone vs cm:0-180: $n points, worst rel $w"; }
  read -r n w <<< "$(worst cm160 h16)"
  awk -v w="$w" 'BEGIN { exit !(w > 1e-3) }' && ok "(size) cm:160-180 with the cut (q <= 40 MeV/c needs theta_cm >= ~147 deg) vs the cut alone: worst rel $w" \
    || bad "cm:160-180 changes nothing: worst rel $w"
fi

echo
if [ "$fail" -eq 0 ]; then echo "PASS: THM spectator-momentum window"; else echo "FAIL"; exit 1; fi
