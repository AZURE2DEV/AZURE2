#!/usr/bin/env bash
#
# Spectator-direction window of a THM experiment: `spectatorAngles=` on an
# `experiment[<name>]` line (ThmDistortion::AngleNodes;
# docs/source/theory/thm_implementation.rst, "Experimental acceptance").  The
# distortion factor R(E) and the DW vertex are averaged over the accepted
# spectator directions, weight d cos(theta_cm) x acceptance: R = <|M|^2> /
# <|M_PW|^2> (dwpw), the DW vertex with the extra |phi~(q)|^2 of the data
# reduction.
#
# Model: tests/18O_p_a_thm, 2H(18O,a15N)n at 54 MeV (neutron spectator, the
# Trojan horse is the target: the quasi-free direction is theta_cm = 180 deg;
# the spectator is slower in the c.m. than the c.m. at the upper energies, so
# a lab angle has two c.m. branches there), and at 100 MeV (one branch).
#
#   (a) a window of zero width is its direction: cm:25-25 and the lab 12-12
#       (one branch) equal spectatorAngle=cm:25 / 12, for R and the DW
#       vertex (1e-8); a 1e-4 deg window around it too;
#   (b) a uniform window equals an independent average of single-direction
#       runs: R's <|M|^2>, <|M_PW|^2> on cm:120-180 against Simpson in
#       cos(theta_cm) over 21 spectatorAngle=cm:<deg> runs, and the DW model
#       against the same runs weighted by |phi~(q)|^2 = (kappa^2 + q^2)^-2
#       (Yukawa: a neutron spectator); quadrature accuracy;
#   (c) the lab window 0-180 (two branches where they exist) equals the c.m.
#       window 0-180, for R and the DW vertex;
#   (d) with the DW vertex a ps window only cuts q: ps=hulthen:0-1000 (no
#       cut) == no ps; a flat acceptance table == the uniform window; a
#       ps cut changes R;
#   (e) refusals.
#
#   ./tests/thm_spectator_angles/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
O18="$HERE/../18O_p_a_thm"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_spectator_angles.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
if command -v timeout >/dev/null 2>&1; then RUN="timeout ${TEST_TIMEOUT:-300}"; else RUN=""; fi
OUT="AZUREOut_aa=1_R=2.out"
KIN="beam=18O target=d spectator=n Ebeam=54"
fail=0

# run NAME BLOCK [EDITS] -- a "calculate with data" run of a copy of
# tests/18O_p_a_thm with BLOCK as its <thm> block, the .azr first passed
# through the awk programs named in EDITS.  Status in NAME/status.
run() {
  local d="$WORK/$1"
  rm -rf "$d"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$O18/data" "$d/"
  cp "$O18/18O_p_a_thm.azr" "$d/run.azr"
  local edit
  for edit in ${3:-}; do awk "${!edit}" "$d/run.azr" > "$d/tmp.azr" && mv "$d/tmp.azr" "$d/run.azr"; done
  [ -z "$2" ] || printf '\n<thm>\n%s\n</thm>\n' "$2" >> "$d/run.azr"
  (cd "$d" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr 2>&1 | head -c 1000000 > log;
   echo "${PIPESTATUS[1]}" > status)
}
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; fail=1; }
ran() { [ "$(cat "$WORK/$1/status")" = 0 ] && [ -f "$WORK/$1/output/$OUT" ] || { bad "run $1 failed"; tail -5 "$WORK/$1/log" | sed 's/^/        /'; return 1; }; }
same() { [ -f "$1" ] && [ -f "$2" ] && [ "$(tr -d '\r' < "$1")" = "$(tr -d '\r' < "$2")" ]; }
# Largest relative difference of the model (column 4) of two runs, point by point.
moddiff() {
  paste -d' ' <(tr -d '\r' < "$WORK/$1/output/$OUT" | awk 'NF > 4 { print $1, $4 }') \
              <(tr -d '\r' < "$WORK/$2/output/$OUT" | awk 'NF > 4 { print $1, $4 }') |
    awk '{ if ($1 != $3) { print "ENERGY"; exit 1 } d = ($4 - $2) / $2; if (d < 0) d = -d; if (d > m) m = d; n++ }
         END { printf "%.2e %d", m, n }'
}
# Largest relative difference of |M|^2, |M_PW|^2 and R of the distortion_point rows.
rowdiff() {
  paste -d' ' <(tr -d '\r' < "$WORK/$1/output/thm_experiments.out" | awk '$1 == "distortion_point" { print $6, $7, $8 }') \
              <(tr -d '\r' < "$WORK/$2/output/thm_experiments.out" | awk '$1 == "distortion_point" { print $6, $7, $8 }') |
    awk '{ for (i = 1; i <= 3; i++) { d = ($(i + 3) - $i) / $i; if (d < 0) d = -d; if (d > m) m = d } n++ }
         END { printf "%.2e %d", m, n }'
}
small() { awk -v v="${1% *}" -v n="${1#* }" -v t="$2" -v k="${3:-3}" 'BEGIN { exit !(n >= k && v < t) }'; }
NOFOLD='/<targetInt>/ { print; T = 1; next } /<\/targetInt>/ { T = 0 } T { next } { print }'
L="experiment[A] segments=1,2 $KIN"

# (a) ----------------------------------------------------------------------
echo "(a) a window of zero width is its direction"
run r_one "$L distortion=coulomb spectatorAngle=cm:25"
run r_win "$L distortion=coulomb spectatorAngles=cm:25-25"
run r_tiny "$L distortion=coulomb spectatorAngles=cm:24.99995-25.00005"
if ran r_one && ran r_win && ran r_tiny; then
  r="$(rowdiff r_one r_win)"; small "$r" 1e-8 && ok "R: cm:25-25 == spectatorAngle=cm:25 (|M|^2, |M_PW|^2, R: ${r% *})" || bad "R cm:25-25: $r"
  r="$(rowdiff r_one r_tiny)"; small "$r" 1e-8 && ok "R: cm:24.99995-25.00005 within ${r% *}" || bad "R tiny window: $r"
  r="$(moddiff r_one r_win)"; small "$r" 1e-8 && ok "R: model ${r% *} at ${r#* } points" || bad "R model: $r"
fi
run d_one "$L distortion=coulomb vertexModel=dw spectatorAngle=cm:25" "NOFOLD"
run d_win "$L distortion=coulomb vertexModel=dw spectatorAngles=cm:25-25" "NOFOLD"
if ran d_one && ran d_win; then
  r="$(moddiff d_one d_win)"; small "$r" 1e-8 && ok "DW: cm:25-25 == spectatorAngle=cm:25, model ${r% *} at ${r#* } points" || bad "DW cm:25-25: $r"
fi
K100="experiment[A] segments=1,2 beam=18O target=d spectator=n Ebeam=100"
run h_one "$K100 distortion=coulomb spectatorAngle=12"
run h_win "$K100 distortion=coulomb spectatorAngles=12-12"
run hd_one "$K100 distortion=coulomb vertexModel=dw spectatorAngle=12" "NOFOLD"
run hd_win "$K100 distortion=coulomb vertexModel=dw spectatorAngles=12-12" "NOFOLD"
if ran h_one && ran h_win && ran hd_one && ran hd_win; then
  r="$(rowdiff h_one h_win)"; small "$r" 1e-8 && ok "R: lab 12-12 == spectatorAngle=12 at 100 MeV (one branch): ${r% *}" || bad "R lab 12-12: $r"
  r="$(moddiff hd_one hd_win)"; small "$r" 1e-8 && ok "DW: lab 12-12 == spectatorAngle=12: model ${r% *}" || bad "DW lab 12-12: $r"
fi

# (b) ----------------------------------------------------------------------
echo "(b) a uniform window == an independent average of single directions"
run r_band "$L distortion=coulomb spectatorAngles=cm:120-180 spectatorAngleNodes=16"
run d_band "$L distortion=coulomb vertexModel=dw spectatorAngles=cm:120-180 spectatorAngleNodes=16" "NOFOLD"
# Simpson in c = cos(theta_cm) on [-1, -0.5], 20 intervals.
NS=20
for j in $(seq 0 $NS); do
  deg="$(awk -v j="$j" -v n="$NS" 'BEGIN { c = -1 + 0.5 * j / n; printf "%.12f", atan2(sqrt(1 - c * c), c) * 180 / atan2(0, -1) }')"
  run r_s$j "$L distortion=coulomb spectatorAngle=cm:$deg"
  run d_s$j "$L distortion=coulomb vertexModel=dw spectatorAngle=cm:$deg" "NOFOLD"
done
if ran r_band && ran d_band; then
  good=1
  for j in $(seq 0 $NS); do ran r_s$j && ran d_s$j || good=0; done
  if [ "$good" = 1 ]; then
    # R: <|M|^2> and <|M_PW|^2> at the three distortion_point energies.
    for j in $(seq 0 $NS); do
      tr -d '\r' < "$WORK/r_s$j/output/thm_experiments.out" | awk -v j="$j" '$1 == "distortion_point" { print j, ++k, $6, $7 }'
    done > "$WORK/r_simpson.txt"
    tr -d '\r' < "$WORK/r_band/output/thm_experiments.out" | awk '$1 == "distortion_point" { print ++k, $6, $7 }' > "$WORK/r_band.txt"
    # Simpson on 20 and on 10 intervals, Richardson: (16 S20 - S10)/15.
    r="$(awk -v n="$NS" 'NR == FNR { c = ($1 == 0 || $1 == n) ? 1 : ($1 % 2 ? 4 : 2); m[$2] += c * $3; p[$2] += c * $4; s[$2] += c
           if ($1 % 2 == 0) { i = $1 / 2; c = (i == 0 || i == n / 2) ? 1 : (i % 2 ? 4 : 2); m2[$2] += c * $3; p2[$2] += c * $4; s2[$2] += c }
           next }
         { a = (16 * m[$1] / s[$1] - m2[$1] / s2[$1]) / 15; b = (16 * p[$1] / s[$1] - p2[$1] / s2[$1]) / 15
           d = ($2 - a) / a; if (d < 0) d = -d; if (d > w) w = d
           d = ($3 - b) / b; if (d < 0) d = -d; if (d > w) w = d; k++ } END { printf "%.2e %d", w, k }' \
         "$WORK/r_simpson.txt" "$WORK/r_band.txt")"
    small "$r" 2e-5 && ok "R: <|M|^2>, <|M_PW|^2> on cm:120-180 == Simpson-Richardson over 21 directions to ${r% *} (${r#* } energies)" \
      || bad "R window against Simpson: $r"
    # DW: the model at the lowest and highest point, weight (kappa^2 + q^2)^-2,
    # q of each direction from its dw_vertex_point rows (MeV/c).
    kappa="$(tr -d '\r' < "$WORK/d_band/output/thm_experiments.out" | awk '/^# k_aA/ { for (i = 1; i <= NF; i++) if ($i == "kappa") { print $(i + 2); exit } }')"
    for j in $(seq 0 $NS); do
      qs="$(tr -d '\r' < "$WORK/d_s$j/output/thm_experiments.out" | awk '$1 == "dw_vertex_point" && $5 == 0 { printf "%s %s ", $2, $3 }')"
      tr -d '\r' < "$WORK/d_s$j/output/$OUT" | awk -v j="$j" -v qs="$qs" 'NF > 4 { print j, $1, $4, qs }'
    done > "$WORK/d_simpson.txt"
    tr -d '\r' < "$WORK/d_band/output/$OUT" | awk 'NF > 4 { print $1, $4 }' > "$WORK/d_band.txt"
    r="$(awk -v n="$NS" -v kap="$kappa" -v hc=197.3269804 '
         NR == FNR { split($0, f, " "); # j E model E1 q1 E2 q2 E3 q3
           q = -1; for (i = 4; i <= 8; i += 2) { d = f[i] - $2; if (d < 0) d = -d; if (d < 1e-9) q = f[i + 1] / hc }
           if (q < 0) next
           c = ($1 == 0 || $1 == n) ? 1 : ($1 % 2 ? 4 : 2); w = c / ((kap * kap + q * q) ^ 2)
           num[$2] += w * $3; den[$2] += w
           if ($1 % 2 == 0) { i = $1 / 2; c = (i == 0 || i == n / 2) ? 1 : (i % 2 ? 4 : 2); w = c / ((kap * kap + q * q) ^ 2)
             num2[$2] += w * $3; den2[$2] += w }
           next }
         { for (e in num) { d = e - $1; if (d < 0) d = -d; if (d < 1e-9) break } if (d >= 1e-9) next
           a = (16 * num[e] / den[e] - num2[e] / den2[e]) / 15; d = ($2 - a) / a; if (d < 0) d = -d; if (d > m) m = d; k++ }
         END { printf "%.2e %d", m, k }' "$WORK/d_simpson.txt" "$WORK/d_band.txt")"
    small "$r" 2e-5 2 && ok "DW: model on cm:120-180 == |phi~|^2-weighted Simpson-Richardson over 21 directions to ${r% *} (${r#* } points, kappa = $kappa fm^-1)" \
      || bad "DW window against Simpson: $r"
  fi
fi

# (c) ----------------------------------------------------------------------
echo "(c) the lab window 0-180 (two branches) == the c.m. window 0-180"
run r_lab "$L distortion=coulomb spectatorAngles=0-180 spectatorAngleNodes=32"
run r_cm "$L distortion=coulomb spectatorAngles=cm:0-180 spectatorAngleNodes=32"
run d_lab "$L distortion=coulomb vertexModel=dw spectatorAngles=0-180 spectatorAngleNodes=32" "NOFOLD"
run d_cm "$L distortion=coulomb vertexModel=dw spectatorAngles=cm:0-180 spectatorAngleNodes=32" "NOFOLD"
if ran r_lab && ran r_cm && ran d_lab && ran d_cm; then
  grep -q "per branch" "$WORK/d_lab/log" && ok "the lab window has two branches here" || bad "no second branch"
  r="$(rowdiff r_lab r_cm)"; small "$r" 1e-7 && ok "R: lab 0-180 == cm:0-180 to ${r% *}" || bad "R lab/cm full: $r"
  r="$(moddiff d_lab d_cm)"; small "$r" 1e-6 && ok "DW: lab 0-180 == cm:0-180, model to ${r% *}" || bad "DW lab/cm full: $r"
fi

# (d) ----------------------------------------------------------------------
echo "(d) the ps window cuts q; acceptance tables"
run d_nocut "$L distortion=coulomb vertexModel=dw spectatorAngles=cm:120-180 ps=hulthen:0-1000" "NOFOLD"
printf '# flat acceptance\n120 1\n150 1\n180 1\n' > "$WORK/flat.dat"
run d_table "$L distortion=coulomb vertexModel=dw spectatorAngles=cm:table:$WORK/flat.dat" "NOFOLD"
run d_plain "$L distortion=coulomb vertexModel=dw spectatorAngles=cm:120-180" "NOFOLD"
run r_cut "$L distortion=coulomb spectatorAngles=cm:120-180 ps=hulthen:0-20"
run r_plain "$L distortion=coulomb spectatorAngles=cm:120-180"
if ran d_nocut && ran d_table && ran d_plain && ran r_cut && ran r_plain; then
  same "$WORK/d_nocut/output/$OUT" "$WORK/d_plain/output/$OUT" && ok "DW: a ps window covering every q changes nothing" \
    || bad "DW: ps=hulthen:0-1000 differs from no ps"
  same "$WORK/d_table/output/$OUT" "$WORK/d_plain/output/$OUT" && ok "DW: a flat acceptance table == the uniform window" \
    || bad "DW: flat table differs"
  r="$(rowdiff r_cut r_plain)"
  awk -v v="${r% *}" 'BEGIN { exit !(v > 1e-6) }' && ok "R: a ps cut at 20 MeV/c changes <|M|^2> by up to ${r% *}" || bad "R ps cut: $r"
  grep -q "|p_s| cut by the ps window" "$WORK/r_cut/output/thm_experiments.out" && ok "the cut is in the distortion line" \
    || bad "no cut in the distortion line"
fi

# (e) ----------------------------------------------------------------------
echo "(e) refused"
refuse() {  # refuse NAME MESSAGE-FRAGMENT BLOCK
  run "$1" "$3"
  if [ "$(cat "$WORK/$1/status")" != 0 ] && grep -q "ERROR: <thm>" "$WORK/$1/log" &&
     grep -q -- "$2" "$WORK/$1/log"; then
    ok "$1: $(grep -m1 'ERROR' "$WORK/$1/log" | tr -d '\r' | cut -c1-150)"
  else
    bad "$1: not refused with '$2' (status $(cat "$WORK/$1/status"))"; tail -3 "$WORK/$1/log" | sed 's/^/        /'
  fi
}
printf '0.4 1.0\n1.0 2.0\n' > "$WORK/w.dat"
printf '10 1\n' > "$WORK/short.dat"
refuse no_distortion "with neither it has nothing to average" "$L spectatorAngles=cm:120-180"
refuse table_distortion "with neither it has nothing to average" "$L distortion=table:$WORK/w.dat spectatorAngles=cm:120-180"
refuse both "exclude each other" "$L distortion=coulomb spectatorAngle=20 spectatorAngles=cm:120-180"
refuse nodes_alone "needs a spectator-direction window" "$L distortion=coulomb spectatorAngleNodes=4"
refuse ps_nodes "the nodes are spectatorAngleNodes=" "$L distortion=coulomb vertexModel=dw ps=hulthen:0-40 psNodes=8 spectatorAngles=cm:120-180"
refuse order "expected thmin-thmax" "$L distortion=coulomb spectatorAngles=cm:50-20"
refuse range "expected thmin-thmax" "$L distortion=coulomb spectatorAngles=10-200"
refuse nodes "1 to 64" "$L distortion=coulomb spectatorAngles=cm:120-180 spectatorAngleNodes=65"
refuse no_file "cannot read the angle table" "$L distortion=coulomb spectatorAngles=cm:table:$WORK/none.dat"
refuse short_table "at least two rows" "$L distortion=coulomb spectatorAngles=table:$WORK/short.dat"
refuse reach_r "no spectator direction of spectatorAngles=cm:0-10 is accepted" "$L distortion=coulomb spectatorAngles=cm:0-10 ps=hulthen:0-40"
refuse reach_dw "no spectator direction of spectatorAngles=cm:0-10 is accepted" "$L distortion=coulomb vertexModel=dw spectatorAngles=cm:0-10 ps=hulthen:0-40"

echo
if [ "$fail" -eq 0 ]; then echo "PASS: THM spectator-direction window"; else echo "FAIL"; exit 1; fi
