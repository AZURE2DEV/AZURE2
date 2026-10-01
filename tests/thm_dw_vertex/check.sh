#!/usr/bin/env bash
#
# Distorted-wave entrance vertex of a THM experiment: `vertexModel=dw` on an
# `experiment[<name>]` line of the <thm> block (ThmDwVertex.h;
# docs/source/theory/thm_implementation.rst, "Distorted-wave entrance
# vertex").  The vertex M_l of every segment of the experiment becomes the
# surface term of the prior-form DWBA built from the experiment's distorted
# waves (distortion=coulomb|optical); R(E) is then not applied.
#
# Model: tests/18O_p_a_thm, whose real reaction 2H(18O,a15N)n at 54 MeV has a
# neutron spectator (entrance p + 18O, l = 0); no folding where the model is
# compared point by point.
#
#   (a) vertexModel=pw is byte-identical to the line without the key, also
#       with distortion=coulomb (R(E) applied as before);
#   (b) plane waves in both channels: the Gram entries of thm_experiments.out
#       are those of the plane-wave vertex at rho = p a, G11 = j_0(rho)^2,
#       G22 = (rho j_0'(rho))^2, G12 = j_0 rho j_0' (1e-8 of the largest), and
#       the model is the plane-wave one up to the kinematic p(E) (within 5 %,
#       the most near the node of M_0 at 0.7 MeV);
#   (c) point Coulomb in d + 18O: R(E) is not applied (no distortion rows),
#       the vertex rows are there, and the model moves; with a ps window the
#       nodes are placed on its reachable part (runs, differs from delta);
#   (d) refusals.
#   The vertex itself is checked against the plane-wave limit (1e-9) and an
#   independent quadrature in tests/reference (ctest thm_dw_vertex), and in
#   the session by tests/pyazr/thm_dw_vertex_test.py.
#
#   ./tests/thm_dw_vertex/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
O18="$HERE/../18O_p_a_thm"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_dw_vertex.XXXXXX")"
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
same() { [ -f "$1" ] && [ -f "$2" ] && [ "$(tr -d '\r' < "$1")" = "$(tr -d '\r' < "$2")" ]; }
ran() { [ "$(cat "$WORK/$1/status")" = 0 ] && [ -f "$WORK/$1/output/$OUT" ] || { bad "run $1 failed"; tail -5 "$WORK/$1/log" | sed 's/^/        /'; return 1; }; }
# "E ratio" per point: model(ON)/model(OFF), energies must agree.
ratios() {
  paste -d' ' <(tr -d '\r' < "$WORK/$1/output/$OUT" | awk 'NF > 4 { print $1, $4 }') \
              <(tr -d '\r' < "$WORK/$2/output/$OUT" | awk 'NF > 4 { print $1, $4 }') |
    awk '{ if ($1 != $3) { print "ENERGY", $1, $3; exit 1 } printf "%.12e %.12e\n", $1, $2 / $4 }'
}
NOFOLD='/<targetInt>/ { print; T = 1; next } /<\/targetInt>/ { T = 0 } T { next } { print }'

# (a) ----------------------------------------------------------------------
echo "(a) vertexModel=pw == no key"
run nokey "experiment[A] segments=1,2 $KIN"
run pw "experiment[A] segments=1,2 $KIN vertexModel=pw"
run r_nokey "experiment[A] segments=1,2 $KIN distortion=coulomb"
run r_pw "experiment[A] segments=1,2 $KIN distortion=coulomb vertexModel=pw"
if ran nokey && ran pw && ran r_nokey && ran r_pw; then
  for f in chiSquared.out normalizations.out parameters.out thm_experiments.out "$OUT"; do
    same "$WORK/nokey/output/$f" "$WORK/pw/output/$f" && ok "$f byte-identical" || bad "$f differs"
    same "$WORK/r_nokey/output/$f" "$WORK/r_pw/output/$f" && ok "$f byte-identical with distortion=coulomb" \
      || bad "$f differs with distortion=coulomb"
  done
fi

# (b) ----------------------------------------------------------------------
echo "(b) plane waves in both channels: the plane-wave vertex"
run b_pw "experiment[A] segments=1,2 $KIN" "NOFOLD"
run b_dw "experiment[A] segments=1,2 $KIN distortion=optical opticalAA=plane opticalSF=plane vertexModel=dw" "NOFOLD"
if ran b_pw && ran b_dw; then
  # dw_vertex_point: E q pa l G11 G22 ReG12 ImG12, l = 0 and 1 (j_1 = sin/x^2 - cos/x).
  r="$(tr -d '\r' < "$WORK/b_dw/output/thm_experiments.out" | awk '$1 == "dw_vertex_point" {
      x = $4; l = $5
      if (l == 0) { j = sin(x) / x; d = cos(x) - sin(x) / x }
      else if (l == 1) { j = sin(x) / (x * x) - cos(x) / x; d = x * (j * (-2) / x + sin(x) / x) }
      else next
      g11 = j * j; g22 = d * d; g12 = j * d; s = g11 > g22 ? g11 : g22
      e = $6 - g11; if (e < 0) e = -e; if (e / s > w) w = e / s
      e = $7 - g22; if (e < 0) e = -e; if (e / s > w) w = e / s
      e = $8 - g12; if (e < 0) e = -e; if (e / s > w) w = e / s
      e = $9; if (e < 0) e = -e; if (e / s > w) w = e / s
      n++ } END { printf "%d %.2e", n, w }')"
  awk -v n="${r% *}" -v w="${r#* }" 'BEGIN { exit !(n >= 3 && w < 1e-8) }' \
    && ok "Gram entries == plane-wave vertex at rho = p a: ${r% *} rows, worst ${r#* } of the largest" \
    || bad "plane-wave Gram entries: $r"
  r="$(ratios b_dw b_pw | awk 'NR == 1 { lo = $2; hi = $2 } { if ($2 < lo) lo = $2; if ($2 > hi) hi = $2; n++ } END { printf "%d %.5f %.5f", n, lo, hi }')"
  set -- $r
  awk -v lo="$2" -v hi="$3" 'BEGIN { exit !(lo > 0.95 && hi < 1.05) }' \
    && ok "model dw(plane)/pw = $2 to $3 at $1 points (p = |k_aA - alpha k_sF| instead of p(q = 0))" \
    || bad "dw(plane)/pw = $r"
  grep -q '^distortion_point' "$WORK/b_dw/output/thm_experiments.out" && bad "R(E) rows with vertexModel=dw" \
    || ok "no R(E) with vertexModel=dw"
fi

# (c) ----------------------------------------------------------------------
echo "(c) point Coulomb in d + 18O"
run c_dw "experiment[A] segments=1,2 $KIN distortion=coulomb vertexModel=dw" "NOFOLD"
run c_ps "experiment[A] segments=1,2 $KIN distortion=coulomb vertexModel=dw ps=hulthen:0-30 psNodes=8" "NOFOLD"
if ran c_dw && ran c_ps; then
  grep -q '^vertex: vertexModel=dw' "$WORK/c_dw/output/thm_experiments.out" && ok "thm_experiments.out: vertex line" \
    || bad "no vertex line"
  n="$(grep -c '^dw_vertex_point' "$WORK/c_dw/output/thm_experiments.out")"
  [ "$n" -ge 3 ] && ok "thm_experiments.out: $n dw_vertex_point rows" || bad "dw_vertex_point rows: $n"
  grep -q '^distortion_point' "$WORK/c_dw/output/thm_experiments.out" && bad "R(E) rows with vertexModel=dw" \
    || ok "R(E) not applied"
  grep -q "The distortion factor R(E) is not applied" "$WORK/c_dw/log" && ok "startup says R(E) is off" \
    || bad "no startup note"
  r="$(ratios c_dw b_dw | awk 'NR == 1 { lo = $2; hi = $2 } { if ($2 < lo) lo = $2; if ($2 > hi) hi = $2 } END { printf "%.4f %.4f", lo, hi }')"
  awk -v lo="${r% *}" -v hi="${r#* }" 'BEGIN { exit !(hi / lo > 1.001) }' \
    && ok "Coulomb/plane model ratio ${r% *} to ${r#* }: energy dependent" || bad "Coulomb/plane: $r"
  same "$WORK/c_dw/output/$OUT" "$WORK/c_ps/output/$OUT" && bad "a ps window changes nothing" \
    || ok "a ps window (nodes on its reachable part) changes the model"
fi

# (d) ----------------------------------------------------------------------
echo "(d) refused"
refuse() {  # refuse NAME MESSAGE-FRAGMENT BLOCK [HEAD]
  local block="$3"
  [ -z "${4:-}" ] || block="$4
$3"
  run "$1" "$block"
  if [ "$(cat "$WORK/$1/status")" != 0 ] && grep -q "ERROR: <thm>" "$WORK/$1/log" &&
     grep -q -- "$2" "$WORK/$1/log"; then
    ok "$1: $(grep -m1 'ERROR' "$WORK/$1/log" | tr -d '\r' | cut -c1-150)"
  else
    bad "$1: not refused with '$2' (status $(cat "$WORK/$1/status"))"; tail -3 "$WORK/$1/log" | sed 's/^/        /'
  fi
}
L="experiment[A] segments=1,2 $KIN"
printf '0.4 1.0\n1.0 2.0\n' > "$WORK/w.dat"
refuse value "expected pw or dw" "$L distortion=coulomb vertexModel=dwba"
refuse no_distortion "needs distortion=coulomb or distortion=optical" "$L vertexModel=dw"
refuse table "needs distortion=coulomb or distortion=optical" "experiment[A] segments=1,2 distortion=table:$WORK/w.dat vertexModel=dw"
refuse ref "belongs to the distortion factor R(E)" "$L distortion=coulomb distortionRef=0.7 vertexModel=dw"
refuse ratio "belongs to the distortion factor R(E)" "$L distortion=coulomb distortionRatio=dw vertexModel=dw"
refuse angle_ps "spectatorAngle= (one direction) and a ps window" "$L distortion=coulomb spectatorAngle=10 ps=hulthen:0-30 vertexModel=dw"
refuse theta "theta= (fixed-angle observable) is not available" "$L distortion=coulomb theta=30-60 vertexModel=dw"
refuse coulomb_integral "cannot be combined with coulombIntegral=1" "$L distortion=coulomb vertexModel=dw" "coulombIntegral=1"
refuse coherent "entranceL=coherent" "$L distortion=coulomb vertexModel=dw" "entranceL=coherent"
refuse spectator_energy "spectatorEnergy for entrance pair" "$L distortion=coulomb vertexModel=dw" "spectatorEnergy=0.3"
refuse reach "do not overlap the ps window" "$L distortion=coulomb vertexModel=dw ps=hulthen:150-200"

echo
if [ "$fail" -eq 0 ]; then echo "PASS: THM distorted-wave vertex"; else echo "FAIL"; exit 1; fi
