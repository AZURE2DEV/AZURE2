#!/usr/bin/env bash
#
# Coulomb consistency of a THM experiment (docs/source/theory/
# thm_implementation.rst, "Coulomb effects: what each option contains";
# CheckThmCoulombConsistency and CheckThmExperiments in src/ThmExperiment.cpp).
# The options that carry Coulomb physics -- the Coulomb term C_l
# (coulombIntegral), the line shape N_C (lineshape=on), the distortion factor
# R(E) (distortion=) and the DW vertex (vertexModel=dw) -- may be combined
# only where they do not count the same interaction twice.
#
# Model: tests/18O_p_a_thm with the made-up charged-spectator reaction
# 18O(3He,a15N)d at 115 MeV (and its real 2H(18O,a15N)n at 54 MeV).
#
#   (a) refused: coulombIntegral=1 with R(E) on a distorted a + A wave
#       (distortion=coulomb; optical with opticalAA coulomb or a potential;
#       the global key before or after the experiment line), with
#       vertexModel=dw, and a ps window with distortionRatio=dw;
#   (b) allowed: coulombIntegral=1 with opticalAA=plane (C_l is then the
#       external term of the plane a + A wave, and is applied), with
#       distortion=table (a WARNING), coulombIntegral=0 written out
#       byte-identical to no key, a ps window with distortionRatio=dwpw;
#   (c) N_C and R(E) are separate factors: without folding,
#       model(R, N_C)/model(R) == model(N_C)/model(none) at every point
#       (1e-9); a neutral spectator (zeta = 0) leaves the DW vertex
#       byte-identical with lineshape=on.
#
#   ./tests/thm_coulomb_consistency/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
O18="$HERE/../18O_p_a_thm"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_coulomb_consistency.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
if command -v timeout >/dev/null 2>&1; then RUN="timeout ${TEST_TIMEOUT:-300}"; else RUN=""; fi
OUT="AZUREOut_aa=1_R=2.out"
KIN="beam=18O target=3He spectator=d Ebeam=115"
KINN="beam=18O target=d spectator=n Ebeam=54"
WS="20,1.2,0.6,5,1.3,0.6,0,0,0,1.3"
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
# "E ratio" per point: model(A)/model(B), energies must agree.
ratios() {
  paste -d' ' <(tr -d '\r' < "$WORK/$1/output/$OUT" | awk 'NF > 4 { print $1, $4 }') \
              <(tr -d '\r' < "$WORK/$2/output/$OUT" | awk 'NF > 4 { print $1, $4 }') |
    awk '{ if ($1 != $3) { print "ENERGY", $1, $3; exit 1 } printf "%.12e %.12e\n", $1, $2 / $4 }'
}
refuse() {  # refuse NAME MESSAGE-FRAGMENT BLOCK
  run "$1" "$3"
  if [ "$(cat "$WORK/$1/status")" != 0 ] && grep -q "ERROR: <thm>" "$WORK/$1/log" &&
     grep -q -- "$2" "$WORK/$1/log"; then
    ok "$1: $(grep -m1 'ERROR' "$WORK/$1/log" | tr -d '\r' | cut -c1-150)"
  else
    bad "$1: not refused with '$2' (status $(cat "$WORK/$1/status"))"; tail -3 "$WORK/$1/log" | sed 's/^/        /'
  fi
}
NOFOLD='/<targetInt>/ { print; T = 1; next } /<\/targetInt>/ { T = 0 } T { next } { print }'
L="experiment[A] segments=1,2 $KIN"

# (a) ----------------------------------------------------------------------
echo "(a) refused combinations"
refuse ci_coulomb "would be counted twice" "coulombIntegral=1
$L distortion=coulomb"
refuse ci_after "would be counted twice" "$L distortion=coulomb
coulombIntegral=1"
refuse ci_optical "would be counted twice" "coulombIntegral=1
$L distortion=optical opticalSF=plane"
refuse ci_optical_ws "would be counted twice" "coulombIntegral=1
$L distortion=optical opticalAA=$WS"
refuse ci_dw "cannot be combined with coulombIntegral=1" "coulombIntegral=1
$L distortion=coulomb vertexModel=dw"
refuse ps_dwratio "Use distortionRatio=dwpw" "$L ps=hulthen:0-30 distortion=coulomb distortionRatio=dw"

# (b) ----------------------------------------------------------------------
echo "(b) allowed combinations"
printf '0.0 1.0\n5.0 2.0\n' > "$WORK/w.dat"
run plain_aa "$L distortion=optical opticalAA=plane"
run plain_aa_ci "coulombIntegral=1
$L distortion=optical opticalAA=plane"
if ran plain_aa && ran plain_aa_ci; then
  ok "coulombIntegral=1 with opticalAA=plane runs"
  same "$WORK/plain_aa/output/$OUT" "$WORK/plain_aa_ci/output/$OUT" && bad "C_l not applied with opticalAA=plane" ||
    ok "C_l applied with opticalAA=plane"
fi
run table_ci "coulombIntegral=1
$L distortion=table:$WORK/w.dat"
if ran table_ci; then
  grep -q "WARNING: <thm> experiment\[A\]: coulombIntegral=1 with distortion=table" "$WORK/table_ci/log" &&
    ok "distortion=table with coulombIntegral=1: runs, warned" || bad "distortion=table with coulombIntegral=1: no warning"
fi
run r_nokey "$L distortion=coulomb"
run r_ci0 "coulombIntegral=0
$L distortion=coulomb"
if ran r_nokey && ran r_ci0; then
  grep -q "WARNING: <thm> experiment\[A\]: coulombIntegral" "$WORK/r_nokey/log" &&
    bad "a Coulomb-consistency warning without the combinations" || ok "no Coulomb-consistency warning without the combinations"
  for f in chiSquared.out normalizations.out thm_experiments.out "$OUT"; do
    same "$WORK/r_nokey/output/$f" "$WORK/r_ci0/output/$f" && ok "coulombIntegral=0: $f byte-identical" ||
      bad "coulombIntegral=0: $f differs"
  done
fi
run ps_dwpw "$L ps=hulthen:0-30 distortion=coulomb distortionRatio=dwpw"
ran ps_dwpw && ok "a ps window with distortionRatio=dwpw runs"

# (c) ----------------------------------------------------------------------
echo "(c) N_C and R(E) are separate factors (no folding)"
run c_none "$L" NOFOLD
run c_ls "$L lineshape=on" NOFOLD
run c_r "$L distortion=coulomb" NOFOLD
run c_r_ls "$L distortion=coulomb lineshape=on" NOFOLD
if ran c_none && ran c_ls && ran c_r && ran c_r_ls; then
  paste -d' ' <(ratios c_r_ls c_r) <(ratios c_ls c_none) > "$WORK/fac.txt"
  worst=$(awk '{ d = $2 / $4 - 1; if (d < 0) d = -d; if (d > m) m = d } END { printf "%.3e", m }' "$WORK/fac.txt")
  spread=$(awk 'NR == 1 { lo = $4; hi = $4 } { if ($4 < lo) lo = $4; if ($4 > hi) hi = $4 } END { printf "%.3f-%.3f", lo, hi }' "$WORK/fac.txt")
  awk -v w="$worst" 'BEGIN { exit !(w + 0 <= 1e-9) }' &&
    ok "model(R,N_C)/model(R) == model(N_C)/model: worst $worst (N_C factor $spread)" ||
    bad "N_C and R do not factorize: worst $worst"
  awk -v s="$spread" 'BEGIN { split(s, a, "-"); exit !(a[2] - a[1] > 0.05) }' &&
    ok "the line shape changes the model ($spread)" || bad "the line shape does nothing here ($spread)"
fi
run n_dw "experiment[A] segments=1,2 $KINN distortion=coulomb vertexModel=dw" NOFOLD
run n_dw_ls "experiment[A] segments=1,2 $KINN distortion=coulomb vertexModel=dw lineshape=on" NOFOLD
if ran n_dw && ran n_dw_ls; then
  same "$WORK/n_dw/output/$OUT" "$WORK/n_dw_ls/output/$OUT" && ok "neutral spectator: DW vertex with lineshape=on byte-identical" ||
    bad "neutral spectator: lineshape=on changes the DW vertex model"
fi

echo
if [ "$fail" -eq 0 ]; then echo "PASS: THM Coulomb consistency"; else echo "FAIL"; exit 1; fi
