#!/usr/bin/env bash
#
# A lab-frame RATIO composite is the ratio of the two terms' LAB cross sections, each term taken to the lab with its
# OWN pair's Jacobian (audit item 15, 2026-10-27).
#
# A lab-frame differential RATIO segment (isDiff 1) is base / denominator in the lab: sigma_lab = sigma_cm / J, so the
# model is (sigma_cm,base / sigma_cm,den) (J_den / J_base).  Up to dev 1b6c346 both factors were computed with the BASE
# segment's entrance and exit pair and the base exit's type (ESegment::CalculateTheoreticalCrossSection): right for a
# denominator of the same pair, and for another particle pair at a general angle without relativistic=1 (the
# non-relativistic Jacobian uses only the point's own two angles), but wrong for another exit pair with relativistic=1
# or at 0 / 180 deg (the mass formula), and for a photon over a particle or the reverse (the wrong kind of Jacobian).
# In mode 3 the lab output (kinflag.dat = 1) of a RATIO test segment was also divided by J_base, because
# EData::MakePoints converted the test points' cross sections while ESegment::Fill skips that for a ratio.
#
# Models from make_models.py (13C + alpha with two neutron exits of the same masses and different Q, the
# tests/sum_lab_frame model; 12C(a,g) external capture, the tests/gamma_lab_frame model), all test segments in
# mode 3 with lab output; the expected ratio is the ratio of the plain segments' lab values:
#   (a) particle, non-relativistic: R1 = n0+/n3- at the same angles, R2 = n0+/(a,a)(60 deg), R3 = n0+/n3-(30 deg),
#       R4 = n0+/n3-(0 deg), R5 = n0+/n0+(60 deg); R6 the c.m. ratio and R7 the angle-integrated ratio = the ratios
#       of the c.m. / angle-integrated outputs (no frame factor);
#   (b) particle_rel (relativistic=1): the same;
#   (c) gamma: RG1 = W(theta)/W(90) (a gamma ratio of one pair), RG2 = (a,g)/(a,a)(60 deg), RG3 = (a,a)/(a,g)(90 deg);
#   (d) the c.m. output (kinflag.dat = 0) of the lab RATIO segments is the same ratio (nothing divides it);
#   (e) mode 1, relativistic=1: a RATIO data segment n0+/(a,a)(60 deg) whose data are the plain segments' lab ratio
#       has chi2 ~ 0;
#   (f) the analytic gradient of (e) with free widths = finite differences (AZURE_GRAD_CHECK).
# dev 1b6c346: (a) R1 3 % at 0 deg and 6.5 % at 180 deg, R4 3.2 %; (b) R1 6.6 %, R2 3.0 %, R3 2.8 %, R4 3.2 %;
# (c) RG2 31 %, RG3 1e-4; (d) every lab RATIO output 17-20 % low (so (a)-(c) fail at the 80 % level there); (e) chi2 143 over
# 48 points at 1 %.
#
#   ./tests/ratio_frame/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/ratio_frame.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
le() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r + 0 <= t + 0) }'; }

# run DIR MODEL MODE KINFLAG [FLAGS] [ENV...] -> $WORK/DIR (a data/ directory placed there beforehand is kept)
run() {
  local dir="$WORK/$1" model="$2" mode="$3" kin="$4" flags="${5:-}"
  shift 5 2> /dev/null || shift $#
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$model.azr" "$dir/t.azr"
  [ "$kin" = 1 ] && echo 1 > "$dir/kinflag.dat"
  # shellcheck disable=SC2086
  (cd "$dir" && for kv in "$@"; do export "$kv"; done
   printf "$mode\n\n\n7\n" | run_limited 600 "$AZURE2_BIN" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}
# block FILE N -> the N-th blank-line-separated block (segments with the same pairs share a file, in file order)
block() { awk -v blk="$2" 'BEGIN { RS = ""; FS = "\n" } NR == blk { print }' "$1"; }
# rows FILE N -> "E_x angle value" of block N (E_x is the same in both output frames; the value precedes the S-factor)
rows() { block "$1" "$2" | awk 'NF >= 5 { printf "%.10f %.10f %.12e\n", $2, $3, $(NF - 1) }'; }
# ratio_dev NUM DEN RATIO -> the largest relative difference between RATIO and NUM / DEN; a denominator with one row
# per energy (a fixed-angle term) is matched by energy, otherwise by energy and angle
ratio_dev() {
  paste "$1" "$3" > "$WORK/nr"
  awk 'FNR == NR { d[$1, $2] = $3; n[$1]++; s[$1] = $3; next }
       { if ($1 != $4 || $2 != $5) { print "grid"; exit }
         den = (n[$1] == 1) ? s[$1] : d[$1, $2]; if (den == "" ) { print "nodenominator"; exit }
         e = $3 / den; x = $6 - e; if (x < 0) x = -x; m = e < 0 ? -e : e; if (m < 1e-300) m = 1e-300
         r = x / m; if (r > w) w = r; c++ }
       END { if (c == 0) print "nan"; else printf "%.3e\n", w }' "$2" "$WORK/nr"
}
# reldiff A B -> the largest relative difference of the third columns
reldiff() {
  paste "$1" "$2" | awk '{ x = $3 - $6; if (x < 0) x = -x; m = $6 < 0 ? -$6 : $6; if (m < 1e-300) m = 1e-300
                           r = x / m; if (r > w) w = r } END { printf "%.3e\n", w }'
}

# ---- (a), (b) particle ratios --------------------------------------------------------------------------------------
# R=3 file: T1, T6 (c.m.), T8 (integrated), R1 .. R7 ; R=5 file: T2, T4 (30 deg), T5 (0 deg), T7 (c.m.), T9 ; R=2: T3
particle_checks() {  # TAG MODEL
  local tag="$1" model="$2" dir="$WORK/$2"
  run "$model" "$model" 3 1
  local a="$dir/output/AZUREOut_aa=2_R=3.extrap" b="$dir/output/AZUREOut_aa=2_R=5.extrap" c="$dir/output/AZUREOut_aa=2_R=2.extrap"
  if [ ! -s "$a" ] || [ ! -s "$b" ] || [ ! -s "$c" ]; then bad "$tag $model: an extrapolation output is missing"; tail -5 "$dir/run.log"; return; fi
  local p="$WORK/$model"
  rows "$a" 1 > "$p.t1"; rows "$b" 1 > "$p.t2"; rows "$c" 1 > "$p.t3"; rows "$b" 2 > "$p.t4"; rows "$b" 3 > "$p.t5"
  rows "$a" 2 > "$p.t6"; rows "$b" 4 > "$p.t7"; rows "$a" 3 > "$p.t8"; rows "$b" 5 > "$p.t9"
  local n; n="$(wc -l < "$p.t1" | tr -d ' ')"
  if [ "$n" -lt 50 ]; then bad "$tag $model: $n points in the n0+ extrapolation"; tail -5 "$dir/run.log"; return; fi
  awk '$2 == 60 { print }' "$p.t1" > "$p.t1_60"
  local i=4 r
  for case in "R1 n0+/n3- at the same angles:t1:t2" "R2 n0+/(a,a)(60 deg):t1:t3" "R3 n0+/n3-(30 deg):t1:t4" \
              "R4 n0+/n3-(0 deg):t1:t5" "R5 n0+/n0+(60 deg), the same pair:t1:t1_60" "R6 c.m. n0+/n3-:t6:t7" \
              "R7 angle-integrated n0+/n3-:t8:t9"; do
    local name="${case%%:*}" rest="${case#*:}" num den
    num="${rest%%:*}"; den="${rest#*:}"
    rows "$a" "$i" > "$p.r$i"
    r="$(ratio_dev "$p.$num" "$p.$den" "$p.r$i")"
    if le "$r" 1e-9; then ok "$tag $name = the ratio of the plain lab values (max rel $r)"; else bad "$tag $name vs the ratio of the plain lab values (max rel $r)"; fi
    i=$((i + 1))
  done
}
particle_checks "(a)" particle
particle_checks "(b)" particle_rel
# not vacuous: the lab ratio differs from the c.m. ratio (R1 vs R6, non-relativistic), and the relativistic
# denominator Jacobian differs from the non-relativistic one (R2 of the two models)
if [ -s "$WORK/particle.r4" ] && [ -s "$WORK/particle.r9" ]; then
  r="$(reldiff "$WORK/particle.r4" "$WORK/particle.r9")"
  if le 1e-3 "$r"; then ok "(a) the frame factors matter: lab ratio R1 differs from the c.m. ratio R6 (max rel $r)"; else bad "(a) vacuous: R1 = R6 (max rel $r)"; fi
fi
if [ -s "$WORK/particle.r5" ] && [ -s "$WORK/particle_rel.r5" ]; then
  r="$(reldiff "$WORK/particle_rel.r5" "$WORK/particle.r5")"
  if le 1e-3 "$r"; then ok "(b) relativistic=1 changes the n0+/(a,a)(60) ratio (max rel $r)"; else bad "(b) vacuous: relativistic R2 = non-relativistic R2 (max rel $r)"; fi
fi

# ---- (c) gamma ratios ----------------------------------------------------------------------------------------------
# R=2 file: G1, G2 (90 deg), RG1, RG2 ; R=1 file: G3 (60 deg), G4, RG3
run gamma gamma 3 1
ga="$WORK/gamma/output/AZUREOut_aa=1_R=2.extrap"; gb="$WORK/gamma/output/AZUREOut_aa=1_R=1.extrap"
if [ -s "$ga" ] && [ -s "$gb" ]; then
  p="$WORK/gamma"
  rows "$ga" 1 > "$p.g1"; rows "$ga" 2 > "$p.g2"; rows "$gb" 1 > "$p.g3"; rows "$gb" 2 > "$p.g4"
  rows "$ga" 3 > "$p.rg1"; rows "$ga" 4 > "$p.rg2"; rows "$gb" 3 > "$p.rg3"
  n="$(wc -l < "$p.g1" | tr -d ' ')"
  [ "$n" -ge 20 ] || { bad "(c) $n points in the (a,g) extrapolation"; tail -5 "$WORK/gamma/run.log"; }
  r="$(ratio_dev "$p.g1" "$p.g2" "$p.rg1")"
  if le "$r" 1e-9; then ok "(c) RG1 W(theta)/W(90) = the ratio of the plain lab values (max rel $r, $n points)"; else bad "(c) RG1 W(theta)/W(90) (max rel $r)"; fi
  r="$(ratio_dev "$p.g1" "$p.g3" "$p.rg2")"
  if le "$r" 1e-9; then ok "(c) RG2 (a,g)/(a,a)(60 deg) = the ratio of the plain lab values (max rel $r)"; else bad "(c) RG2 (a,g)/(a,a)(60 deg) (max rel $r)"; fi
  r="$(ratio_dev "$p.g4" "$p.g2" "$p.rg3")"
  if le "$r" 1e-9; then ok "(c) RG3 (a,a)/(a,g)(90 deg) = the ratio of the plain lab values (max rel $r)"; else bad "(c) RG3 (a,a)/(a,g)(90 deg) (max rel $r)"; fi
else
  bad "(c) the gamma extrapolation wrote no output"; tail -5 "$WORK/gamma/run.log"
fi

# ---- (d) the c.m. output of a lab RATIO test segment is the same ratio ---------------------------------------------
run particle_cm particle 3 0
pa="$WORK/particle_cm/output/AZUREOut_aa=2_R=3.extrap"
if [ -s "$pa" ] && [ -s "$WORK/particle.r5" ]; then
  rows "$pa" 5 > "$WORK/particle_cm.r5"
  r="$(reldiff "$WORK/particle_cm.r5" "$WORK/particle.r5")"
  if le "$r" 1e-12; then ok "(d) kinflag 0 and 1 print the same lab ratio for R2 (max rel $r)"; else bad "(d) the lab output of the RATIO test segment is scaled (kinflag 0 vs 1: max rel $r)"; fi
else
  bad "(d) the c.m. run wrote no output"; tail -5 "$WORK/particle_cm/run.log"
fi

# ---- (e) a RATIO data segment, (f) its gradient ---------------------------------------------------------------------
if [ -s "$WORK/particle_rel.t1" ] && [ -s "$WORK/particle_rel.t3" ]; then
  # data/ratio.dat: E_lab, angle (15 .. 165 step 30), n0+ / (a,a)(60) from the plain lab values, 1 % errors
  block "$WORK/particle_rel/output/AZUREOut_aa=2_R=3.extrap" 1 | awk 'NF >= 5 { printf "%.10f %.10f %.12e\n", $1, $3, $4 }' > "$WORK/t1lab"
  block "$WORK/particle_rel/output/AZUREOut_aa=2_R=2.extrap" 1 | awk 'NF >= 5 { printf "%.10f %.12e\n", $1, $4 }' > "$WORK/t3lab"
  awk 'FNR == NR { d[$1] = $2; next } ($2 - 15) % 30 == 0 && $2 >= 15 && $2 <= 165 { v = $3 / d[$1]; printf "%.6f  %.1f  %.12e  %.12e\n", $1, $2, v, 0.01 * v }' \
      "$WORK/t3lab" "$WORK/t1lab" > "$WORK/ratio.dat"
  nd="$(wc -l < "$WORK/ratio.dat" | tr -d ' ')"
  mkdir -p "$WORK/rdata/data" && cp "$WORK/ratio.dat" "$WORK/rdata/data/ratio.dat"
  run rdata ratio_data 1 1
  chi="$(awk -F, 'NR > 1 && $1 + 0 == 1 { gsub(/ /, "", $2); print $2 }' "$WORK/rdata/output/chiSquared.out")"
  if [ -z "$chi" ]; then
    bad "(e) no chi2 for the RATIO data segment"; tail -5 "$WORK/rdata/run.log"
  elif le "$chi" 1e-6; then
    ok "(e) RATIO data segment n0+/(a,a)(60 deg), relativistic=1, against the plain segments' lab ratio: chi2 $chi ($nd points)"
  else
    bad "(e) RATIO data segment vs the plain segments' lab ratio: chi2 $chi (expected ~0, $nd points)"
  fi
  # (f): data off the model (a few % per point) so chi2 and its gradient are not at the rounding floor
  mkdir -p "$WORK/rfree/data"
  awk '{ printf "%s  %s  %.12e  %s\n", $1, $2, $3 * (1 + 0.02 * (NR % 5 - 2)), $4 }' "$WORK/ratio.dat" > "$WORK/rfree/data/ratio.dat"
  run rfree ratio_free 2 0 "--use-gradient" AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_H=1e-7
  awk '/^\[grad-check\] idx / { a = f = ""; for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1) }
         print a, f }' "$WORK/rfree/run.log" > "$WORK/rows"
  nrows="$(wc -l < "$WORK/rows" | tr -d ' ')"
  if [ "$nrows" -lt 5 ]; then
    bad "(f) gradient check listed $nrows parameters"; grep -i "grad\|error" "$WORK/rfree/run.log" | head -5
  else
    worst="$(awk '{ a = $1 < 0 ? -$1 : $1; if (a > m) m = a; d = $1 - $2; if (d < 0) d = -d; if (d > w) w = d }
               END { printf "%.2e", w / m }' "$WORK/rows")"
    if le "$worst" 1e-5; then ok "(f) analytic gradient of the RATIO segment = finite differences, $nrows parameters, worst $worst of the largest"
    else bad "(f) gradient vs finite differences, worst $worst of the largest (> 1e-5)"; fi
  fi
else
  bad "(e) no relativistic plain outputs to build the ratio data from"
fi

if [ "$failures" -eq 0 ]; then echo "ratio_frame: PASS"; else echo "ratio_frame: $failures check(s) FAILED"; fi
[ "$failures" -eq 0 ]
