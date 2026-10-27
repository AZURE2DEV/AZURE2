#!/usr/bin/env bash
#
# A lab-frame differential SUM composite adds its terms in the lab frame (audit item 7, 2026-10-26).
#
# A lab-frame differential segment is compared in the c.m. frame of its base pair (sigma_cm = J_base sigma_lab).
# A component to another exit pair is evaluated in its own c.m. frame, with its own Jacobian, so the model of the
# segment is sum_i s_i (J_base / J_i) sigma_cm,i (ESegment::TermFrameFactor).  Before the fix the c.m. terms were
# added as they were: 13C(a,n0+ + n3-) at 6 MeV differed by -3 % (0 deg) to +6 % (150 deg) from the lab sum.
# Models from make_models.py (13C + alpha, exits n + 16O(6.05, 0+) and n + 16O(6.13, 3-), same masses, different Q):
#
#   (a) mode 3 with lab output (kinflag.dat = 1): SUM(T1 + T2) = T1 + T2, SUM of the same pair twice = 2 T1,
#       T1 + 0.5 T2 with the scaling token; the c.m. SUM = the c.m. terms; and the two frames really differ
#       (the c.m. SUM is not the sum of the lab outputs);
#   (b) mode 1: a SUM data segment with a fixed 0.05 MeV energy shift, whose data are the fitted lab values of the
#       two plain segments under the same shift, has chi2 ~ 0 (the frame factors follow the shift);
#   (c) the analytic gradient (mapped-components path and the per-term path) = finite differences on the SUM
#       segment with free widths (AZURE_GRAD_CHECK).
#
#   ./tests/sum_lab_frame/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/sum_lab_frame.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
le() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r + 0 <= t + 0) }'; }

# run DIR MODEL MODE KINFLAG [FLAGS] [ENV...] -> $WORK/DIR
run() {
  local dir="$WORK/$1" model="$2" mode="$3" kin="$4" flags="${5:-}"
  shift 5 2> /dev/null || shift $#
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$model.azr" "$dir/t.azr"
  cp -r "$HERE/models/data" "$dir/"
  [ "$kin" = 1 ] && echo 1 > "$dir/kinflag.dat"
  # shellcheck disable=SC2086
  (cd "$dir" && for kv in "$@"; do export "$kv"; done
   printf "$mode\n\n\n7\n" | run_limited 300 "$AZURE2_BIN" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}
# block FILE N -> the N-th blank-line-separated block of an output file (segments with the same pairs share a file)
block() { awk -v blk="$2" 'BEGIN { RS = ""; FS = "\n" } NR == blk { print }' "$1"; }
# xs FILE N -> energy, angle and cross section (the column before the S-factor) of block N
xs() { block "$1" "$2" | awk 'NF >= 5 { printf "%.12e %.12e %.12e\n", $1, $3, $(NF - 1) }'; }
# maxrel A B -> largest relative difference of the third columns (and a check that energies and angles agree)
maxrel() {
  paste "$1" "$2" | awk '
    { if ($1 != $4 || $2 != $5) { print "grid"; exit }
      d = $3 - $6; if (d < 0) d = -d; m = $6 < 0 ? -$6 : $6; if (m < 1e-300) m = 1e-300
      r = d / m; if (r > w) w = r; n++ }
    END { if (n == 0) print "nan"; else printf "%.3e\n", w }'
}

# ---- (a) extrapolations ------------------------------------------------------------------------------------------
run lab extrap 3 1
run cm extrap 3 0
f3="$WORK/lab/output/AZUREOut_aa=2_R=3.extrap"; f5="$WORK/lab/output/AZUREOut_aa=2_R=5.extrap"
c3="$WORK/cm/output/AZUREOut_aa=2_R=3.extrap"; c5="$WORK/cm/output/AZUREOut_aa=2_R=5.extrap"
if [ -s "$f3" ] && [ -s "$f5" ] && [ -s "$c3" ] && [ -s "$c5" ]; then
  xs "$f3" 1 > "$WORK/t1"; xs "$f5" 1 > "$WORK/t2"; xs "$f3" 2 > "$WORK/t3"; xs "$f3" 3 > "$WORK/t4"
  xs "$f3" 4 > "$WORK/t5"; xs "$f5" 2 > "$WORK/t6"; xs "$f3" 5 > "$WORK/t7"; xs "$f3" 6 > "$WORK/t8"
  n="$(wc -l < "$WORK/t1" | tr -d ' ')"
  if [ "$n" -lt 50 ]; then bad "(a) $n points in the lab extrapolation"; tail -5 "$WORK/lab/run.log"; fi
  paste "$WORK/t1" "$WORK/t2" | awk '{ printf "%.12e %.12e %.12e\n", $1, $2, $3 + $6 }' > "$WORK/t12"
  r="$(maxrel "$WORK/t3" "$WORK/t12")"
  if le "$r" 1e-9; then ok "(a) lab SUM(n0+ + n3-) = n0+ + n3- in the lab (max rel $r, $n points)"; else bad "(a) lab SUM vs the lab sum of the terms (max rel $r)"; fi
  awk '{ printf "%.12e %.12e %.12e\n", $1, $2, 2 * $3 }' "$WORK/t1" > "$WORK/t11"
  r="$(maxrel "$WORK/t4" "$WORK/t11")"
  if le "$r" 1e-9; then ok "(a) SUM of the same pair twice = 2 x the plain segment (max rel $r)"; else bad "(a) SUM of the same pair twice vs 2 x plain (max rel $r)"; fi
  paste "$WORK/t1" "$WORK/t2" | awk '{ printf "%.12e %.12e %.12e\n", $1, $2, $3 + 0.5 * $6 }' > "$WORK/t12h"
  r="$(maxrel "$WORK/t8" "$WORK/t12h")"
  if le "$r" 1e-9; then ok "(a) scaling token: SUM = n0+ + 0.5 n3- in the lab (max rel $r)"; else bad "(a) scaled SUM vs n0+ + 0.5 n3- (max rel $r)"; fi
  # the c.m. segments (isDiff 5): no frame factor, the c.m. sum
  xs "$c3" 4 > "$WORK/c5"; xs "$c5" 2 > "$WORK/c6"; xs "$c3" 5 > "$WORK/c7"
  paste "$WORK/c5" "$WORK/c6" | awk '{ printf "%.12e %.12e %.12e\n", $1, $2, $3 + $6 }' > "$WORK/c56"
  r="$(maxrel "$WORK/c7" "$WORK/c56")"
  if le "$r" 1e-9; then ok "(a) c.m. SUM = the sum of the c.m. terms (max rel $r)"; else bad "(a) c.m. SUM vs the c.m. sum (max rel $r)"; fi
  # the frames differ: the lab terms added in the c.m. output of the lab SUM (block 2 of the c.m. run) are not the
  # c.m. sum of the lab terms' own c.m. values (blocks 1 of each file) -- the Jacobians of the two pairs differ
  xs "$c3" 1 > "$WORK/l1cm"; xs "$c5" 1 > "$WORK/l2cm"; xs "$c3" 2 > "$WORK/l3cm"
  paste "$WORK/l1cm" "$WORK/l2cm" | awk '{ printf "%.12e %.12e %.12e\n", $1, $2, $3 + $6 }' > "$WORK/l12cm"
  r="$(paste "$WORK/l3cm" "$WORK/l12cm" | awk '{ d = $3 - $6; if (d < 0) d = -d; m = $6 < 0 ? -$6 : $6; if (m < 1e-300) m = 1e-300; r = d / m; if (r > w) w = r } END { printf "%.3e", w }')"
  if le 1e-3 "$r"; then ok "(a) the Jacobians of the two pairs differ: the lab SUM is not the plain c.m. sum (max rel $r)"; else bad "(a) the test is vacuous: lab SUM = plain c.m. sum (max rel $r)"; fi
else
  bad "(a) an extrapolation is missing"; tail -5 "$WORK/lab/run.log"; tail -5 "$WORK/cm/run.log"
fi

# ---- (b) a SUM data segment under an energy shift ---------------------------------------------------------------
run plain plain_data 1 1
p3="$WORK/plain/output/AZUREOut_aa=2_R=3.out"; p5="$WORK/plain/output/AZUREOut_aa=2_R=5.out"
if [ -s "$p3" ] && [ -s "$p5" ]; then
  # data/sum.dat: the data grid (unshifted lab energies and angles) with the sum of the fitted lab values, 1 % errors
  paste "$HERE/models/data/d1.dat" <(awk 'NF >= 5 { print $4 }' "$p3") <(awk 'NF >= 5 { print $4 }' "$p5") |
    awk '{ v = $5 + $6; printf "%.10g  %.10g  %.12e  %.12e\n", $1, $2, v, 0.01 * v }' > "$WORK/sum.dat"
  nd="$(wc -l < "$WORK/sum.dat" | tr -d ' ')"
  mkdir -p "$WORK/sumdata/data" && cp "$WORK/sum.dat" "$WORK/sumdata/data/sum.dat"   # run() merges models/data in
  run sumdata sum_data 1 1
  chi="$(awk -F, 'NR > 1 && $1 + 0 == 1 { gsub(/ /, "", $2); print $2 }' "$WORK/sumdata/output/chiSquared.out")"
  if [ -z "$chi" ]; then
    bad "(b) no chi2 for the SUM data segment"; tail -5 "$WORK/sumdata/run.log"
  elif le "$chi" 1e-6; then
    ok "(b) SUM data segment with a 0.05 MeV shift against the shifted plain segments' lab sum: chi2 $chi ($nd points)"
  else
    bad "(b) SUM data segment vs the shifted plain segments' lab sum: chi2 $chi (expected ~0)"
  fi
else
  bad "(b) the plain data run wrote no output"; tail -5 "$WORK/plain/run.log"
fi

# ---- (c) gradients -----------------------------------------------------------------------------------------------
gradcheck() {  # NAME ENV...
  local name="$1"; shift
  # data off the model (a few % per point), so chi2 and its gradient are not at the rounding floor
  mkdir -p "$WORK/$name/data"
  awk '{ printf "%s  %s  %.12e  %s\n", $1, $2, $3 * (1 + 0.02 * (NR % 5 - 2)), $4 }' "$WORK/sum.dat" > "$WORK/$name/data/sum.dat"
  run "$name" sum_free 2 0 "--use-gradient" AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_H=1e-7 "$@"
  awk '/^\[grad-check\] idx / { idx = $3; a = f = ""; j = "-"
         for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
         print idx, a, f, j }' "$WORK/$name/run.log" | sort -n > "$WORK/$name.rows"
  local nrows worst
  nrows="$(wc -l < "$WORK/$name.rows" | tr -d ' ')"
  if [ "$nrows" -lt 5 ]; then
    bad "(c) $name: gradient check listed $nrows parameters"; grep -i "grad\|error" "$WORK/$name/run.log" | head -5; return
  fi
  worst="$(awk '{ a = $2 < 0 ? -$2 : $2; if (a > m) m = a; g[NR] = $2; f[NR] = $3; j[NR] = $4 }
             END { for (i = 1; i <= NR; i++) { d = g[i] - f[i]; if (d < 0) d = -d; e = j[i] - f[i]; if (e < 0) e = -e
                                               if (d > w) w = d; if (j[i] != "-" && e > w) w = e }
                   printf "%.2e", w / m }' "$WORK/$name.rows")"
  if le "$worst" 1e-5; then
    ok "(c) $name: adjoint and Jacobian-path gradients = finite differences, $nrows parameters, worst $worst of the largest"
  else
    bad "(c) $name: gradient vs finite differences, worst $worst of the largest (> 1e-5)"
  fi
}
if [ -s "$WORK/sum.dat" ]; then
  gradcheck mapped
  gradcheck perterm AZURE_NO_COMPONENT_MAPPING=1
else
  bad "(c) no data for the gradient check"
fi

echo
if [ "$failures" -eq 0 ]; then echo "sum_lab_frame: PASS"; else echo "sum_lab_frame: FAIL ($failures)"; exit 1; fi
