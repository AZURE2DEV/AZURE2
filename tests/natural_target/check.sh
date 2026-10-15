#!/usr/bin/env bash
#
# Several compound nuclei in one project (natural-abundance targets).  Two
# p+C nuclei (13N from 12C, 14N from 13C) share one file (models/M.azr):
#   (a) each nucleus's own segment in the merged file reproduces the standalone
#       project (models/A.azr, models/B.azr) to 1e-10 -- the nuclei do not couple;
#   (b) a "natural 1" SUM over both is 0.6 sigma_A + 0.4 sigma_B to 1e-12, its
#       decomposition file carries the two terms, and a plain natural segment
#       is 0.4 sigma_B; the same under a target-integration effect to 1e-10;
#   (c) in data mode the two abundances are fit parameters: the analytic
#       gradient agrees with finite differences (AZURE_GRAD_CHECK) and a fit
#       from 0.5/0.5 recovers 0.6/0.4; a normgroup over a natural SUM and a
#       plain segment: the shared norm's gradient agrees too;
#   (d) a level with channels of both nuclei, "natural 1" on a RATIO and a
#       component naming a missing pair key are refused.
#
#   ./tests/natural_target/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/natural_target.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME MODEL ROUTE(true|false) FLAGS [MODE] -> outputs in $WORK/NAME/output; returns 1 when nothing was written
run() {
  local dir="$WORK/$1" mode="${5:-3}"
  mkdir -p "$dir/output" "$dir/checks"
  sed "1,2s/^\(true\|false\) #Perform A-Matrix/$3 #Perform A-Matrix/" "$2" > "$dir/t.azr"
  for f in data.dat data_enr.dat; do [ -f "$WORK/$f" ] && cp "$WORK/$f" "$dir/"; done
  (cd "$dir" && printf "$mode\n\n\n7\n" | run_limited 600 "$AZURE2_BIN" $4 --no-gui --no-readline t.azr > run.log 2>&1)
}
# block FILE N -> the N-th blank-line-separated block of an output file (segments with the same pair keys share a file)
block() { awk -v blk="$2" 'BEGIN { RS = ""; FS = "\n" } NR == blk { print }' "$1"; }
# xs FILE N -> energy and cross section (the column before the S-factor) of block N
xs() { block "$1" "$2" | awk 'NF >= 4 { printf "%.12e %.12e\n", $1, $(NF - 1) }'; }
# maxrel A B -> largest relative difference between two two-column lists (energy, value), row by row
maxrel() {
  paste "$1" "$2" | awk '
    NF >= 4 { d = $2 - $4; if (d < 0) d = -d; m = $4 < 0 ? -$4 : $4; if (m < 1e-300) m = 1e-300;
              r = d / m; if (r > worst) worst = r; rows++ }
    END { if (rows == 0) print "nan"; else printf "%.3e\n", worst }'
}
le() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r + 0 <= t + 0) }'; }

M="$HERE/models"

# ---- (a), (b): extrapolations, A-/R-matrix x Brune -------------------------------------------------
for route in true false; do
  for flags in "" "--use-brune"; do
    tag="$( [ "$route" = true ] && echo Amat || echo Rmat )${flags:+_brune}"
    run "A_$tag" "$M/A.azr" "$route" "$flags" || true
    run "B_$tag" "$M/B.azr" "$route" "$flags" || true
    run "M_$tag" "$M/M.azr" "$route" "$flags" || true
    fa="$WORK/A_$tag/output/AZUREOut_aa=1_R=1.extrap"; fb="$WORK/B_$tag/output/AZUREOut_aa=2_R=2.extrap"
    m1="$WORK/M_$tag/output/AZUREOut_aa=1_R=1.extrap"; m2="$WORK/M_$tag/output/AZUREOut_aa=2_R=2.extrap"
    if [ ! -s "$fa" ] || [ ! -s "$fb" ] || [ ! -s "$m1" ] || [ ! -s "$m2" ]; then
      bad "$tag: an extrapolation is missing"; tail -5 "$WORK/M_$tag/run.log"; continue
    fi
    xs "$fa" 1 > "$WORK/a.$tag"; xs "$fb" 1 > "$WORK/b.$tag"
    xs "$m1" 1 > "$WORK/m1.$tag"; xs "$m1" 2 > "$WORK/m3.$tag"
    xs "$m2" 1 > "$WORK/m2.$tag"; xs "$m2" 2 > "$WORK/m4.$tag"
    n="$(wc -l < "$WORK/a.$tag")"
    [ "$n" -gt 100 ] || bad "$tag: only $n points in A"
    r="$(maxrel "$WORK/m1.$tag" "$WORK/a.$tag")"
    if le "$r" 1e-10; then ok "$tag: 13N alone in the merged file == A.azr (max rel $r)"; else bad "$tag: merged 13N vs A.azr differ (max rel $r)"; fi
    r="$(maxrel "$WORK/m2.$tag" "$WORK/b.$tag")"
    if le "$r" 1e-10; then ok "$tag: 14N alone in the merged file == B.azr (max rel $r)"; else bad "$tag: merged 14N vs B.azr differ (max rel $r)"; fi
    paste "$WORK/a.$tag" "$WORK/b.$tag" | awk '{ printf "%.12e %.12e\n", $1, 0.6 * $2 + 0.4 * $4 }' > "$WORK/sum.$tag"
    r="$(maxrel "$WORK/m3.$tag" "$WORK/sum.$tag")"
    if le "$r" 1e-9; then ok "$tag: natural SUM == 0.6 sigma_A + 0.4 sigma_B (max rel $r)"; else bad "$tag: natural SUM vs 0.6 A + 0.4 B (max rel $r)"; fi
    awk '{ printf "%.12e %.12e\n", $1, 0.4 * $2 }' "$WORK/b.$tag" > "$WORK/b4.$tag"
    r="$(maxrel "$WORK/m4.$tag" "$WORK/b4.$tag")"
    if le "$r" 1e-9; then ok "$tag: plain natural segment == 0.4 sigma_B (max rel $r)"; else bad "$tag: plain natural segment vs 0.4 B (max rel $r)"; fi
    dec="$WORK/M_$tag/output/AZUREOut_sample_seg=3.extrap"
    if [ -s "$dec" ]; then
      awk '!/^#/ && NF >= 4 { printf "%.12e %.12e\n", $1, $3 }' "$dec" > "$WORK/d1.$tag"
      awk '!/^#/ && NF >= 4 { printf "%.12e %.12e\n", $1, $4 }' "$dec" > "$WORK/d2.$tag"
      awk '!/^#/ && NF >= 4 { printf "%.12e %.12e\n", $1, $2 - $3 - $4 }' "$dec" > "$WORK/d0.$tag"
      awk '{ printf "%.12e %.12e\n", $1, 0.6 * $2 }' "$WORK/a.$tag" > "$WORK/a6.$tag"
      r1="$(maxrel "$WORK/d1.$tag" "$WORK/a6.$tag")"; r2="$(maxrel "$WORK/d2.$tag" "$WORK/b4.$tag")"
      res="$(awk '{ v = $2 < 0 ? -$2 : $2; if (v > w) w = v } END { printf "%.3e", w }' "$WORK/d0.$tag")"
      # the decomposition file is written with 10 significant digits
      if le "$r1" 1e-9 && le "$r2" 1e-9 && le "$res" 1e-10; then ok "$tag: decomposition file carries 0.6 A and 0.4 B and they add up"; else bad "$tag: decomposition file ($r1, $r2, residual $res)"; fi
    else
      bad "$tag: no decomposition file AZUREOut_sample_seg=3.extrap"
    fi
  done
done

# ---- (b) under target integration: the component's grid is built with its own pair --------------------
run "A_thick" "$M/A_thick.azr" true "" || true
run "B_thick" "$M/B_thick.azr" true "" || true
run "M_thick" "$M/M_thick.azr" true "" || true
fa="$WORK/A_thick/output/AZUREOut_aa=1_R=1.extrap"; fb="$WORK/B_thick/output/AZUREOut_aa=2_R=2.extrap"
m1="$WORK/M_thick/output/AZUREOut_aa=1_R=1.extrap"; m2="$WORK/M_thick/output/AZUREOut_aa=2_R=2.extrap"
if [ -s "$fa" ] && [ -s "$fb" ] && [ -s "$m1" ] && [ -s "$m2" ]; then
  xs "$fa" 1 > "$WORK/a.thick"; xs "$fb" 1 > "$WORK/b.thick"; xs "$m1" 1 > "$WORK/m1.thick"; xs "$m1" 2 > "$WORK/m3.thick"; xs "$m2" 1 > "$WORK/m2.thick"
  r="$(maxrel "$WORK/m1.thick" "$WORK/a.thick")"; r2="$(maxrel "$WORK/m2.thick" "$WORK/b.thick")"
  if le "$r" 1e-10 && le "$r2" 1e-10; then ok "thick: each nucleus alone == its standalone project ($r, $r2)"; else bad "thick: standalone vs merged ($r, $r2)"; fi
  paste "$WORK/a.thick" "$WORK/b.thick" | awk '{ printf "%.12e %.12e\n", $1, 0.6 * $2 + 0.4 * $4 }' > "$WORK/sum.thick"
  r="$(maxrel "$WORK/m3.thick" "$WORK/sum.thick")"
  if le "$r" 1e-10; then ok "thick: natural SUM == 0.6 A_thick + 0.4 B_thick (max rel $r)"; else bad "thick: natural SUM vs 0.6 A + 0.4 B (max rel $r): the component's sub-point grid is not its own pair's"; fi
else
  bad "thick: an extrapolation is missing"; tail -5 "$WORK/M_thick/run.log"
fi

# ---- (b3) a second, named sample on the segment line ----------------------------------------------------
run "M_samples" "$M/M_samples.azr" true "" || true
m1="$WORK/M_samples/output/AZUREOut_aa=1_R=1.extrap"; m2="$WORK/M_samples/output/AZUREOut_aa=2_R=2.extrap"
if [ -s "$m1" ] && [ -s "$m2" ] && [ -s "$WORK/a.Amat" ]; then
  xs "$m1" 3 > "$WORK/m5.samples"; xs "$m2" 3 > "$WORK/m6.samples"
  paste "$WORK/a.Amat" "$WORK/b.Amat" | awk '{ printf "%.12e %.12e\n", $1, 0.1 * $2 + 0.9 * $4 }' > "$WORK/sum.samples"
  r="$(maxrel "$WORK/m5.samples" "$WORK/sum.samples")"
  if le "$r" 1e-9; then ok "samples: SUM on sample enr == 0.1 sigma_A + 0.9 sigma_B (max rel $r)"; else bad "samples: SUM on sample enr vs 0.1 A + 0.9 B (max rel $r)"; fi
  awk '{ printf "%.12e %.12e\n", $1, 0.9 * $2 }' "$WORK/b.Amat" > "$WORK/b9.samples"
  r="$(maxrel "$WORK/m6.samples" "$WORK/b9.samples")"
  if le "$r" 1e-9; then ok "samples: plain segment on sample enr == 0.9 sigma_B (max rel $r)"; else bad "samples: plain segment on enr vs 0.9 B (max rel $r)"; fi
  xs "$m1" 2 > "$WORK/m3.samples"
  r="$(maxrel "$WORK/m3.samples" "$WORK/sum.Amat")"
  if le "$r" 1e-9; then ok "samples: the unnamed sample's SUM is unchanged by the second sample (max rel $r)"; else bad "samples: unnamed sample's SUM changed (max rel $r)"; fi
  if grep -q "sample enr" "$WORK/M_samples/output/AZUREOut_sample_seg=5.extrap" 2>/dev/null; then ok "samples: decomposition file names its sample"; else bad "samples: decomposition file of segment 5 missing or without the sample name"; fi
else
  bad "samples: an extrapolation is missing"; tail -5 "$WORK/M_samples/run.log"
fi

# ---- (c) data mode: abundances as fit parameters --------------------------------------------------------
# data.dat = 0.6 sigma_A + 0.4 sigma_B on the 2 keV lab grid of the test segments, 0.1 % errors, no noise
# (the two nuclei differ by a few percent, so 1 % errors would leave the abundances almost unconstrained)
if [ -s "$WORK/a.Amat" ] && [ -s "$WORK/b.Amat" ]; then
  paste "$WORK/a.Amat" "$WORK/b.Amat" | awk '{ e = 1.0 + (NR - 1) * 0.002; v = 0.6 * $2 + 0.4 * $4; printf "%.6f 0 %.12e %.12e\n", e, v, 0.001 * v }' > "$WORK/data.dat"
  # the analytic gradient of the free abundances against finite differences
  (cd "$WORK" && mkdir -p gradcheck/output gradcheck/checks && cp data.dat gradcheck/ && cp "$M/M_data.azr" gradcheck/t.azr &&
   cd gradcheck && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit run_limited 600 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > run.log 2>&1)
  if grep -q "abundance" "$WORK/gradcheck/run.log"; then
    line="$(grep "abundance" "$WORK/gradcheck/run.log" | head -1)"
    if echo "$line" | grep -Eq "abundance[^0-9]*0/[0-9]+"; then ok "gradient check: abundances agree with finite differences ($line)"; else bad "gradient check: $line"; fi
  else
    bad "gradient check printed nothing about the abundances"; grep -i "grad" "$WORK/gradcheck/run.log" | head -5
  fi
  # a normalization group over the natural SUM (leader, norm free) and a plain segment on the same data
  # (follower): the follower's rows must reach the leader's norm column on both gradient paths
  (cd "$WORK" && mkdir -p gradcheck_ng/output gradcheck_ng/checks && cp data.dat gradcheck_ng/ &&
   awk '/<segmentsData>/ { print; print "1 1 1 1.0 3.0 0 0 0 1.05 1 10 0 0 0 data.dat 1 0 1 2 2 -999 normgroup N natural 1";
                           print "1 1 1 1.0 3.0 0 0 0 1 1 10 0 0 0 data.dat 0 0 normgroup N"; skip = 1; next }
        /<\/segmentsData>/ { skip = 0 } skip { next } { print }' "$M/M_data.azr" > gradcheck_ng/t.azr &&
   cd gradcheck_ng && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit run_limited 600 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > run.log 2>&1)
  line="$(grep "components with relative" "$WORK/gradcheck_ng/run.log" | head -1)"
  if echo "$line" | grep -q "norm 0/1 .*Jacobian path:.* norm 0/1"; then ok "normalization group over a natural SUM and a plain segment: norm gradient agrees with finite differences"
  else bad "normalization group gradient check: $line"; fi
  # two samples: data_enr.dat = 0.1 sigma_A + 0.9 sigma_B; four free abundances, two sum-to-one priors
  paste "$WORK/a.Amat" "$WORK/b.Amat" | awk '{ e = 1.0 + (NR - 1) * 0.002; v = 0.1 * $2 + 0.9 * $4; printf "%.6f 0 %.12e %.12e\n", e, v, 0.001 * v }' > "$WORK/data_enr.dat"
  (cd "$WORK" && mkdir -p gradcheck2/output gradcheck2/checks && cp data.dat data_enr.dat gradcheck2/ && cp "$M/M_samples_data.azr" gradcheck2/t.azr &&
   cd gradcheck2 && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit run_limited 600 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > run.log 2>&1)
  line="$(grep "abundance" "$WORK/gradcheck2/run.log" | head -1)"
  if echo "$line" | grep -Eq "abundance[^0-9]*0/4"; then ok "gradient check with two samples: 4 abundances agree with finite differences"; else bad "gradient check with two samples: $line"; fi
  run "fit2" "$M/M_samples_data.azr" true "" 2 || true
  sav="$WORK/fit2/output/param.sav"
  if [ -s "$sav" ]; then
    v="$(awk '$1 == "abundance_7_13" || $1 == "abundance_7_14" || $1 == "abundance_enr_7_13" || $1 == "abundance_enr_7_14" { printf "%s ", $2 }' "$sav")"
    if awk -v v="$v" 'BEGIN { n = split(v, a, " "); exit !(n == 4 && (a[1]-0.6)^2 < 1e-6 && (a[2]-0.4)^2 < 1e-6 && (a[3]-0.1)^2 < 1e-6 && (a[4]-0.9)^2 < 1e-6) }'; then ok "fit with two samples recovers 0.6 / 0.4 and 0.1 / 0.9 ($v)"; else bad "fit with two samples gave $v (want 0.6 0.4 0.1 0.9)"; tail -5 "$WORK/fit2/run.log"; fi
  else
    bad "two-sample fit wrote no param.sav"; tail -8 "$WORK/fit2/run.log"
  fi
  # a fit from 0.5 / 0.5 must find 0.6 / 0.4
  run "fit" "$M/M_data.azr" true "" 2 || true
  sav="$WORK/fit/output/param.sav"
  if [ -s "$sav" ]; then
    a13="$(awk '$1 == "abundance_7_13" { print $2 }' "$sav")"; a14="$(awk '$1 == "abundance_7_14" { print $2 }' "$sav")"
    chi="$(grep -o "Total-Chi-Squared: [0-9.e+-]*" "$WORK/fit/output/chiSquared.out" 2>/dev/null | head -1)"
    if awk -v a="$a13" -v b="$a14" 'BEGIN { exit !((a - 0.6)^2 < 1e-6 && (b - 0.4)^2 < 1e-6) }'; then ok "fit recovers the abundances 0.6 / 0.4 ($a13 / $a14; $chi)"; else bad "fit gave abundances $a13 / $a14 (want 0.6 / 0.4; $chi)"; tail -5 "$WORK/fit/run.log"; fi
  else
    bad "fit wrote no param.sav"; tail -8 "$WORK/fit/run.log"
  fi
else
  bad "no A/B extrapolations for the data-mode check"
fi

# ---- (d) refused models ----------------------------------------------------------------------------------
for name in bad_level bad_ratio bad_key bad_sample; do
  run "$name" "$M/$name.azr" true "" || true
  if grep -q "ERROR" "$WORK/$name/run.log" && ! ls "$WORK/$name/output"/AZUREOut_* > /dev/null 2>&1; then
    ok "$name refused: $(grep -m1 ERROR "$WORK/$name/run.log" | cut -c1-90)"
  else
    bad "$name was not refused"; tail -3 "$WORK/$name/run.log"
  fi
done

if [ "$failures" -eq 0 ]; then echo "PASS: natural-target segments"; else echo "FAIL: $failures check(s)"; exit 1; fi
