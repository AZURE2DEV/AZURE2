#!/usr/bin/env bash
#
# Prior parameter covariance ("priorCovariance=<file>" in <config>): a multivariate Gaussian prior
# (p - mu)^T C^-1 (p - mu) on the parameters the file lists (names, means, covariance), on their values as the
# parameter vector holds them (the param.sav space).  Model: the p+p model of tests/norm_groups (five levels, the
# exact cross sections at 40 and 90 degrees times 0.8 / 0.9, 1 % errors); every project is built here.
#
#   (a) a DIAGONAL prior on two norms reproduces the norms' own normError priors (switched off): the same data
#       chi2, Total-Prior-Covariance-Chi-Squared = the normError terms exactly (1e-10; Total-Norm-Chi-Squared is
#       written with 6 digits), at moved norms (mode 1) and after an --use-lm fit (same norms, same chi2);
#   (b) a CORRELATED 2x2 prior on a norm and a reduced width: the term at moved parameters = the hand-computed
#       quadratic form d^T C^-1 d (1e-10);
#   (c) AZURE_GRAD_CHECK (--use-gradient) with every level parameter free and a correlated 3x3 prior on a norm, a
#       width and an energy shift (the shift's gradient takes the segment-local finite-difference path): analytic
#       = finite differences (1e-6 of the largest), the Jacobian path = the adjoint gradient (1e-10 of the largest,
#       the finite-differenced shift column to 1e-8 of its own), and the gradient with the prior minus the gradient
#       without it = 2 C^-1 (p - mu) on the three parameters (1e-7) and 0 elsewhere (adjoint and Jacobian paths);
#   (d) a norm the data do not constrain (its segment's errors x 1e8): a fit returns the prior mean; with a
#       correlated prior (rho 0.95) and the other norm pinned by the data, it lands on the conditional mean
#       mu_B + C_AB/C_AA (A - mu_A), 0.024 away from mu_B -- to 1e-4, MIGRAD (--use-gradient) and --use-lm;
#   (e) refused with an ERROR and no output: an unknown name, a non-symmetric matrix, a non-positive-definite one,
#       names/means and row-count and row-length mismatches, a #parametrization tag of another basis, a missing
#       file; a matching tag is accepted; a listed parameter fixed in the run is marginalized with a Note (the term
#       is the 1-parameter marginal (n - mu)^2 / C_nn, 1e-10);
#   (f) pyazr.prior_from_fit: an --use-lm fit with --covariance-band (covariance.dat with its #columns header) ->
#       the prior file (names = the header, means = param.sav, the basis tag) -> a second fit from the .azr's start
#       with the prior and the same data returns the first fit's widths to 1e-6 (--use-lm; MIGRAD to 1e-2 of each
#       prior sigma), the same data chi2 and a prior term ~ 0;
#   (g) MCMC (RWA): the logprior column of samples.mcmc is the Gaussian log density of the unconstrained norm
#       (1e-10), which REPLACES its normError auto prior; the chain's mean and spread are the prior's; a prior on a
#       level parameter is refused in a physical-parameter chain.
#
#   ./tests/prior_covariance/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PY="$(command -v python3 || command -v python)"
BASE="$ROOT/tests/norm_groups/norm_groups.azr"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/prior_covariance.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
# Every comparison insists on numbers: an empty field (a crashed run, a missing line) or nan is a FAIL, never 0.
num='^[-+]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$'
isnum() { awk -v a="$1" -v re="$num" 'BEGIN { exit !(a ~ re) }'; }
le() { awk -v a="$1" -v b="$2" -v re="$num" 'BEGIN { exit !(a ~ re && a + 0 <= b + 0) }'; }
gt() { awk -v a="$1" -v b="$2" -v re="$num" 'BEGIN { exit !(a ~ re && a + 0 > b + 0) }'; }
# relclose A B TOL: both numbers and |A-B| <= TOL*max(|B|,1e-300)
relclose() { awk -v a="$1" -v b="$2" -v t="$3" -v re="$num" 'BEGIN { if (!(a ~ re && b ~ re)) exit 1
  d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; if (m == 0) m = 1e-300; exit !(d <= t * m) }'; }
absclose() { awk -v a="$1" -v b="$2" -v t="$3" -v re="$num" 'BEGIN { if (!(a ~ re && b ~ re)) exit 1
  d = a - b; if (d < 0) d = -d; exit !(d <= t) }'; }

# mk NAME LEVELS PRIOR SEGMENT_LINE... -> $WORK/NAME/t.azr (+ data/ with xs_90_big.dat, the 90-degree file with
# its errors x 1e8).  LEVELS: fixall (every energy and width fixed), fixE (energies fixed), free.  PRIOR: the
# priorCovariance value, or "" for none.
mk() {
  local name="$1" levels="$2" prior="$3"; shift 3
  local dir="$WORK/$name"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$ROOT/tests/norm_groups/data" "$dir/data"
  awk '{ $4 = $4 * 1e8; print }' "$dir/data/xs_90_s09.dat" > "$dir/data/xs_90_big.dat"
  printf '%s\n' "$@" > "$dir/segs.txt"
  awk -v lv="$levels" -v pc="$prior" -v segs="$dir/segs.txt" '
    /<levels>/ { inlev = 1; print; next }
    /<\/levels>/ { inlev = 0 }
    inlev && NF > 12 && lv == "fixall" { $4 = 1; $11 = 1; print; next }
    inlev && NF > 12 && lv == "fixE" { $4 = 1; print; next }
    /<\/config>/ && pc != "" { print "priorCovariance=" pc }
    /<segmentsData>/ { print; while ((getline l < segs) > 0) print l; skip = 1; next }
    /<\/segmentsData>/ { skip = 0 }
    skip { next }
    { print }' "$BASE" > "$dir/t.azr"
}
# run NAME INPUT FLAGS: the CLI in $WORK/NAME, log in run.log; the exit status is kept in $WORK/NAME/status
run() {
  (cd "$WORK/$1" && printf "$2" | run_limited 600 "$AZURE2_BIN" $3 --no-gui --no-readline t.azr > run.log 2>&1
   echo $? > status)
}
need() {  # NAME: chiSquared.out must exist
  if [ ! -s "$WORK/$1/output/chiSquared.out" ]; then bad "$1: no chiSquared.out"; tail -8 "$WORK/$1/run.log"; return 1; fi
}
chi() { awk -v k="$2:" '{ for (i = 1; i < NF; i++) if ($i == k) print $(i + 1) }' "$WORK/$1/output/chiSquared.out" 2>/dev/null; }
par() { awk -v n="$2" '$1 == n { print $2 }' "$WORK/$1/output/param.sav" 2>/dev/null; }
# setpar IN OUT name=value...: a parameter file with those values (the others as IN has them)
setpar() {
  local in="$1" out="$2"; shift 2
  awk -v kv="$*" 'BEGIN { n = split(kv, a, " "); for (i = 1; i <= n; i++) { split(a[i], b, "="); v[b[1]] = b[2] } }
       ($1 in v) { $2 = v[$1] } { print }' "$in" > "$out"
}

SEG1='1 1 1 0 100 40 40 4 1 1 10 0 0 0 data/xs_40_s08.dat 0 0'
SEG2='1 1 1 0 100 90 90 4 1 1 5 0 0 0 data/xs_90_s09.dat 0 0'
SEG1_0='1 1 1 0 100 40 40 4 1 1 0 0 0 0 data/xs_40_s08.dat 0 0'
SEG2_0='1 1 1 0 100 90 90 4 1 1 0 0 0 0 data/xs_90_s09.dat 0 0'

# ---- (a) a diagonal prior = the normError priors -------------------------------------------------------------------
mk a_norm fixall "" "$SEG1" "$SEG2"
mk a_prior fixall prior.txt "$SEG1_0" "$SEG2_0"
printf '# the normError priors of the two segments: 10 %% and 5 %% of 1\nsegment_1_norm segment_2_norm\n1 1\n0.01 0\n0 0.0025\n' \
  > "$WORK/a_prior/prior.txt"
run a_norm '1\n\n\n7\n' ""
if need a_norm; then
  setpar "$WORK/a_norm/output/param.par" "$WORK/moved.par" segment_1_norm=1.07 segment_2_norm=0.96
  for d in a_norm a_prior; do cp "$WORK/moved.par" "$WORK/$d/moved.par"; run $d '1\nmoved.par\n\n7\n' ""; done
  if need a_norm && need a_prior; then
    want=1.13   # ((1.07 - 1)/0.1)^2 + ((0.96 - 1)/0.05)^2 = 0.49 + 0.64
    c1="$(chi a_norm Total-Chi-Squared)"; c2="$(chi a_prior Total-Chi-Squared)"
    n1="$(chi a_norm Total-Norm-Chi-Squared)"; p2="$(chi a_prior Total-Prior-Covariance-Chi-Squared)"
    if isnum "$c1" && [ "$c1" = "$c2" ]; then ok "(a) same data chi2 ($c1)"; else bad "(a) data chi2 $c1 vs $c2"; fi
    if relclose "$p2" "$want" 1e-10 && relclose "$n1" "$p2" 1e-5; then
      ok "(a) Total-Prior-Covariance-Chi-Squared $p2 = the normError terms $want (1e-10; Total-Norm-Chi-Squared $n1)"
    else bad "(a) prior term $p2, normError terms $n1, exact $want"; fi
    if grep -q "Note: prior covariance prior.txt: a multivariate Gaussian prior on 2 parameters" "$WORK/a_prior/run.log"
    then ok "(a) the prior is noted in the log"; else bad "(a) no prior Note in the log"; fi
    if [ -z "$(chi a_norm Total-Prior-Covariance-Chi-Squared)" ]; then ok "(a) no prior field without the key"
    else bad "(a) a prior field without the key"; fi
  fi
  # the same through a fit (both norms free, levels fixed): one minimum
  mk a_norm_fit fixall "" "$SEG1" "$SEG2"; mk a_prior_fit fixall prior.txt "$SEG1_0" "$SEG2_0"
  cp "$WORK/a_prior/prior.txt" "$WORK/a_prior_fit/"
  run a_norm_fit '2\n\n\n\n7\n' "--use-lm"; run a_prior_fit '2\n\n\n\n7\n' "--use-lm"
  if need a_norm_fit && need a_prior_fit; then
    f1="$(par a_norm_fit segment_1_norm) $(par a_norm_fit segment_2_norm)"; f2="$(par a_prior_fit segment_1_norm) $(par a_prior_fit segment_2_norm)"
    t1="$(chi a_norm_fit Total-Chi-Squared)"; t2="$(chi a_prior_fit Total-Chi-Squared)"
    n1="$(chi a_norm_fit Total-Norm-Chi-Squared)"; p2="$(chi a_prior_fit Total-Prior-Covariance-Chi-Squared)"
    set -- $f1 $f2
    if [ $# -eq 4 ] && relclose "$1" "$3" 1e-6 && relclose "$2" "$4" 1e-6 && [ "$t1" = "$t2" ] && relclose "$n1" "$p2" 1e-5
    then ok "(a) --use-lm fit: the same norms ($f1 / $f2), data chi2 $t1, penalty $n1 = prior term $p2"
    else bad "(a) --use-lm fit: norms $f1 vs $f2, chi2 $t1 vs $t2, penalty $n1 vs prior $p2"; fi
  fi
fi

# ---- (b) a correlated 2x2 prior: the quadratic form -----------------------------------------------------------------
mk b free prior.txt "$SEG1_0" "$SEG2_0"
printf 'segment_1_norm width_1_1\n1.02 0.6\n0.0004 0.00012\n0.00012 0.0009\n' > "$WORK/b/prior.txt"
run b '1\n\n\n7\n' ""
if need b; then
  setpar "$WORK/b/output/param.par" "$WORK/b/moved.par" segment_1_norm=1.07 width_1_1=0.63
  run b '1\nmoved.par\n\n7\n' ""
  if need b; then
    p="$(chi b Total-Prior-Covariance-Chi-Squared)"
    want="$(awk 'BEGIN { a = 0.0004; c = 0.00012; e = 0.0009; x = 1.07 - 1.02; y = 0.63 - 0.6; det = a * e - c * c
                         printf "%.17g", (e * x * x - 2 * c * x * y + a * y * y) / det }')"
    if relclose "$p" "$want" 1e-10; then ok "(b) correlated 2x2: term $p = d^T C^-1 d $want (1e-10)"
    else bad "(b) correlated 2x2: term $p, expected $want"; fi
  fi
fi

# ---- (c) gradient check --------------------------------------------------------------------------------------------
SEG1_SHIFT='1 1 1 0 100 40 40 4 1 1 10 0.002 0.01 1 data/xs_40_s08.dat 0 0'
SEG2_10='1 1 1 0 100 90 90 4 1 1 10 0 0 0 data/xs_90_s09.dat 0 0'
mk c_prior free prior.txt "$SEG1_SHIFT" "$SEG2_10"
mk c_none free "" "$SEG1_SHIFT" "$SEG2_10"
cat > "$WORK/c_prior/prior.txt" <<'EOF'
segment_1_norm width_2_1 segment_1_energy_shift
1.05 0.7 0.004
1e-4 2e-5 1e-7
2e-5 1e-3 -2e-6
1e-7 -2e-6 1e-5
EOF
run c_none '1\n\n\n7\n' ""
if [ -s "$WORK/c_none/output/param.par" ]; then
  setpar "$WORK/c_none/output/param.par" "$WORK/c_start.par" segment_1_norm=1 width_2_1=0.75 segment_1_energy_shift=0.002
  for d in c_prior c_none; do
    cp "$WORK/c_start.par" "$WORK/$d/start.par"
    (cd "$WORK/$d" && printf '2\n\nstart.par\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
       run_limited 600 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > run.log 2>&1)
    awk '{ sub(/\r$/, "") }
         /^\[grad-check\] idx / { a = f = j = ""
           for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
           print $3, a, f, j, $4 }' "$WORK/$d/run.log" | sort -n > "$WORK/$d/rows"
  done
  nrow="$(awk 'END { print NR }' "$WORK/c_prior/rows")"
  if [ "$nrow" -ne 14 ] || [ "$(awk 'END { print NR }' "$WORK/c_none/rows")" -ne 14 ]; then
    bad "(c) gradient check listed $nrow parameters (expected 14: 5 energies, 6 widths, 2 norms, 1 shift)"; tail -5 "$WORK/c_prior/run.log"
  else
    read -r fd jac jsh <<< "$(awk '{ a = $2 < 0 ? -$2 : $2; if (a > m) m = a
                                     d = $2 - $3; if (d < 0) d = -d; if (d > w) w = d
                                     e = $2 - $4; if (e < 0) e = -e
                                     if ($5 == "shift") { s = e / a } else if (e > v) v = e }
                                   END { printf "%.2e %.2e %.2e", w / m, v / m, s }' "$WORK/c_prior/rows")"
    le "$fd" 1e-6 && ok "(c) with the prior: gradient = finite differences for 14 parameters ($fd of the largest)" \
      || bad "(c) with the prior: gradient vs finite differences $fd of the largest"
    if le "$jac" 1e-10 && le "$jsh" 1e-8; then ok "(c) Jacobian path = adjoint ($jac of the largest; the shift column $jsh of its own)"
    else bad "(c) Jacobian path vs adjoint: $jac of the largest, shift column $jsh"; fi
    # the prior's own gradient: (with - without) = 2 C^-1 (p - mu) at the start (indices: norm 11, width_2_1 3, shift 13)
    read -r dg djac drest <<< "$("$PY" -I - "$WORK/c_prior/rows" "$WORK/c_none/rows" <<'PYEOF'
import sys
import numpy as np
w = {int(r.split()[0]): [float(v) for v in r.split()[1:4]] for r in open(sys.argv[1])}
o = {int(r.split()[0]): [float(v) for v in r.split()[1:4]] for r in open(sys.argv[2])}
C = np.array([[1e-4, 2e-5, 1e-7], [2e-5, 1e-3, -2e-6], [1e-7, -2e-6, 1e-5]])
d = np.array([1.0 - 1.05, 0.75 - 0.7, 0.002 - 0.004])
g = 2 * np.linalg.solve(C, d)
idx = [11, 3, 13]
big = max(abs(v[0]) for v in w.values())
dg = max(abs((w[i][0] - o[i][0]) - g[k]) / abs(g).max() for k, i in enumerate(idx))
dj = max(abs((w[i][2] - o[i][2]) - g[k]) / abs(g).max() for k, i in enumerate(idx))
rest = max(max(abs(w[i][0] - o[i][0]), abs(w[i][2] - o[i][2])) for i in w if i not in idx) / big
print(f"{dg:.2e} {dj:.2e} {rest:.2e}")
PYEOF
)"
    if le "$dg" 1e-7 && le "$djac" 1e-7 && le "$drest" 1e-9; then
      ok "(c) gradient with minus without the prior = 2 C^-1 (p - mu): adjoint $dg, Jacobian path $djac; elsewhere 0 ($drest)"
    else bad "(c) prior gradient: adjoint $dg, Jacobian path $djac (of the largest prior component), elsewhere $drest"; fi
  fi
else
  bad "(c) no param.par for the start"
fi

# ---- (d) an unconstrained parameter: the prior mean / the conditional mean -----------------------------------------
SEGA='1 1 1 0 100 40 40 4 1 1 0 0 0 0 data/xs_40_s08.dat 0 0'
SEGB='1 1 1 0 100 90 90 4 1 1 0 0 0 0 data/xs_90_big.dat 0 0'
for m in "--use-gradient" "--use-lm"; do
  tag="${m#--use-}"
  mk "d1_$tag" fixall prior.txt "$SEGA" "$SEGB"
  printf 'segment_2_norm\n0.9\n0.0004\n' > "$WORK/d1_$tag/prior.txt"
  run "d1_$tag" '2\n\n\n\n7\n' "$m"
  if need "d1_$tag"; then
    b="$(par "d1_$tag" segment_2_norm)"; a="$(par "d1_$tag" segment_1_norm)"
    if absclose "$b" 0.9 1e-4 && absclose "$a" 1.25 1e-3; then ok "(d) $m: the unconstrained norm at the prior mean ($b; 0.9), the other at the data's ($a)"
    else bad "(d) $m: unconstrained norm $b (prior mean 0.9), other $a (data 1.25)"; fi
  fi
  mk "d2_$tag" fixall prior.txt "$SEGA" "$SEGB"
  printf 'segment_1_norm segment_2_norm\n1.2 0.9\n0.01 0.00475\n0.00475 0.0025\n' > "$WORK/d2_$tag/prior.txt"
  run "d2_$tag" '2\n\n\n\n7\n' "$m"
  if need "d2_$tag"; then
    b="$(par "d2_$tag" segment_2_norm)"; a="$(par "d2_$tag" segment_1_norm)"
    want="$(awk -v a="$a" 'BEGIN { printf "%.10g", 0.9 + 0.00475 / 0.01 * (a - 1.2) }')"
    if absclose "$b" "$want" 1e-4 && gt "$(awk -v b="$b" 'BEGIN { print b - 0.9 }')" 0.01; then
      ok "(d) $m: correlated prior moves it to mu_B + C_AB/C_AA (A - mu_A) = $want ($b; A = $a)"
    else bad "(d) $m: correlated prior: $b, expected $want (A = $a)"; fi
  fi
done

# ---- (e) refusals; a fixed parameter is marginalized ---------------------------------------------------------------
refuse() {  # NAME WHAT PRIOR-TEXT...: refused with an ERROR, nothing written
  local name="$1" what="$2"; shift 2
  mk "$name" fixall prior.txt "$SEG1_0" "$SEG2_0"
  printf '%s\n' "$@" > "$WORK/$name/prior.txt"
  run "$name" '1\n\n\n7\n' ""
  local out; out="$(ls "$WORK/$name/output")"
  if grep -q '^\*\*ERROR: priorCovariance' "$WORK/$name/run.log" && [ -z "$out" ] && [ "$(cat "$WORK/$name/status")" != 0 ]; then
    ok "(e) $what: refused ($(grep -m1 '^\*\*ERROR: priorCovariance' "$WORK/$name/run.log" | cut -c1-110)...), no output"
  else bad "(e) $what: not refused cleanly (output: '${out}', status $(cat "$WORK/$name/status"))"; tail -3 "$WORK/$name/run.log"; fi
}
refuse e_unknown "unknown name" 'segment_1_norm segment_9_norm' '1 1' '0.01 0' '0 0.01'
refuse e_nonsym "non-symmetric" 'segment_1_norm segment_2_norm' '1 1' '0.01 0.001' '0.0011 0.01'
refuse e_nonpd "not positive definite" 'segment_1_norm segment_2_norm' '1 1' '0.01 0.02' '0.02 0.01'
refuse e_means "names/means mismatch" 'segment_1_norm segment_2_norm' '1' '0.01 0' '0 0.01'
refuse e_rows "too many rows" 'segment_1_norm segment_2_norm' '1 1' '0.01 0' '0 0.01' '0 0'
refuse e_rowlen "row length" 'segment_1_norm segment_2_norm' '1 1' '0.01 0 0' '0 0.01'
refuse e_tag "#parametrization 2 in a Brune run" '#parametrization 2' 'segment_1_norm segment_2_norm' '1 1' '0.01 0' '0 0.01'
mk e_missing fixall nothere.txt "$SEG1_0" "$SEG2_0"
run e_missing '1\n\n\n7\n' ""
if grep -q '^\*\*ERROR: priorCovariance: cannot open' "$WORK/e_missing/run.log" && [ -z "$(ls "$WORK/e_missing/output")" ]; then
  ok "(e) missing file: refused, no output"; else bad "(e) missing file not refused cleanly"; fi
mk e_tagok fixall prior.txt "$SEG1_0" "$SEG2_0"
printf '#parametrization 1.000000e+00 0.000000e+00\nsegment_1_norm segment_2_norm\n1 1\n0.01 0\n0 0.01\n' > "$WORK/e_tagok/prior.txt"
run e_tagok '1\n\n\n7\n' ""
if need e_tagok && isnum "$(chi e_tagok Total-Prior-Covariance-Chi-Squared)"; then ok "(e) the run's own tag (param.sav form) is accepted"
else bad "(e) the run's own tag was not accepted"; fi
# a listed parameter fixed in the run (energy_2: every level fixed) -> marginalized
mk e_fixed fixall prior.txt "$SEG1_0" "$SEG2_0"
printf 'segment_1_norm energy_2\n1.02 0.81\n0.0004 0.00015\n0.00015 0.0001\n' > "$WORK/e_fixed/prior.txt"
run e_fixed '1\n\n\n7\n' ""
if need e_fixed; then
  setpar "$WORK/e_fixed/output/param.par" "$WORK/e_fixed/moved.par" segment_1_norm=1.07
  run e_fixed '1\nmoved.par\n\n7\n' ""
  p="$(chi e_fixed Total-Prior-Covariance-Chi-Squared)"
  want="$(awk 'BEGIN { printf "%.17g", (1.07 - 1.02)^2 / 0.0004 }')"
  if relclose "$p" "$want" 1e-10 && grep -q "Note: prior covariance prior.txt: energy_2 is fixed in this run and marginalized" "$WORK/e_fixed/run.log"; then
    ok "(e) a fixed listed parameter: Note, marginalized; term $p = (n - mu)^2/C_nn (1e-10)"
  else bad "(e) fixed listed parameter: term $p (marginal $want), Note: $(grep -c 'marginalized' "$WORK/e_fixed/run.log")"; fi
fi

# ---- (f) prior_from_fit round trip ---------------------------------------------------------------------------------
SEGF1='1 1 1 0 100 40 40 4 1.25 0 0 0 0 0 data/xs_40_s08.dat 0 0'
SEGF2='1 1 1 0 100 90 90 4 1.25 0 0 0 0 0 data/xs_90_s09.dat 0 0'
mk f_a fixE "" "$SEGF1" "$SEGF2"
run f_a '2\n\n\n7\n' "--use-lm --covariance-band"
if need f_a && [ -s "$WORK/f_a/output/covariance.dat" ]; then
  "$PY" -I - "$ROOT/pyazr/priors.py" "$WORK/f_a/output" "$WORK/f_prior.txt" > "$WORK/f_helper.log" 2>&1 <<'PYEOF'
import importlib.util, sys
spec = importlib.util.spec_from_file_location("priors", sys.argv[1])
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
m.prior_from_fit(sys.argv[2] + "/param.sav", sys.argv[2] + "/covariance.dat", sys.argv[3])
PYEOF
  cols="$(awk '$1 == "#columns" { $1 = ""; sub(/^ /, ""); print }' "$WORK/f_a/output/covariance.dat")"
  names="$(grep -v '^#' "$WORK/f_prior.txt" | head -1)"
  means_ok="$(awk -v sav="$WORK/f_a/output/param.sav" 'BEGIN { while ((getline l < sav) > 0) { split(l, t, " "); v[t[1]] = t[2] } }
     /^#/ { next } { n++ } n == 1 { for (i = 1; i <= NF; i++) nm[i] = $i; k = NF }
     n == 2 { bad = 0; for (i = 1; i <= k; i++) { d = $i - v[nm[i]]; if (d < 0) d = -d; if (d > 1e-15 * (v[nm[i]] < 0 ? -v[nm[i]] : v[nm[i]])) bad = 1 } print bad ? "no" : "yes" }' "$WORK/f_prior.txt")"
  if [ "$cols" = "width_1_1 width_2_1 width_3_1 width_4_1 width_4_2 width_5_1" ] && [ "$names" = "$cols" ] && [ "$means_ok" = yes ] &&
     [ "$(head -1 "$WORK/f_prior.txt")" = "#parametrization 1" ]; then
    ok "(f) prior_from_fit: names = covariance.dat's #columns ($cols), means = param.sav, tag #parametrization 1"
  else bad "(f) prior_from_fit file: columns '$cols', names '$names', means $means_ok, first line '$(head -1 "$WORK/f_prior.txt")'"; cat "$WORK/f_helper.log"; fi
  for m in "--use-lm" "--use-gradient"; do
    tag="${m#--use-}"
    mk "f_b_$tag" fixE prior.txt "$SEGF1" "$SEGF2"
    cp "$WORK/f_prior.txt" "$WORK/f_b_$tag/prior.txt"
    run "f_b_$tag" '2\n\n\n\n7\n' "$m"
    need "f_b_$tag" || continue
    # worst |w_B - w_A| relative to w_A, and in units of the prior sigma
    read -r wrel wsig <<< "$("$PY" -I - "$WORK/f_a/output/param.sav" "$WORK/f_b_$tag/output/param.sav" "$WORK/f_prior.txt" <<'PYEOF'
import sys
import numpy as np
def sav(p):
    return {t[0]: float(t[1]) for t in (l.split() for l in open(p)) if len(t) > 1 and not t[0].startswith("#")}
a, b = sav(sys.argv[1]), sav(sys.argv[2])
rows = [l.split() for l in open(sys.argv[3]) if l.strip() and not l.startswith("#")]
names, C = rows[0], np.array([[float(v) for v in r] for r in rows[2:]])
rel = max(abs(b[n] - a[n]) / abs(a[n]) for n in names)
sig = max(abs(b[n] - a[n]) / np.sqrt(C[i, i]) for i, n in enumerate(names))
print(f"{rel:.2e} {sig:.2e}")
PYEOF
)"
    ca="$(chi f_a Total-Chi-Squared)"; cb="$(chi "f_b_$tag" Total-Chi-Squared)"; pb="$(chi "f_b_$tag" Total-Prior-Covariance-Chi-Squared)"
    if [ "$m" = "--use-lm" ]; then good=0; le "$wrel" 1e-6 && good=1
    else good=0; le "$wsig" 1e-2 && good=1; fi
    if [ $good = 1 ] && isnum "$ca" && [ "$ca" = "$cb" ] && le "$pb" 1e-3; then
      ok "(f) $m from the .azr start with the prior: the first fit's widths (worst $wrel relative, $wsig sigma), chi2 $cb = $ca, prior term $pb"
    else bad "(f) $m: widths off by $wrel relative / $wsig sigma, chi2 $cb vs $ca, prior term $pb"; fi
  done
else
  bad "(f) the first fit wrote no covariance.dat"
fi

# ---- (g) MCMC --------------------------------------------------------------------------------------------------------
SEGB_10='1 1 1 0 100 90 90 4 1 1 10 0 0 0 data/xs_90_big.dat 0 0'
mk g_mcmc fixall prior.txt "$SEGA" "$SEGB_10"
printf 'segment_2_norm\n0.9\n0.0004\n' > "$WORK/g_mcmc/prior.txt"
run g_mcmc '6\n\n\n16\n400\n1\n\n2\nyes\nyes\n' ""
if [ -s "$WORK/g_mcmc/output/samples.mcmc" ]; then
  read -r lpd mean sd <<< "$("$PY" -I - "$WORK/g_mcmc/output/samples.mcmc" <<'PYEOF'
import sys
import numpy as np
d = np.genfromtxt(sys.argv[1], delimiter=",", skip_header=1)
x = d[:, 6]
lp = -0.5 * ((x - 0.9) / 0.02) ** 2 - np.log(0.02 * np.sqrt(2 * np.pi))
s = d[d[:, 0] > 150, 6]
print(f"{np.max(np.abs(lp - d[:, 4])):.2e} {s.mean():.6f} {s.std():.6f}")
PYEOF
)"
  if le "$lpd" 1e-10 && grep -q "replaces the diagonal prior of segment_2_norm" "$WORK/g_mcmc/run.log"; then
    ok "(g) MCMC logprior = the prior's Gaussian log density ($lpd), replacing the normError auto prior (noted)"
  else bad "(g) MCMC logprior off by $lpd, replacement note $(grep -c 'replaces the diagonal' "$WORK/g_mcmc/run.log")"; fi
  if absclose "$mean" 0.9 0.006 && absclose "$sd" 0.02 0.005; then ok "(g) the unconstrained norm's chain: mean $mean, sd $sd (prior 0.9 +- 0.02)"
  else bad "(g) chain mean $mean, sd $sd (prior 0.9 +- 0.02)"; fi
else bad "(g) no samples.mcmc"; tail -5 "$WORK/g_mcmc/run.log"; fi
mk g_phys free prior.txt "$SEGA" "$SEGB_10"
printf 'width_1_1\n0.6\n0.0004\n' > "$WORK/g_phys/prior.txt"
run g_phys '6\n\n\n8\n100\n1\n\n1\nno\nyes\n' ""
if grep -q '^\*\*ERROR: priorCovariance: width_1_1 is a level parameter' "$WORK/g_phys/run.log" && [ ! -s "$WORK/g_phys/output/samples.mcmc" ]; then
  ok "(g) a level-parameter prior is refused in a physical-parameter chain"
else bad "(g) physical-parameter chain with a width prior not refused"; fi

echo
if [ "$failures" -ne 0 ]; then echo "prior_covariance: $failures failure(s)"; exit 1; fi
echo "prior_covariance: all checks passed"
