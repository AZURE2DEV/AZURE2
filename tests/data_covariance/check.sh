#!/usr/bin/env bash
#
# Data covariance matrices ("covariance <file>" on segment lines), on the p+12C/p+13C project of
# tests/angle_mapping (models written by make_models.py; see README.md):
#
#   (a) a diagonal covariance diag(e_i^2) gives the chi2 of no covariance, segment by segment;
#   (b) dense covariance (blocks across files, absolute vectors, rows of unused records and of
#       a file no segment fits): the engine's r^T M^-1 r equals an independent numpy evaluation
#       from the output files, and the per-segment partition adds up (needs python3 + numpy);
#   (c) the same with model-scaled relvector terms, the model taken from covarianceModel.out;
#   (d) the analytic gradient (adjoint and whitened residual-Jacobian paths) agrees with central
#       finite differences of chi2 for every free parameter (levels, norms, a shift-group shift,
#       a sqrt(E) shift, a shift function, abundances), dense and relvector;
#   (e) MIGRAD (--use-gradient), --use-lm and --use-gsl-lm fit the dense model to one minimum;
#   (f) broken inputs stop the run with an ERROR: wrong record count, a data file missing from
#       the header, overlapping windows on one file, a matrix that is not positive definite, an
#       asymmetric block, a covariance block on a test segment.
#
#   ./tests/data_covariance/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PYTHON="$(command -v python3 || command -v python || true)"
if [ -n "$PYTHON" ] && ! "$PYTHON" -c "import numpy" > /dev/null 2>&1; then PYTHON=""; fi

WORK="$(mktemp -d "${TMPDIR:-/tmp}/data_covariance.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME MODEL MODE "FLAGS" [VAR=VALUE ...] -> a fresh copy of MODEL.azr (as t.azr) in $WORK/NAME
run() {
  local name="$1" model="$2" mode="$3" flags="$4"
  shift 4
  local dir="$WORK/$name"
  rm -rf "$dir"
  mkdir -p "$dir/output" "$dir/checks"
  cp -r "$HERE/../angle_mapping/data" "$dir/"
  cp "$HERE"/cov_*.txt "$dir/"
  cp "$HERE/$model.azr" "$dir/t.azr"
  (cd "$dir" && for v in "$@"; do export "$v"; done &&
   printf "$mode\n\n\n7\n" | run_limited 900 "$AZURE2_BIN" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}
# the segment table of chiSquared.out ("segment chi2" per line)
segchi() { awk -F, '$1 ~ /^[0-9]+$/ { print $1, $2 }' "$WORK/$1/output/chiSquared.out"; }

# ---- (a) diagonal covariance = no covariance ----------------------------------------------------
run plain plain 1 ""
run diag diag 1 ""
if [ ! -s "$WORK/diag/output/covariance.out" ]; then
  bad "diagonal covariance: no covariance.out"; grep -i "error" "$WORK/diag/run.log" | head -3
else
  worst="$(join <(segchi plain) <(segchi diag) | awk '{ d = $2 - $3; if (d < 0) d = -d; s = $2 < 0 ? -$2 : $2; if (s < 1) s = 1
                                                       if (d / s > w) w = d / s; n++ } END { printf "%.1e %d", w, n }')"
  eng="$(awk -F, '!/^#/ { d = $5 - $6; if (d < 0) d = -d; printf "%.1e", d / $6 }' "$WORK/diag/output/covariance.out")"
  if awk -v w="$worst" -v e="$eng" 'BEGIN { split(w, x, " "); exit !(x[2] >= 14 && x[1] + 0 <= 1e-6 && e + 0 <= 1e-10) }'; then
    ok "diag(e^2) covariance: per-segment chi2 equal to no covariance (${worst#* } segments), r^T M^-1 r = sum r^2 to $eng"
  else
    bad "diag(e^2) covariance vs none: segments $worst, r^T M^-1 r vs sum r^2 $eng"
  fi
fi

# ---- (b), (c) dense and relvector against numpy ---------------------------------------------------
for m in dense rel; do
  run "$m" "$m" 1 ""
  if [ ! -s "$WORK/$m/output/covariance.out" ]; then
    bad "$m: no covariance.out"; grep -i "error" "$WORK/$m/run.log" | head -3
    continue
  fi
  if [ -z "$PYTHON" ]; then
    echo "  skip  $m: independent chi2 (no python3 with numpy)"
    continue
  fi
  res="$("$PYTHON" "$HERE/verify.py" "$WORK/$m")"
  if awk -v r="$res" 'BEGIN { split(r, x, " "); d = x[1] - x[2]; if (d < 0) d = -d
                              exit !(x[1] > 0 && d / x[1] <= 1e-8 && x[3] + 0 <= 1e-5) }'; then
    ok "$m: engine r^T M^-1 r = numpy v^T C^-1 v (${res% *}), segment partition = chiSquared.out (${res##* })"
  else
    bad "$m: numpy chi2, engine chi2, segment partition: $res"
  fi
done
if [ -s "$WORK/rel/output/covarianceModel.out" ] && [ ! -e "$WORK/dense/output/covarianceModel.out" ]; then
  ok "relvector run writes the frozen model to covarianceModel.out (the absolute one does not)"
else
  bad "covarianceModel.out: rel $(ls "$WORK/rel/output/covarianceModel.out" 2>&1 | head -1), dense $(ls "$WORK/dense/output/covarianceModel.out" 2>&1 | head -1)"
fi

# ---- (d) gradients ----------------------------------------------------------------------------------
rows() {
  awk '/^\[grad-check\] idx / { idx = $3; a = f = ""; j = "-"
         for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
         print idx, a, f, j }' "$WORK/$1/run.log" | sort -n
}
for m in dense rel; do
  run "g_$m" "$m" 2 "--use-gradient" AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_H=1e-8
  rows "g_$m" > "$WORK/g_$m.rows"
  nrows="$(wc -l < "$WORK/g_$m.rows" | tr -d ' ')"
  kinds="$(awk '/^\[grad-check\] idx / { print $4 }' "$WORK/g_$m/run.log" | sort -u | tr '\n' ' ')"
  if [ "$nrows" -lt 20 ]; then
    bad "$m: gradient check listed $nrows parameters"; grep -i "grad\|error" "$WORK/g_$m/run.log" | head -5
    continue
  fi
  worst="$(awk '{ a = $2 < 0 ? -$2 : $2; if (a > m) m = a; g[NR] = $2; f[NR] = $3; j[NR] = $4 }
             END { for (i = 1; i <= NR; i++) { d = g[i] - f[i]; if (d < 0) d = -d; e = j[i] - f[i]; if (e < 0) e = -e
                                               if (d > w) w = d; if (j[i] != "-" && e > w) w = e }
                   printf "%.2e", w / m }' "$WORK/g_$m.rows")"
  if awk -v w="$worst" 'BEGIN { exit !(w + 0 <= 1e-6) }'; then
    ok "$m: adjoint and Jacobian-path gradients = finite differences, $nrows parameters ($kinds), worst $worst of the largest"
  else
    bad "$m: gradient vs finite differences, worst $worst of the largest (> 1e-6)"
  fi
done

# ---- (e) three fitters --------------------------------------------------------------------------------
objective() {  # data + norm + shift + abundance chi2 of a finished fit
  awk '/^Total-Chi-Squared/ { for (i = 1; i < NF; i++) if ($i ~ /Chi-Squared:$/) s += $(i + 1); printf "%.4f", s }' \
    "$WORK/$1/output/chiSquared.out"
}
declare -a objs=()
for v in "migrad:--use-gradient" "lm:--use-lm" "gsl:--use-gsl-lm"; do
  tag="${v%%:*}"
  run "f_$tag" dense 2 "${v#*:}"
  objs+=("$tag=$(objective "f_$tag")")
done
start="$(objective dense)"
spread="$(printf '%s\n' "${objs[@]}" | awk -F= -v s="$start" '{ x = $2 + 0; if (NR == 1 || x < lo) lo = x; if (NR == 1 || x > hi) hi = x; if (x >= s) up++ }
                                                                END { printf "%.3f %d", (hi - lo) / lo, up + 0 }')"
if awk -v s="$spread" 'BEGIN { split(s, x, " "); exit !(x[1] + 0 <= 2e-3 && x[2] == 0) }'; then
  ok "MIGRAD, LM and GSL fits of the dense model: ${objs[*]} (start $start), within ${spread% *} of each other"
else
  bad "fits of the dense model: ${objs[*]} (start $start), spread/rises $spread"
fi

# ---- (f) broken inputs ---------------------------------------------------------------------------------
broken() {  # NAME MODE "description" EXPECTED-ERROR-REGEX : edit the copy with the rest of the arguments (a command)
  local name="$1" mode="$2" what="$3" pat="$4"
  shift 4
  run "$name" dense 1 ""
  local dir="$WORK/$name"
  rm -rf "$dir/output"/*
  (cd "$dir" && "$@")
  (cd "$dir" && printf "$mode\n\n\n7\n" | run_limited 300 "$AZURE2_BIN" --no-gui --no-readline t.azr > run2.log 2>&1)
  if grep -Eq "ERROR.*($pat)" "$dir/run2.log" && [ -z "$(ls "$dir/output")" ]; then
    ok "refused: $what"
  else
    bad "not refused: $what ($(grep -m1 -i "error" "$dir/run2.log"))"
  fi
}
broken b_count 1 "a header with a wrong record count" "has 61 data rows, the header says 60" \
  sed -i.bak 's|^file data/p12c_90.dat 61|file data/p12c_90.dat 60|' cov_dense.txt
broken b_missing 1 "a segment's data file missing from the header" "not among the" \
  sed -i.bak 's|data/p12c_90s.dat 0 0|data/p12c_90s.dat 0 0 covariance cov_dense.txt|' t.azr
broken b_overlap 1 "two segments fitting one record" "fitted twice" \
  sed -i.bak 's|0.3 2.3 145 155 1|0.3 2.3 115 125 1|' t.azr
broken b_notpd 1 "a matrix that is not positive definite" "not positive definite" \
  sh -c "printf 'block 2\\n1 2\\n0 10 10 0\\n' >> cov_dense.txt"
broken b_asym 1 "an asymmetric block" "not symmetric" \
  sh -c "printf 'block 2\\n1 2\\n0 1e-3 2e-3 0\\n' >> cov_dense.txt"
broken b_test 3 "a covariance block on a test segment" "carries a covariance block" \
  sed -i.bak 's|^</segmentsTest>|1 1 1 0.5 2.0 0.1 60 60 0 1 covariance cov_dense.txt\n</segmentsTest>|' t.azr

echo
if [ "$failures" -eq 0 ]; then echo "data_covariance: all checks passed"; else echo "data_covariance: $failures check(s) FAILED"; fi
[ "$failures" -eq 0 ]
