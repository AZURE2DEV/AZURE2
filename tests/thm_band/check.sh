#!/usr/bin/env bash
#
# The CLI cross-section band (--covariance-band) of a THM (HOES) point is
# sqrt(g Sigma g^T) with g the derivative of the THM model *as it lies against
# the data*: q(p) = m(p) n*(p0)/n*(p), m the HOES model and n* the profiled
# norm of the segment.  The band used to take the T-matrix adjoint row for
# HOES points, which differentiates a different observable.
#
# Model: tests/18O_p_a_thm (two 1/2+ levels, two THM segments with profiled
# norms) with all six R-matrix parameters freed.  A hand-made 6x6 covariance
# (output/covariance.dat, correlated) is put in place and a plain calculation
# with --covariance-band reads it.  The reference g is taken here by central
# differences of the CLI's own output: each parameter of param.par is moved by
# +-h, the model column m and the scaled-data column d n* of
# AZUREOut_aa=1_R=2.out are read, and q_i = m_i D_i(p0)/D_i(p) (D = d n*, the
# conversion to the CM frame cancels in the ratio).  Every point whose band is
# not negligible must agree to 2e-3 relative.
#
# A second pass groups the two segments into one THM experiment with a linear
# background (<thm> experiment[A] segments=1,2 background=linear): the curve
# written is then m + b(E), with the norm and background profiled together,
# and q = (m + b) n*(p0)/n*(p) -- the same finite differences of the output.
# A third adds a free coherent background (cbackground=, two cbkg values):
# it enters the HOES amplitude, so its columns belong in the band (8x8).
#
#   ./tests/thm_band/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../18O_p_a_thm"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TOP="$(mktemp -d "${TMPDIR:-/tmp}/thm_band.XXXXXX")"
trap 'rm -rf "$TOP"' EXIT
if command -v timeout >/dev/null 2>&1; then RUN="timeout ${TEST_TIMEOUT:-300}"; else RUN=""; fi
OUT="AZUREOut_aa=1_R=2.out"

# band_case LABEL BLOCK [NAMES SIGMAS] -- the whole check on the project with
# BLOCK (lines) as its <thm> block ("" = none); NAMES/SIGMAS: free parameters
# beyond the six R-matrix ones and their standard deviations.
band_case() {
WORK="$TOP/$1"
echo "== $1"
# Project with every level energy and channel width free (fields 4 and 11).
mkdir -p "$WORK/proj"
cp -r "$SRC/data" "$WORK/proj/"
awk '/<levels>/ { L = 1; print; next } /<\/levels>/ { L = 0 }
     L && NF > 0 { $4 = 0; $11 = 0 } { print }' "$SRC/18O_p_a_thm.azr" > "$WORK/proj/run.azr"
[ -z "$2" ] || printf '\n<thm>\n%s\n</thm>\n' "$2" >> "$WORK/proj/run.azr"

# run DIR PARFILE [FLAGS] -- a "calculate with data" run from a copy of proj.
run() {
  local d="$WORK/$1"
  rm -rf "$d"
  cp -r "$WORK/proj" "$d"
  mkdir -p "$d/output" "$d/checks"
  [ -z "${4:-}" ] || cp "$4" "$d/output/covariance.dat"
  # The parameter file is handed over relative to the run directory: the path
  # is read from stdin, where MSYS does not translate a POSIX path, and the
  # native Windows binary cannot open /tmp/...
  local par=""
  [ -z "$2" ] || { cp "$2" "$d/start.par"; par="start.par"; }
  (cd "$d" && printf '1\n%s\n\n7\n' "$par" | $RUN "$AZURE2_BIN" --no-gui --no-readline $3 run.azr 2>&1 \
     | head -c 1000000 > log)
  [ -f "$d/output/$OUT" ] || { echo "FAIL: run $1 produced no $OUT"; tail -5 "$d/log"; exit 1; }
}

# Base run: writes param.par (the internal parameter values).  Re-write it at
# full precision so every run below starts from the identical point.
run base "" ""
awk '{ printf "%s %.15e %.15e\n", $1, $2, $3 }' "$WORK/base/output/param.par" > "$WORK/p0.par"
names=(energy_1 width_1_1 width_1_2 energy_2 width_2_1 width_2_2 ${3:-})
sigmas="2e-3 0.02 0.01 3e-3 0.02 0.005 ${4:-}"
np=${#names[@]}
for n in "${names[@]}"; do
  grep -q "^$n " "$WORK/p0.par" || { echo "FAIL: $n not in param.par"; exit 1; }
done

# Covariance, in the order of the free parameters (= the names above).
# sigma = (2 keV, 0.02, 0.01, 3 keV, 0.02, 0.005, extra), correlations -0.4
# (E1,E2), 0.5 (w11,w12), 0.3 (w21,w22), -0.2 (E1,w21); with two extra
# parameters 0.3 between them and -0.25 (w11, first extra).
awk -v n="$np" -v sig="$sigmas" 'BEGIN {
  split(sig, s, " ")
  for (i = 1; i <= n; i++) for (j = 1; j <= n; j++) c[i, j] = (i == j)
  c[1, 4] = c[4, 1] = -0.4; c[2, 3] = c[3, 2] = 0.5; c[5, 6] = c[6, 5] = 0.3
  c[1, 5] = c[5, 1] = -0.2
  if (n >= 8) { c[7, 8] = c[8, 7] = 0.3; c[2, 7] = c[7, 2] = -0.25 }
  for (i = 1; i <= n; i++) { line = ""
    for (j = 1; j <= n; j++) line = line sprintf("%.17e%s", c[i, j] * s[i] * s[j], j < n ? " " : "")
    print line }
}' > "$WORK/cov.dat"

# The band run.
run band "$WORK/p0.par" "--covariance-band" "$WORK/cov.dat"
[ -f "$WORK/band/output/$OUT.band" ] || { echo "FAIL: no $OUT.band"; tail -5 "$WORK/band/log"; exit 1; }
grep -q "Writing cross-section uncertainty bands" "$WORK/band/log" \
  || { echo "FAIL: band not written"; grep -i band "$WORK/band/log"; exit 1; }

# Central differences.  Energies h = 1e-5 MeV, widths (and cbkg values)
# h = 1e-4 |x|.  Files are numbered: a cbkg name holds '/'.
k=0
for n in "${names[@]}"; do
  k=$((k + 1))
  for sgn in p m; do
    awk -v n="$n" -v sgn="$sgn" '
      $1 == n { h = (n ~ /^energy/) ? 1e-5 : 1e-4 * ($2 < 0 ? -$2 : $2)
                $2 = (sgn == "p") ? $2 + h : $2 - h
                printf "%s %.15e %.15e\n", $1, $2, $3; print h > "/dev/stderr"; next }
      { print }' "$WORK/p0.par" > "$WORK/p$k.$sgn.par" 2> "$WORK/p$k.h"
    run "fd_${k}_$sgn" "$WORK/p$k.$sgn.par" ""
  done
done

# Compare.  Output lines: E Ex angle m S D ... (blank lines between segments).
awk -v np="$np" -v W="$WORK" -v OUT="$OUT" '
  function readout(file, M, D,    i) {
    i = 0
    while ((getline line < file) > 0) {
      nf = split(line, f, " "); if (nf < 6) continue
      i++; M[i] = f[4]; D[i] = f[6]
    }
    close(file); return i
  }
  BEGIN {
    n0 = readout(W "/band/output/" OUT, M0, D0)
    # band file: E Ex angle xs dxs S dS
    nb = 0
    while ((getline line < (W "/band/output/" OUT ".band")) > 0) {
      nf = split(line, f, " "); if (nf < 7) continue
      nb++; XS[nb] = f[4]; DX[nb] = f[5]
    }
    if (nb != n0 || n0 == 0) { printf "FAIL: %d band lines vs %d output lines\n", nb, n0; exit 1 }
    for (k = 1; k <= np; k++) {
      getline h < (W "/p" k ".h")
      readout(W "/fd_" k "_p/output/" OUT, Mp, Dp)
      readout(W "/fd_" k "_m/output/" OUT, Mm, Dm)
      for (i = 1; i <= n0; i++)
        g[i, k] = (Mp[i] * D0[i] / Dp[i] - Mm[i] * D0[i] / Dm[i]) / (2 * h)
    }
    # covariance
    r = 0
    while ((getline line < (W "/cov.dat")) > 0) { r++; split(line, f, " "); for (c = 1; c <= np; c++) S[r, c] = f[c] }
    maxb = 0
    for (i = 1; i <= n0; i++) {
      v = 0
      for (a = 1; a <= np; a++) for (b = 1; b <= np; b++) v += g[i, a] * S[a, b] * g[i, b]
      ref[i] = (v > 0) ? sqrt(v) : 0
      if (ref[i] > maxb) maxb = ref[i]
      if (XS[i] != M0[i]) { printf "FAIL: band xs %s != model %s at line %d\n", XS[i], M0[i], i; exit 1 }
    }
    bad = 0; worst = 0; used = 0
    for (i = 1; i <= n0; i++) {
      if (ref[i] < 1e-3 * maxb) continue
      used++
      e = (DX[i] - ref[i]) / ref[i]; if (e < 0) e = -e
      if (e > worst) worst = e
      if (e > 2e-3) { bad++; if (bad <= 5) printf "  point %d: band %.6e, finite differences %.6e (rel %.2e)\n", i, DX[i], ref[i], e }
    }
    if (used < 0.9 * n0) { printf "FAIL: only %d of %d points compared\n", used, n0; exit 1 }
    if (bad) { printf "FAIL: %d of %d THM band points off the finite-difference J Sigma J^T\n", bad, used; exit 1 }
    printf "PASS: THM band = sqrt(J Sigma J^T) of the profiled model, %d parameters, at %d points (worst rel %.2e)\n", np, used, worst
  }'
}

band_case per_segment "" || exit 1
band_case experiment_linear "experiment[A] segments=1,2 background=linear" || exit 1
band_case coherent_background "experiment[A] segments=1,2 cbackground=1/2+:2=0.3,-0.2" \
  "cbkg_A_1/2+_2_1/2,0,1/2,1_re0 cbkg_A_1/2+_2_1/2,0,1/2,1_im0" "0.05 0.04" || exit 1
