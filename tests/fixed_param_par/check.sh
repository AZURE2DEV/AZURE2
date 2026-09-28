#!/usr/bin/env bash
#
# A parameter fixed in an external param.par stays fixed under every
# minimizer.
#
# MIGRAD minimizes the AZUREParams it is handed, whose fixed flags include an
# external param.par's ("fixed" in the fourth column).  The Levenberg-Marquardt
# (--use-lm) and GSL trust-region (--use-gsl-lm) solvers built their free set
# from a fresh mask of the .azr alone, so they fitted a parameter the user had
# fixed there -- and the band covariance of such a fit had a column too many.
#
# Model: the one-level p+7Li -> a+a project of tests/identical_entrance_reaction
# (2+, E, Gamma_p, Gamma_a all free in the .azr) against nine synthetic points
# made with E = 17.762 MeV, Gamma_p = 0.9, Gamma_a = 0.15.  param.par starts
# from the .azr values and fixes width_1_1 (Gamma_p) at 0.857457.  Each of
# MIGRAD, --use-lm and --use-gsl-lm must leave width_1_1 exactly there, move
# energy_1 and width_1_2, and reach the same chi2 (to 1e-5); --use-lm
# --covariance-band must save a 2x2 covariance.  Seven fits of ~0.1 s.
#
#   ./tests/fixed_param_par/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/fixed_param_par.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
if command -v timeout >/dev/null 2>&1; then RUN="timeout ${TEST_TIMEOUT:-120}"; else RUN=""; fi
fail=0

mkdir -p "$WORK/proj/data"
sed -n '1,/<\/levels>/p' "$HERE/../identical_entrance_reaction/projects/be8.azr" > "$WORK/proj/run.azr"
cat >> "$WORK/proj/run.azr" <<'EOF'
<segmentsData>
1 1 2 0 100 0 180 0 1 0 0 0 0 0 data/p7li_aa_int.dat 0 0
</segmentsData>
<segmentsTest>
</segmentsTest>
<targetInt>
</targetInt>
EOF
cat > "$WORK/proj/data/p7li_aa_int.dat" <<'EOF'
0.3 0 3.315412e-02 9.946235e-04
0.4 0 1.237896e-01 3.713687e-03
0.5 0 4.449601e-01 1.334880e-02
0.55 0 6.902651e-01 2.070795e-02
0.6 0 7.241651e-01 2.172495e-02
0.65 0 5.533478e-01 1.660043e-02
0.7 0 3.866714e-01 1.160014e-02
0.8 0 2.048842e-01 6.146525e-03
1.0 0 8.804675e-02 2.641402e-03
EOF
# In the project, so every run names it relative to its own directory: the
# path is read from stdin, where MSYS does not translate a POSIX path, and the
# native Windows binary cannot open /tmp/... (it re-prompted, took the next
# blank line for "new file" and fitted without the fixed flag).
cat > "$WORK/proj/fix.par" <<'EOF'
energy_1 1.7755100e+01 1.7755100e+00
width_1_1 8.5745701e-01 8.5745701e-02 fixed
width_1_2 1.2520344e-01 1.2520344e-02
segment_1_energy_shift 0.0000000e+00 5.0000000e-04
EOF

# fit NAME STDIN [flags] -- a fit (menu 2) from a copy of proj
fit() {
  local d="$WORK/$1" in="$2"
  shift 2
  cp -r "$WORK/proj" "$d"
  mkdir -p "$d/output" "$d/checks"
  (cd "$d" && printf "$in" | $RUN "$AZURE2_BIN" --no-gui --no-readline "$@" run.azr 2>&1 | head -c 1000000 > log)
  [ -f "$d/output/param.sav" ] || { echo "  FAIL  $1: no param.sav"; tail -3 "$d/log" | sed 's/^/        /'; fail=1; return 1; }
}
value() { awk -v n="$2" '$1 == n { print $2 }' "$1"; }

chis=()
for m in migrad lm gsl-lm; do
  flag=""; [ "$m" = migrad ] || flag="--use-$m"
  # Menu 2 asks about the band first (n), then the parameter file.
  fit "$m" "2\nn\nfix.par\n\n7\n" $flag || continue
  sav="$WORK/$m/output/param.sav"
  w11="$(value "$sav" width_1_1)"; e1="$(value "$sav" energy_1)"; w12="$(value "$sav" width_1_2)"
  chi="$(grep -oE 'Total Chi-Squared: [0-9.eE+-]+' "$WORK/$m/log" | awk '{ print $3 }')"
  chis+=("$chi")
  echo "  $m: width_1_1 $w11  energy_1 $e1  width_1_2 $w12  chi2 $chi"
  if awk -v w="$w11" 'BEGIN { d = (w - 0.85745701) / 0.85745701; exit !(d < 1e-12 && d > -1e-12) }'; then
    echo "  ok    $m keeps the param.par-fixed width_1_1"
  else
    echo "  FAIL  $m moved width_1_1, fixed in param.par"; fail=1
  fi
  if awk -v e="$e1" -v w="$w12" 'BEGIN { de = (e - 17.7551) / 17.7551; dw = (w - 0.12520344) / 0.12520344
        exit !((de > 1e-5 || de < -1e-5) && (dw > 1e-3 || dw < -1e-3)) }'; then
    echo "  ok    $m fits the free parameters"
  else
    echo "  FAIL  $m did not move energy_1 and width_1_2"; fail=1
  fi
done
if [ "${#chis[@]}" -eq 3 ] && awk -v a="${chis[0]}" -v b="${chis[1]}" -v c="${chis[2]}" 'BEGIN {
      d1 = (b - a) / a; d2 = (c - a) / a; if (d1 < 0) d1 = -d1; if (d2 < 0) d2 = -d2
      exit !(d1 < 1e-5 && d2 < 1e-5) }'; then
  echo "  ok    MIGRAD, LM and GSL-LM reach the same minimum"
else
  echo "  FAIL  the minimizers disagree: ${chis[*]}"; fail=1
fi

# The band covariance of an LM fit spans the free R-matrix parameters only.
if fit lmband "2\nfix.par\n\n7\n" --use-lm --covariance-band; then
  cov="$WORK/lmband/output/covariance.dat"
  if [ -f "$cov" ] && awk 'NF { r++; if (NF != 2) bad = 1 } END { exit !(r == 2 && !bad) }' "$cov"; then
    echo "  ok    LM band covariance is 2x2 (energy_1, width_1_2)"
  else
    echo "  FAIL  LM band covariance is not 2x2"; cat "$cov" 2>/dev/null | head -4; fail=1
  fi
fi

if [ "$fail" -ne 0 ]; then echo "FAIL"; exit 1; fi
echo "PASS: a param.par-fixed parameter stays fixed under MIGRAD, LM and GSL-LM"
