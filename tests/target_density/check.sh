#!/usr/bin/env bash
#
# Fitted target density ("densityfit vary sigma [max]" and "densitygroup <label>" on a <targetInt> line).
# Models from make_models.py: tests/13N, 12C(p,p) at 84.3 (and 114.5) deg lab, 0.3-0.9 MeV, levels fixed, a thin
# target with a constant stopping cross section (10 keV at 1e18 atoms/cm^2).
#
#   1. data made at N = 1.1e18 (true.azr, 1 % errors): MIGRAD with the analytic gradient, --use-lm and --use-gsl-lm
#      recover density_1 = 1.1e18 from 1e18 to 1e-4, with a final chi2 below 0.01 (122 at the start);
#   2. the gradient: d chi2 / d N from finite differences agrees with the Jacobian path (AZURE_GRAD_CHECK);
#   3. a 2 % prior about 1e18 pulls the fit between the two, Total-Density-Chi-Squared = ((N - 1e18)/2e16)^2;
#   4. a density group: two segments, one parameter (density_1; no density_2), recovered from both;
#   5. densityfit on a line without target integration is refused.
#
#   ./tests/target_density/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/target_density.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME AZR MODE FLAGS [env...] -> $WORK/NAME (mode 1 calculate, 2 fit)
run() {
  local name="$1" azr="$2" mode="$3" flags="$4"; shift 4
  local dir="$WORK/$name"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$SRC/data" "$dir/data"
  [ -d "$WORK/fitdata" ] && ln -s "$WORK/fitdata" "$dir/fitdata"
  cp "$azr" "$dir/t.azr"
  local input='1\n\n\n7\n'
  [ "$mode" = 2 ] && input='2\nn\n\n\n'
  (cd "$dir" && for kv in "$@"; do export "$kv"; done
   printf "$input" | run_limited 900 "$AZURE2_BIN" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}
out() { echo "$WORK/$1/output/AZUREOut_aa=1_R=1.out"; }
par() { awk -v n="$2" '$1 == n { print $2; exit }' "$WORK/$1/output/param.sav"; }
chi() { sed -n "s/.*$2: \([^ ]*\).*/\1/p" "$WORK/$1/output/chiSquared.out" | head -1; }
near() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a / b - 1; if (d < 0) d = -d; exit !(d < t) }'; }

# ---- data from the truth, N = 1.1e18 --------------------------------------------------------
run true "$HERE/true.azr" 1 ""
if [ -s "$(out true)" ]; then
  mkdir -p "$WORK/fitdata"
  awk 'NF >= 7' "$(out true)" > "$WORK/true_rows"
  # lab value = model_cm x (data_lab / data_cm): the frame factor of each point, as the engine applies it
  i=0
  for f in meyer_84 meyer_114; do
    awk -v skip="$i" 'NR == FNR { if (FNR > skip) { m[FNR - skip] = $4; d[FNR - skip] = $6 } next }
                      NF >= 4 && $1 >= 0.3 && $1 <= 0.9 { k++; v = m[k] * $3 / d[k]; printf "%s %s %.12e %.12e\n", $1, $2, v, 0.01 * v }' \
        "$WORK/true_rows" "$SRC/data/$f.dat" > "$WORK/fitdata/$f.dat"
    i=$((i + $(awk '$1 >= 0.3 && $1 <= 0.9' "$SRC/data/$f.dat" | grep -c .)))
  done
  ok "data at N = 1.1e18: $(grep -c . "$WORK/fitdata/meyer_84.dat") + $(grep -c . "$WORK/fitdata/meyer_114.dat") points"
else
  bad "true.azr: no output"; tail -5 "$WORK/true/run.log"
fi

# ---- 1. fits recover the density --------------------------------------------------------------
for v in "fit:--use-gradient" "fit_lm:--use-lm" "fit_gsl:--use-gsl-lm"; do
  name="${v%%:*}"; flags="${v#*:}"
  run "$name" "$HERE/fit.azr" 2 "$flags"
  n="$(par "$name" density_1)"; c="$(chi "$name" Total-Chi-Squared)"
  if [ -n "$n" ] && near "$n" 1.1e18 1e-4 && awk -v c="${c:-1e9}" 'BEGIN { exit !(c < 1e-2) }'; then
    ok "$name ($flags): density_1 $n from 1e18 (truth 1.1e18), chi2 $c"
  else bad "$name ($flags): density_1 ${n:-missing}, chi2 ${c:-missing}"; tail -5 "$WORK/$name/run.log"; fi
done

# ---- 2. the gradient ------------------------------------------------------------------------
run grad "$HERE/fit.azr" 2 "--use-gradient" AZURE_GRAD_CHECK=exit
line="$(grep "components with relative disagreement" "$WORK/grad/run.log" | head -1)"
if echo "$line" | grep -q "density 0/1;.*density 0/1"; then ok "d chi2/d N: finite differences = the Jacobian path ($line)"
else bad "gradient check: ${line:-no [grad-check] line}"; tail -5 "$WORK/grad/run.log"; fi

# ---- 3. a prior -------------------------------------------------------------------------------
run prior "$HERE/fit_prior.azr" 2 "--use-gradient"
n="$(par prior density_1)"; dc="$(chi prior Total-Density-Chi-Squared)"
if [ -n "$n" ] && [ -n "$dc" ] && awk -v n="$n" -v s="$dc" 'BEGIN { p = ((n - 1e18) / 2e16)^2; d = s - p; if (d < 0) d = -d
       exit !(n > 1.0e18 && n < 1.1e18 && d <= 1e-4 * (p + 1)) }'; then
  ok "a 2 % prior pulls N to $n; Total-Density-Chi-Squared $dc = ((N - 1e18)/2e16)^2"
else bad "prior fit: density_1 ${n:-missing}, Total-Density-Chi-Squared ${dc:-missing}"; tail -5 "$WORK/prior/run.log"; fi

# ---- 4. a density group -----------------------------------------------------------------------
run group "$HERE/fit_group.azr" 2 "--use-gradient"
n="$(par group density_1)"; n2="$(par group density_2)"; c="$(chi group Total-Chi-Squared)"
if [ -n "$n" ] && [ -z "$n2" ] && near "$n" 1.1e18 1e-4 && awk -v c="${c:-1e9}" 'BEGIN { exit !(c < 1e-2) }'; then
  ok "a density group: one parameter, density_1 $n from both segments, chi2 $c"
else bad "group fit: density_1 ${n:-missing}, density_2 '${n2}', chi2 ${c:-missing}"; tail -5 "$WORK/group/run.log"; fi
if grep -q "the group's density is its first line's" "$WORK/group/run.log"; then ok "the second line's own density is noted and replaced"
else bad "no note about the group member's own density"; fi

# ---- 5. refused -------------------------------------------------------------------------------
run bad "$HERE/bad.azr" 1 ""
if grep -q "densityfit / densitygroup need target integration" "$WORK/bad/run.log" && [ ! -s "$(out bad)" ]; then
  ok "densityfit on a line without target integration is refused"
else bad "densityfit without target integration was not refused"; tail -5 "$WORK/bad/run.log"; fi

if [ "$failures" -eq 0 ]; then echo "PASS target_density"; exit 0; fi
echo "FAIL target_density ($failures)"; exit 1
