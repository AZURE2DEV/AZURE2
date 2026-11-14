#!/usr/bin/env bash
#
# Fitted target density (areal_density ... fit [prior ...] on a target record; one target record named by several data
# sets is a density group; the old "densityfit" and "densitygroup" blocks of a <targetInt> line are the converter's).
# Models from make_models.py (deleted; README.md: how the models were made): tests/13N, 12C(p,p) at 84.3 (and 114.5) deg lab, 0.3-0.9 MeV, levels fixed, a thin
# target with a constant stopping cross section (10 keV at 1e18 atoms/cm^2).
#
#   1. data made at N = 1.1e18 (true.azr, 1 % errors): MIGRAD with the analytic gradient, --use-lm and --use-gsl-lm
#      recover the areal density N = 1.1e18 from 1e18 to 1e-4, with a final chi2 below 0.01 (122 at the start);
#   2. the gradient: d chi2 / d N from finite differences agrees with the Jacobian path (AZURE_GRAD_CHECK);
#   3. a 2 % prior about 1e18 pulls the fit between the two, Total-Density-Chi-Squared = ((N - 1e18)/2e16)^2;
#   4. a density group: two segments, one parameter (t.areal_density, no other), recovered from both;
#   5. densityfit on a line without target integration is refused (by the converter);
#   6. a density beyond the target's grid_max (from a parameter file): one warning with the parameter's name.
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
  local name="$1" azr="$2" command="$3" flags="$4"; shift 4
  local dir="$WORK/$name"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$SRC/data" "$dir/data"
  [ -d "$WORK/fitdata" ] && ln -s "$WORK/fitdata" "$dir/fitdata"
  cp "$azr" "$dir/t.azr"
  (cd "$dir" && for kv in "$@"; do export "$kv"; done
   run_limited 900 "$AZURE2_BIN" "$command" $flags --no-gui t.azr > run.log 2>&1)
}
out() { echo "$WORK/$1/output/AZUREOut_aa=1_R=1.out"; }
par() { awk -v n="$2" '$1 == n { print $2; exit }' "$WORK/$1/output/param.sav"; }
chi() { sed -n "s/.*$2: \([^ ]*\).*/\1/p" "$WORK/$1/output/chiSquared.out" | sed -n 1p; }
near() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a / b - 1; if (d < 0) d = -d; exit !(d < t) }'; }

# ---- data from the truth, N = 1.1e18 --------------------------------------------------------
run true "$HERE/true.azr" calculate ""
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
  run "$name" "$HERE/fit.azr" fit "$flags"
  n="$(par "$name" tgt_meyer_84.areal_density)"; c="$(chi "$name" Total-Chi-Squared)"
  if [ -n "$n" ] && near "$n" 1.1e18 1e-4 && awk -v c="${c:-1e9}" 'BEGIN { exit !(c ~ /^ *[-+]?\.?[0-9]/ && c < 1e-2) }'; then
    ok "$name ($flags): N $n from 1e18 (truth 1.1e18), chi2 $c"
  else bad "$name ($flags): N ${n:-missing}, chi2 ${c:-missing}"; tail -5 "$WORK/$name/run.log"; fi
done

# ---- 2. the gradient ------------------------------------------------------------------------
run grad "$HERE/fit.azr" fit "--use-gradient" AZURE_GRAD_CHECK=exit
line="$(grep "components with relative disagreement" "$WORK/grad/run.log" | sed -n 1p)"
if grep -q "density 0/1;.*density 0/1" <<< "$line"; then ok "d chi2/d N: finite differences = the Jacobian path ($line)"
else bad "gradient check: ${line:-no [grad-check] line}"; tail -5 "$WORK/grad/run.log"; fi

# ---- 3. a prior -------------------------------------------------------------------------------
run prior "$HERE/fit_prior.azr" fit "--use-gradient"
n="$(par prior tgt_meyer_84.areal_density)"; dc="$(chi prior Total-Density-Chi-Squared)"
if [ -n "$n" ] && [ -n "$dc" ] && awk -v n="$n" -v s="$dc" 'BEGIN { p = ((n - 1e18) / 2e16)^2; d = s - p; if (d < 0) d = -d
       exit !(n ~ /^ *[-+]?\.?[0-9]/ && n > 1.0e18 && n < 1.1e18 && d <= 1e-4 * (p + 1)) }'; then
  ok "a 2 % prior pulls N to $n; Total-Density-Chi-Squared $dc = ((N - 1e18)/2e16)^2"
else bad "prior fit: N ${n:-missing}, Total-Density-Chi-Squared ${dc:-missing}"; tail -5 "$WORK/prior/run.log"; fi

# ---- 4. a density group -----------------------------------------------------------------------
run group "$HERE/fit_group.azr" fit "--use-gradient"
n="$(par group t.areal_density)"; n2="$(awk '$1 ~ /areal_density/ && $1 != "t.areal_density" { print $1 }' "$WORK/group/output/param.sav")"; c="$(chi group Total-Chi-Squared)"
if [ -n "$n" ] && [ -z "$n2" ] && near "$n" 1.1e18 1e-4 && awk -v c="${c:-1e9}" 'BEGIN { exit !(c ~ /^ *[-+]?\.?[0-9]/ && c < 1e-2) }'; then
  ok "a density group: one parameter, N $n from both segments, chi2 $c"
else bad "group fit: N ${n:-missing}, density_2 '${n2}', chi2 ${c:-missing}"; tail -5 "$WORK/group/run.log"; fi
# the old group's second line carried a density of its own, which the old engine replaced by the group's: the
# converter drops it with a note (the format-3 group is one target record)
OLDG="$HERE/../convert/fixtures/tests/target_density"
(cd "$HERE" && "$AZURE2_BIN" --convert "$OLDG/fit_group.azr" "$WORK/group_conv.azr" --to 3 > "$WORK/group_conv.log" 2>&1)
if grep -q "own density 2e+18 was not used by the old engine (the group's is its first line's): dropped" "$WORK/group_conv.log"; then
  ok "the old group member's own density: dropped with a note by the converter"
else bad "no converter note about the group member's own density"; fi

# ---- 5. refused -------------------------------------------------------------------------------
# the old file (tests/convert/fixtures): a densityfit on a line without target integration, refused by the old engine;
# the converter refuses it with that reason (format 3 has a fitted areal density only on a target record)
(cd "$HERE" && "$AZURE2_BIN" --convert "$OLDG/bad.azr" "$WORK/bad_conv.azr" --to 3 > "$WORK/bad_conv.log" 2>&1)
if grep -q "densityfit / densitygroup need target integration" "$WORK/bad_conv.log" && [ ! -f "$WORK/bad_conv.azr" ]; then
  ok "densityfit on a line without target integration is refused (by the converter)"
else bad "densityfit without target integration was not refused"; tail -3 "$WORK/bad_conv.log"; fi

# ---- 6. a value beyond the grids ----------------------------------------------------------------
# a parameter file that puts the group's density above its grid_max (1.5e18): the calculation uses the grid's end and
# says so once, naming the parameter by its format-3 name and the field to widen
run over "$HERE/fit_group.azr" calculate ""
printf '    t.areal_density   2.0000000e+18   1.0000000e+16\n' > "$WORK/over/over.sav"
(cd "$WORK/over" && run_limited 900 "$AZURE2_BIN" calculate --parameter-file=over.sav --no-gui t.azr > run.log 2>&1)
said="$(grep -c "WARNING: t.areal_density = 2e+18 is outside \[1e+16, 1.5e+18\], the range its sub-point grids span; the calculation uses the nearest end.  Give the target's areal_density a larger grid_max." "$WORK/over/run.log")"
if [ "$said" = 1 ] && ! grep -q "density_[0-9]\|densityfit" "$WORK/over/run.log"; then
  ok "a density beyond grid_max: one warning naming t.areal_density and grid_max"
else bad "a density beyond grid_max: warning lines $said: $(grep -m1 'is outside' "$WORK/over/run.log")"; fi

if [ "$failures" -eq 0 ]; then echo "PASS target_density"; exit 0; fi
echo "FAIL target_density ($failures)"; exit 1
