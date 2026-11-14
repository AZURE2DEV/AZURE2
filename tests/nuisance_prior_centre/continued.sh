#!/usr/bin/env bash
#
# A fit continued from its fitted values computes the fit it continues (docs/upgrade/parameter_kinds.md section 4,
# "The layout of a fitted value"): the grids and limits of a fitted areal density, Doppler temperature, stopping scale
# and kernel parameter are laid out from the prior's centre, so a project with the fitted value written as the start
# (the prior unchanged) is laid out as the fit was.  For each of the four kinds, on projects of the tree:
#
#   (1) the fit (MIGRAD, --use-gradient), its param.sav; the fitted value written into the record (the value token
#       only, as save_fit and the web UI's "accept" write it);
#   (2) the project and the written project calculated at the fit's param.sav: chiSquared.out identical (the total
#       and every data set, to the digits printed; the start differs, nothing else);
#   (3) a second fit started from the written project stays at the first's minimum: chi^2 plus the priors equal to
#       the six digits printed, the fitted value within 1e-2 of its uncertainty (MIGRAD's second start refines the
#       first's EDM of up to 1e-6 along the flat direction density-prior, so the value is not equal digit for digit);
#   (4) a start outside the range laid out from the centre (D.doppler 400 K under prior 300 20: [75, 363] K) is laid
#       out from the start, and one line of output says so;
#   (5) with PYAZR_ENGINE_DIR: save_fit(verify=True) of each fit (it checks the written project's chi^2 against the
#       fit's to 1e-6), the layout value of the written project equal to the fit's (density_max, doppler_max, the
#       apparatus layout), and in (4) the start's layout (484 K = 1.21 x 400 K).
#
#   density   tests/target_density/fit_prior.azr (data at N = 1.1e18, prior 1e18 2 %)
#   doppler   tests/doppler_fgm fit_prior (the dop450 values, 20 % errors, prior 300 20)
#   scale     tests/stopping_yield/models/fit.azr, the densities fixed at the truth, scale 1 fit prior 1 4 %
#   kernel    tests/rectangular_resolution/models/fit.azr, width 14 keV fit prior 14 3 (data at 20 keV)
#
# Usage: tests/nuisance_prior_centre/continued.sh path/to/AZURE2     (about a minute at OMP_NUM_THREADS=2)

export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TESTS="$HERE/.."
AZURE2_BIN="${1:?usage: continued.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PY="$(command -v python3 || command -v python)"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/nuisance_continued.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# project NAME SRC_DIR PROJECT SED: $WORK/NAME with SRC_DIR's files linked and PROJECT edited by SED as t.azr
project() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  for f in "$2"/*; do
    case "$(basename "$f")" in output | checks | *.azr) ;; *) ln -s "$f" "$dir/$(basename "$f")" ;; esac
  done
  sed -e "$4" "$2/$3" > "$dir/t.azr"
}
# azure2 DIR PROJECT COMMAND [FLAGS]: the program in DIR
azure2() { (cd "$1" && run_limited 900 "$AZURE2_BIN" $3 $4 --no-gui "$2" > "${2%.azr}.$5log" 2>&1); }
par() { awk -v n="$2" '$1 == n { print $2; exit }' "$1"; }
total() { sed -n 's/.*Total-Chi-Squared: \([^ ]*\).*/\1/p' "$1" | tr -d '\r' | head -1; }
# the objective of chiSquared.out: the data's chi^2 and every prior term (the Total-...-Chi-Squared fields), 6 digits
objective() { tr -d '\r' < "$1" | awk '$1 == "Total-Chi-Squared:" { for (i = 1; i < NF; i++) if ($i ~ /^Total-.*Chi-Squared:$/) s += $(i + 1)
                                       printf "%.6g\n", s; exit }'; }

# continued NAME PARAM FIELD_TEXT UNIT_TEXT: (1)-(3) for the project $WORK/NAME/t.azr, whose parameter PARAM is
# written on the line containing FIELD_TEXT followed by its value and UNIT_TEXT
continued() {
  local dir="$WORK/$1" name="$2" field="$3" unit="$4" v v2
  azure2 "$dir" t.azr fit --use-gradient fit.
  cp "$dir/output/param.sav" "$dir/fit.sav" 2> /dev/null; cp "$dir/output/chiSquared.out" "$dir/fit_chi.out" 2> /dev/null
  v="$(par "$dir/fit.sav" "$name")"
  if [ -z "$v" ]; then bad "($1) the fit wrote no $name"; tail -5 "$dir/t.fit.log"; return; fi
  # (1) the value token replaced, the prior untouched
  awk -v f="$field" -v u="$unit" -v v="$v" 'index($0, f) {
      i = index($0, f) + length(f); rest = substr($0, i); sub(/^ *[^ ]+/, "", rest)
      if (substr(rest, 1, length(u) + 1) == " " u) $0 = substr($0, 1, i - 1) " " v rest } { print }' \
    "$dir/t.azr" > "$dir/written.azr"
  if ! grep -qF "$field $v $unit" "$dir/written.azr"; then bad "($1) $name = $v not written"; return; fi
  # (2) both projects at the fit's parameters
  azure2 "$dir" t.azr 'calculate --parameter-file=fit.sav' "" at.; cp "$dir/output/chiSquared.out" "$dir/t_chi.out"
  azure2 "$dir" written.azr 'calculate --parameter-file=fit.sav' "" at.; cp "$dir/output/chiSquared.out" "$dir/w_chi.out"
  if [ -s "$dir/t_chi.out" ] && [ "$(tr -d '\r' < "$dir/t_chi.out")" = "$(tr -d '\r' < "$dir/w_chi.out")" ]; then
    ok "($1) the written project ($name = $v) computes the fit's chi^2, total and per data set: $(total "$dir/w_chi.out")"
  else
    bad "($1) the written project's chi^2 differs from the fit's:"
    diff <(tr -d '\r' < "$dir/t_chi.out") <(tr -d '\r' < "$dir/w_chi.out") | head -6
  fi
  # (3) a second fit from the written project
  azure2 "$dir" written.azr fit --use-gradient fit.
  v2="$(par "$dir/output/param.sav" "$name")"
  local e o1 o2
  e="$(awk -v n="$name" '$1 == n { print $3; exit }' "$dir/fit.sav")"
  o1="$(objective "$dir/fit_chi.out")"; o2="$(objective "$dir/output/chiSquared.out")"
  if awk -v a="$o1" -v b="$o2" -v p="$v" -v q="$v2" -v e="$e" 'BEGIN { d = a - b; if (d < 0) d = -d; f = p - q; if (f < 0) f = -f
       exit !(q != "" && b != "" && d <= 1e-5 * a && f <= 1e-2 * e) }'; then
    ok "($1) a second fit from the written project stays at the minimum: chi^2 + priors $o2 (first $o1), $name $v2 (first $v, uncertainty $e)"
  else
    bad "($1) the second fit moved: chi^2 + priors '$o2' (first $o1), $name '$v2' (first $v, uncertainty $e)"
  fi
}

# ---- the projects
SRC13N="$TESTS/13N"
project density "$TESTS/target_density" fit_prior.azr ''
ln -s "$SRC13N/data" "$WORK/density/data"
mkdir -p "$WORK/true/output" "$WORK/true/checks"; ln -s "$SRC13N/data" "$WORK/true/data"
cp "$TESTS/target_density/true.azr" "$WORK/true/t.azr"; azure2 "$WORK/true" t.azr calculate
truth="$WORK/true/output/AZUREOut_aa=1_R=1.out"
if [ -s "$truth" ]; then  # the data at N = 1.1e18 as tests/target_density/check.sh makes them
  mkdir -p "$WORK/density/fitdata"; awk 'NF >= 7' "$truth" > "$WORK/true_rows"; i=0
  for f in meyer_84 meyer_114; do
    awk -v skip="$i" 'NR == FNR { if (FNR > skip) { m[FNR - skip] = $4; d[FNR - skip] = $6 } next }
                      NF >= 4 && $1 >= 0.3 && $1 <= 0.9 { k++; v = m[k] * $3 / d[k]; printf "%s %s %.12e %.12e\n", $1, $2, v, 0.01 * v }' \
      "$WORK/true_rows" "$SRC13N/data/$f.dat" > "$WORK/density/fitdata/$f.dat"
    i=$((i + $(awk '$1 >= 0.3 && $1 <= 0.9' "$SRC13N/data/$f.dat" | grep -c .)))
  done
  continued density tgt_meyer_84.areal_density "areal_density " "atoms/cm2 fit prior"
else bad "(density) true.azr gave no output"; fi

if "$PY" -I "$TESTS/doppler_fgm/make_models.py" "$TESTS/reaction_rate/models/n12C.azr" "$WORK/dop_models" > /dev/null 2>&1; then
  mkdir -p "$WORK/dop450/output" "$WORK/dop450/checks"; ln -s "$WORK/dop_models/data" "$WORK/dop450/data"
  cp "$WORK/dop_models/dop450.azr" "$WORK/dop450/t.azr"; azure2 "$WORK/dop450" t.azr calculate
  if "$PY" -I "$TESTS/doppler_fgm/make_models.py" fitdata "$TESTS/reaction_rate/models/n12C.azr" "$WORK/dop450" \
       "$WORK/dop_models" > /dev/null 2>&1; then
    project doppler "$WORK/dop_models" fit_prior.azr ''
    continued doppler D.doppler "free_gas" "K fit prior"
    # (4) a start outside the centre's range
    project outside "$WORK/dop_models" fit_prior.azr 's/free_gas 300 K fit prior 300 20/free_gas 400 K fit prior 300 20/'
    azure2 "$WORK/outside" t.azr calculate
    want="Note: D.doppler starts at 400 K, outside [75, 363] K, the range laid out from its prior's centre 300 K;"
    if grep -qF "$want" "$WORK/outside/t.log"; then ok "(4) a start outside [75, 363] K: laid out from the start, said once"
    else bad "(4) no line '$want'"; fi
  else bad "(doppler) the fit data were not made"; fi
else bad "(doppler) the models were not made"; fi

project scale "$TESTS/stopping_yield/models" fit.azr \
  's#areal_density  5e+17 atoms/cm2 fit$#areal_density  6e+17 atoms/cm2 fixed#; s#areal_density  1e+18 atoms/cm2 fit$#areal_density  9e+17 atoms/cm2 fixed#'
continued scale stopping_scale.TiN.1H "scale" "fit prior"

project kernel "$TESTS/rectangular_resolution/models" fit.azr 's/rectangular width 14 keV fit$/rectangular width 0.014 MeV fit prior 0.014 0.003/'
continued kernel det_fit.energy_resolution.width "rectangular width" "MeV fit prior"

# (5) the pyazr side (needs the compiled module)
if [ -n "${PYAZR_ENGINE_DIR:-}" ]; then
  out="$("$PY" "$HERE/continued.py" "$WORK" 2>&1)"; rc=$?
  echo "$out" | grep -E '^  (ok|FAIL) ' || true
  [ $rc -eq 0 ] || bad "(5) pyazr: $(echo "$out" | grep -E "failure|Error" | tail -1)"
else
  echo "  skip  (5) save_fit and the layouts: PYAZR_ENGINE_DIR is not set"
fi

if [ "$failures" -eq 0 ]; then echo "PASS: a continued fit computes the fit it continues"; else echo "FAIL: $failures check(s)"; exit 1; fi
