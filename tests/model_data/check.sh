#!/usr/bin/env bash
#
# Model data (<modelData>, src/ModelData.cpp, docs/upgrade/project_format3.md section 3.10): T2's level 5
# (tests/upgrade/T2/T2_true.azr, the generator's parameters) against values computed by hand from the widths written
# there (README.md), in the three bases; the term's chi^2 in chiSquared.out; a project with only model data; the
# gradient with the term against finite differences; a fit pulled by a displaced resonance strength; the refusals.
# One function per quantity of the table (test_<quantity>: repo_consistency check 27).
#
#   ./tests/model_data/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
T2="$HERE/../upgrade/T2"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/model_data.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
close() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; exit !(b ~ /^ *[-+]?\.?[0-9]/ && a != "" && d <= t * m) }'; }
check() {  # check LABEL VALUE REFERENCE TOL
  if close "$2" "$3" "$4"; then ok "$1: $2 (reference $3)"; else bad "$1: ${2:-missing}, reference $3 (rel. tol. $4)"; fi
}

# project NAME FILE: $WORK/NAME/m.azr, T2_true.azr (or FILE) with the <modelData> section of $WORK/NAME.md appended
project() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$T2/data" "$dir/data"
  cp "$T2/T2_true.azr" "$dir/m.azr"
  [ -f "$dir.md" ] && cat "$dir.md" >> "$dir/m.azr"
  return 0
}
# run NAME MODE [FLAGS...]: the engine on $WORK/NAME/m.azr (output/, run.log); env vars pass through
run() {
  local name="$1" command="$2"
  shift 2
  (cd "$WORK/$name" && run_limited 600 "$AZURE2_BIN" "$command" --no-gui "$@" m.azr > run.log 2>&1)
}
# model NAME DATUM [column]: the model value (or another column) of a datum in modelData.out
model() {
  awk -v n="$2" -v c="${3:-3}" '$1 == n { print $c; exit }' "$WORK/$1/output/modelData.out"
}
total() { sed -n "s/.*$2: \([^ ]*\).*/\1/p" "$WORK/$1/output/chiSquared.out"; }

# The data: each value is the hand value of README.md plus one error, so that every residual is -1.
cat > "$WORK/t2.md" <<'MD'
<modelData>
datum wg_L5
    quantity   resonance_strength
    level      5
    J          3-
    reaction   15N(p,a)12C
    value      2.7318028 keV
    error      0.1 keV
    source     "hand value 2.6318028 keV + 1 sigma (README.md)"
datum Ex_L5
    quantity   resonance_energy
    level      5
    J          3-
    frame      E_x
    value      13.2702 MeV
    error      1 keV
datum Ecm_L5
    quantity   resonance_energy
    level      5
    J          3-
    frame      E_cm
    pair       p+15N
    value      1142.8 keV
    error      1 keV
datum Elab_L5
    quantity   resonance_energy
    level      5
    J          3-
    frame      E_lab
    pair       p+15N
    value      1.2195150815 MeV
    error      1 keV
datum G_L5
    quantity   total_width
    level      5
    J          3-
    value      15.622 keV
    error      1 keV
datum Gp_L5
    quantity   partial_width
    level      5
    J          3-
    pair       p+15N
    value      1802 eV
    error      100 eV
datum BR_L5
    quantity   branching_ratio
    level      5
    J          3-
    reaction   15N(p,a)12C
    value      89.360005 %
    error      1 %
datum a_p15N
    quantity   scattering_length
    pair       p+15N
    s          1
    value      0 fm
    error      100 fm
</modelData>
MD
project t2
run t2 calculate
project formal; cp "$WORK/t2/m.azr" "$WORK/formal/m.azr"; run formal calculate --no-brune
project park; cp "$WORK/t2/m.azr" "$WORK/park/m.azr"; run park calculate --use-park

# Each quantity in the three bases: the hand value (Brune), and the same value in the formal and Park bases, to the
# tolerance of the formal -> Brune iteration (1e-6 MeV in the energies, level_parameters.md section 4).
basis_check() {  # basis_check LABEL DATUM HAND TOL
  check "$1 (Brune), hand value" "$(model t2 "$2")" "$3" "$4"
  check "$1, formal basis = Brune" "$(model formal "$2")" "$(model t2 "$2")" 2e-6
  check "$1, Park basis = Brune" "$(model park "$2")" "$(model t2 "$2")" 2e-6
}
test_resonance_strength() {
  # omega gamma = (7/4) 1.702 x 12.92 / 14.622 keV (omega = 7 / (2 x 2) for p + 15N)
  basis_check "omega gamma of level 5 (keV)" wg_L5 2.6318028 5e-7
  grep -q "W_BROAD_LEVEL.*wg_L5" "$WORK/t2/run.log" && ok "W_BROAD_LEVEL for wg_L5 (Gamma / E_r = 1.3 %)" || bad "no W_BROAD_LEVEL for wg_L5"
}
test_resonance_energy() {
  basis_check "E_x of level 5 (MeV)" Ex_L5 13.2692 1e-9
  basis_check "E_cm of level 5 (keV): 13.2692 - 12.1274 MeV" Ecm_L5 1141.8 1e-9
  basis_check "E_lab of level 5 (MeV): E_cm (m_p + m_15N)/m_15N" Elab_L5 1.2185150815 1e-9
}
test_total_width() { basis_check "Gamma of level 5 (keV): 1.702 + 12.92 + 0" G_L5 14.622 1e-9; }
test_partial_width() { basis_check "Gamma_p of level 5 (eV)" Gp_L5 1702 1e-9; }
test_branching_ratio() { basis_check "Gamma_a0 / Gamma of level 5 (%)" BR_L5 88.360005 1e-7; }
test_scattering_length() {
  # The same s-wave length as scatteringLengths.out (p + 15N, 1-); its chi^2 (value 0, error 100 fm) is (a0/100)^2.
  local a0
  a0="$(awk '$1 == 1 && $3 == "1-" { print $5; exit }' "$WORK/t2/output/scatteringLengths.out")"
  check "scattering length p+15N s = 1 = scatteringLengths.out" "$(model t2 a_p15N)" "$a0" 1e-6
  check "its chi^2 in modelData.out = (a0 / 100 fm)^2" "$(awk '$1 == "a_p15N" { print $NF; exit }' "$WORK/t2/output/modelData.out")" \
    "$(awk -v a="$a0" 'BEGIN { printf "%.10g", (a / 100) ^ 2 }')" 1e-4
}
test_anc() {
  # Level 1 moved below the p + 15N threshold (12.0 MeV): its p channel (l = 0, s = 1) is closed and the value of
  # <levels> there is its ANC, 2.0 fm^-1/2, the value the datum must give back in each basis.  Its alpha widths are
  # made small (100 eV, 0): with T2's 80.9 keV and 50 eV no real reduced width reproduces them at 12.0 MeV (the
  # engine warns "Denominator less than zero" and reports another ANC).
  cat > "$WORK/anc.md" <<'MD'
<modelData>
datum C_L1
    quantity   anc
    level      1
    J          1-
    pair       p+15N
    l          0
    s          1
    value      2.04 fm^-1/2
    error      0.04 fm^-1/2
</modelData>
MD
  local b
  for b in "" --no-brune --use-park; do
    project "anc$b"
    sed 's/^1 -1 12.420800 1 1 1 2 0 1 1 1 496 /1 -1 12.420800 1 1 1 2 0 1 1 1 2.0 /
         s/^1 -1 12.420800 1 1 2 0 2 1 1 1 80900 /1 -1 12.420800 1 1 2 0 2 1 1 1 100 /
         s/^1 -1 12.420800 1 1 3 4 2 1 1 1 50 /1 -1 12.420800 1 1 3 4 2 1 1 1 0 /
         s/^1 -1 12.420800 /1 -1 12.000000 /' "$T2/T2_true.azr" > "$WORK/anc$b/m.azr"
    cat "$WORK/anc.md" >> "$WORK/anc$b/m.azr"
    run "anc$b" calculate $b
  done
  check "ANC of level 1 (fm^-1/2), Brune, as written" "$(model anc C_L1)" 2.0 1e-9
  check "ANC, formal basis = Brune" "$(model anc--no-brune C_L1)" "$(model anc C_L1)" 2e-6
  check "ANC, Park basis = Brune" "$(model anc--use-park C_L1)" "$(model anc C_L1)" 2e-6
}

echo "-- every quantity, three bases"
test_resonance_strength
test_resonance_energy
test_total_width
test_partial_width
test_branching_ratio
test_scattering_length
test_anc

echo "-- the term in the totals"
# Seven level data with residual -1 and the scattering length, (a0/100)^2 of it.
want="$(awk -v a="$(model t2 a_p15N)" 'BEGIN { printf "%.10g", 7 + (a / 100) ^ 2 }')"
check "Total-Model-Data-Chi-Squared = 7 + (a0/100 fm)^2" "$(total t2 Total-Model-Data-Chi-Squared)" "$want" 2e-6
check "modelData.out total = chiSquared.out's" "$(awk '/^Total model-data chi2/ { print $NF }' "$WORK/t2/output/modelData.out")" \
  "$(total t2 Total-Model-Data-Chi-Squared)" 1e-6
grep -q "^Model data (<modelData>)" "$WORK/t2/run.log" && ok "the table is in the run log" || bad "the table is not in the run log"

echo "-- a project with only model data (no segments)"
project only
# No data sets, and nothing that refers to them (the THM experiments, the targets, beams and detectors, the prior
# covariance).
awk '/^<(data|thm|targets)>/ { print; skip = 1; next } /^<\/(data|thm|targets)>/ { skip = 0 } !skip' \
  "$T2/T2_true.azr" | grep -v priorCovariance > "$WORK/only/m.azr"
cat "$WORK/t2.md" >> "$WORK/only/m.azr"
run only calculate
check "model data alone: omega gamma as with segments" "$(model only wg_L5)" "$(model t2 wg_L5)" 1e-12
check "model data alone: Total-Model-Data-Chi-Squared" "$(total only Total-Model-Data-Chi-Squared)" "$(total t2 Total-Model-Data-Chi-Squared)" 1e-12

echo "-- gradient with the term (AZURE_GRAD_CHECK) and a fit pulled by a displaced strength"
# The fit project: data set 5 (the implanted target across level 5) alone, level 5's energy and alpha0 width free.
# The THM experiment goes with its data sets (an experiment whose data sets are all inactive would be ignored with a
# warning).
fitproj() {  # fitproj NAME WG_VALUE: the fit project with a strength datum of value WG_VALUE (keV), or none
  project "$1"
  awk '$1 == "dataset" { n++; ds = 1 }
       ds && (/^$/ || /^<\/data>/) { if (n != 5) print "    active         no"; ds = 0 }
       /^experiment\[/ { next }
       /^3 -1 13.269200 / { $4 = 0; if ($6 == 2) $11 = 0 }
       { print }' "$T2/T2_true.azr" > "$WORK/$1/m.azr"
  if [ -n "$2" ]; then
    printf '<modelData>\ndatum wg\n    quantity resonance_strength\n    level 5\n    J 3-\n    reaction 15N(p,a)12C\n    value %s keV\n    error 0.02 keV\n</modelData>\n' "$2" >> "$WORK/$1/m.azr"
  fi
}
# The term's share of the gradient: AZURE_GRAD_CHECK's analytic and finite-difference gradients with the datum minus
# those without it (the data part's own finite-difference error, 5e-5 relative in E here, cancels).
# The tolerance is 5e-6: the finite-difference side is the difference of two finite differences and moves with the
# machine.  --no-brune, parameter 23, analytic / finite differences: -3138.73784458 / -3138.73798362 (4e-8) with the
# Coulomb table file, -3138.73840766 / -3138.7372511 (3.7e-7) without it, -3138.73552009 / -3138.7389307 (1.1e-6) on
# the CI's macOS machine (Intel, Apple's clang); the other five parameters agree to 5e-8 or better on all three.
gradline() { awk -v i="$2" '/\[grad-check\] idx/ && $3 == i { for (k = 1; k <= NF; k++) if ($k == "analytic") a = $(k + 1); else if ($k == "fd") f = $(k + 1); print a, f; exit }' "$WORK/$1/run.log"; }
for b in "" --no-brune --use-park; do
  fitproj "gradw$b" 2.6918028
  fitproj "grad0$b" ""
  AZURE_GRAD_CHECK=exit run "gradw$b" fit --use-gradient $b
  AZURE_GRAD_CHECK=exit run "grad0$b" fit --use-gradient $b
  n="$(grep -c '\[grad-check\] idx' "$WORK/gradw$b/run.log")"
  [ "$n" -ge 2 ] || bad "gradient ${b:---use-brune}: $n parameters checked"
  for i in $(awk '/\[grad-check\] idx/ { print $3 }' "$WORK/gradw$b/run.log"); do
    read -r aw fw <<< "$(gradline "gradw$b" "$i")"
    read -r a0 f0 <<< "$(gradline "grad0$b" "$i")"
    check "gradient ${b:---use-brune}, parameter $i: the term's share, analytic = finite differences" \
      "$(awk -v x="$aw" -v y="$a0" 'BEGIN { printf "%.12g", x - y }')" "$(awk -v x="$fw" -v y="$f0" 'BEGIN { printf "%.12g", x - y }')" 5e-6
  done
done
fitproj free ""
fitproj pulled 2.6918028  # + 3 sigma
run free fit
run pulled fit
GA() {  # the fitted alpha0 width of level 5 (eV), from parameters.out
  awk '/^J = 3.0-/ { lv = ($6 + 0 > 13.2) } lv && $1 == "R" && $3 == 2 {
         f = ($13 == "keV") ? 1e3 : (($13 == "MeV") ? 1e6 : (($13 == "meV") ? 1e-3 : 1)); print $12 * f; exit }' \
    "$WORK/$1/output/parameters.out"
}
a_free="$(GA free)"; a_pulled="$(GA pulled)"
if awk -v a="$a_free" -v b="$a_pulled" 'BEGIN { exit !(a != "" && b != "" && b > a) }'; then
  ok "Gamma_a0 of level 5 pulled up by omega gamma + 3 sigma: $a_free -> $a_pulled eV (model omega gamma $(model pulled wg) keV)"
else bad "Gamma_a0 of level 5 not pulled up: ${a_free:-missing} -> ${a_pulled:-missing}"; fi
t="$(total pulled Total-Model-Data-Chi-Squared)"
if [ -n "$t" ]; then ok "the pulled fit reports Total-Model-Data-Chi-Squared $t"; else bad "no Total-Model-Data-Chi-Squared in the pulled fit"; fi
# The same fit by Levenberg-Marquardt (the term's rows and their Jacobian).
fitproj lm 2.6918028
run lm fit --use-lm
a_lm="$(GA lm)"
if awk -v a="$a_free" -v b="$a_lm" 'BEGIN { exit !(a != "" && b != "" && b > a) }'; then ok "--use-lm: Gamma_a0 pulled up as well: $a_free -> $a_lm eV"
else bad "--use-lm: Gamma_a0 ${a_free:-missing} -> ${a_lm:-missing}"; fi

echo "-- refusals"
refuse() {  # refuse LABEL CODE MD [SED]: the run stops with the code in an ERROR line and writes no chiSquared.out
  local dir="$WORK/refuse_$2"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$T2/data" "$dir/data"
  sed "${4:-}" "$T2/T2_true.azr" > "$dir/m.azr"
  printf '%b' "$3" >> "$dir/m.azr"
  (cd "$dir" && run_limited 120 "$AZURE2_BIN" calculate --no-gui m.azr > run.log 2>&1)
  if grep -q "^ERROR: m.azr:[0-9]*: wg[.A-Za-z_]*: error $2: datum 'wg'" "$dir/run.log" && [ ! -s "$dir/output/chiSquared.out" ]; then ok "refused: $1 ($2)"
  else bad "not refused: $1 ($2)"; grep -i "error" "$dir/run.log" | sed -n '1,3p'; fi
}
wg() { printf '<modelData>\\ndatum wg\\n    quantity resonance_strength\\n    level %s\\n    J %s\\n    reaction %s\\n    value 2.6 keV\\n    error 0.1 keV\\n%s</modelData>\\n' "$1" "$2" "$3" "${4:-}"; }
refuse "wrong J^pi for the level ID" E_LEVEL_JPI "$(wg 5 2+ '15N(p,a)12C')"
refuse "unknown level ID" E_UNKNOWN_LEVEL "$(wg 9 3- '15N(p,a)12C')"
refuse "strength of a level below the entrance threshold" E_OUT_OF_RANGE "$(wg 1 1- '15N(p,a)12C')" 's/^1 -1 12.420800 /1 -1 12.000000 /'
refuse "a reaction no pair has" E_UNKNOWN_PAIR "$(wg 5 3- '15N(p,d)14N')"
# A fourth pair, d + 14N, given only to level 1 (1-): no channel of level 5 (3-).
refuse "an entrance pair that is no channel of the level" E_NOT_A_CHANNEL "$(wg 5 3- '14N(d,a)12C')" \
  '/^1 -1 12.420800 1 1 3 /a\
1 -1 12.420800 1 1 4 2 2 1 1 1 0.1 1.0 1 1.0 1 0.0 2.014102 14.003074 1 7 12.1274 20.7354 0 0 0.0 0 5.0 0 0 0'
refuse "two data of one name" E_DUPLICATE_NAME "$(wg 5 3- '15N(p,a)12C' 'datum wg\n    quantity total_width\n    level 5\n    J 3-\n    value 14 keV\n    error 1 keV\n')"

if [ "$failures" -eq 0 ]; then echo "PASS model_data"; exit 0; fi
echo "FAIL model_data ($failures)"; exit 1
