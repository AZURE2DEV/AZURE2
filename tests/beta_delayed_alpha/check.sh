#!/usr/bin/env bash
#
# What a beta-delayed spectrum (entrance pair of type 20) may carry (docs/upgrade/data_model.md S2).  Its energies are
# the light decay particle's (lab), converted with the decay pair's masses, so
#   (a) a data set without its decay pair as the exit (16N(b-,g)16C[total], the old exit -1 with data type 3) is
#       refused at load: its energies would be converted with the beta pair's masses, as a beam;
#   (b) a target (a target integration), (c) a beam (a beam profile) and (d) a target with electron screening are
#       refused (E_FIELD_NOT_APPLICABLE): a decay has no beam and the beta pair's charges are not a reaction's;
#   (e) a Gaussian resolution (old grammar only, see below) and (f) an energy shift of the measured particle energy run (a detector
#       resolution and a calibration; the sub-points move by the shift times the decay pair's c.m. factor).
#
#   ./tests/beta_delayed_alpha/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/beta_delayed_alpha.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
# variant NAME FIELDS [RECORDS]: the model with these fields added to its data set bda (printf escapes) and these
# records in <targets>, menu 1.
variant() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -r "$HERE/data" "$WORK/$1/"
  awk -v fields="$(printf "$2")" -v records="$(printf "${3:-}")" '
    /^dataset bda$/ { inbda = 1 }
    inbda && (/^$/ || /^<\/data>/) { if (fields != "") print fields; inbda = 0 }
    { print }
    /^<targets>$/ && records != "" { print records }' "$HERE/16N.azr" > "$WORK/$1/m.azr"
  (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
}
refused() { grep -q "$2" "$WORK/$1/run.log" && [ ! -s "$WORK/$1/output/chiSquared.out" ]; }
ran() { grep -q '^Total' "$WORK/$1/output/chiSquared.out" 2> /dev/null && ! grep -qi nan "$WORK/$1/output/chiSquared.out"; }

# (a) the old segment line with exit -1 and data type 3 converts to the total-capture reaction
mkdir -p "$WORK/a/output" "$WORK/a/checks"
cp -r "$HERE/data" "$WORK/a/"
sed 's/^    reaction       16N(b-,a)12C$/    reaction       16N(b-,g)16C[total]/' "$HERE/16N.azr" > "$WORK/a/m.azr"
(cd "$WORK/a" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
if refused a "needs its decay pair as the exit pair"; then ok "(a) exit -1: refused"
else bad "(a) exit -1: $(grep -m1 ERROR "$WORK/a/run.log")"; fi

NA="E_FIELD_NOT_APPLICABLE: dataset 'bda', field"
variant b '    target         layer' 'target layer\n    kind           layer\n    active         16N\n    areal_density  1e18 atoms/cm2 fixed\n    stopping       formula "10" eV*cm2/1e15atoms'
if refused b "$NA 'target' does not belong here: quantity 'beta_delayed_spectrum' does not use it"; then ok "(b) target integration: refused"
else bad "(b) target integration: $(grep -m1 ERROR "$WORK/b/run.log")"; fi

variant c '    beam           profile\n    quadrature     points 10' 'beam profile\n    profile        skewed_gaussian xi 1.0 MeV omega 0.1 MeV alpha 0 weight 1'
if refused c "$NA 'beam' does not belong here: quantity 'beta_delayed_spectrum' does not use it"; then ok "(c) beam profile: refused"
else bad "(c) beam profile: $(grep -m1 ERROR "$WORK/c/run.log")"; fi

variant d '    target         screened' 'target screened\n    kind           none\n    screening      1000 eV fixed'
if refused d "$NA 'target' does not belong here: quantity 'beta_delayed_spectrum' does not use it"; then ok "(d) electron screening: refused"
else bad "(d) electron screening: $(grep -m1 ERROR "$WORK/d/run.log")"; fi

# (e) is not stated in format 3: a detector's energy_resolution is built only with a beam profile (E_NOT_BUILT), and a
# beam spread does not apply to a decay (E_FIELD_NOT_APPLICABLE); the converter refuses the old line for that reason.
# It runs on the old original (tests/convert/fixtures) with the old <targetInt> line until the engine's old reading
# path goes (docs/upgrade/format3_tree.md section 4).
mkdir -p "$WORK/e/output" "$WORK/e/checks"
cp -r "$HERE/data" "$WORK/e/"
awk '{ print } /^<targetInt>$/ { print "1 \"1\" 20 1 0.02 0 0 \"\" 0 0 0 0 \"\" 0 0 0.04 20 50" }' \
  "$HERE/../convert/fixtures/tests/beta_delayed_alpha/16N.azr" > "$WORK/e/m.azr"
(cd "$WORK/e" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
if ran e; then ok "(e) Gaussian resolution: runs ($(grep '^Total' "$WORK/e/output/chiSquared.out"))"
else bad "(e) Gaussian resolution: $(grep -m1 ERROR "$WORK/e/run.log")"; fi

variant f '    shift          offset 0.01 MeV fixed'
if ran f; then ok "(f) energy shift: runs ($(grep '^Total' "$WORK/f/output/chiSquared.out"))"
else bad "(f) energy shift: $(grep -m1 ERROR "$WORK/f/run.log")"; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
