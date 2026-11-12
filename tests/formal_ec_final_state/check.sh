#!/usr/bin/env bash
#
# External capture to a final state that is not the first level of its J^pi group, in the formal
# parametrization (--no-brune).
#
# The external-capture terms need the final state's OBSERVED energy and its on-resonance (Brune) reduced
# width amplitude, which is what the ANC is made of.  In the formal parametrization the working energy
# and amplitude of a level are those only for the level the boundary conditions were taken from, the first
# of its J^pi group.  Before the fix the formal values were used for every final state, so the capture
# cross section depended on the ORDER of the <levels> lines: with the 1/2- ground state of 13N listed
# after the 1/2- background pole (models/A_swap.azr) the --no-brune 12C(p,g0) cross section was up to
# 9.7 times the one of the same levels in the other order (models/A_orig.azr), while --use-brune and the
# elastic scattering of both runs agreed exactly.
#
# Checked here, for a two-level (A) and a three-level (B) 1/2- group:
#   - capture and elastic scattering do not depend on the order of the levels, in either parametrization;
#   - the formal and the Brune parametrization give the same cross sections for both orders;
#   - --no-transform (formal input parameters, no observed values to go by) says that it cannot correct;
#   - with capture pseudo-data: chi^2 is the same for both orders and both parametrizations, the analytic
#     gradient steps aside (finite differences) for the level order that needs the correction and is still
#     used and right for the other, and a formal fit from the corrected order comes back to the data.
#
#   ./tests/formal_ec_final_state/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/formal_ec_final_state.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

run() {  # NAME MODEL FLAGS
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$2" "$dir/t.azr"
  (cd "$dir" && printf '3\nn\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" $3 --no-gui --no-readline t.azr > run.log 2>&1)
}
xs() { awk 'NF >= 4 { printf "%.12e %.12e\n", $1, $(NF - 1) }' "$1"; }
maxrel() {
  paste "$1" "$2" | awk '
    NF >= 4 { d = $2 - $4; if (d < 0) d = -d; m = $4 < 0 ? -$4 : $4; if (m < 1e-300) m = 1e-300;
              r = d / m; if (r != r || r > worst) worst = r; rows++ }
    END { if (rows == 0) print "nan"; else printf "%.3e\n", worst }'
}
le() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r == r && r + 0 <= t + 0) }'; }
same() {  # LABEL FILE1 FILE2 TOLERANCE
  local r; r="$(maxrel "$2" "$3")"
  if le "$r" "$4"; then ok "$1 (max rel $r)"; else bad "$1: max rel $r > $4"; fi
}

CAP="output/AZUREOut_aa=1_R=2.extrap"; ELA="output/AZUREOut_aa=1_R=1.extrap"
M="$HERE/models"
# The formal run reaches the Brune one through an eigenvalue problem (TransformIn) and its inverse
# (the final state's on-resonance amplitude): agreement to rounding of those, not bit for bit.
TOL=1e-7
for pair in "A_orig A_swap" "B_orig B_mid"; do
  set -- $pair; first="$1"; second="$2"
  missing=0
  for m in "$first" "$second"; do
    for f in use-brune no-brune; do
      run "${m}_$f" "$M/$m.azr" "--$f" || true
      for o in "$CAP" "$ELA"; do
        if [ ! -s "$WORK/${m}_$f/$o" ]; then
          bad "$m --$f: no $(basename "$o") ($(tail -1 "$WORK/${m}_$f/run.log" | cut -c1-80))"; missing=1
        else
          xs "$WORK/${m}_$f/$o" > "$WORK/${m}_$f.$( [ "$o" = "$CAP" ] && echo cap || echo ela )"
        fi
      done
    done
  done
  [ "$missing" -eq 0 ] || continue
  n="$(wc -l < "$WORK/${first}_use-brune.cap")"
  [ "$n" -gt 100 ] || bad "$first: only $n capture points"
  if grep -q "nan" "$WORK/${second}_no-brune.cap"; then bad "$second --no-brune: NaN in the capture cross section"; fi
  same "$second vs $first, --use-brune, capture" "$WORK/${second}_use-brune.cap" "$WORK/${first}_use-brune.cap" 1e-10
  same "$second vs $first, --use-brune, elastic" "$WORK/${second}_use-brune.ela" "$WORK/${first}_use-brune.ela" 1e-10
  same "$second vs $first, --no-brune, capture" "$WORK/${second}_no-brune.cap" "$WORK/${first}_no-brune.cap" "$TOL"
  same "$second vs $first, --no-brune, elastic" "$WORK/${second}_no-brune.ela" "$WORK/${first}_no-brune.ela" "$TOL"
  same "$first, formal vs Brune, capture" "$WORK/${first}_no-brune.cap" "$WORK/${first}_use-brune.cap" "$TOL"
  same "$second, formal vs Brune, capture" "$WORK/${second}_no-brune.cap" "$WORK/${second}_use-brune.cap" "$TOL"
  same "$second, formal vs Brune, elastic" "$WORK/${second}_no-brune.ela" "$WORK/${second}_use-brune.ela" "$TOL"
  # the corrected final state must not be announced as displaced: nothing was fitted
  if grep -q "external-capture final state in the formal parametrization" "$WORK/${second}_no-brune/run.log"; then
    bad "$second --no-brune: a pole-offset warning in an unfitted run"
  else ok "$second --no-brune: no pole-offset warning"; fi
done

# Formal input parameters (--no-transform): no observed values, so no correction; the run must say so
# for the level order that needs one and stay quiet for the other.
run "A_swap_nt" "$M/A_swap.azr" "--no-brune --no-transform" || true
run "A_orig_nt" "$M/A_orig.azr" "--no-brune --no-transform" || true
W="external capture in the formal parametrization without the input parameter transformation"
if grep -q "$W" "$WORK/A_swap_nt/run.log"; then ok "--no-transform, ground state second: warned"; else bad "--no-transform, ground state second: no warning"; fi
if grep -q "$W" "$WORK/A_orig_nt/run.log"; then bad "--no-transform, ground state first: warned without cause"; else ok "--no-transform, ground state first: no warning"; fi
if grep -q "$W" "$WORK/A_swap_no-brune/run.log" 2>/dev/null; then bad "transformed run carries the --no-transform warning"; fi

# Data mode.  Pseudo-data: the Brune capture cross section of A_orig (every fourth point, 5 %), entered in
# the laboratory frame; then the 1/2+ width is moved from 34 to 30 keV so that chi^2 is far from zero.
# chiSquared.out carries six digits, hence the 2e-5.
SRC="$WORK/A_orig_use-brune/$CAP"
if [ ! -s "$SRC" ]; then
  bad "no Brune capture cross section to make the pseudo-data from"
else
  mkdir -p "$WORK/fitdata"
  awk 'NF >= 5 && (++n % 4) == 1 { printf "%.10e 0 %.10e %.10e\n", $1 * (1.00728 + 12) / 12, $(NF - 1), 0.05 * $(NF - 1) }' \
    "$SRC" > "$WORK/fitdata/pseudo.dat"
  fitmodel() {  # NAME MODEL
    local dir="$WORK/$1"
    mkdir -p "$dir/output" "$dir/checks" "$dir/data"
    cp "$WORK/fitdata/pseudo.dat" "$dir/data/"
    awk -v seg="1 1 2 0 5 0 180 0 1 0 0 0 0 0 data/pseudo.dat 0 0" '
      /^<segmentsData>/ { print; print seg; next }
      !moved && $3 == "2.3682" && sub(/ 34000 /, " 30000 ") { moved = 1 }
      { print }' "$2" > "$dir/t.azr"
  }
  total() { awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$1" 2>/dev/null; }
  calc() {  # NAME FLAGS -> total chi^2 of a calculation
    rm -f "$WORK/$1/output/chiSquared.out"
    (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" $2 --no-gui --no-readline t.azr > calc.log 2>&1)
    total "$WORK/$1/output/chiSquared.out"
  }
  close() { awk -v a="$1" -v b="$2" 'BEGIN { if (a == "" || b == "" || b + 0 == 0) exit 1; d = (a - b) / b; if (d < 0) d = -d; exit !(d <= 2e-5) }'; }
  fitmodel fit_orig "$M/A_orig.azr"; fitmodel fit_swap "$M/A_swap.azr"
  ref="$(calc fit_orig --use-brune)"
  if [ -n "$ref" ] && awk -v c="$ref" 'BEGIN { exit !(c > 100) }'; then ok "pseudo-data: chi^2 $ref at the moved start (Brune, ground state first)"
  else bad "pseudo-data: chi^2 '$ref' at the moved start, expected > 100"; fi
  for c in "fit_orig --no-brune" "fit_swap --use-brune" "fit_swap --no-brune"; do
    set -- $c; got="$(calc "$1" "$2")"
    if close "$got" "$ref"; then ok "chi^2 $1 $2 = $got"; else bad "chi^2 $1 $2 = '$got', expected $ref"; fi
  done

  gradcheck() {  # NAME -> the [grad-check] summary line of a formal fit start
    (cd "$WORK/$1" && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
       run_limited 600 "$AZURE2_BIN" --no-brune --use-gradient --no-gui --no-readline t.azr > grad.log 2>&1)
    awk '{ sub(/\r$/, "") } /^\[grad-check\] chi2 / { print; exit }' "$WORK/$1/grad.log"
  }
  line="$(gradcheck fit_swap)"
  case "$line" in
    *"BAILED (all finite differences)"*) ok "ground state second: the analytic gradient steps aside" ;;
    *) bad "ground state second: the analytic gradient was not refused ('$(echo "$line" | cut -c1-90)')" ;;
  esac
  line="$(gradcheck fit_orig)"
  case "$line" in
    *"analytic block used"*"energy 0/"*"gamma 0/"*) ok "ground state first: analytic gradient used, no disagreement with finite differences" ;;
    *) bad "ground state first: analytic gradient ('$(echo "$line" | cut -c1-140)')" ;;
  esac

  # a formal fit from the corrected order: back to the data (the true parameters give chi^2 = 0)
  rm -f "$WORK/fit_swap/output/chiSquared.out"
  (cd "$WORK/fit_swap" && printf '2\n\n\n7\n' | run_limited 900 "$AZURE2_BIN" --no-brune --use-gradient --no-gui --no-readline t.azr > fit.log 2>&1)
  got="$(total "$WORK/fit_swap/output/chiSquared.out")"
  if [ -n "$got" ] && awk -v c="$got" 'BEGIN { exit !(c == c && c < 1) }'; then ok "formal fit, ground state second: chi^2 $ref -> $got"
  else bad "formal fit, ground state second: chi^2 '$got' ($(tail -1 "$WORK/fit_swap/fit.log" | cut -c1-80))"; fi
fi

if [ "$failures" -eq 0 ]; then echo "PASS: formal-parametrization external-capture final states"; else echo "FAIL: $failures check(s)"; exit 1; fi
