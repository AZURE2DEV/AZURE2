#!/usr/bin/env bash
#
# The spin-orbit term of a hybrid pair (docs/upgrade/spin_orbit_potential.md section 4.2; README.md here):
#   1. the p1/2, p3/2 and s1/2 phase shifts of models/phase.azr (levels without widths: the potential's phase) equal
#      -phi of physics/tests/data/spin_orbit.csv (make_spin_orbit_reference.py, independent of the engine) to 1e-6 rad;
#   2. unitarity: n + 16O with the term, sigma_tot (optical theorem) = sigma_el (unitarity.py);
#   3. spin_orbit_potential.azr, models/capture.azr (12C(p,g)13N with external capture) and models/phase.azr at 1
#      and 3 threads: every output file identical;
#   4. the analytic gradient of spin_orbit_potential.azr equals finite differences (AZURE_GRAD_CHECK), and its chi^2
#      at the start equals expected/chiSquared.out;
#   5. the term reaches the capture path: chi^2 of models/capture.azr differs from the same model without it.
# The capture integrals per J against an independent integration: tests/pyazr/spin_orbit_potential_test.py.
#
#   ./tests/spin_orbit_potential/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/spin_orbit_potential.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME MODEL MENU THREADS [env ...] [-- flags]: menu 1 (calculation), 2 (fit) or 3 (extrapolation) in $WORK/NAME
run() {
  local name="$1" model="$2" menu="$3" threads="$4"; shift 4
  local envs=() flags=()
  while [ $# -gt 0 ] && [ "$1" != "--" ]; do envs+=("$1"); shift; done
  [ $# -gt 0 ] && shift && flags=("$@")
  mkdir -p "$WORK/$name/output" "$WORK/$name/checks"
  [ -e "$WORK/$name/data" ] || ln -s "$(dirname "$model")/data" "$WORK/$name/data"
  cp "$model" "$WORK/$name/m.azr"
  (cd "$WORK/$name" && rm -rf output/* checks/* && printf '%s\n\n\n7\n' "$menu" |
     env OMP_NUM_THREADS="$threads" "${envs[@]}" $(guard_command 900) "$AZURE2_BIN" "${flags[@]}" --no-gui --no-readline m.azr > run.log 2>&1)
}

# ---- 1. phase shifts against the independent integration
run phase "$HERE/models/phase.azr" 3 1
if res="$(python3 -I "$HERE/phase_check.py" "$WORK/phase/output/AZUREOut_aa=1_R=1.extrap" \
            "$HERE/../../physics/tests/data/spin_orbit.csv" 2>&1)"; then
  ok "phase shifts = -phi of make_spin_orbit_reference.py ($(echo "$res" | tail -1))"
else
  bad "phase shifts against make_spin_orbit_reference.py:"; echo "$res"
fi

# ---- 2. unitarity
if res="$(python3 -I "$HERE/unitarity.py" "$AZURE2_BIN" "$WORK/unitarity" 2>&1)"; then
  ok "unitarity with the term: $(echo "$res" | tail -1)"
else
  bad "unitarity with the term:"; echo "$res" | tail -5
fi

# ---- 3. 1 against 3 threads
for m in main:"$HERE/spin_orbit_potential.azr":1 capture:"$HERE/models/capture.azr":1 phase:"$HERE/models/phase.azr":3; do
  IFS=: read -r name model menu <<< "$m"
  run "${name}_1" "$model" "$menu" 1
  run "${name}_3" "$model" "$menu" 3
  if [ -z "$(ls "$WORK/${name}_1/output")" ]; then
    bad "$name: no output ($(grep -m1 ERROR "$WORK/${name}_1/run.log"))"
  elif diff -r "$WORK/${name}_1/output" "$WORK/${name}_3/output" > /dev/null; then
    ok "$name: every output file identical at 1 and 3 threads ($(ls "$WORK/${name}_1/output" | wc -l | tr -d ' ') files)"
  else
    bad "$name: outputs differ between 1 and 3 threads"
  fi
done

# ---- 4. gradient against finite differences
run grad "$HERE/spin_orbit_potential.azr" 2 2 AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 -- --use-gradient
chi="$(awk '/^\[grad-check\] chi2 / { sub(",", "", $3); print $3 }' "$WORK/grad/run.log")"
want="$(awk '/^Total-Chi-Squared:/ { print $2 }' "$HERE/expected/chiSquared.out")"
if [ -n "$chi" ] && awk -v a="$chi" -v b="$want" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(d <= 1e-5) }'; then
  ok "chi^2 at the start $chi = expected $want"
else
  bad "chi^2 at the start '$chi', expected $want"
fi
awk '/^\[grad-check\] idx / { a = f = ""
       for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1) }
       print a, f }' "$WORK/grad/run.log" > "$WORK/rows"
n="$(wc -l < "$WORK/rows" | tr -d ' ')"
fd="$(awk '{ a = $1 < 0 ? -$1 : $1; if (a > m) m = a; d = $1 - $2; if (d < 0) d = -d; if (d > w) w = d }
           END { if (m > 0) printf "%.2e", w / m; else print "nan" }' "$WORK/rows")"
if [ "$n" -ge 5 ] && awk -v w="$fd" 'BEGIN { exit !(w ~ /^ *[-+]?\.?[0-9]/ && w + 0 <= 1e-5) }'; then
  ok "gradient = finite differences for $n parameters ($fd of the largest)"
else
  bad "gradient vs finite differences: $n parameters, $fd of the largest"
fi

# ---- 5. the term reaches the capture path
mkdir -p "$WORK/noso" && ln -s "$HERE/models/data" "$WORK/noso/data"
sed '/^    spin_orbit /d' "$HERE/models/capture.azr" > "$WORK/noso/capture.azr"
run capture_noso "$WORK/noso/capture.azr" 1 1
a="$(awk '/^Total-Chi-Squared:/ { print $2 }' "$WORK/capture_1/output/chiSquared.out")"
b="$(awk '/^Total-Chi-Squared:/ { print $2 }' "$WORK/capture_noso/output/chiSquared.out")"
if [ -n "$a" ] && [ -n "$b" ] && [ "$a" != "$b" ]; then
  ok "capture: chi^2 $a with the term, $b without"
else
  bad "capture: chi^2 '$a' with the term, '$b' without"
fi

if [ "$failures" -eq 0 ]; then echo "spin_orbit_potential: all checks passed"; else echo "spin_orbit_potential: $failures check(s) FAILED"; fi
exit $(( failures > 0 ))
