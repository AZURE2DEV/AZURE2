#!/usr/bin/env bash
#
# The quantity spin_observable in a reaction, checked by time reversal, its refusals, and its analytic gradient.
#
# The 4He toy has p+t (pair 1) and n+3He (pair 2); the same nine spin observables are computed for t(p,n)3He (data
# sets grid1_2 ... grid1_10) and 3He(n,p)t (grid2_2 ... grid2_10) on matching grids (equal excitation energy, equal
# c.m. angle).  Time reversal (Ohlsen 1972) requires
#   forward b:y      = inverse e1:y      forward t:y      = inverse e2:y
#   forward e1:y     = inverse b:y       forward e2:y     = inverse t:y
#   forward b:y>e1:y = inverse b:y>e1:y  forward b:y>e2:y = inverse t:y>e1:y
#   forward t:y>e1:y = inverse b:y>e2:y  forward t:y>e2:y = inverse t:y>e2:y
# to round-off (the engine reaches 1e-9; the check allows 1e-6).
#
# Refused, each with its message: an operator that does not follow the grammar (the validator, E_SPIN_OPERATOR), a
# Cartesian tensor on a spin-1/2 beam, a rank above 2j (E_REACTION_NOT_ALLOWED: the pairs' spins decide), a particle
# named twice and an exit particle on the entrance side of '>' (the validator).
#
# The analytic gradient (--use-gradient, AZURE_GRAD_CHECK) against finite differences on this project and on its two
# neighbours (spin_observables_p3He, spin_observables_dd): vector and tensor analysing powers, spin correlations,
# polarizations and transfers; worst difference 1e-5 of the largest component.
#
#   ./tests/spin_observables_reaction/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/spin_observables_reaction.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run DIR COMMAND AZR DATA [VAR=VALUE ...]: the project AZR of DIR with the data directory DATA -> DIR/run.log
run() {
  local d="$1" command="$2" azr="$3" data="$4"
  shift 4
  mkdir -p "$d/output" "$d/checks"
  [ -e "$d/data" ] || cp -r "$data" "$d/data"
  (cd "$d" && env "$@" $(guard_command 900) "$AZURE2_BIN" $command "$azr" < /dev/null > run.log 2>&1)
}

# --- time reversal ---
mkdir -p "$WORK/tr"; cp "$HERE/spin_observables_reaction.azr" "$WORK/tr/"
run "$WORK/tr" calculate spin_observables_reaction.azr "$HERE/data"
F="$WORK/tr/output/AZUREOut_aa=1_R=2.out"; R="$WORK/tr/output/AZUREOut_aa=2_R=1.out"
if [ ! -s "$F" ] || [ ! -s "$R" ]; then
  bad "no output (see run.log)"; sed 's/^/    /' "$WORK/tr/run.log" | tail -5
else
  # block b of a file: blocks are separated by blank lines; block 1 is the cross section, 2..10 the operators
  # b:y t:y e1:y e2:y b:y>e1:y b:y>e2:y t:y>e1:y t:y>e2:y b:y*t:y
  pairs="2:4 3:5 4:2 5:3 6:6 7:8 8:7 9:9"
  for p in $pairs; do
    fb=${p%%:*}; rb=${p##*:}
    res=$(awk -v fb="$fb" -v rb="$rb" '
      FNR == 1 { blk = 1; file++ }
      NF == 0 { if (prev) blk++; prev = 0; next }
      { prev = 1 }
      file == 1 && blk == fb { f[++nf] = $4 }
      file == 2 && blk == rb { r[++nr] = $4 }
      END { if (nf == 0 || nf != nr) { print "count " nf " " nr; exit }
            m = 0; a = 0
            for (i = 1; i <= nf; i++) { d = f[i] - r[i]; if (d < 0) d = -d; if (d > m) m = d
                                        x = f[i]; if (x < 0) x = -x; if (x > a) a = x }
            printf "%.3e %.3f\n", m, a }' "$F" "$R")
    set -- $res
    if [ "$1" = count ]; then bad "forward block $fb / inverse block $rb: point counts $2 $3"
    elif awk -v m="$1" -v a="$2" 'BEGIN { exit !(m ~ /^ *[-+]?\.?[0-9]/ && a ~ /^ *[-+]?\.?[0-9]/ && m + 0 < 1e-6 && a + 0 > 0.05) }'; then
      ok "forward block $fb = inverse block $rb (max |diff| $1, max |value| $2)"
    else bad "forward block $fb vs inverse block $rb: max |diff| $1, max |value| $2"
    fi
  done
fi

# --- refusals ---
# refuse LABEL OPERATOR MESSAGE: the first operator of the project (b:y of grid1_2) replaced by OPERATOR
refuse() {
  local d="$WORK/bad_$1"; mkdir -p "$d"
  awk -v new="$2" '!done && $1 == "operator" { sub(/b:y$/, new); done = 1 } { print }' \
    "$HERE/spin_observables_reaction.azr" > "$d/m.azr"
  grep -qF "operator       $2" "$d/m.azr" || bad "operator $2 was not written into the project"
  run "$d" calculate m.azr "$HERE/data"
  if grep -qF "$3" "$d/run.log" && [ ! -s "$d/output/chiSquared.out" ]; then ok "refused $2 ($3)"
  else bad "operator $2 was not refused with \"$3\""; grep -i "error\|E_" "$d/run.log" | sed -n '1,3p' | sed 's/^/    /'
  fi
}
refuse token   'b:w'        'E_SPIN_OPERATOR'
refuse token2  'b:w'        'unknown operator'
refuse tensor  'b:yy'       'Cartesian tensor operators are defined for spin 1 only'
refuse rank    'b:T20'      'rank 2 exceeds 2j'
refuse twice   'b:y*b:x'    'appears twice'
refuse side    'e1:y>b:y'   'exit particle'

# --- the gradient ---
gradient() {  # NAME: the project tests/NAME/NAME.azr
  local d="$WORK/grad_$1" src="$HERE/../$1"
  mkdir -p "$d"; cp "$src/$1.azr" "$d/"
  run "$d" "fit --use-gradient" "$1.azr" "$src/data" AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_H=1e-7
  awk '/^\[grad-check\] idx / { a = f = ""
         for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1) }
         print a, f }' "$d/run.log" > "$d/grad.rows"
  local nrows worst
  nrows="$(wc -l < "$d/grad.rows" | tr -d ' ')"
  if grep -q "adjoint bailed" "$d/run.log"; then
    bad "$1: the adjoint bailed"; grep "adjoint bailed" "$d/run.log" | sed -n '1,2p'
  elif [ "$nrows" -lt 5 ]; then
    bad "$1: the gradient check listed $nrows parameters"; grep -i "grad\|error" "$d/run.log" | sed -n '1,5p'
  else
    worst="$(awk '{ a = $1 < 0 ? -$1 : $1; if (a > m) m = a; d = $1 - $2; if (d < 0) d = -d; if (d > w) w = d }
               END { printf "%.2e", w / m }' "$d/grad.rows")"
    if awk -v w="$worst" 'BEGIN { exit !(w ~ /^ *[-+]?\.?[0-9]/ && w + 0 <= 1e-5) }'; then
      ok "$1: analytic gradient = finite differences, $nrows parameters, worst $worst of the largest"
    else
      bad "$1: gradient vs finite differences, worst $worst of the largest (> 1e-5)"
    fi
  fi
}
gradient spin_observables_reaction
gradient spin_observables_p3He
gradient spin_observables_dd

[ "$failures" -eq 0 ] && echo "spin_observables_reaction: all checks passed" || { echo "spin_observables_reaction: $failures failure(s)"; exit 1; }
