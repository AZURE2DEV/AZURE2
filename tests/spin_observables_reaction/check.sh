#!/usr/bin/env bash
#
# isDiff 9 in a reaction, checked by time reversal.  The 4He toy has p+t (pair 1) and n+3He (pair 2);
# the same nine spin observables are computed for t(p,n)3He (segments 2-10) and 3He(n,p)t (segments 12-20)
# on matching grids (equal excitation energy, equal c.m. angle).  Time reversal (Ohlsen 1972) requires
#   forward b:y      = inverse e1:y      forward t:y      = inverse e2:y
#   forward e1:y     = inverse b:y       forward e2:y     = inverse t:y
#   forward b:y>e1:y = inverse b:y>e1:y  forward b:y>e2:y = inverse t:y>e1:y
#   forward t:y>e1:y = inverse b:y>e2:y  forward t:y>e2:y = inverse t:y>e2:y
# to round-off (the engine reaches 1e-9; the check allows 1e-6).  It also checks that start-up refuses
# tokens the pairs cannot carry, with a message: a bad token, a target operator on a spin-0 target, a
# Cartesian tensor on a spin-1/2 beam, an exit operator of the wrong particle.
#
#   ./tests/spin_observables_reaction/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/spin_observables_reaction.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

run() {  # dir azr -> log
  mkdir -p "$1/output" "$1/checks"; ln -sf "$HERE/data" "$1/data"
  (cd "$1" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline "$2" > run.log 2>&1)
}

# --- time reversal ---
mkdir -p "$WORK/tr"; cp "$HERE/spin_observables_reaction.azr" "$WORK/tr/"
run "$WORK/tr" spin_observables_reaction.azr
F="$WORK/tr/output/AZUREOut_aa=1_R=2.out"; R="$WORK/tr/output/AZUREOut_aa=2_R=1.out"
if [ ! -s "$F" ] || [ ! -s "$R" ]; then
  bad "no output (see run.log)"; sed 's/^/    /' "$WORK/tr/run.log" | tail -5
else
  # block b of a file: blocks are separated by blank lines; block 1 is the cross section, 2..10 the tokens
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
    elif awk -v m="$1" -v a="$2" 'BEGIN { exit !(m < 1e-6 && a > 0.05) }'; then
      ok "forward block $fb = inverse block $rb (max |diff| $1, max |value| $2)"
    else bad "forward block $fb vs inverse block $rb: max |diff| $1, max |value| $2"
    fi
  done
fi

# --- start-up refusals ---
refuse() {  # label, sed expression applied to the first isDiff 9 line, expected message fragment
  local d="$WORK/bad_$1"; mkdir -p "$d"
  sed "0,/ 9 b:y /s// 9 $2 /" "$HERE/spin_observables_reaction.azr" > "$d/m.azr"
  run "$d" m.azr
  if grep -q "ERROR: Data segment #2: .*$3" "$d/run.log"; then ok "refused $2 ($3)"
  else bad "token $2 was not refused with \"$3\""; sed 's/^/    /' "$d/run.log" | grep -i error | head -3
  fi
}
refuse token   'b:w'        'unknown operator'
refuse tensor  'b:yy'       'Cartesian tensor operators are defined for spin 1 only'
refuse rank    'b:T20'      'rank 2 exceeds 2j'
refuse twice   'b:y*b:x'    'appears twice'
refuse side    'e1:y>b:y'   "exit particle"

[ "$failures" -eq 0 ] && echo "spin_observables_reaction: all checks passed" || { echo "spin_observables_reaction: $failures failure(s)"; exit 1; }
