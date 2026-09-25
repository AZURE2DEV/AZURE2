#!/usr/bin/env bash
#
# Every key of the optional <thm> block (Config::ReadThmBlock) reaches the
# engine and moves the THM observable by the pinned amount, the defaults are
# what a file without the block gets, and a malformed block stops the run.
#
# Model: tests/7Li_p_a (THM segment, entrance pair 5 = 7Li+p with two entrance
# l in one channel spin, 30 keV Gaussian).  Each case copies the project into a
# temporary directory, appends "\n<thm>\n...\n</thm>\n" (the .azr has no
# trailing newline) and reads the total chi2 from output/chiSquared.out.
#
#   ./tests/thm_options/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../7Li_p_a"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_options.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
if command -v timeout >/dev/null 2>&1; then RUN="timeout ${TEST_TIMEOUT:-600}"; else RUN=""; fi

fail=0

# run NAME BLOCK -> prints total chi2 (empty if none); AZURE2's exit status is
# left in $WORK/NAME/status (run is called in a subshell).
# BLOCK is the literal text appended after the last line of the .azr ("" = none).
run() {
  local d="$WORK/$1"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$SRC/data" "$d/"
  cp "$SRC/7Li_p_a.azr" "$d/run.azr"
  [ -z "$2" ] || printf '%b' "$2" >> "$d/run.azr"
  [ "${3:-}" != sorted ] || sort_levels "$d/run.azr"
  (cd "$d" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr > log 2>&1)
  echo $? > "$d/status"
  [ -f "$d/output/chiSquared.out" ] || return 0
  grep -oE 'Total-Chi-Squared: [0-9.eE+-]+' "$d/output/chiSquared.out" | awk '{print $2}'
}
# Rewrite <levels> sorted by (J, parity, E), the order the GUI writes.
sort_levels() {
  python3 - "$1" <<'PY'
import sys
f = sys.argv[1]; t = open(f).read()
a = t.index('<levels>') + len('<levels>\n'); b = t.index('</levels>')
bl = [x for x in t[a:b].split('\n\n') if x.strip()]
bl.sort(key=lambda x: (float(x.split()[0]), int(x.split()[1]), float(x.split()[2])))
open(f, 'w').write(t[:a] + '\n\n'.join(bl) + '\n' + t[b:])
PY
}
block() { printf '\\n<thm>\\n%s\\n</thm>\\n' "$1"; }

# expect NAME ACTUAL EXPECTED [RELTOL]; RELTOL 0 = string-identical output
expect() {
  local tol="${4:-1e-3}"
  if [ -z "$2" ]; then
    echo "  FAIL  $1: no chi2 produced"; fail=1; return
  fi
  if [ "$tol" = 0 ]; then
    if [ "$2" = "$3" ]; then echo "  ok    $1: $2"; else
      echo "  FAIL  $1: $2, expected identical $3"; fail=1; fi
    return
  fi
  if awk -v a="$2" -v b="$3" -v t="$tol" 'BEGIN{d=(a-b)/b; if(d<0)d=-d; exit !(d<=t)}'; then
    echo "  ok    $1: $2 (pin $3)"
  else
    echo "  FAIL  $1: $2, expected $3 (rel $tol)"; fail=1
  fi
}

# (a) defaults: no block == every key written with its default value
none="$(run none "")"
expect "no block (vertex=constant pin)" "$none" 2137.83
defaults="$(run defaults "$(block "entranceL=incoherent
vertex=constant
kinematics=lacognata
coulombIntegral=0
spectatorEnergy=0")")"
expect "explicit defaults == no block" "$defaults" "$none" 0
empty="$(run empty "$(block "# comment only")")"
expect "empty block == no block" "$empty" "$none" 0

# (b) earlier defaults.  vertex=perlevel (alias real) was the default until
# 2026-09-25 (pin 2180.69); with the coherent l sum as well it is the AZURE2
# of before Sep 2026 (legacy pin 2740.48; measured 2740.55, rel 3e-5).
perlevel="$(run perlevel "$(block vertex=perlevel)")"
expect "vertex=perlevel (previous default pin)" "$perlevel" 2180.69
expect "vertex=real == vertex=perlevel" "$(run real "$(block vertex=real)")" "$perlevel" 0
expect "entranceL=coherent" "$(run coherent "$(block entranceL=coherent)")" 3195.37
expect "entranceL=coherent + vertex=perlevel (legacy pin)" \
  "$(run legacy "$(block "entranceL=coherent
vertex=perlevel")")" 2740.48

# (c) kinematic factor conventions
expect "kinematics=lacognata == no block" "$(run lacognata "$(block kinematics=lacognata)")" "$none" 0
expect "kinematics=triple"    "$(run triple    "$(block kinematics=triple)")"    2380.33
expect "kinematics=kf3body"   "$(run kf3body   "$(block kinematics=kf3body)")"   2624.63
expect "kinematics=lambda32"  "$(run lambda32  "$(block kinematics=lambda32)")"  5499.45

# (d) on-shell vertex and the external Coulomb term
expect "vertex=onshell"    "$(run onshell "$(block vertex=onshell)")"    2070.17
expect "coulombIntegral=1" "$(run coulomb "$(block coulombIntegral=1)")" 2111.79

# (e) spectator energy: global and per pair (pair 5 is the entrance pair)
expect "spectatorEnergy=0 == no block" "$(run sp0 "$(block spectatorEnergy=0)")" "$none" 0
spg="$(run spg "$(block spectatorEnergy=0.5)")"
expect "spectatorEnergy=0.5" "$spg" 2195.06
expect "spectatorEnergy[5]=0.5 == global 0.5" "$(run sp5 "$(block "spectatorEnergy[5]=0.5")")" "$spg" 0
expect "spectatorEnergy[3]=0.5 (not a THM pair) == no block" \
  "$(run sp3 "$(block "spectatorEnergy[3]=0.5")")" "$none" 0

# the result does not depend on the order of the levels in the file (the
# constant vertex once took B_c from the first level read)
if command -v python3 > /dev/null; then
  expect "levels sorted as the GUI writes them == file order" "$(run sorted "" sorted)" "$none" 0
  expect "vertex=perlevel, sorted == file order" \
    "$(run plsorted "$(block vertex=perlevel)" sorted)" "$(run plfile "$(block vertex=perlevel)")" 0
fi

# (f)/(g) malformed blocks: AZURE2 prints an ERROR, exits non-zero, writes nothing
refuse() {   # refuse NAME TEXT PATTERN
  local c; c="$(run "$1" "$2")"
  local st; st="$(cat "$WORK/$1/status")"
  if [ -z "$c" ] && [ "$st" -ne 0 ] && grep -q "$3" "$WORK/$1/log"; then
    echo "  ok    $1 refused (exit $st): $(grep -m1 'ERROR' "$WORK/$1/log")"
  else
    echo "  FAIL  $1: exit $st, chi2 '${c}', log:"; sed 's/^/        /' "$WORK/$1/log" | tail -5
    fail=1
  fi
}
refuse misspelt     "$(block vertx=onshell)"             "ERROR: <thm> line not understood: 'vertx=onshell'"
refuse badvalue     "$(block kinematics=kf2body)"        "ERROR: <thm> line not understood"
refuse negative     "$(block spectatorEnergy=-0.1)"      "ERROR: <thm> line not understood"
refuse unterminated "\n<thm>\nvertex=onshell\n"          "ERROR: <thm> block is not terminated"

if [ "$fail" -eq 0 ]; then
  echo "  ok    all <thm> options behave as pinned"; exit 0
else
  echo "  FAIL  <thm> options"; exit 1
fi
