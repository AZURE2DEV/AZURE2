#!/usr/bin/env bash
#
# The published THM analysis of 18O(p,alpha)15N is reproduced, and only with
# the right entrance vertex.
#
# La Cognata, Spitaleri & Mukhamedzhanov, ApJ 723 (2010) 1512: two interfering
# 1/2+ levels (660/799 keV), formal parameters of Table 3 with B_c = S_c(E_1),
# 17 keV Gaussian resolution.  The project in this directory holds those
# parameters (as Brune parameters; Brune is on, the CLI default); its segment 2
# is the digitized band mid-line of their Fig. 4 on a 5 keV grid.
#
# The check takes AZURE2's folded HOES curve of segment 2 (fit column of
# output/AZUREOut_aa=1_R=2.out), turns it into the S(E) shape of the figure
# (x P_0(E) exp(2 pi eta), from lc723_band.txt), fits ONE normalization N to
# the band mid-line minus the published linear background, S = N*s + bg, and
# compares with the mid-line over 0.505-0.895 MeV:
#
#   default (vertex=constant): rms <= 12 % overall, <= 7 % in 0.56-0.84 MeV
#                              (the reproduction measured 10.1 % / 5.1 %)
#   vertex=perlevel:           peak-region rms >= 15 %
#                              (the reproduction measured 22.9-23.1 %)
#
# i.e. it encodes the verdict of the reproduction: the formal vertex with the
# channel constant B_c reproduces the band, the per-level S_c(E_lambda) vertex
# does not.  It also checks that the project's shared "1,2" <targetInt> line
# folds both segments as one line per segment would.  Pure bash + awk.
#
#   ./tests/18O_p_a_thm/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/18O_p_a_thm.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
. "$HERE/../lib/guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-600}")"

fail=0

# run NAME BLOCK -> leaves $WORK/NAME/output/AZUREOut_aa=1_R=2.out
run() {
  local d="$WORK/$1"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$HERE/data" "$d/"
  cp "$HERE/18O_p_a_thm.azr" "$d/run.azr"
  [ -z "$2" ] || printf '<thm>\n%s\n</thm>\n' "$2" >> "$d/run.azr"
  (cd "$d" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr > log 2>&1)
}

# compare NAME -> prints "rms pkrms max inband N" (empty on a malformed output)
compare() {
  local out="$WORK/$1/output/AZUREOut_aa=1_R=2.out"
  [ -f "$out" ] || return 0
  # segment 2 is the second blank-line-separated block of the output file
  awk '
    FNR == NR { if ($1 !~ /^#/ && NF >= 5) { n++; E[n] = $1; mid[n] = $2; half[n] = $3; bg[n] = $4; cf[n] = $5 } ; next }
    NF == 0 { if (inblk) { blk++; inblk = 0 } ; next }
    { inblk = 1; if (blk == 1) { m++; Eo[m] = $1; fit[m] = $4 } }
    END {
      if (m != n || n == 0) { exit 1 }
      for (i = 1; i <= n; i++) {
        d = Eo[i] - E[i]; if (d < 0) d = -d
        if (d > 1e-3) exit 1
        s[i] = fit[i] * cf[i]; sxy += s[i] * (mid[i] - bg[i]); sxx += s[i] * s[i]
      }
      N = sxy / sxx
      for (i = 1; i <= n; i++) {
        tot = N * s[i] + bg[i]; r = (tot - mid[i]) / mid[i]
        a2 += r * r; ar = (r < 0 ? -r : r); if (ar > mx) mx = ar
        if (E[i] > 0.56 && E[i] < 0.84) { p2 += r * r; np++ }
        d = tot - mid[i]; if (d < 0) d = -d; if (d <= half[i]) inb++
      }
      printf "%.2f %.2f %.2f %.2f %.4g\n", 100 * sqrt(a2 / n), 100 * sqrt(p2 / np), 100 * mx, inb / n, N
    }' "$HERE/lc723_band.txt" "$out"
}

le() { awk -v a="$1" -v b="$2" 'BEGIN{exit !(a <= b)}'; }

run default ""
read -r rms pk mx inb N <<< "$(compare default)"
if [ -z "${rms:-}" ]; then
  echo "  FAIL  default: no usable segment-2 curve"; sed 's/^/        /' "$WORK/default/log" | tail -5; fail=1
else
  echo "  default (vertex=constant): rms $rms %, 0.56-0.84 MeV rms $pk %, max $mx %, in band $inb"
  if le "$rms" 12 && le "$pk" 7; then
    echo "  ok    reproduces the published band (rms <= 12 %, peak rms <= 7 %)"
  else
    echo "  FAIL  published band not reproduced (want rms <= 12 %, peak rms <= 7 %)"; fail=1
  fi
fi

# One <targetInt> line for both segments ("1,2", as the project has it) is the
# same as one line per segment.  The shared effect's sigma was converted to
# the c.m. once per listed segment, so segment 2 was folded with 17 x 0.947
# keV instead of 17 keV.
run split ""
# Byte-identical files, without cmp (not in every MSYS2 image); \r stripped so a
# CRLF-writing Windows build compares like the others.
same() { [ "$(tr -d '\r' < "$1")" = "$(tr -d '\r' < "$2")" ]; }
awk '$2 == "\"1,2\"" { l = $0; sub(/"1,2"/, "\"1\"", l); print l; sub(/"1,2"/, "\"2\""); print; next }
     { print }' "$WORK/split/run.azr" > "$WORK/split/run.azr.tmp" && mv "$WORK/split/run.azr.tmp" "$WORK/split/run.azr"
(cd "$WORK/split" && rm -rf output && mkdir output &&
   printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr > log 2>&1)
if [ "$(grep -c '^1  *"[12]"' "$WORK/split/run.azr")" -eq 2 ] && [ -s "$WORK/split/output/chiSquared.out" ] &&
   same "$WORK/default/output/chiSquared.out" "$WORK/split/output/chiSquared.out" &&
   same "$WORK/default/output/AZUREOut_aa=1_R=2.out" "$WORK/split/output/AZUREOut_aa=1_R=2.out"; then
  echo "  ok    targetInt \"1,2\" == one line per segment (identical output)"
else
  echo "  FAIL  targetInt \"1,2\" differs from one line per segment:"
  paste "$WORK/default/output/chiSquared.out" "$WORK/split/output/chiSquared.out" 2>/dev/null | sed 's/^/        /'
  fail=1
fi

run perlevel "vertex=perlevel"
read -r rms pk mx inb N <<< "$(compare perlevel)"
if [ -z "${rms:-}" ]; then
  echo "  FAIL  perlevel: no usable segment-2 curve"; sed 's/^/        /' "$WORK/perlevel/log" | tail -5; fail=1
else
  echo "  vertex=perlevel:           rms $rms %, 0.56-0.84 MeV rms $pk %, max $mx %, in band $inb"
  if le 15 "$pk"; then
    echo "  ok    vertex=perlevel misses the band (peak rms >= 15 %)"
  else
    echo "  FAIL  vertex=perlevel is no longer clearly worse (peak rms $pk %, want >= 15 %)"; fail=1
  fi
fi

if [ "$fail" -eq 0 ]; then
  echo "  ok    18O(p,alpha) THM: La Cognata et al. 2010 reproduced with vertex=constant only"; exit 0
else
  echo "  FAIL  18O(p,alpha) THM reproduction"; exit 1
fi
