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
# does not.  It also checks that the project's shared beam record folds both
# data sets as one beam record per data set would, that a target on a THM data
# set is refused (and the old <targetInt> line by the converter), and that an
# experiment that only inactive data sets name is ignored.  Pure bash + awk.
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
  (cd "$d" && $RUN "$AZURE2_BIN" calculate --no-gui run.azr > log 2>&1)
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

le() { awk -v a="$1" -v b="$2" 'BEGIN{exit !(a ~ /^ *[-+]?\.?[0-9]/ && a <= b)}'; }

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

# One beam record for both data sets (beam_targetInt1, as the project has it) is the same as one beam record per data
# set.  The shared effect's sigma was converted to the c.m. once per member, so data set 2 was folded with 17 x 0.947
# keV instead of 17 keV.
run split ""
# Byte-identical files, without cmp (not in every MSYS2 image); \r stripped so a
# CRLF-writing Windows build compares like the others.
same() { [ "$(tr -d '\r' < "$1")" = "$(tr -d '\r' < "$2")" ]; }
awk '$1 == "dataset" { ds = $2 } /^<\/data>/ { ds = "" }
     ds == "lc723_band_mid" && $1 == "beam" { sub(/beam_targetInt1/, "beam_band") }
     /^<\/targets>/ { print "beam beam_band"; print "    spread         0.017952 MeV lab" }
     { print }' "$WORK/split/run.azr" > "$WORK/split/run.azr.tmp" && mv "$WORK/split/run.azr.tmp" "$WORK/split/run.azr"
(cd "$WORK/split" && rm -rf output && mkdir output &&
   $RUN "$AZURE2_BIN" calculate --no-gui run.azr > log 2>&1)
if [ "$(grep -c '^beam beam_' "$WORK/split/run.azr")" -eq 2 ] && [ -s "$WORK/split/output/chiSquared.out" ] &&
   same "$WORK/default/output/chiSquared.out" "$WORK/split/output/chiSquared.out" &&
   same "$WORK/default/output/AZUREOut_aa=1_R=2.out" "$WORK/split/output/AZUREOut_aa=1_R=2.out"; then
  echo "  ok    one beam record for both data sets == one per data set (identical output)"
else
  echo "  FAIL  one beam record for both data sets differs from one per data set:"
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

# Target integration does not apply to a THM data set: its energies are the
# reconstructed x + A energies, and its sub-points go below the entrance
# threshold, where the stopping power was evaluated at E < 0 (chi2 = nan).
# Format 3 refuses a target on a THM data set when the file is read; the old
# line (a <targetInt> line on a THM segment) is refused by the converter with
# the old engine's reason.
d="$WORK/thick"
mkdir -p "$d/output" "$d/checks"
cp -r "$HERE/data" "$d/"
awk '$1 == "dataset" { ds = $2 }
     ds == "lc723_thm_points" && $1 == "thm" { print "    target         tgt_thick" }
     /^<\/targets>/ { print "target tgt_thick"; print "    kind           layer"; print "    active         18O"
                      print "    areal_density  1e+21 atoms/cm2 fixed"; print "    stopping       formula \"a0/sqrt(x)\" a0=1e-14 MeV*cm2/atom" }
     { print }' "$HERE/18O_p_a_thm.azr" > "$d/run.azr"
(cd "$d" && $RUN "$AZURE2_BIN" calculate --no-gui run.azr 2>&1 | head -c 1000000 > log)
if grep -q "E_FIELD_NOT_APPLICABLE: dataset 'lc723_thm_points', field 'target' does not belong here: .* only when thm is not given" "$d/log" &&
   [ ! -f "$d/output/chiSquared.out" ]; then
  echo "  ok    a target on a THM data set refused"
else
  echo "  FAIL  a target on a THM data set not refused"; tail -3 "$d/log" | sed 's/^/        /'; fail=1
fi
ORIG="$HERE/../convert/fixtures/tests/18O_p_a_thm/18O_p_a_thm.azr"
awk '/<targetInt>/ { print; print "1  \"1\"  20  0  0.0001  1  1e+21  \"a0/sqrt(x)\"  1  1e-14  0  0  0  \"\"  0  0  0.8  5  50"; T = 1; next }
     /<\/targetInt>/ { T = 0 } !T { print }' "$ORIG" > "$d/old.azr"
(cd "$d" && "$AZURE2_BIN" --convert old.azr conv.azr --to 3 > conv.log 2>&1; echo $? > conv.status)
if [ "$(cat "$d/conv.status")" != 0 ] && [ ! -f "$d/conv.azr" ] &&
   grep -q "data line 1 is a THM segment with target integration, which the old engine refused" "$d/conv.log"; then
  echo "  ok    the old <targetInt> line on a THM segment refused by the converter"
else
  echo "  FAIL  the old <targetInt> line on a THM segment not refused by the converter"; tail -3 "$d/conv.log" | sed 's/^/        /'; fail=1
fi

# An experiment that only inactive data sets name has no member: it is ignored with a warning, as the old engine
# ignored an experiment with no segment in use, and the run equals the run without it -- also when no <thm> section
# defines it (idle_undefined; an active data set naming an undefined experiment is refused, tests/thm_experiment).
for v in idle idle_undefined plain; do
  d="$WORK/$v"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$HERE/data" "$d/"
  awk '$1 == "dataset" { ds = $2 }
       ds == "lc723_thm_points" && $1 == "thm" { print "    active         no"; if (V != "plain") { print "    thm            A"; next } }
       { print }' V="$v" "$HERE/18O_p_a_thm.azr" > "$d/run.azr"
  [ "$v" != idle ] || printf '<thm>\nexperiment[A] background=linear\n</thm>\n' >> "$d/run.azr"
  (cd "$d" && $RUN "$AZURE2_BIN" calculate --no-gui run.azr > log 2>&1)
done
if grep -q "WARNING: <thm> experiment\[A\] has no data set in use; ignored." "$WORK/idle/log" &&
   [ -s "$WORK/plain/output/chiSquared.out" ] &&
   same "$WORK/idle/output/chiSquared.out" "$WORK/plain/output/chiSquared.out" &&
   same "$WORK/idle/output/AZUREOut_aa=1_R=2.out" "$WORK/plain/output/AZUREOut_aa=1_R=2.out" &&
   same "$WORK/idle_undefined/output/chiSquared.out" "$WORK/plain/output/chiSquared.out" &&
   same "$WORK/idle_undefined/output/AZUREOut_aa=1_R=2.out" "$WORK/plain/output/AZUREOut_aa=1_R=2.out"; then
  echo "  ok    an experiment of inactive data sets only: ignored with a warning, the run as without it"
else
  echo "  FAIL  an experiment of inactive data sets only: not ignored as without it"; tail -3 "$WORK/idle/log" | sed 's/^/        /'; fail=1
fi

if [ "$fail" -eq 0 ]; then
  echo "  ok    18O(p,alpha) THM: La Cognata et al. 2010 reproduced with vertex=constant only"; exit 0
else
  echo "  FAIL  18O(p,alpha) THM reproduction"; exit 1
fi
