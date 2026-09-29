#!/usr/bin/env bash
#
# THM experiments: `experiment[<name>] key=value ...` lines of the <thm> block
# (ThmExperiment.h; docs/source/theory/thm_implementation.rst, "THM
# experiments").  The segments of an experiment share one profiled norm and an
# optional background b0 + b1 E + b2 E^2, both eliminated by closed-form
# linear least squares.
#
# Model: tests/18O_p_a_thm (two THM segments with free norms, 39 + 79 points,
# all parameters fixed; 0.1 s a run).  Each case appends a <thm> block to a
# temporary copy of the .azr.
#
#   (a) one-segment experiments without background change nothing: every
#       output file is byte-identical to the run without the block;
#   (b) two segments sharing one norm: n* and chi2 equal a hand computation
#       from the per-segment sums S_mm, S_md, S_dd (models from run (a), data
#       from the data files) to 1e-8, and both segments carry that n*;
#   (c) background: data made as s m + a0 + a1 E (m the model of run (a),
#       E the c.m. energy) are fitted with n = 1/s, b_k = a_k/s recovered to
#       1e-8 and chi2 ~ 0, and the fitted curve written is the model plus b(E);
#   (d) what AZURE2 refuses (ERROR line, non-zero exit): the reserved keys,
#       an unknown key or nuclide, a malformed value, partial kinematics, a
#       segment in two experiments, a segment that is not THM or has a fixed
#       norm, kinematics that do not make the entrance pair.
#
#   ./tests/thm_experiment/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../18O_p_a_thm"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_experiment.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
if command -v timeout >/dev/null 2>&1; then RUN="timeout ${TEST_TIMEOUT:-300}"; else RUN=""; fi
OUT="AZUREOut_aa=1_R=2.out"
fail=0

# run NAME BLOCK [AWK] -- a "calculate with data" run of a copy of the project
# with BLOCK (literal lines) as its <thm> block ("" = none), the .azr first
# passed through the awk program AWK if given.  Exit status in NAME/status.
run() {
  local d="$WORK/$1"
  rm -rf "$d"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$SRC/data" "$d/"
  if [ -n "${3:-}" ]; then awk "$3" "$SRC/18O_p_a_thm.azr" > "$d/run.azr"; else cp "$SRC/18O_p_a_thm.azr" "$d/run.azr"; fi
  [ -z "$2" ] || printf '\n<thm>\n%s\n</thm>\n' "$2" >> "$d/run.azr"
  (cd "$d" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr 2>&1 | head -c 1000000 > log;
   echo "${PIPESTATUS[1]}" > status)
}
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; fail=1; }
# Byte-identical files, without cmp (not in every MSYS2 image); \r stripped.
same() { [ -f "$1" ] && [ -f "$2" ] && [ "$(tr -d '\r' < "$1")" = "$(tr -d '\r' < "$2")" ]; }
# field FILE KEY -> the value after "KEY" in thm_experiments.out (first match)
field() { tr -d '\r' < "$1" | awk -v k="$2" '$1 == k || $1 == k ":" { print $2; exit }'; }

# (a) ----------------------------------------------------------------------
echo "(a) one-segment experiments without background == no experiment"
run base ""
[ "$(cat "$WORK/base/status")" = 0 ] && [ -f "$WORK/base/output/$OUT" ] || { echo "FAIL: base run"; tail -5 "$WORK/base/log"; exit 1; }
run one "experiment[A] segments=1
experiment[B] segments=2 beam=18O target=d spectator=n Ebeam=54   # kinematics: printed only"
for f in chiSquared.out normalizations.out "$OUT"; do
  if same "$WORK/base/output/$f" "$WORK/one/output/$f"; then ok "$f byte-identical"; else bad "$f differs"; fi
done
[ -f "$WORK/one/output/thm_experiments.out" ] && ok "thm_experiments.out written" || bad "no thm_experiments.out"
grep -q "Trojan horse d = x + n, B(x+s) = 2.2245" "$WORK/one/log" && ok "kinematics summary printed (B(d) = 2.2246 MeV)" \
  || bad "no kinematics summary"

# The models m (column 4 of the output, 10 digits), the c.m. energies, and
# the data: one output line per data-file line (all points are in range).
# model.txt: segment E_cm m d e
tr -d '\r' < "$WORK/base/output/$OUT" | awk -v D1="$SRC/data/lc723_thm_points.dat" -v D2="$SRC/data/lc723_band_mid.dat" '
  BEGIN { seg = 1; blank = 0 }
  NF == 0 { if (n[seg] > 0 && !blank) seg++; blank = 1; next }
  { blank = 0; n[seg]++; E[seg, n[seg]] = $1; M[seg, n[seg]] = $4 }
  END {
    for (s = 1; s <= 2; s++) {
      file = (s == 1) ? D1 : D2; k = 0
      while ((getline line < file) > 0) { split(line, f, " "); if (f[1] == "") continue; k++; d[k] = f[3]; e[k] = f[4] }
      close(file)
      if (k != n[s]) { print "MISMATCH", s, k, n[s]; exit 1 }
      for (i = 1; i <= k; i++) printf "%d %s %s %s %s\n", s, E[s, i], M[s, i], d[i], e[i]
    }
  }' > "$WORK/model.txt" || { echo "FAIL: output and data files do not line up"; cat "$WORK/model.txt"; exit 1; }

# (b) ----------------------------------------------------------------------
echo "(b) two segments sharing one profiled norm"
run shared "experiment[A] segments=1,2"
X="$WORK/shared/output/thm_experiments.out"
read -r hn hchi <<< "$(awk '$5 != 0 { w = 1 / ($5 * $5); Smm += $3 * $3 * w; Smd += $3 * $4 * w; Sdd += $4 * $4 * w }
  END { printf "%.12e %.12e\n", Smm / Smd, Sdd - Smd * Smd / Smm }' "$WORK/model.txt")"
en="$(field "$X" norm)"; echi="$(field "$X" chi2)"
near() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(a != "" && d <= t) }'; }
near "$en" "$hn" 1e-8 && ok "n* = sum S_mm / sum S_md: $en (hand $hn)" || bad "n* $en, hand $hn"
near "$echi" "$hchi" 1e-8 && ok "chi2 = sum S_dd - (sum S_md)^2 / sum S_mm: $echi (hand $hchi)" || bad "chi2 $echi, hand $hchi"
read -r n1 n2 c1 c2 tot <<< "$(tr -d '\r' < "$WORK/shared/output/chiSquared.out" | awk -F, '
  NR == 2 { n1 = $4; c1 = $2 } NR == 3 { n2 = $4; c2 = $2 }
  /Total-Chi-Squared/ { split($0, f, " "); t = f[2] } END { print n1, n2, c1, c2, t }')"
[ "$n1" = "$n2" ] && near "$n1" "$hn" 1e-5 && ok "chiSquared.out: both segments carry n* ($n1)" || bad "segment norms $n1 / $n2 (n* $hn)"
near "$tot" "$hchi" 1e-5 && ok "chiSquared.out total = the experiment's chi2 ($tot)" || bad "total $tot vs $hchi"
near "$(awk -v a="$c1" -v b="$c2" 'BEGIN { print a + b }')" "$tot" 1e-5 && ok "segment chi2 add up" || bad "segments $c1 + $c2 vs $tot"
nn="$(tr -d '\r' < "$WORK/shared/output/normalizations.out" | awk -F, 'NR > 1 { print $9 }' | sort -u | wc -l | tr -d ' ')"
[ "$nn" = 1 ] && ok "normalizations.out: one shared norm" || bad "normalizations.out has $nn different norms"

# (c) ----------------------------------------------------------------------
echo "(c) a linear background is recovered"
S_TRUE=33333.25; A0=0.3125; A1=-0.21
mkdir -p "$WORK/syn"
for s in 1 2; do
  f=$([ $s = 1 ] && echo lc723_thm_points.dat || echo lc723_band_mid.dat)
  # Keep the energies of the data file; replace sigma and its error (5 %;
  # zero where the original error is zero).
  awk -v s="$s" -v S="$S_TRUE" -v a0="$A0" -v a1="$A1" -v F="$SRC/data/$f" '
    $1 == s { k++; y[k] = S * $3 + a0 + a1 * $2; z[k] = ($5 == 0) ? 0 : 0.05 * y[k] }
    END { i = 0; while ((getline line < F) > 0) { split(line, t, " "); if (t[1] == "") continue; i++
            printf "%s %s %.17e %.17e\n", t[1], t[2], y[i], z[i] } }' "$WORK/model.txt" > "$WORK/syn/$f"
done
run bkg "experiment[A] segments=1-2 background=linear"
cp "$WORK/syn/"*.dat "$WORK/bkg/data/"
(cd "$WORK/bkg" && rm -rf output && mkdir output && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr > log 2>&1)
X="$WORK/bkg/output/thm_experiments.out"
en="$(field "$X" norm)"; b0="$(field "$X" b0)"; b1="$(field "$X" b1)"; echi="$(field "$X" chi2)"
near "$en" "$(awk -v s="$S_TRUE" 'BEGIN { printf "%.17e", 1 / s }')" 1e-8 && ok "n = 1/s recovered: $en" || bad "n $en, want 1/$S_TRUE"
near "$b0" "$(awk -v s="$S_TRUE" -v a="$A0" 'BEGIN { printf "%.17e", a / s }')" 1e-8 && ok "b0 = a0/s recovered: $b0" || bad "b0 $b0"
near "$b1" "$(awk -v s="$S_TRUE" -v a="$A1" 'BEGIN { printf "%.17e", a / s }')" 1e-8 && ok "b1 = a1/s recovered: $b1" || bad "b1 $b1"
awk -v c="$echi" 'BEGIN { exit !(c != "" && c < 1e-8) }' && ok "chi2 ~ 0 ($echi)" || bad "chi2 $echi"
grep -q "^status: profiled" "$X" && ok "status: profiled" || bad "status: $(grep status "$X")"
# The fitted curve written is m + b(E) = the data scaled by n (chi2 ~ 0).
worst="$(tr -d '\r' < "$WORK/bkg/output/$OUT" | awk 'NF >= 6 { d = ($4 - $6) / $6; if (d < 0) d = -d; if (d > w) w = d } END { printf "%.2e", w }')"
awk -v w="$worst" 'BEGIN { exit !(w < 1e-8) }' && ok "output: model + b(E) = data x n (worst rel $worst)" || bad "output curve off the data by $worst"
# The background is not a shift of the data: without it the same data do not fit.
run nobkg "experiment[A] segments=1-2"
cp "$WORK/syn/"*.dat "$WORK/nobkg/data/"
(cd "$WORK/nobkg" && rm -rf output && mkdir output && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr > log 2>&1)
awk -v c="$(field "$WORK/nobkg/output/thm_experiments.out" chi2)" 'BEGIN { exit !(c > 1) }' \
  && ok "without the background the same data give chi2 > 1" || bad "no-background chi2 is small"

# (d) ----------------------------------------------------------------------
echo "(d) refused"
refuse() {  # refuse NAME MESSAGE-FRAGMENT BLOCK [AWK]
  run "$1" "$3" "${4:-}"
  if [ "$(cat "$WORK/$1/status")" != 0 ] && grep -q "ERROR: <thm> experiment\[" "$WORK/$1/log" &&
     grep -q -- "$2" "$WORK/$1/log"; then
    ok "$1: $(grep -m1 'ERROR' "$WORK/$1/log" | tr -d '\r')"
  else
    bad "$1: not refused with '$2' (status $(cat "$WORK/$1/status"))"; tail -3 "$WORK/$1/log" | sed 's/^/        /'
  fi
}
for k in ps theta distortion; do
  refuse "reserved_$k" "not implemented yet" "experiment[A] segments=1,2 $k=1"
done
refuse unknown_key "unknown key 'foo'" "experiment[A] segments=1 foo=1"
refuse bad_nuclide "unknown nuclide '8Be'" "experiment[A] segments=1-2 beam=8Be target=d spectator=n Ebeam=54"
refuse bad_explicit "expected Z,A,mass" "experiment[A] segments=1-2 beam=18O target=2,1,2.0135 spectator=n Ebeam=54"
refuse partial_kinematics "all four or none" "experiment[A] segments=1-2 beam=18O target=d"
refuse bad_background "expected none, const, linear or quadratic" "experiment[A] segments=1 background=cubic"
refuse bad_segments "expected segment numbers" "experiment[A] segments=1,x"
refuse bad_ebeam "Ebeam='-3'" "experiment[A] segments=1-2 beam=18O target=d spectator=n Ebeam=-3"
refuse no_segments "segments= is required" "experiment[A] background=linear"
refuse repeated_key "is given twice" "experiment[A] segments=1
experiment[A] segments=2"
refuse two_experiments "already in experiment\[A\]" "experiment[A] segments=1,2
experiment[B] segments=2"
refuse no_such_segment "has only 2 line" "experiment[A] segments=1-3"
refuse not_thm "is not a THM segment" "experiment[A] segments=1,2" \
  '/<segmentsData>/ { S = 1; print; next } S && NF > 8 { S++; if (S == 3) $8 = 0 } /<\/segmentsData>/ { S = 0 } { print }'
refuse fixed_norm "has a fixed norm" "experiment[A] segments=1,2" \
  '/<segmentsData>/ { S = 1; print; next } S && NF > 8 { S++; if (S == 3) $10 = 0 } /<\/segmentsData>/ { S = 0 } { print }'
refuse wrong_pair "does not give the entrance pair" "experiment[A] segments=1-2 beam=18O target=t spectator=n Ebeam=54"

echo
if [ "$fail" -eq 0 ]; then echo "PASS: THM experiments (shared norm, background, refusals)"; else echo "FAIL"; exit 1; fi
