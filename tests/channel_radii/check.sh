#!/usr/bin/env bash
#
# Per-channel radii: the <channelRadii> block (2026-11-10).  A channel takes its pair's radius (column 28) unless a
# <channelRadii> line "<pair key> <l> <radius> [s=<s>] [J=<J>]" gives it its own.
#
# Truth without the feature: a partial wave's phase shift depends only on its own channel's radius (p + 12C: one
# channel per J^pi, so every J-group is one partial wave), so the l = 1 waves of a model whose l = 1 radius is
# overridden must equal those of the same model with the PAIR radius changed, and every other wave must equal the
# unchanged model's -- at the level of the observed parameters in the file, on every route.
#
#   (a) p + 12C (tests/hard_sphere/models/no_dummies.azr, a = 5.0 fm), override l = 1 -> 4.0 fm: the phase shifts of
#       1/2- and 3/2- equal the a = 4.0 model's, those of 1/2+, 3/2+, 5/2+ the a = 5.0 model's (1e-9), on the
#       A- and R-matrix routes, with and without --use-brune; and the two references differ (teeth);
#   (b) the J qualifier: "1 1 4.0 J=3/2" moves only 3/2- (1/2- stays at a = 5.0);
#   (c) a block giving every l the pair's radius is identical to no block; a block giving every l 4.0 fm is identical
#       to the a = 4.0 model (cross sections, both test segments);
#   (d) 13N (capture into the 1/2- ground state, E1 + E2, external capture on): every l of pair 1 at 4.0 fm by
#       <channelRadii> gives every point (AZUREOut) of the model with pair 1 at 4.0 fm (capture integrals, ANC
#       conversion, target integration); l = 4 alone at 4.0 fm (no capture transition) runs and changes the result;
#   (e) refused with an ERROR and no output: l = 0 alone at 4.0 fm in 13N (E1 from l = 0 into the l = 1 final
#       channel at 3.6 fm), an unknown pair key, a gamma pair, a non-positive radius, a malformed s=, a repeated line,
#       an unterminated block, a channel matched by an s line and a J line;  warned about: a line matching no channel;
#   (f) the compound-nucleus check prints the per-channel radii.
#
#   ./tests/channel_radii/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/../hard_sphere/models/no_dummies.azr"
N13="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/channel_radii.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# radius FILE R: every <levels> line of pair 1 with column 28 (the channel radius) set to R
radius() {
  awk -v r="$2" '/^<levels>/ { inl = 1; print; next } /^<\/levels>/ { inl = 0 }
                 inl && NF >= 28 && $6 == 1 { $28 = r } { print }' "$1"
}
# block FILE LINES...: a <channelRadii> block with the given lines right after </levels>
block() {
  local f="$1"; shift
  printf '%s\n' "$@" > "$WORK/.blk"
  awk -v b="$WORK/.blk" '{ print } /^<\/levels>/ { print "<channelRadii>"; while ((getline l < b) > 0) print l; print "</channelRadii>" }' "$f"
}
# segs FILE TESTLINES...: replace the <segmentsTest> lines
segs() {
  local f="$1"; shift
  printf '%s\n' "$@" > "$WORK/.seg"
  awk -v s="$WORK/.seg" '/^<segmentsTest>/ { print; while ((getline l < s) > 0) print l; skip = 1; next }
                         /^<\/segmentsTest>/ { skip = 0 } !skip { print }' "$f"
}
# route FILE true|false
route() { sed "1,2s/^\(true\|false\) #Perform A-Matrix/$2 #Perform A-Matrix/" "$1"; }

# extrap NAME AZR FLAGS -> mode 3 into $WORK/NAME, 0 if an .extrap was written
extrap() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$2" "$dir/t.azr"
  (cd "$dir" && printf '3\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" $3 --no-gui --no-readline t.azr > run.log 2>&1)
  [ -s "$dir/output/AZUREOut_aa=1_R=1.extrap" ]
}
# maxdiff A B -> largest |a - b| / max(|a|, |b|, 1e-12) over every numeric column of two output files (cross sections
# in barns are small numbers, so the comparison is relative; the floor keeps exact zeros and phases at zero sane)
maxdiff() {
  paste "$1" "$2" | awk '
    NF >= 2 && NF % 2 == 0 { n = NF / 2; for (i = 1; i <= n; i++) { a = $i + 0; b = $(i + n) + 0; d = a - b; if (d < 0) d = -d;
      m = 1e-12; if (a > m) m = a; if (-a > m) m = -a; if (b > m) m = b; if (-b > m) m = -b; if (d / m > w) w = d / m }; rows++ }
    END { if (rows == 0) print "nan"; else printf "%.3e\n", w }'
}
# outs DIR -> the 13N AZUREOut files of a run, concatenated
outs() { cat "$WORK/$1/output/AZUREOut_aa=1_R=1.out" "$WORK/$1/output/AZUREOut_aa=1_R=2.out"; }
le() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r != "nan" && r + 0 <= t + 0) }'; }
gt() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r != "nan" && r + 0 > t + 0) }'; }

radius "$BASE" 4.0 > "$WORK/b4.azr"
cp "$BASE" "$WORK/a5.azr"
grep -q " 5.0 0 0 0$" "$BASE" || { bad "fixture: pair radius 5.0 not found in $BASE"; exit 1; }

# ---- (a), (b): phase shifts wave by wave ---------------------------------------------------------------------
# wave J l  (the p + 12C waves the three levels and the hard-sphere groups give)
WAVES="0.5:0 0.5:1 1.5:1 1.5:2 2.5:2"
for rt in true false; do
  for flags in "" "--use-brune"; do
    tag="$( [ "$rt" = true ] && echo A || echo R )${flags:+_brune}"
    for w in $WAVES; do
      J="${w%%:*}"; L="${w##*:}"; seg="1 1 1 0.5 3.5 0.1 0 0 0 2 $J $L"
      for m in a5 b4; do
        route "$WORK/$m.azr" "$rt" | segs /dev/stdin "$seg" > "$WORK/${m}_p.azr"
        extrap "${tag}_${m}_${J}_${L}" "$WORK/${m}_p.azr" "$flags" || bad "$tag $m J=$J l=$L: no output"
      done
      block "$WORK/a5.azr" "1 1 4.0" | route /dev/stdin "$rt" | segs /dev/stdin "$seg" > "$WORK/m_p.azr"
      extrap "${tag}_m_${J}_${L}" "$WORK/m_p.azr" "$flags" || { bad "$tag override l=1, J=$J l=$L: no output"; tail -5 "$WORK/${tag}_m_${J}_${L}/run.log"; continue; }
      block "$WORK/a5.azr" "1 1 4.0 J=3/2" | route /dev/stdin "$rt" | segs /dev/stdin "$seg" > "$WORK/j_p.azr"
      extrap "${tag}_j_${J}_${L}" "$WORK/j_p.azr" "$flags" || { bad "$tag override l=1 J=3/2, J=$J l=$L: no output"; continue; }
      fa="$WORK/${tag}_a5_${J}_${L}/output/AZUREOut_aa=1_R=1.extrap"
      fb="$WORK/${tag}_b4_${J}_${L}/output/AZUREOut_aa=1_R=1.extrap"
      fm="$WORK/${tag}_m_${J}_${L}/output/AZUREOut_aa=1_R=1.extrap"
      fj="$WORK/${tag}_j_${J}_${L}/output/AZUREOut_aa=1_R=1.extrap"
      teeth="$(maxdiff "$fa" "$fb")"
      gt "$teeth" 1e-6 || bad "$tag J=$J l=$L: a = 4.0 and 5.0 give the same phase shift ($teeth): no teeth"
      if [ "$L" = 1 ]; then want="$fb"; wname="a=4.0"; else want="$fa"; wname="a=5.0"; fi
      r="$(maxdiff "$fm" "$want")"
      if le "$r" 1e-9; then ok "$tag l=1 override: J=$J l=$L == $wname model ($r; references differ by $teeth)"; else bad "$tag l=1 override: J=$J l=$L differs from the $wname model ($r)"; fi
      if [ "$J" = 1.5 ] && [ "$L" = 1 ]; then jwant="$fb"; jname="a=4.0"; else jwant="$fa"; jname="a=5.0"; fi
      r="$(maxdiff "$fj" "$jwant")"
      if le "$r" 1e-9; then ok "$tag l=1 J=3/2 override: J=$J l=$L == $jname model ($r)"; else bad "$tag l=1 J=3/2 override: J=$J l=$L differs from the $jname model ($r)"; fi
    done
  done
done

# ---- (c): all-l blocks against the plain models (cross sections) ----------------------------------------------
same5="$(block "$WORK/a5.azr" "1 0 5.0" "1 1 5.0" "1 2 5.0" "1 3 5.0" "1 4 5.0")"
all4="$(block "$WORK/a5.azr" "1 0 4.0" "1 1 4.0" "1 2 4.0" "1 3 4.0" "1 4 4.0")"
printf '%s\n' "$same5" > "$WORK/same5.azr"; printf '%s\n' "$all4" > "$WORK/all4.azr"
for rt in true false; do
  for flags in "" "--use-brune"; do
    tag="$( [ "$rt" = true ] && echo A || echo R )${flags:+_brune}"
    for m in a5 b4 same5 all4; do route "$WORK/$m.azr" "$rt" > "$WORK/${m}_x.azr"; extrap "${tag}_x_$m" "$WORK/${m}_x.azr" "$flags" || bad "$tag $m: no output"; done
    f() { echo "$WORK/${tag}_x_$1/output/AZUREOut_aa=1_R=1.extrap"; }
    r="$(maxdiff "$(f same5)" "$(f a5)")"; if le "$r" 0; then ok "$tag every l at the pair radius == no block (identical)"; else bad "$tag every l at the pair radius differs from no block ($r)"; fi
    r="$(maxdiff "$(f all4)" "$(f b4)")"; if le "$r" 1e-12; then ok "$tag every l at 4.0 fm == pair radius 4.0 ($r)"; else bad "$tag every l at 4.0 fm differs from pair radius 4.0 ($r)"; fi
  done
done

# ---- (d), (e, capture): 13N -------------------------------------------------------------------------------------
# calc NAME AZR -> mode 1 in a copy of tests/13N (data/ linked), 0 if chiSquared.out was written
calc() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  ln -sf "$N13/data" "$dir/data"
  cp "$2" "$dir/13N.azr"
  (cd "$dir" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline 13N.azr > run.log 2>&1)
  [ -s "$dir/output/chiSquared.out" ]
}
radius "$N13/13N.azr" 4.0 > "$WORK/n_b4.azr"
block "$N13/13N.azr" "1 0 4.0" "1 1 4.0" "1 2 4.0" "1 3 4.0" "1 4 4.0" > "$WORK/n_all4.azr"
block "$N13/13N.azr" "1 4 4.0" > "$WORK/n_l4.azr"
block "$N13/13N.azr" "1 0 4.0" > "$WORK/n_l0.azr"
calc n_a "$N13/13N.azr" || bad "13N: no chiSquared.out"
calc n_b4 "$WORK/n_b4.azr" || bad "13N pair radius 4.0: no chiSquared.out"
if calc n_all4 "$WORK/n_all4.azr"; then
  r="$(maxdiff <(outs n_all4) <(outs n_b4))"
  t="$(maxdiff <(outs n_a) <(outs n_b4))"
  if le "$r" 1e-9 && gt "$t" 1e-3; then ok "13N: every l at 4.0 fm == pair radius 4.0, every point of every segment ($r; 3.6 vs 4.0 differ by $t)"; else bad "13N: every l at 4.0 fm vs pair radius 4.0: $r (teeth $t)"; fi
else bad "13N every l at 4.0 fm: no chiSquared.out"; tail -5 "$WORK/n_all4/run.log"; fi
if calc n_l4 "$WORK/n_l4.azr"; then
  r="$(maxdiff <(outs n_l4) <(outs n_a))"
  if gt "$r" 1e-12; then ok "13N: l = 4 alone at 4.0 fm (no capture transition) runs and changes the calculation ($r)"; else bad "13N: l = 4 at 4.0 fm left the calculation unchanged ($r)"; fi
else bad "13N l = 4 at 4.0 fm: refused"; tail -5 "$WORK/n_l4/run.log"; fi
if calc n_l0 "$WORK/n_l0.azr"; then bad "13N: l = 0 alone at 4.0 fm (E1 into the l = 1 final channel at 3.6 fm) was not refused";
elif grep -q "ERROR: <channelRadii>: external capture" "$WORK/n_l0/run.log"; then ok "13N: l = 0 alone at 4.0 fm refused (capture transition across two radii)";
else bad "13N: l = 0 alone at 4.0 fm failed without the capture ERROR"; tail -5 "$WORK/n_l0/run.log"; fi

# ---- (e): malformed blocks ---------------------------------------------------------------------------------------
refused() {  # NAME PATTERN LINES...
  local name="$1" pat="$2"; shift 2
  block "$WORK/a5.azr" "$@" > "$WORK/$name.azr"
  if extrap "neg_$name" "$WORK/$name.azr" ""; then bad "$name: not refused"
  elif grep -q "$pat" "$WORK/neg_$name/run.log"; then ok "$name refused"
  else bad "$name: no '$pat' in the log"; tail -5 "$WORK/neg_$name/run.log"; fi
}
refused unknown_key "names pair key 9" "9 1 4.0"
refused nonpositive "the radius positive" "1 1 0"
refused bad_spin "is not s=<s> or J=<J>" "1 1 4.0 s=1/3"
refused repeated "repeats an earlier line" "1 1 4.0" "1 1 4.5"
refused ambiguous "matches one line naming s and another naming J" "1 1 4.5 s=1/2" "1 1 4.2 J=3/2"
refused junk "is not \"<pair key> <l> <radius" "1 one 4.0"
awk '{ print } /^<\/levels>/ { print "<channelRadii>"; print "1 1 4.0" }' "$WORK/a5.azr" | grep -v "^</channelRadii>" > "$WORK/open.azr"
if extrap neg_open "$WORK/open.azr" ""; then bad "unterminated block: not refused"; elif grep -q "has no </channelRadii>" "$WORK/neg_open/run.log"; then ok "unterminated block refused"; else bad "unterminated block: no message"; tail -3 "$WORK/neg_open/run.log"; fi
block "$N13/13N.azr" "2 1 4.0" > "$WORK/gamma.azr"
if calc neg_gamma "$WORK/gamma.azr"; then bad "gamma pair: not refused"; elif grep -q "a gamma pair; only particle pairs" "$WORK/neg_gamma/run.log"; then ok "gamma pair refused"; else bad "gamma pair: no message"; tail -3 "$WORK/neg_gamma/run.log"; fi
block "$WORK/a5.azr" "1 9 4.0" > "$WORK/unused.azr"
if extrap unused "$WORK/unused.azr" "" && grep -q "WARNING: <channelRadii> entry for pair 1, l = 9 matches no channel" "$WORK/unused/run.log"; then
  r="$(maxdiff "$WORK/unused/output/AZUREOut_aa=1_R=1.extrap" "$WORK/A_x_a5/output/AZUREOut_aa=1_R=1.extrap")"
  if le "$r" 0; then ok "a line matching no channel is warned about and changes nothing"; else bad "unused line changed the result ($r)"; fi
else bad "unused line: no warning or no output"; tail -3 "$WORK/unused/run.log"; fi

# ---- (f): the compound-nucleus check lists the channel radii ---------------------------------------------------
block "$WORK/a5.azr" "1 1 4.0 J=3/2" | sed '5s/^none #Compond/screen #Compond/' > "$WORK/report.azr"
if extrap report "$WORK/report.azr" "" && grep -q "Channel Radius  l = 1, J = 1.5: 4" "$WORK/report/run.log"; then ok "compound-nucleus check lists the per-channel radius"; else bad "compound-nucleus check: no per-channel radius line"; grep -i "radius" "$WORK/report/run.log" | head -5; fi

if [ "$failures" -eq 0 ]; then echo "PASS: per-channel radii"; else echo "FAIL: $failures check(s)"; exit 1; fi
