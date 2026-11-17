#!/usr/bin/env bash
#
# Per-channel radii: the `channel_radius` records of <channelRadii>.  A channel takes its pair's radius (column 28 of
# <levels>) unless a record (pair, l, radius, optionally s and J) gives it its own.
#
# Truth without the feature: a partial wave's phase shift depends only on its own channel's radius (p + 12C: one
# channel per J^pi, so every J-group is one partial wave), so the l = 1 waves of a model whose l = 1 radius is given by
# a record must equal those of the same model with the PAIR radius changed, and every other wave must equal the
# unchanged model's.
#
#   (a) p + 12C (tests/hard_sphere/models/no_dummies.azr, a = 5.0 fm), a record l = 1 -> 4.0 fm: the phase shifts of
#       1/2- and 3/2- equal the a = 4.0 model's, those of 1/2+, 3/2+, 5/2+ the a = 5.0 model's (1e-9), on the A- and
#       R-matrix routes, with and without --no-brune; and the two references differ (teeth);
#   (b) the J field: l = 1, J = 3/2 moves only 3/2- (1/2- stays at a = 5.0);
#   (c) records giving every l the pair's radius are identical to no record; records giving every l 4.0 fm are
#       identical to the a = 4.0 model (both extrapolations);
#   (d) 13N (capture into the 1/2- ground state, external capture on): every l of p+12C at 4.0 fm by records gives
#       every point of the model with the pair at 4.0 fm (capture integrals, ANC conversion, target integration);
#       l = 0 alone at 4.0 fm (E1 from l = 0 into the l = 1 final channel at 3.6 fm) is refused;
#   (e) refused: a photon pair, a selection of channels given twice, a channel named by a record with s and one with
#       J; warned about, and changing nothing: a record that names no channel.
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
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

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
# records FILE SPEC...: a <channelRadii> section after </levels>; SPEC is l:radius[:s=<s>][:J=<J>] for pair p+12C
records() {
  local f="$1" n=0 spec; shift
  : > "$WORK/.blk"
  for spec in "$@"; do
    n=$((n + 1))
    echo "$spec" | awk -F: -v n="$n" '{ print "channel_radius r" n; print "    pair p+12C"; print "    l " $1; print "    radius " $2 " fm";
      for (i = 3; i <= NF; i++) { split($i, kv, "="); print "    " kv[1] " " kv[2] } }' >> "$WORK/.blk"
  done
  awk -v b="$WORK/.blk" '{ print } /^<\/levels>/ { print "<channelRadii>"; while ((getline l < b) > 0) print l; close(b); print "</channelRadii>" }' "$f"
}
# waves FILE: the extrapolations replaced by the phase shifts of the five p + 12C waves (J, l; s = 1/2), seven
# energies each, in the order of WAVES: rows 1-7 of the output are the first wave, 8-14 the second, ...
WAVES="1/2:0 1/2:1 3/2:1 3/2:2 5/2:2"
waves() {
  awk -v W="$WAVES" '/^<extrapolations>/ { print; n = split(W, w, " ");
      for (i = 1; i <= n; i++) { split(w[i], jl, ":"); print "extrapolation ph" i; print "    reaction 12C(p,p)12C";
        print "    quantity phase_shift"; print "    unit deg"; print "    J " jl[1]; print "    l " jl[2]; print "    s 1/2";
        print "    energies list 0.5 1 1.5 2 2.5 3 3.5 MeV lab" }
      skip = 1; next }
    /^<\/extrapolations>/ { skip = 0 } !skip { print }' "$1"
}
# rows FILE K: the seven rows of wave K (1-based)
rows() { awk -v k="$2" 'NF > 0 { n++; if (n > 7 * (k - 1) && n <= 7 * k) print }' "$1"; }
# route FILE true|false: the A-matrix line of <config>
route() { awk -v v="$2" 'NR <= 2 && /#Perform A-Matrix/ { sub(/^(true|false)/, v) } { print }' "$1"; }

# extrap NAME AZR FLAGS -> extrapolate in $WORK/NAME; 0 if the first extrapolation was written
extrap() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$2" "$dir/t.azr"
  (cd "$dir" && run_limited 300 "$AZURE2_BIN" extrapolate $3 --no-gui t.azr < /dev/null > run.log 2>&1)
  [ -s "$dir/output/AZUREOut_aa=1_R=1.extrap" ]
}
# maxdiff A B -> largest |a - b| / max(|a|, |b|, 1e-12) over every numeric column of two output files
maxdiff() {
  paste "$1" "$2" | awk '
    NF >= 2 && NF % 2 == 0 { n = NF / 2; for (i = 1; i <= n; i++) { a = $i + 0; b = $(i + n) + 0; d = a - b; if (d < 0) d = -d;
      m = 1e-12; if (a > m) m = a; if (-a > m) m = -a; if (b > m) m = b; if (-b > m) m = -b; if (d / m > w) w = d / m }; rows++ }
    END { if (rows == 0) print "nan"; else printf "%.3e\n", w }'
}
le() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r ~ /^ *[-+]?\.?[0-9]/ && r + 0 <= t + 0) }'; }
gt() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r ~ /^ *[-+]?\.?[0-9]/ && r + 0 > t + 0) }'; }
out1() { echo "$WORK/$1/output/AZUREOut_aa=1_R=1.extrap"; }

awk '/^<levels>/ { inl = 1; next } /^<\/levels>/ { inl = 0 } inl && NF >= 28 && $6 == 1 && $28 + 0 != 5.0 { bad = 1 } END { exit bad }' "$BASE" ||
  { bad "fixture: the pair radius of $BASE is not 5.0 fm"; exit 1; }
cp "$BASE" "$WORK/a5.azr"
radius "$BASE" 4.0 > "$WORK/b4.azr"
records "$BASE" 1:4.0 > "$WORK/m.azr"
records "$BASE" 1:4.0:J=3/2 > "$WORK/j.azr"

# ---- (a), (b): phase shifts wave by wave ---------------------------------------------------------------------
for rt in true false; do
  for flags in "" "--no-brune"; do
    tag="$( [ "$rt" = true ] && echo A || echo R )${flags:+_formal}"
    good=1
    for m in a5 b4 m j; do
      route "$WORK/$m.azr" "$rt" > "$WORK/.r.azr"; waves "$WORK/.r.azr" > "$WORK/.w.azr"
      extrap "${tag}_w_$m" "$WORK/.w.azr" "$flags" || { bad "$tag $m: no phase shifts"; tail -5 "$WORK/${tag}_w_$m/run.log"; good=0; }
    done
    [ "$good" = 1 ] || continue
    n="$(awk 'NF > 0 { n++ } END { print n + 0 }' "$(out1 "${tag}_w_a5")")"
    [ "$n" = 35 ] || { bad "$tag: $n rows of phase shifts, expected 35 (five waves, seven energies)"; continue; }
    k=0
    for w in $WAVES; do
      k=$((k + 1)); J="${w%%:*}"; L="${w##*:}"
      for m in a5 b4 m j; do rows "$(out1 "${tag}_w_$m")" "$k" > "$WORK/.$m"; done
      teeth="$(maxdiff "$WORK/.a5" "$WORK/.b4")"
      gt "$teeth" 1e-6 || bad "$tag J=$J l=$L: a = 4.0 and 5.0 give the same phase shift ($teeth): no teeth"
      if [ "$L" = 1 ]; then want=b4; wname="a=4.0"; else want=a5; wname="a=5.0"; fi
      r="$(maxdiff "$WORK/.m" "$WORK/.$want")"
      if le "$r" 1e-9; then ok "$tag l=1 record: J=$J l=$L == $wname model ($r; references differ by $teeth)"; else bad "$tag l=1 record: J=$J l=$L differs from the $wname model ($r)"; fi
      if [ "$J" = 3/2 ] && [ "$L" = 1 ]; then want=b4; wname="a=4.0"; else want=a5; wname="a=5.0"; fi
      r="$(maxdiff "$WORK/.j" "$WORK/.$want")"
      if le "$r" 1e-9; then ok "$tag l=1 J=3/2 record: J=$J l=$L == $wname model ($r)"; else bad "$tag l=1 J=3/2 record: J=$J l=$L differs from the $wname model ($r)"; fi
    done
  done
done

# ---- (c): every l against the plain models (cross sections) ----------------------------------------------------
records "$BASE" 0:5.0 1:5.0 2:5.0 3:5.0 4:5.0 > "$WORK/same5.azr"
records "$BASE" 0:4.0 1:4.0 2:4.0 3:4.0 4:4.0 > "$WORK/all4.azr"
for rt in true false; do
  for flags in "" "--no-brune"; do
    tag="$( [ "$rt" = true ] && echo A || echo R )${flags:+_formal}"
    for m in a5 b4 same5 all4; do route "$WORK/$m.azr" "$rt" > "$WORK/.r.azr"; extrap "${tag}_x_$m" "$WORK/.r.azr" "$flags" || bad "$tag $m: no output"; done
    r="$(maxdiff "$(out1 "${tag}_x_same5")" "$(out1 "${tag}_x_a5")")"; if le "$r" 0; then ok "$tag every l at the pair radius == no record (identical)"; else bad "$tag every l at the pair radius differs from no record ($r)"; fi
    r="$(maxdiff "$(out1 "${tag}_x_all4")" "$(out1 "${tag}_x_b4")")"; if le "$r" 1e-12; then ok "$tag every l at 4.0 fm == pair radius 4.0 ($r)"; else bad "$tag every l at 4.0 fm differs from pair radius 4.0 ($r)"; fi
    r="$(maxdiff "$(out1 "${tag}_x_a5")" "$(out1 "${tag}_x_b4")")"; gt "$r" 1e-3 || bad "$tag: the a = 4.0 and 5.0 cross sections do not differ ($r): no teeth"
  done
done

# ---- (d): 13N, external capture -----------------------------------------------------------------------------------
# calc NAME AZR -> calculate in a copy of tests/13N (its data copied), 0 if chiSquared.out was written
calc() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -r "$N13/data" "$dir/data"
  cp "$2" "$dir/13N.azr"
  (cd "$dir" && run_limited 600 "$AZURE2_BIN" calculate --no-gui 13N.azr < /dev/null > run.log 2>&1)
  [ -s "$dir/output/chiSquared.out" ]
}
outs() { cat "$WORK/$1"/output/AZUREOut_aa=1_R=*.out; }
radius "$N13/13N.azr" 4.0 > "$WORK/n_b4.azr"
records "$N13/13N.azr" 0:4.0 1:4.0 2:4.0 3:4.0 4:4.0 > "$WORK/n_all4.azr"
records "$N13/13N.azr" 0:4.0 > "$WORK/n_l0.azr"
calc n_a "$N13/13N.azr" || { bad "13N: no chiSquared.out"; tail -5 "$WORK/n_a/run.log"; }
calc n_b4 "$WORK/n_b4.azr" || { bad "13N pair radius 4.0: no chiSquared.out"; tail -5 "$WORK/n_b4/run.log"; }
if calc n_all4 "$WORK/n_all4.azr"; then
  outs n_all4 > "$WORK/o_all4"; outs n_b4 > "$WORK/o_b4"; outs n_a > "$WORK/o_a"
  r="$(maxdiff "$WORK/o_all4" "$WORK/o_b4")"; t="$(maxdiff "$WORK/o_a" "$WORK/o_b4")"
  if le "$r" 1e-9 && gt "$t" 1e-3; then ok "13N: every l at 4.0 fm == pair radius 4.0, every point of every data set ($r; 3.6 vs 4.0 differ by $t)"; else bad "13N: every l at 4.0 fm vs pair radius 4.0: $r (teeth $t)"; fi
else bad "13N every l at 4.0 fm: no chiSquared.out"; tail -5 "$WORK/n_all4/run.log"; fi
if calc n_l0 "$WORK/n_l0.azr"; then bad "13N: l = 0 alone at 4.0 fm (E1 into the l = 1 final channel at 3.6 fm) was not refused";
elif grep -q "external capture across two channel radii" "$WORK/n_l0/run.log"; then ok "13N: l = 0 alone at 4.0 fm refused (a capture transition across two radii)";
else bad "13N: l = 0 alone at 4.0 fm failed without the capture ERROR"; tail -5 "$WORK/n_l0/run.log"; fi

# ---- (e): refusals and the warning ------------------------------------------------------------------------------
refused() {  # NAME PATTERN AZR
  if extrap "neg_$1" "$3" ""; then bad "$1: not refused"
  elif grep -q "$2" "$WORK/neg_$1/run.log"; then ok "$1 refused ($2)"
  else bad "$1: no '$2' in the log"; tail -5 "$WORK/neg_$1/run.log"; fi
}
records "$BASE" 1:4.0 1:4.5 > "$WORK/twice.azr"
refused "a selection given twice" "E_CHANNEL_RADIUS_TWICE" "$WORK/twice.azr"
records "$BASE" 1:4.5:s=1/2 1:4.2:J=3/2 > "$WORK/ambiguous.azr"
refused "a channel named by an s record and a J record" "its radius is not defined" "$WORK/ambiguous.azr"
records "$BASE" 1:0 > "$WORK/zero.azr"
refused "a radius of zero" "radius" "$WORK/zero.azr"
records "$N13/13N.azr" 1:4.0 | awk '{ sub(/pair p\+12C/, "pair g+13N"); print }' > "$WORK/gamma.azr"
if calc neg_gamma "$WORK/gamma.azr"; then bad "a photon pair: not refused"; elif grep -q "E_CHANNEL_RADIUS_PAIR_KIND" "$WORK/neg_gamma/run.log"; then ok "a photon pair refused (E_CHANNEL_RADIUS_PAIR_KIND)"; else bad "a photon pair: no E_CHANNEL_RADIUS_PAIR_KIND"; tail -5 "$WORK/neg_gamma/run.log"; fi
records "$BASE" 9:4.0 > "$WORK/unused.azr"
if extrap unused "$WORK/unused.azr" "" && grep -q "WARNING: channel_radius r1 .*gives no channel of the model its radius" "$WORK/unused/run.log"; then
  r="$(maxdiff "$(out1 unused)" "$(out1 A_x_a5)")"
  if le "$r" 0; then ok "a record that names no channel is warned about and changes nothing"; else bad "the unused record changed the result ($r)"; fi
else bad "the unused record: no warning or no output"; tail -3 "$WORK/unused/run.log"; fi

if [ "$failures" -eq 0 ]; then echo "PASS: per-channel radii"; else echo "FAIL: $failures check(s)"; exit 1; fi
