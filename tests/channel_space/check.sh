#!/usr/bin/env bash
#
# The channel-space solution of the level matrix (ChannelSpace.h,
# AMatrixFunc::SolveChannelSpace) against the level-space one, on a synthetic
# p+12C project with many levels per J-group and few channels
# (channel_space.azr): J = 1/2+ with eight levels, 3/2- with seven, each with
# four channels (p0, two p1 channel spins of a synthetic 2+ state, gamma), and
# a one-level 5/2+ group that stays in level space.  The elastic data include
# a point within 1e-9 MeV of a level (a formal pole of the Lane-Thomas level
# matrix), where the channel route must hand over to the level route.
#
#   (a) every output file of "calculate with data" agrees with
#       AZURE_LEVEL_SPACE=1 to 1e-10 of the largest value of its column --
#       Brune (default), Lane-Thomas (--no-brune) and Park (--use-park);
#   (b) the bilinear forms G^T A G of every channel-space solve agree with
#       the level matrix's to 1e-11 of the largest (AZURE_CHANNEL_VERIFY); the
#       J-groups with more levels than channels are solved in channel space,
#       and the near-pole point falls back to the level route (AZURE_PROFILE
#       counters);
#   (c) the analytic gradient (adjoint and residual-Jacobian paths) agrees with
#       central finite differences of chi^2 in the three formalisms, and with
#       the level-space gradient to 1e-9 of its largest component.
#
#   ./tests/channel_space/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/channel_space.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME MODE "FLAGS" [VAR=VALUE ...] -> a fresh copy of the project in $WORK/NAME, run there
run() {
  local name="$1" mode="$2" flags="$3"
  shift 3
  local dir="$WORK/$name"
  rm -rf "$dir"
  mkdir -p "$dir/output" "$dir/checks"
  cp -r "$HERE/data" "$dir/"
  cp "$HERE/channel_space.azr" "$dir/t.azr"
  (cd "$dir" && for v in "$@"; do export "$v"; done &&
   printf "$mode\n\n\n7\n" | run_limited 900 "$AZURE2_BIN" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}

# maxdiff A B -> largest |a - b| over the numeric fields of two files, relative to
# the largest |value| of the field's column; "lines" if the files differ in shape.
maxdiff() {
  awk 'function abs(x) { return x < 0 ? -x : x }
       FNR == NR { for (i = 1; i <= NF; i++) a[FNR, i] = $i; na[FNR] = NF; nla = FNR; next }
       { if (na[FNR] != NF) { shape = 1; next }
         for (i = 1; i <= NF; i++) {
           x = a[FNR, i]; y = $i
           if (x == y) continue
           if (x !~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/ || y !~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/) { shape = 1; continue }
           d = abs(x - y); if (d > dmax[i]) dmax[i] = d
         }
         for (i = 1; i <= NF; i++) if ($i ~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/) {
           v = abs($i); if (v > cmax[i]) cmax[i] = v }
         nlb = FNR }
       END { if (shape || nla != nlb) { print "lines"; exit }
             w = 0; for (i in dmax) { r = (cmax[i] > 0) ? dmax[i] / cmax[i] : dmax[i]; if (r > w) w = r }
             printf "%.2e\n", w }' "$1" "$2"
}

profile_count() {  # NAME "counter label" -> the count AZURE_PROFILE reported last for it (0 if none)
  awk -v k="$2" 'index($0, k) { n = $NF } END { print n + 0 }' "$WORK/$1/run.log"
}

# ---- (a), (b) outputs, channel route against level route ---------------------------------------
for mode in brune:"" lanethomas:"--no-brune" park:"--use-park"; do
  tag="${mode%%:*}"; flags="${mode#*:}"
  run "c_$tag" 1 "$flags" AZURE_PROFILE=1 AZURE_CHANNEL_VERIFY=1
  run "l_$tag" 1 "$flags" AZURE_LEVEL_SPACE=1
  n="$(ls "$WORK/c_$tag/output" | wc -l | tr -d ' ')"
  if [ "$n" -lt 4 ] || [ ! -s "$WORK/c_$tag/output/chiSquared.out" ]; then
    bad "$tag: calculation wrote no output ($n files)"; tail -5 "$WORK/c_$tag/run.log"; continue
  fi
  worst=0; worstfile=""
  for f in "$WORK/c_$tag/output"/*; do
    g="$WORK/l_$tag/output/$(basename "$f")"
    if [ ! -f "$g" ]; then bad "$tag: $(basename "$f") missing on the level route"; continue; fi
    d="$(maxdiff "$f" "$g")"
    if [ "$d" = lines ]; then bad "$tag: $(basename "$f") differs in shape between the routes"; continue; fi
    if awk -v d="$d" -v w="$worst" 'BEGIN { exit !(d + 0 > w + 0) }'; then worst="$d"; worstfile="$(basename "$f")"; fi
  done
  if awk -v w="$worst" 'BEGIN { exit !(w + 0 <= 1e-10) }'; then
    ok "$tag: outputs ($n files) agree with the level route (worst $worst of the column scale${worstfile:+, $worstfile})"
  else
    bad "$tag: outputs differ from the level route by $worst of the column scale in $worstfile"
  fi
  # the bilinear forms themselves, at every channel-space solve (AZURE_CHANNEL_VERIFY)
  v="$(awk '/^\[channel-verify\]/ { for (i = 1; i <= NF; i++) if ($i == "route") v = $(i + 1) } END { print v + 0 }' "$WORK/c_$tag/run.log")"
  if awk -v v="$v" 'BEGIN { exit !(v + 0 <= 1e-11) }'; then
    ok "$tag: bilinear forms G^T A G agree with the level matrix to $v of the largest"
  else
    bad "$tag: bilinear forms differ from the level matrix by $v of the largest"
  fi
  nc="$(profile_count "c_$tag" "channel-space solves")"
  [ "$nc" -gt 0 ] && ok "$tag: $nc J-group solves in channel space" || bad "$tag: no channel-space solve"
  if [ "$tag" = lanethomas ]; then
    np="$(profile_count "c_$tag" "near pole")"
    [ "$np" -gt 0 ] && ok "$tag: the point at a level energy falls back to the level route ($np solves)" \
      || bad "$tag: no near-pole fallback"
  fi
done

# ---- (c) gradients -----------------------------------------------------------------------------
rows() {
  awk '/^\[grad-check\] idx / { idx = $3; a = f = ""; j = "-"
         for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
         print idx, a, f, j }' "$WORK/$1/run.log" | sort -n
}
gradcheck() {  # NAME "FLAGS" [VAR=VALUE ...]
  local name="$1" flags="$2"
  shift 2
  run "$name" 2 "--use-gradient $flags" AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
    AZURE_GRAD_CHECK_H=1e-7 "$@"
  rows "$name" > "$WORK/$name.rows"
}
for mode in brune:"" lanethomas:"--no-brune" park:"--use-park"; do
  tag="${mode%%:*}"; flags="${mode#*:}"
  gradcheck "g_$tag" "$flags"
  gradcheck "gl_$tag" "$flags" AZURE_LEVEL_SPACE=1
  nrows="$(wc -l < "$WORK/g_$tag.rows" | tr -d ' ')"
  if [ "$nrows" -lt 10 ]; then
    bad "$tag: gradient check listed $nrows parameters"; grep -i "grad" "$WORK/g_$tag/run.log" | head -5
    continue
  fi
  # relative to the largest component: 1e-5 against finite differences (Park: 1e-4,
  # its level matrix carries dS/dE differentiated numerically)
  tol=1e-5; [ "$tag" = park ] && tol=1e-4
  worst="$(awk '{ a = $2 < 0 ? -$2 : $2; if (a > m) m = a; g[NR] = $2; f[NR] = $3; j[NR] = $4 }
             END { for (i = 1; i <= NR; i++) { d = g[i] - f[i]; if (d < 0) d = -d; e = j[i] - f[i]; if (e < 0) e = -e
                                               if (d > w) w = d; if (j[i] != "-" && e > w) w = e }
                   printf "%.2e", w / m }' "$WORK/g_$tag.rows")"
  if awk -v w="$worst" -v t="$tol" 'BEGIN { exit !(w + 0 <= t + 0) }'; then
    ok "$tag: analytic and Jacobian-path gradients agree with finite differences for $nrows parameters (worst $worst of the largest)"
  else
    bad "$tag: gradient vs finite differences, worst $worst of the largest (> $tol)"
  fi
  same="$(join "$WORK/g_$tag.rows" "$WORK/gl_$tag.rows" | awk '{ a = $2 < 0 ? -$2 : $2; if (a > m) m = a
              d = $2 - $5; if (d < 0) d = -d; if (d > w) w = d
              if ($4 != "-" && $7 != "-") { e = $4 - $7; if (e < 0) e = -e; if (e > w) w = e }; n++ }
              END { if (n == 0 || m == 0) print "nan"; else printf "%.2e %d", w / m, n }')"
  if awk -v s="$same" 'BEGIN { split(s, x, " "); exit !(x[2] > 0 && x[1] + 0 <= 1e-9) }'; then
    ok "$tag: channel-space gradient = level-space gradient (${same% *} of the largest, ${same#* } parameters)"
  else
    bad "$tag: channel-space vs level-space gradient: $same"
  fi
done

if [ "$failures" -eq 0 ]; then echo "channel_space: all checks passed"; else echo "channel_space: $failures check(s) FAILED"; fi
[ "$failures" -eq 0 ]
