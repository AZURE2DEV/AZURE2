#!/usr/bin/env bash
#
# The light particle of a pair in an excited state: the optional 32nd column e1 (MeV) of a <levels> line
# (2026-11-05).  PPair keeps e1 and the heavy partner's e2 apart; GetExE() is e1 + e2, the channel's excitation
# above the ground-state threshold.  Two particles are identical when Z, mass, spin, parity, g-factor AND
# excitation are equal, so 12C(2+)+12C(2+) is entered as e1 = e2 = 4.43891 MeV, j1 = j2 = 2.
# Models from make_models.py (see there); comparisons in analyze.py.
#
#   (a) byte identity: an explicit "0" in column 32 on every line gives byte-identical output files (and run
#       log, progress bars aside) to the same project without it -- alpha + 12C in modes 1 and 3, and tests/13N (capture, external
#       capture, many segments) in mode 1;
#   (b) moving the excitation: alpha' + 12C(2+) entered with the excitation on the heavy partner (e2 = E*,
#       J2 = 2+) and on the light one (e1 = E*, J1 = 2+, J2 = 0+) give the same (a,a'), (a',a), elastic and
#       a' + 12C(2+) elastic cross sections, angle-integrated, c.m. and lab differential, extrapolated and at
#       lab-frame data points, to 1e-10; without column 32 the light version differs (the check sees e1);
#   (c) mutual excitation: 12C(2+)+12C(2+) entered as e1 = e2 is identical -- the hard-sphere J^pi groups
#       (includeAllJPi, l <= 4) hold exactly the l + s even channels, its elastic c.m. cross section is
#       symmetric about 90 deg, and a level listing an l + s odd channel is refused; the same pair entered the
#       old way (e1 = 0, e2 = 2 E*) is not identical: all channels, and an asymmetric cross section;
#   (d) refused with an ERROR and no output: e1 on a gamma pair, a negative e1, two e1 for one pair key, and a
#       UPOS (secondary-gamma) segment whose exit pair has e1 (the correlation is that of the heavy partner's
#       decay); the same segment with the excitation on the heavy partner runs.
#
# The symmetry in (c) is checked to 1e-6, not 1e-10: the Coulomb amplitude uses pi = 3.141592650
# (include/Constants.h), which breaks sigma(theta) = sigma(180 - theta) by ~1e-7 at 10/170 deg for any
# identical pair -- 1.4e-7 for ground-state 12C+12C, a model untouched by this change.
#
#   ./tests/light_excitation/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PY="$(command -v python3 || command -v python)"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/light_excitation.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
# a must be a number: an empty string (analyze.py failed) or nan/inf is a failure, never 0
num='^[-+]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$'
le() { awk -v a="$1" -v b="$2" -v re="$num" 'BEGIN { exit !(a ~ re && a + 0 <= b + 0) }'; }
gt() { awk -v a="$1" -v b="$2" -v re="$num" 'BEGIN { exit !(a ~ re && a + 0 > b + 0) }'; }
an() { "$PY" -I "$HERE/analyze.py" "$@"; }
# two run logs alike, apart from the threaded per-segment progress bars.  Several threads draw their bars at
# once with unsynchronized writes, so the bar text (and the \n that ends each bar) interleaves differently run
# to run: a fragment such as "[*****   ] 80%" or "11  " can land on a line of its own, after another thread's
# \n, without "Segment #" (macOS and Windows CI, 016bcfe; it passed by luck on 6868d44).  Ignore every line
# that holds a \r or "Segment #", is blank, or is made only of characters a bar can contain; compare the rest.
# Byte comparisons in Python: the Windows CI runner has no diff or cmp.  A mismatch prints the first
# differing line (samelog) or file (same_tree), so the next CI log says what moved.
samelog() {
  "$PY" -I -c 'import sys
bar = set(b"Segment#[]*% \t0123456789")
def f(p):
    out = []
    for l in open(p, "rb").read().split(b"\n"):
        if b"\r" in l or b"Segment #" in l or set(l) <= bar:
            continue
        out.append(l)
    return out
a, b = f(sys.argv[1]), f(sys.argv[2])
if a != b:
    i = next((k for k in range(min(len(a), len(b))) if a[k] != b[k]), min(len(a), len(b)))
    print("    samelog: %d vs %d lines, first difference at %d:" % (len(a), len(b), i))
    print("      %r\n      %r" % (a[i] if i < len(a) else None, b[i] if i < len(b) else None))
    sys.exit(1)' "$1" "$2"
}
# same_tree <dir1> <dir2>: the same file names with byte-identical contents (intEC.dat: numbers to 1e-12)
same_tree() {
  "$PY" -I -c 'import os, sys
def files(top):
    return sorted(os.path.relpath(os.path.join(d, f), top) for d, _, fs in os.walk(top) for f in fs)
a, b = sys.argv[1:3]
na, nb = files(a), files(b)
rd = lambda p: open(p, "rb").read()
if na != nb:
    print("    same_tree: file lists differ: %s" % sorted(set(na) ^ set(nb)))
    sys.exit(1)
def numeric_same(pa, pb, tol=1e-12):
    # intEC.dat holds the EC integrals to 17 digits; macOS CI wrote last-digit differences between two runs
    # that should agree (2026-11-06), every other file byte-identical.  Compare it number by number.
    ta, tb = rd(pa).split(), rd(pb).split()
    if len(ta) != len(tb):
        return False, "token counts %d vs %d" % (len(ta), len(tb))
    worst = 0.0
    for x, y in zip(ta, tb):
        if x == y:
            continue
        try:
            fx, fy = float(x), float(y)
        except ValueError:
            return False, "text differs: %r vs %r" % (x, y)
        d = abs(fx - fy) / max(abs(fx), abs(fy), 1e-300)
        worst = max(worst, d)
    return worst <= tol, "largest relative difference %.3g" % worst
bad = []
for n in na:
    pa, pb = os.path.join(a, n), os.path.join(b, n)
    if rd(pa) == rd(pb):
        continue
    if os.path.basename(n) == "intEC.dat":
        ok, why = numeric_same(pa, pb)
        print("    same_tree: %s not byte-identical, %s (tolerance 1e-12)" % (n, why))
        if ok:
            continue
    bad.append(n)
if bad:
    print("    same_tree: contents differ: %s" % bad)
    sys.exit(1)' "$1" "$2"
}

"$PY" -I "$HERE/make_models.py" "$WORK/models" || { bad "models not built"; exit 1; }

# run <dir> <model .azr> <mode>
run() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  [ -e "$dir/data" ] || ln -s "$(dirname "$2")/data" "$dir/data"
  cp "$2" "$dir/t.azr"
  (cd "$dir" && printf "$3\n\n\n7\n" | run_limited 900 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
}
M="$WORK/models"

# (a) byte identity
for mode in 1 3; do
  run a_plain_$mode "$M/aC_heavy.azr" $mode
  run a_zero_$mode "$M/aC_heavy0.azr" $mode
  n="$(ls "$WORK/a_plain_$mode/output" | wc -l)"
  if [ "$n" -lt 3 ]; then bad "(a) mode $mode: only $n output files"; continue; fi
  if same_tree "$WORK/a_plain_$mode/output" "$WORK/a_zero_$mode/output" &&
     samelog "$WORK/a_plain_$mode/run.log" "$WORK/a_zero_$mode/run.log"; then
    ok "(a) alpha + 12C, mode $mode: column 32 = 0 gives byte-identical output ($n files) and log"
  else bad "(a) alpha + 12C, mode $mode: output differs with an explicit column 32 = 0"; fi
done
mkdir -p "$WORK/n13/data"
cp "$HERE/../13N/data/"* "$WORK/n13/data/"
cp "$HERE/../13N/13N.azr" "$WORK/n13/13N.azr"
awk '/^<\/levels>/ { inl = 0 } { if (inl && NF >= 31) print $0 " 0"; else print } /^<levels>/ { inl = 1 }' \
  "$WORK/n13/13N.azr" > "$WORK/n13/13N_zero.azr"
if [ "$(awk 'NF == 32' "$WORK/n13/13N_zero.azr" | wc -l)" -lt 5 ]; then bad "(a) 13N: no 32-column lines written"; fi
run n13_plain "$WORK/n13/13N.azr" 1
run n13_zero "$WORK/n13/13N_zero.azr" 1
n="$(ls "$WORK/n13_plain/output" | wc -l)"
if [ "$n" -ge 3 ] && [ -s "$WORK/n13_plain/output/chiSquared.out" ] &&
   same_tree "$WORK/n13_plain/output" "$WORK/n13_zero/output" &&
   samelog "$WORK/n13_plain/run.log" "$WORK/n13_zero/run.log"; then
  ok "(a) tests/13N, mode 1: column 32 = 0 gives byte-identical output ($n files) and log"
else bad "(a) tests/13N, mode 1: output differs with an explicit column 32 = 0 (or no output)"; fi

# (b) moving the excitation between the partners
for mode in 1 3; do
  run b_heavy_$mode "$M/aC_heavy.azr" $mode
  run b_light_$mode "$M/aC_light.azr" $mode
  read -r w nf nn <<< "$(an compare "$WORK/b_heavy_$mode/output" "$WORK/b_light_$mode/output" 'AZUREOut_*' chiSquared.out)"
  if le "$w" 1e-10 && [ "${nf:-0}" -ge 1 ]; then
    ok "(b) mode $mode: excitation on the light vs the heavy partner, max rel difference $w ($nf files, $nn numbers)"
  else bad "(b) mode $mode: light vs heavy partner excitation differ, max rel '$w'"; fi
done
for seg in 'aa=1_R=2' 'aa=2_R=1' 'aa=1_R=1' 'aa=2_R=2'; do
  [ -s "$WORK/b_heavy_3/output/AZUREOut_$seg.extrap" ] || bad "(b) no AZUREOut_$seg.extrap"
done
run b_noex "$M/aC_light_noex.azr" 3
read -r w rest <<< "$(an compare "$WORK/b_heavy_3/output" "$WORK/b_noex/output" 'AZUREOut_aa=1_R=2.extrap')"
if gt "$w" 0.01; then ok "(b) without column 32 the light version differs by up to $w: e1 is what moves the threshold"
else bad "(b) without column 32 the light version differs by only '$w': the check does not discriminate"; fi

# (c) mutual excitation, identical bosons
run c_new "$M/cc_new.azr" 3
run c_old "$M/cc_old.azr" 3
read -r nbad ng nc <<< "$(an channels "$WORK/c_new/run.log" 1 2 4 4 1)"
if [ "$nbad" = 0 ] && [ "${ng:-0}" -ge 10 ]; then
  ok "(c) e1 = e2: the $ng hard-sphere J^pi groups hold exactly the $nc channels with l + s even"
else bad "(c) e1 = e2: hard-sphere channels differ from the l + s even set in '$nbad' groups"; fi
read -r nbad ng nc_old <<< "$(an channels "$WORK/c_old/run.log" 0 2 4 4 1)"
if [ "$nbad" = 0 ] && [ "${nc_old:-0}" -gt "${nc:-0}" ]; then
  ok "(c) e1 = 0, e2 = 2E* (old way): not identical, all $nc_old channels"
else bad "(c) old way: '$nbad' groups differ from the full channel set ($nc_old vs $nc channels)"; fi
read -r w np <<< "$(an symmetry "$WORK/c_new/output/AZUREOut_aa=1_R=1.extrap")"
if le "$w" 1e-6 && [ "${np:-0}" -ge 50 ]; then ok "(c) e1 = e2: elastic c.m. cross section symmetric about 90 deg ($w, $np angle pairs)"
else bad "(c) e1 = e2: elastic cross section not symmetric about 90 deg: '$w'"; fi
read -r w np <<< "$(an symmetry "$WORK/c_old/output/AZUREOut_aa=1_R=1.extrap")"
if gt "$w" 0.1; then ok "(c) old way: asymmetric by up to $w"
else bad "(c) old way: asymmetry only '$w': the check does not discriminate"; fi
run c_bad "$M/cc_bad.azr" 3
if grep -q "^ERROR: Identical-particle pair (Z=6, A=12, j=2) has a channel with J=4, L=2, S=3" "$WORK/c_bad/run.log" &&
   [ -z "$(ls "$WORK/c_bad/output")" ]; then ok "(c) e1 = e2: an l + s odd channel on a level is refused"
else bad "(c) e1 = e2: the l + s odd channel was not refused"; fi

# (d) refusals
refused() {   # <name> <model> <message start> [mode, default 3]
  run "$1" "$M/$2.azr" "${4:-3}"
  if grep -q "^ERROR: $3" "$WORK/$1/run.log" && [ -z "$(ls "$WORK/$1/output")" ]; then ok "(d) $2: ERROR, no output"
  else bad "(d) $2: not refused ('ERROR: $3' missing or output written)"; sed 's/^/        /' "$WORK/$1/run.log" | tail -5; fi
}
refused d_gamma bad_gamma "pair 3 is a gamma pair (type 10) but has a light-particle excitation energy"
refused d_negative bad_negative "pair 2 has a negative light-particle excitation energy"
refused d_mixed bad_mixed "pair 2 has two light-particle excitation energies"
refused d_upos bad_upos "segment #3 (UPOS): exit pair 2 has an excited light particle" 1
run d_upos_heavy "$M/upos_heavy.azr" 1
if [ -s "$WORK/d_upos_heavy/output/chiSquared.out" ] && ! grep -q "^ERROR" "$WORK/d_upos_heavy/run.log"; then
  ok "(d) the same UPOS segment with the excitation on the heavy partner runs"
else bad "(d) the UPOS segment with the excitation on the heavy partner did not run"; fi

if [ "$failures" -eq 0 ]; then echo "light_excitation: all checks passed"; exit 0; fi
echo "light_excitation: $failures check(s) FAILED"
exit 1
