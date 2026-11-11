#!/usr/bin/env bash
#
# A <targetInt> line shared by segments of DIFFERENT entrance pairs converts its lab energies with each point's own
# pair (audit 2026-10-25, tier 3 entry 3 / checklist t3-22; 2026-10-27).
#
# Every segment listed on a <targetInt> line, and every component of a SUM/RATIO or natural-target segment, points to
# the same TargetEffect.  Up to dev 1599725 the Gaussian sigma (and a `sigmalab` expression, the beam profile and, since
# 5ef20e5, the Bohr straggling width) was converted lab -> c.m. ONCE, with the factor m2/(m1+m2) of the first segment
# on the line, so every other pair sharing it got a kernel scaled by the wrong mass ratio (n + 20Ne vs n + 22Ne: 0.43 %;
# a beam-profile centroid: tens of keV).  Now the line keeps its lab values and each point converts them with its own
# factor (EPoint::SetTargetKinematics).
#
# Model: tests/natural_neon (n + 20Ne pair 2, n + 22Ne pair 4) with two plain angle-integrated elastic data segments at
# the same 41 lab energies (0.85 - 1.25 MeV, across the 21Ne / 23Ne resonances) and a 30 keV lab beam sigma:
#   shared: one beam record for both data sets;  split: one beam record per data set;  wrong: split with data set 2's
#   spread scaled by f(20Ne)/f(22Ne), which is what the old build computed for the shared line.
#   (a) shared == split for both segments to 1e-9 relative;
#   (b) wrong differs from split on segment 2 by more than 1e-5 somewhere (4e-5 here: so (a) has teeth), segment 1 identical.
# dev 1599725: (a) fails on segment 2 by the (b) amount (up to ~1e-3 near the resonances).
#
#   ./tests/shared_effect_pairs/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/../natural_neon"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/shared_effect_pairs.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

SIG=0.030
RATIO="$(awk 'BEGIN { f20 = 19.99244018 / (1.008664916 + 19.99244018); f22 = 21.991385109 / (1.00866491588 + 21.991385109); printf "%.10f", f20 / f22 }')"
SIG2="$(awk -v s="$SIG" -v r="$RATIO" 'BEGIN { printf "%.8f", s * r }')"

# make_model NAME BEAM1 SPREAD1 BEAM2 SPREAD2: the natural-neon project with two plain data sets (n + 20Ne, n + 22Ne) at the
# same lab energies, each folded with its beam record (a Gaussian lab spread, MeV; the same name twice: one record)
make_model() {
  local name="$1" b1="$2" s1="$3" b2="$4" s2="$5"
  mkdir -p "$WORK/$name/output" "$WORK/$name/checks" "$WORK/$name/data"
  awk 'BEGIN { for (i = 0; i <= 40; i++) printf "%.4f 0 2.0 0.1\n", 0.85 + i * 0.01 }' > "$WORK/$name/data/shared.dat"
  {
    for k in 1 2; do
      if [ $k = 1 ]; then r="20Ne(n,n)20Ne"; b="$b1"; else r="22Ne(n,n)22Ne"; b="$b2"; fi
      printf 'dataset s%s\n    reaction       %s\n    quantity       sigma\n    unit           b\n' "$k" "$r"
      printf '    file           data/shared.dat\n    columns        E_lab[MeV] _ value error\n'
      printf '    energy_window  0.5 1.5 MeV\n    norm           1 fixed\n    beam           %s\n    quadrature     points 30\n\n' "$b"
    done
  } > "$WORK/$name/data.txt"
  { printf 'beam %s\n    spread         %s MeV lab\n' "$b1" "$s1"
    [ "$b2" = "$b1" ] || printf '\nbeam %s\n    spread         %s MeV lab\n' "$b2" "$s2"; } > "$WORK/$name/targets.txt"
  awk -v d="$WORK/$name/data.txt" -v t="$WORK/$name/targets.txt" '
    { sub(/\r$/, "") }
    /^<data>$/ { print; while ((getline l < d) > 0) print l; skip = 1; next }
    /^<extrapolations>$/ { print; skip = 1; next }
    /^<targets>$/ { print; while ((getline l < t) > 0) print l; skip = 1; next }
    /^<\/data>$|^<\/extrapolations>$|^<\/targets>$/ { skip = 0 }
    !skip { print }' "$BASE/natural_neon.azr" > "$WORK/$name/t.azr"
}
run() { (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 900 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1); }
# yields NAME PAIR -> "E fit" rows of the data set with entrance = exit = PAIR
yields() { awk '{ sub(/\r$/, "") } NF >= 5 { printf "%.10f %.12e\n", $1, $4 }' "$WORK/$1/output/AZUREOut_aa=$2_R=$2.out" 2> /dev/null; }
maxrel() { paste "$1" "$2" | awk 'BEGIN { m = 0 } { if ($4 == 0) { if ($2 != 0) m = 1 } else { r = $2 / $4 - 1; if (r < 0) r = -r; if (r > m) m = r } } END { printf "%.3e", m }'; }

make_model shared B "$SIG" B "$SIG"
make_model split  B1 "$SIG" B2 "$SIG"
make_model wrong  B1 "$SIG" B2 "$SIG2"
for m in shared split wrong; do run "$m"; yields "$m" 2 > "$WORK/$m.y2"; yields "$m" 4 > "$WORK/$m.y4"; done
n2="$(wc -l < "$WORK/shared.y2" | tr -d ' ')"; n4="$(wc -l < "$WORK/shared.y4" | tr -d ' ')"
if [ "$n2" -eq 41 ] && [ "$n4" -eq 41 ]; then ok "the runs wrote 41 yields per segment (n + 20Ne, n + 22Ne)"; else bad "shared run wrote $n2 / $n4 yields, expected 41 / 41"; tail -3 "$WORK/shared/run.log"; fi

r1="$(maxrel "$WORK/shared.y2" "$WORK/split.y2")"; r2="$(maxrel "$WORK/shared.y4" "$WORK/split.y4")"
if awk -v a="$r1" -v b="$r2" 'BEGIN { exit !(a <= 1e-9 && b <= 1e-9) }'; then ok "(a) one line for both pairs == one line per pair: n + 20Ne to $r1, n + 22Ne to $r2"; else bad "(a) the shared line differs from one line per pair: n + 20Ne by $r1, n + 22Ne by $r2 (expected <= 1e-9: the 22Ne kernel was converted with 20Ne's mass ratio)"; fi
w1="$(maxrel "$WORK/wrong.y2" "$WORK/split.y2")"; w2="$(maxrel "$WORK/wrong.y4" "$WORK/split.y4")"
if awk -v a="$w1" -v b="$w2" 'BEGIN { exit !(a <= 1e-9 && b > 1e-5) }'; then ok "(b) the 20Ne-ratio kernel on the 22Ne segment moves its yields by up to $w2 (segment 1 untouched: $w1), so (a) has teeth"; else bad "(b) the wrong-ratio kernel moves the 22Ne yields by only $w2 (segment 1: $w1); expected > 1e-5 and 0"; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
