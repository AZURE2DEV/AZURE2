#!/usr/bin/env bash
#
# A phase-shift segment is one (J, l, s) partial wave: the "spin <s>" keyword, and an ambiguous (J, l) is refused
# (audit 2026-10-25, tier 3 / checklist t3-23; 2026-10-27).
#
# A phase-shift segment is specified by (J, l).  Up to dev 62c7ed2 the diagonal U-matrix elements of EVERY channel spin
# with that (J, l) were accumulated in one TempTMatrix and half the argument of the sum reported: for a target with spin
# (p + 3He 1P1 / 3P1, p + 7Be J = 1 l = 1 with s = 1 and 2) that is not the phase shift of either wave.  Now the segment
# line may carry "spin <s>", the phase branch keeps only that channel spin, and EData::Initialize refuses a segment whose
# (J, l) has several channel spins in the level scheme without the keyword, or whose spin does not exist.
#
# Model: tests/scattering_length/be7p_fit.azr (p + 7Be, s = 1 or 2) plus two 1+ levels with l = 1 channels: level A
# (Ex 1.2) couples only to s = 1, level B (Ex 1.7) only to s = 2 (the other spin's width is 0), so the s = 1 wave is
# unaffected by B and the s = 2 one by A.  Truth for each spin: the same model with the other spin's channel lines
# removed, which is unambiguous and the same on every build.  Phase-shift data segment J = 1, l = 1 at lab 0.3 - 2.0 MeV.
#   (a) "spin 1" on the full model == the s = 1-only model at every energy (1e-9);  (b) the same for "spin 2";
#   (c) the two phases differ from each other by more than 1 degree somewhere (teeth);
#   (d) the full model without "spin" stops with an ERROR naming s = 1 and 2 and writes no chiSquared.out;
#   (e) "spin 1.5" (no such channel) stops with an ERROR.
# dev 62c7ed2: (a) and (b) fail (the keyword is ignored, both spins summed), (d) and (e) fail (it runs).
#
#   ./tests/phase_channel_spin/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/../scattering_length"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/phase_channel_spin.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# level lines (the file stores 2s and 2l): 1+ level A at Ex 1.2 (s = 1 l = 1 width 60 keV, s = 2 width 0),
# level B at Ex 1.7 (s = 1 width 0, s = 2 width 90 keV); pair columns copied from the fixture's p + 7Be lines
TAIL="0.5 1 1.5 -1 0 1.00783 7.01693 1 4 0.1375 0.1375 0 0 0.0 0 4.3 0 0 0"
A1="1 1 1.200000 0 1 1 2 2 4 1 0 60000 $TAIL"
A2="1 1 1.200000 0 1 1 4 2 4 1 0 0 $TAIL"
B1="1 1 1.700000 0 1 1 2 2 5 1 0 0 $TAIL"
B2="1 1 1.700000 0 1 1 4 2 5 1 0 90000 $TAIL"

# make_model NAME KEYWORD-TAIL LEVELS...: the fixture with its data segments replaced by one J = 1, l = 1 phase-shift
# segment (tail appended), its test segments removed, and the given level lines added
make_model() {
  local name="$1" tail="$2"; shift 2
  mkdir -p "$WORK/$name/output" "$WORK/$name/checks" "$WORK/$name/data"
  awk 'BEGIN { for (i = 0; i <= 34; i++) printf "%.3f 0.0 10.0 0.5\n", 0.3 + i * 0.05 }' > "$WORK/$name/data/p11.dat"
  printf '%s\n' "$@" > "$WORK/$name/levels"
  awk -v levelsfile="$WORK/$name/levels" -v seg="1 1 1 0 10 0 180 2 1 1 1 0 0 0 0 0 data/p11.dat 0 0$tail" '
    { sub(/\r$/, "") }
    /^<\/levels>$/ { while ((getline l < levelsfile) > 0) print l }
    /^<segmentsData>$/ { print; print seg; skip = 1; next }
    /^<segmentsTest>$/ { print; skip = 1; next }
    /^<\/segmentsData>$|^<\/segmentsTest>$/ { skip = 0 }
    !skip { print }' "$BASE/be7p_fit.azr" > "$WORK/$name/t.azr"
}
run() { (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1); }
phases() { awk '{ sub(/\r$/, "") } NF >= 5 { printf "%.10f %.10f\n", $1, $4 }' "$WORK/$1/output/AZUREOut_aa=1_R=1.out" 2> /dev/null; }
maxdiff() { paste "$1" "$2" | awk 'BEGIN { m = 0 } { d = $2 - $4; if (d < 0) d = -d; if (d > m) m = d } END { printf "%.3e", m }'; }

make_model full1 " spin 1" "$A1" "$A2" "$B1" "$B2"
make_model full2 " spin 2" "$A1" "$A2" "$B1" "$B2"
make_model only1 "" "$A1" "$B1"
make_model only2 "" "$A2" "$B2"
make_model nospin "" "$A1" "$A2" "$B1" "$B2"
make_model badspin " spin 1.5" "$A1" "$A2" "$B1" "$B2"
for m in full1 full2 only1 only2 nospin badspin; do run "$m"; phases "$m" > "$WORK/$m.ph"; done
n="$(wc -l < "$WORK/only1.ph" | tr -d ' ')"
if [ "$n" -eq 35 ]; then ok "the single-spin models wrote the 35 phase shifts"; else bad "only1 wrote $n phases, expected 35"; tail -3 "$WORK/only1/run.log"; fi

d1="$(maxdiff "$WORK/full1.ph" "$WORK/only1.ph")"; d2="$(maxdiff "$WORK/full2.ph" "$WORK/only2.ph")"
n1="$(wc -l < "$WORK/full1.ph" | tr -d ' ')"; n2="$(wc -l < "$WORK/full2.ph" | tr -d ' ')"
if [ "$n1" -eq 35 ] && awk -v d="$d1" 'BEGIN { exit !(d <= 1e-9) }'; then ok "(a) \"spin 1\" on the full model == the s = 1 wave alone ($d1 deg)"; else bad "(a) \"spin 1\" differs from the s = 1 wave alone by $d1 deg ($n1 phases; expected <= 1e-9)"; tail -2 "$WORK/full1/run.log"; fi
if [ "$n2" -eq 35 ] && awk -v d="$d2" 'BEGIN { exit !(d <= 1e-9) }'; then ok "(b) \"spin 2\" on the full model == the s = 2 wave alone ($d2 deg)"; else bad "(b) \"spin 2\" differs from the s = 2 wave alone by $d2 deg ($n2 phases; expected <= 1e-9)"; fi
d12="$(maxdiff "$WORK/only1.ph" "$WORK/only2.ph")"
if awk -v d="$d12" 'BEGIN { exit !(d > 1.0) }'; then ok "(c) the two waves differ by up to $d12 deg, so (a) and (b) have teeth"; else bad "(c) the two waves differ by only $d12 deg"; fi
if [ ! -f "$WORK/nospin/output/chiSquared.out" ] && grep -q "ERROR: phase-shift segment #1 (J = 1, l = 1) has 2 channel spins" "$WORK/nospin/run.log"; then ok "(d) the ambiguous segment without \"spin\" is refused: $(grep -o 's = 1 2' "$WORK/nospin/run.log" | head -1)"; else bad "(d) the ambiguous segment without \"spin\" was not refused"; grep -i "error" "$WORK/nospin/run.log" | head -2; fi
if [ ! -f "$WORK/badspin/output/chiSquared.out" ] && grep -q "ERROR: phase-shift segment #1 (J = 1, l = 1) asks for channel spin s = 1.5" "$WORK/badspin/run.log"; then ok "(e) \"spin 1.5\" (no such channel) is refused"; else bad "(e) \"spin 1.5\" was not refused"; grep -i "error" "$WORK/badspin/run.log" | head -2; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
