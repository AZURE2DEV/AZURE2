#!/usr/bin/env bash
#
# Progress bars of the external-capture amplitude loops (2026-11-06).
#
# EData::CalculateECAmplitudes and EData::InitializeComponentSegments compute the amplitudes of a segment's points
# in an OpenMP loop and redraw a 25-cell bar ("\r\tSegment #<key> [****     ] 36%").  Up to dev 0119771 every thread
# incremented the shared point counter without synchronization and wrote its bar to the output stream piece by
# piece, without a lock: the bars of several threads interleaved character by character ("Segment #\r\tSegment #",
# "Segment # 11   [****"), and a bar's fragment could land on a line of its own after another thread's newline.
# Besides the garbled log (and the flaky tests/light_excitation log compare on the macOS and Windows runners), the
# GUI's output stream is not thread-safe.  Now the counter is atomic and one thread at a time draws a whole bar.
#
# Model: tests/13N (external capture on every capture segment), and the same with segment 1 a SUM of its own term
# and one identical capture component (the component-segment loop), mode 1, OMP_NUM_THREADS=8 (THREADS), RUNS runs
# each (default 3).
#   (a) every progress line of the log, split at \r, is a sequence of whole bars
#       "\tSegment #<key> [<25 cells>] <n>%", the last one 100 %;
#   (b) no line outside them is made only of bar characters (a stray fragment);
#   (c) one finished bar per EC segment, in segment order, and the same chi2 as the single-thread run.
#
#   ./tests/progress_bars/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
THREADS="${THREADS:-8}"
RUNS="${RUNS:-3}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PY="$(command -v python3 || command -v python)"
SRC="$HERE/../13N"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/progress_bars.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# models: plain 13N, and segment 1 (artemov, 12C(p,g)) a SUM of itself and one identical capture component
cp "$SRC/13N.azr" "$WORK/plain.azr"
"$PY" -I - "$SRC/13N.azr" "$WORK/sum.azr" <<'PYEOF' || { bad "models not built"; exit 1; }
import sys
lines = open(sys.argv[1]).read().split("\n")
i0 = lines.index("<segmentsData>")
f = lines[i0 + 1].split()
if f[14] != "data/artemov.dat" or f[15:17] != ["0", "0"]:
    sys.exit("segment 1 is not the plain artemov line")
lines[i0 + 1] = " ".join(f[:15]) + " 1 0 1 1 2 -999 0"   # SUM, one component: entrance 1, exit 2, same angles
open(sys.argv[2], "w").write("\n".join(lines))
PYEOF

# run NAME MODEL THREADS
run() {
  local d="$WORK/$1"
  mkdir -p "$d/output" "$d/checks"; ln -s "$SRC/data" "$d/data"; cp "$WORK/$2.azr" "$d/t.azr"
  (cd "$d" && printf '1\n\n\n7\n' | OMP_NUM_THREADS="$3" run_limited 900 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
}
chi2() { grep -h 'Total Chi-Squared' "$WORK/$1/run.log" | tail -1 | awk '{ print $NF }'; }

# analyze LOG -> "nbars nbad nstray keys" ; prints the offending fragments to stderr
analyze() {
  "$PY" -I -c '
import re, sys
# thm: the segment number is right-aligned in its setw(12) field (dev: left), hence " *".
bar = re.compile(rb"\tSegment # *(-?[0-9]+(?: \(Total: [0-9]+\))?) +\[[* ]{25}\] ([0-9]+)%")
chars = set(b"Segment#[]*% \t0123456789():Tola")
done, bad, stray = [], [], []
for line in open(sys.argv[1], "rb").read().split(b"\n"):
    if b"\r" in line or line.startswith(b"\tSegment #"):
        frags = line.split(b"\r")
        ms = [bar.fullmatch(f) for f in frags]
        if not all(ms) or ms[-1].group(2) != b"100" or len({m.group(1) for m in ms}) != 1:
            bad.append(line)
        else:
            done.append(ms[-1].group(1).split()[0].decode())
    elif line.strip() and set(line) <= chars and not line.startswith(b"Segment #"):
        stray.append(line)
for l in bad[:3] + stray[:3]:
    sys.stderr.write("        %r\n" % l[:160])
print(len(done), len(bad), len(stray), ",".join(done) or "-")' "$1"
}

for model in plain sum; do
  run "${model}_1" "$model" 1
  ref="$(chi2 "${model}_1")"
  read -r n1 b1 s1 keys1 <<< "$(analyze "$WORK/${model}_1/run.log" 2> /dev/null)"
  if [ -z "$ref" ] || [ "${n1:-0}" -lt 1 ] || [ "$b1" -ne 0 ] || [ "$s1" -ne 0 ]; then
    bad "($model) single-thread reference: chi2 '$ref', $n1 bars, $b1 bad, $s1 stray"; continue
  fi
  for r in $(seq 1 "$RUNS"); do
    name="${model}_${THREADS}t_$r"
    run "$name" "$model" "$THREADS"
    read -r n b s keys <<< "$(analyze "$WORK/$name/run.log")"
    c="$(chi2 "$name")"
    if [ "$b" -eq 0 ]; then ok "(a) $model, $THREADS threads, run $r: $n progress lines, every fragment a whole bar ending at 100 %"
    else bad "(a) $model, $THREADS threads, run $r: $b of $((n + b)) progress lines garbled"; fi
    if [ "$s" -eq 0 ]; then ok "(b) $model, run $r: no stray bar fragments"
    else bad "(b) $model, run $r: $s stray bar fragments on lines of their own"; fi
    if [ "$keys" = "$keys1" ] && [ "$c" = "$ref" ]; then ok "(c) $model, run $r: bars for segments $keys, chi2 $c as with 1 thread"
    else bad "(c) $model, run $r: bars for '$keys' (1 thread: '$keys1'), chi2 '$c' (1 thread: '$ref')"; fi
  done
done
# the SUM model must have exercised the component loop
if grep -q "Initializing Component Segments (1 components)" "$WORK/sum_1/run.log" 2> /dev/null ||
   grep -qE "Component Segments \([1-9]" "$WORK/sum_1/run.log" 2> /dev/null; then ok "the SUM model computes component EC amplitudes"
else bad "the SUM model has no component segment with external capture"; grep -m3 -i component "$WORK/sum_1/run.log"; fi

[ "$failures" -eq 0 ] && echo "progress_bars: all checks passed" || echo "progress_bars: $failures check(s) FAILED"
exit $((failures > 0))
