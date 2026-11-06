#!/usr/bin/env bash
#
# Separation energies against masses (2026-10-29, CNuc::CheckMassSeparationConsistency).
#
# Two particle pairs of one compound nucleus must satisfy S_b - S_a = (m1b + m2b - m1a - m2a) c^2.  AZURE2 takes the
# thresholds from the separation energies and the kinematics from the masses, so a mismatch changes nothing in the
# calculation, but it is how typos and mixed mass conventions show up: a scan of 6,326 archive projects found 879
# with mismatches above 100 keV (one 17O+a project at 73 MeV).  The engine now says so when it reads the pairs.
#
# Model: tests/relativistic/models/li7pn_ang_rel.azr (p + 7Li and n + 7Be in 8Be, atomic masses, S_n = S_p - Q from the
# same masses: consistent to < 1 eV), mode 1.
#   (a) as it is: no warning, and the run completes;
#   (b) S of pair 2 lowered by 5 keV (raising it would move the n + 7Be threshold past a data point): one WARNING
#       naming pairs 2 and 1 and "-5.0 keV"; the run still completes;
#   (c) the proton mass of pair 1 replaced by the bare proton (1.00727646688 u): a WARNING of -511.0 keV with the
#       electron-mass hint;
#   (d) every mass replaced by its mass number: no WARNING, one NOTE that the pairs use mass numbers;
#   (e) the model of (b) with mass numbers, filled by pyazr's AzrModel.fill_masses (AME2020): no WARNING and no
#       NOTE, and the run completes (skipped without python3).
#
#   ./tests/mass_separation_check/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../relativistic/models"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/mass_separation_check.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# model NAME AWK-PROGRAM: li7pn_ang_rel.azr with the <levels> lines rewritten by the program (fields as in NucLine:
# $6 pair key, $18 m1, $19 m2, $22 entrance S, $23 S), run in mode 1
model() {
  local d="$WORK/$1"; mkdir -p "$d/output" "$d/checks"; cp -r "$SRC/data" "$d/"
  awk -v OFS=' ' "BEGIN{inl=0} /<levels>/{inl=1; print; next} /<\\/levels>/{inl=0} inl && NF > 25 { $2 } {print}" \
    "$SRC/li7pn_ang_rel.azr" > "$d/t.azr"
  (cd "$d" && printf '1\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
}
warnings() { grep -c "separation energies and masses disagree" "$WORK/$1/run.log"; }
completed() { [ -s "$WORK/$1/output/chiSquared.out" ]; }

model plain ""
if [ "$(warnings plain)" -eq 0 ] && completed plain; then ok "(a) consistent masses and separation energies: no warning, run completes"
else bad "(a) consistent model: $(warnings plain) warning(s), chiSquared.out $(ls "$WORK/plain/output" | grep -c chiSquared)"; fi

model shifted 'if ($6 == 2) { $22 = sprintf("%.10f", $22 - 0.005); $23 = sprintf("%.10f", $23 - 0.005) }'
w="$(grep -m1 "separation energies and masses disagree" "$WORK/shifted/run.log")"
if [ "$(warnings shifted)" -eq 1 ] && echo "$w" | grep -q "pair 2 .*pair 1 .*(-5.0 keV)" && completed shifted; then
  ok "(b) S of pair 2 lowered by 5 keV: one warning (-5.0 keV), run completes"
else bad "(b) 5 keV shift: $(warnings shifted) warning(s): ${w:0:200}"; fi

model bare 'if ($6 == 1) { $18 = "1.00727646688" }'
w="$(grep -m1 "separation energies and masses disagree" "$WORK/bare/run.log")"
if echo "$w" | grep -q "(-511.0 keV).*1 electron mass"; then ok "(c) bare-proton mass among atomic masses: -511.0 keV with the electron-mass hint"
else bad "(c) bare proton: ${w:0:300}"; fi

model massnum '$18 = sprintf("%d", $18 + 0.5); $19 = sprintf("%d", $19 + 0.5)'
if [ "$(warnings massnum)" -eq 0 ] && grep -q "use mass numbers as masses" "$WORK/massnum/run.log"; then
  ok "(d) mass numbers as masses: no warning, one note"
else bad "(d) mass numbers: $(warnings massnum) warning(s), note $(grep -c 'mass numbers as masses' "$WORK/massnum/run.log")"; fi

PY="$(command -v python3 || true)"
if [ -n "$PY" ]; then
  d="$WORK/filled"; mkdir -p "$d/output" "$d/checks"; cp -r "$SRC/data" "$d/"
  # mass numbers and a 5 keV S error, then the fill from the table
  awk -v OFS=' ' 'BEGIN{inl=0} /<levels>/{inl=1; print; next} /<\/levels>/{inl=0}
       inl && NF > 25 { $18 = sprintf("%d", $18 + 0.5); $19 = sprintf("%d", $19 + 0.5);
                        if ($6 == 2) $23 = sprintf("%.10f", $23 - 0.005) } {print}' \
    "$SRC/li7pn_ang_rel.azr" > "$d/numbers.azr"
  if "$PY" - "$HERE/../.." "$d/numbers.azr" "$d/t.azr" > "$d/fill.log" 2>&1 <<'PYEOF'
# azrfile.py by path: no package import, so no numpy or compiled module is needed
import importlib.util, os, sys
spec = importlib.util.spec_from_file_location("azrfile", os.path.join(sys.argv[1], "pyazr", "azrfile.py"))
azrfile = importlib.util.module_from_spec(spec)
spec.loader.exec_module(azrfile)
m = azrfile.AzrModel.from_file(sys.argv[2])
print(m.fill_masses())
m.write(sys.argv[3])
PYEOF
  then
    (cd "$d" && printf '1\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
    if [ "$(warnings filled)" -eq 0 ] && ! grep -q "mass numbers as masses" "$d/run.log" && completed filled; then
      ok "(e) mass numbers filled by pyazr from AME2020: no warning, no note, run completes"
    else bad "(e) pyazr fill: $(warnings filled) warning(s), note $(grep -c 'mass numbers as masses' "$d/run.log")"; fi
  else bad "(e) pyazr fill_masses failed: $(tail -3 "$d/fill.log")"; fi
else
  echo "  skip  (e) no python3"
fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
