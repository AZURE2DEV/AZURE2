#!/usr/bin/env bash
#
# The engine runs a format-3 project (docs/upgrade/project_format3.md step 3b, G3): Config::ReadConfigFile reads it
# through records -> Validate -> <levels> tables -> Schema::Build, and the set-up consumes the specs.
#   1. T2's old original (tests/convert/fixtures/tests/upgrade/T2), converted with --exact-grids, gives the old file's
#      outputs byte for byte, param.par after the converter's parameter-name map (calculation, menu 1, Brune basis);
#      the conversion is the committed tests/upgrade/T2/T2.azr;
#   2. T2_true with the <modelData> section below, in both formats: equal chiSquared.out and modelData.out;
#   3. a data set named with 60 characters: its parameter names whole in param.par (no fixed-size name buffer);
#   4. --describe of the converted T2: the pair keys, the resolved data file, the parameters, the rows read and kept
#      (counted again here from the files and windows); --json parses;
#   5. a format-3 file with an error is refused with the reader's diagnostic;
#   6. the converted T2 at 1 and 2 threads: identical outputs;
#   7. a run without data (menu 3) does not need the data files: with the data directory away, T2's extrapolations give
#      the outputs of the run with data present, byte for byte (W_DATA_FILE_NOT_READ); a data run refuses (E_DATA_FILE).
#
#   ./tests/format3_engine/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
T2="$HERE/../upgrade/T2"
T2OLD="$HERE/../convert/fixtures/tests/upgrade/T2"   # the old originals of T2.azr and T2_true.azr (format3_tree.md)
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/format3_engine.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# setup NAME FILE: $WORK/NAME with T2's data, covariance and prior files and FILE as m.azr
setup() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  for f in data cov_diff.txt prior_cov.txt; do ln -s "$T2/$f" "$WORK/$1/$f"; done
  cp "$2" "$WORK/$1/m.azr"
  for f in "$(dirname "$2")"/*.f3.txt; do [ -f "$f" ] && cp "$f" "$WORK/$1/"; done  # the converter's side files
  return 0
}
# run NAME THREADS: a calculation (menu 1) of $WORK/NAME/m.azr
run() {
  (cd "$WORK/$1" && rm -rf output/* checks/* && printf '1\n\n\n7\n' |
     OMP_NUM_THREADS="$2" run_limited 600 "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
}
# convert NAME: $WORK/NAME/m.azr -> m3.azr (exact grids) with its report
convert() {
  (cd "$WORK/$1" && "$AZURE2_BIN" --convert m.azr m3.azr --to 3 --exact-grids --json-report report.json > convert.log 2>&1)
}

# ---- 1. T2 in both formats
setup old "$T2OLD/T2.azr"
run old 1
setup new "$T2OLD/T2.azr"
if convert new; then
  cmp -s "$WORK/new/m3.azr" "$T2/T2.azr" && ok "the committed T2.azr is the conversion of its original, byte for byte" ||
    bad "the committed T2.azr differs from the conversion of its original"
  mv "$WORK/new/m3.azr" "$WORK/new/m.azr"
  run new 1
  grep -q '^ERROR' "$WORK/new/run.log" && bad "T2 format 3: the engine refused it: $(grep -m1 '^ERROR' "$WORK/new/run.log")"
  same=1
  [ -s "$WORK/old/output/param.par" ] || bad "T2 (old file): no output ($(grep -m1 '^ERROR' "$WORK/old/run.log"))"
  for f in $(cd "$WORK/old/output" && ls); do
    [ "$f" = param.par ] || [ "$f" = normalizations.out ] && continue
    cmp -s "$WORK/old/output/$f" "$WORK/new/output/$f" || { same=0; bad "T2 format 3: output/$f differs"; }
  done
  [ "$same" = 1 ] && ok "T2 format 3: every output file but param.par and normalizations.out byte for byte the old file's"
  # param.par: the old names mapped by the converter's report.  One recorded difference (project_format3.md 8.1, the
  # director's decision on G2 item 4): a member of a varied shared normalization has a fixed parameter of its own,
  # <record>.value[<data set>], also where its old line did not vary its norm (T2's s03, whose own '1 0 0' the
  # converter drops with a note); that data set is then listed in normalizations.out too, and a follower's fixed
  # parameters take their start and step from the record, not from its own old line.  Nothing else may differ.
  if python3 - "$WORK/old/output" "$WORK/new/output" "$WORK/new/report.json" <<'PY'
import json, sys
old, new, report = sys.argv[1:4]
names = json.load(open(report))["parameters"]
a = [l.split() for l in open(old + "/param.par")]
b = [l.split() for l in open(new + "/param.par")]
a = [[names.get(t[0], t[0])] + t[1:] for t in a]
extra = [t for t in b if t and "[" in t[0] and t[0] not in [x[0] for x in a if x]]
followers = {t[0][t[0].index("[") + 1:-1] for t in extra}
b = [t for t in b if t not in extra]
same_rows = [x == y or (x and y and x[0] == y[0] and "[" in x[0]) for x, y in zip(a, b)]
na, nb = open(old + "/normalizations.out").read().splitlines(), open(new + "/normalizations.out").read().splitlines()
nb = [l for l in nb if not any(f"data/{f}." in l for f in followers if l not in na)]
print("    followers with a parameter of their own:", sorted(followers))
sys.exit(0 if len(a) == len(b) and all(same_rows) and na == nb and followers == {"s03_diff150"} else 1)
PY
  then ok "T2 format 3: param.par and normalizations.out equal after the name map, but for the follower s03's row"
  else bad "T2 format 3: param.par or normalizations.out differ after the name map"; fi
else
  bad "T2: the conversion failed: $(head -3 "$WORK/new/convert.log")"
fi

# ---- 2. model data in both formats
MD='<modelData>
datum wg
    quantity resonance_strength
    level 5
    J 3-
    reaction 15N(p,a)12C
    value 2.6 keV
    error 0.1 keV
datum gamma5
    quantity total_width
    level 5
    J 3-
    value 14 keV
    error 1 keV
</modelData>'
setup mdold "$T2OLD/T2_true.azr"
printf '%s\n' "$MD" >> "$WORK/mdold/m.azr"
setup mdnew "$WORK/mdold/m.azr"
run mdold 1
if convert mdnew; then
  mv "$WORK/mdnew/m3.azr" "$WORK/mdnew/m.azr"
  run mdnew 1
  for f in chiSquared.out modelData.out; do
    if [ -s "$WORK/mdold/output/$f" ] && cmp -s "$WORK/mdold/output/$f" "$WORK/mdnew/output/$f"; then
      ok "model data: $f equal in both formats"
    else
      bad "model data: $f differs (or is missing)"
    fi
  done
else
  bad "model data: the conversion failed: $(head -3 "$WORK/mdnew/convert.log")"
fi

# ---- 3. a 60-character data set name
long="d_$(printf 'x%.0s' $(seq 1 58))"
if [ -f "$WORK/new/m.azr" ]; then
  first="$(awk '$1 == "dataset" { print $2; exit }' "$WORK/new/m.azr")"
  setup long "$WORK/new/m.azr"
  sed -i "s/^dataset $first\\b/dataset $long/" "$WORK/long/m.azr"
  run long 1
  if grep -q "^$long.shift.offset " "$WORK/long/output/param.par" && grep -q "^$long.shift.sqrt " "$WORK/long/output/param.par"; then
    ok "a 60-character data set name: its parameters named whole in param.par"
  else
    bad "a 60-character data set name: $long.shift.offset not in param.par ($(grep -m1 '^ERROR' "$WORK/long/run.log"))"
  fi
fi

# ---- 4. --describe
if [ -f "$WORK/new/m.azr" ]; then
  d="$(cd "$WORK/new" && "$AZURE2_BIN" --describe m.azr)"
  if echo "$d" | grep -q "entrance pair 1 (p+15N), exit pair 2 (a+12C)" && echo "$d" | grep -q "data file $WORK/new/data/" &&
     echo "$d" | grep -q ".shift.offset (energy shift, fixed)"; then
    ok "--describe: pair keys with their pairs, the resolved data file, the parameters"
  else
    bad "--describe: $(echo "$d" | head -8)"
  fi
  if (cd "$WORK/new" && "$AZURE2_BIN" --describe m.azr --json) | python3 -c 'import json,sys; d=json.load(sys.stdin); sys.exit(0 if d["datasets"] and d["datasets"][0]["terms"][0]["entrance_pair"] > 0 else 1)'; then
    ok "--describe --json: a JSON document with the data sets and their pair keys"
  else
    bad "--describe --json: not a JSON document with the data sets"
  fi
  # The rows of each data set's measurement: read and kept as counted here from the file and the window.
  if (cd "$WORK/new" && "$AZURE2_BIN" --describe m.azr --json) | python3 -c '
import json, sys
d = json.load(sys.stdin)
bad = []
for s in d["datasets"]:
    if not s["active"]:
        continue
    lo, hi = s["E_lab_window_MeV"]
    alo, ahi = s.get("angle_window_deg", [None, None])
    read = kept = 0
    es = []
    for line in open(s["data_file_resolved"]):
        t = line.split()
        if not t or t[0].startswith("#"):
            continue
        read += 1
        e, a = float(t[0]), float(t[1])
        if lo <= e <= hi and (alo is None or alo <= a <= ahi):
            kept += 1
            es.append(e)
    want = (read, kept, [min(es), max(es)] if es else None)
    got = (s["rows_read"], s["rows_kept"], s.get("kept_E_lab_MeV"))
    if want != got:
        bad.append((s["name"], want, got))
if bad:
    print(bad)
sys.exit(1 if bad else 0)'; then
    ok "--describe --json: rows read and kept and the kept E_lab range of every data set, as counted from its file"
  else
    bad "--describe --json: rows read or kept differ from the file and the window"
  fi
  if echo "$d" | grep -Eq "^  rows: [0-9]+ read, [0-9]+ kept, E_lab "; then
    ok "--describe: the rows read and kept and their E_lab range"
  else
    bad "--describe: no rows line"
  fi
fi

# ---- 5. a format-3 file with an error
if [ -f "$WORK/new/m.azr" ]; then
  setup broken "$WORK/new/m.azr"
  sed -i '0,/^ *reaction /s/^\( *reaction \).*/\115N(p,d)14N/' "$WORK/broken/m.azr"
  run broken 1
  if grep -q 'E_UNKNOWN_PAIR' "$WORK/broken/run.log" && grep -q 'is refused' "$WORK/broken/run.log"; then
    ok "a format-3 file with an unknown pair: refused with E_UNKNOWN_PAIR"
  else
    bad "a format-3 file with an unknown pair: $(grep -m2 'ERROR' "$WORK/broken/run.log")"
  fi
fi

# ---- 6. 1 against 2 threads
if [ -f "$WORK/new/m.azr" ]; then
  setup two "$WORK/new/m.azr"
  run two 2
  if diff -rq "$WORK/new/output" "$WORK/two/output" > /dev/null; then
    ok "T2 format 3: 1 and 2 threads give identical outputs"
  else
    bad "T2 format 3: 1 and 2 threads differ: $(diff -rq "$WORK/new/output" "$WORK/two/output" | head -2)"
  fi
fi

# ---- 7. extrapolations without the data files
run3() {
  (cd "$WORK/$1" && rm -rf output/* checks/* && printf '3\nn\n\n\n7\n' |
     OMP_NUM_THREADS=1 run_limited 600 "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
}
if [ -f "$WORK/new/m.azr" ]; then
  # T2 has no extrapolation: one is added (15N(p,a)12C on a list of lab energies)
  awk '/^<\/extrapolations>/ { print "extrapolation xa\n    reaction       15N(p,a)12C\n    quantity       sigma\n    unit           b\n    energies       list 0.3 0.5 0.7 0.9 MeV lab" } { print }' \
    "$WORK/new/m.azr" > "$WORK/x3.azr"
  grep -q '^extrapolation xa' "$WORK/x3.azr" || bad "the extrapolation was not added (no <extrapolations> section)"
  setup x3data "$WORK/x3.azr"
  for f in "$WORK/new"/*.f3.txt; do [ -f "$f" ] && cp "$f" "$WORK/x3data/"; done
  run3 x3data
  setup x3nodata "$WORK/x3data/m.azr"
  rm "$WORK/x3nodata/data"
  run3 x3nodata
  n=$(ls "$WORK/x3data/output" | wc -l)
  if [ "$n" -gt 0 ] && diff -r "$WORK/x3data/output" "$WORK/x3nodata/output" > /dev/null &&
     grep -q 'W_DATA_FILE_NOT_READ' "$WORK/x3nodata/run.log" && ! grep -q '^ERROR' "$WORK/x3nodata/run.log"; then
    ok "extrapolations without the data files: the $n outputs of the run with data, byte for byte (a warning per file)"
  else
    bad "extrapolations without the data files: $(diff -rq "$WORK/x3data/output" "$WORK/x3nodata/output" | head -2) $(grep -m1 '^ERROR' "$WORK/x3nodata/run.log")"
  fi
  run x3nodata 1
  if grep -q 'E_DATA_FILE' "$WORK/x3nodata/run.log" && [ ! -s "$WORK/x3nodata/output/chiSquared.out" ]; then
    ok "a data run without the data files: refused with E_DATA_FILE"
  else
    bad "a data run without the data files: not refused with E_DATA_FILE"
  fi
fi

if [ "$failures" -eq 0 ]; then echo "PASS format3_engine"; exit 0; fi
echo "FAIL format3_engine ($failures)"; exit 1
