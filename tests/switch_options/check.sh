#!/usr/bin/env bash
#
# The switches that can move a number are options of the program (docs/upgrade/evaluation.md section 9), and an
# option has exactly the effect of its environment variable.  For each of the fourteen, on a project where it acts:
#
#   1. the run with the variable and the run with the option write the same output files, byte for byte, and
#      --describe names the value with its source (environment, option);
#   2. the run without either differs, or the line says "the same files as without" (a route equal to rounding, or a
#      switch whose effect this project cannot show);
#   3. precedence: the option over the options file over the variable (--refine-quadrature, K = 4 / 3 / 2), the
#      --options-file option over AZURE_OPTIONS_FILE;
#   4. refusals alike: --refine-quadrature=65 and AZURE_REFINE_QUADRATURE=65 exit 2 with the same message under their
#      own names (also from the options file); a value option without its value exits 2; an empty SRIM directory
#      refuses a project with a target alike by option and by variable;
#   5. --threads=N: N >= 1 (0, x, a bare --threads exit 2); the run's header and --describe say the count and its
#      source, the option over OMP_NUM_THREADS; --threads=1 computes what OMP_NUM_THREADS=1 does.
#
#   ./tests/switch_options/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
TESTS="$(cd "$HERE/.." && pwd)"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TABLE="${AZURE2_COULOMB_TABLE:-}"
unset AZURE2_COULOMB_TABLE AZURE2_SRIM_DATA_DIR AZURE_OPTIONS_FILE AZURE_REFINE_QUADRATURE AZURE2_EC_TABLE_NODES \
  AZURE2_SWAP_COUPLING_ORDER AZURE_BOUND_SHIFT_NUMERICAL AZURE_NO_NEUTRAL_CLOSED_FORM AZURE_SHIFT_DE_1EV \
  AZURE_LEVEL_SPACE AZURE2_THM_FD AZURE_GRAD_ADJOINT AZURE_GRAD_DENSE AZURE_TE_WEIGHTS_FD

WORK="$(mktemp -d "${TMPDIR:-/tmp}/switch_options.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "${WORK:?}"' EXIT
export HOME="$WORK/home"  # no ~/.azure2/coulomb_table.bin: without a table file the runs use the exact functions
mkdir -p "$HOME"

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

CALC=calculate
FIT=fit
EXTRAP=extrapolate
# run NAME FOLDER AZR MENU [VAR=VALUE ...] -- [FLAGS ...]: the project FOLDER/AZR (a folder under tests/) in a fresh
# copy, the command COMMAND, the variables and flags given; output in $WORK/NAME/output, the log in $WORK/NAME/log.
run() {
  local name="$1" folder="$2" azr="$3" command="$4"; shift 4
  local vars=()
  while [ $# -gt 0 ] && [ "$1" != -- ]; do vars+=("$1"); shift; done
  [ $# -gt 0 ] && shift
  local dir="${WORK:?}/${name:?}"
  rm -rf "$dir"; mkdir -p "$dir"
  cp -R "$TESTS/$folder/." "$dir/"
  rm -rf "$dir/output" "$dir/checks"; mkdir -p "$dir/output" "$dir/checks"
  (cd "$dir" && env "${vars[@]}" "$AZURE2_BIN" $command --no-gui "$@" "$azr" > log 2>&1)
  echo $? > "$dir/status"
}
same() { diff -r "$WORK/$1/output" "$WORK/$2/output" > /dev/null 2>&1; }
filled() { [ "$(cat "$WORK/$1/status")" = 0 ] && [ -n "$(ls "$WORK/$1/output")" ]; }

# switch VARIABLE VALUE OPTION FOLDER AZR MENU [FLAGS ...]: cases 1 and 2 for one switch.
switch() {
  local var="$1" value="$2" option="$3" folder="$4" azr="$5" command="$6"; shift 6
  local tag; tag="$(echo "$var" | tr 'A-Z' 'a-z')"
  run "$tag.base" "$folder" "$azr" "$command" -- "$@"
  run "$tag.var" "$folder" "$azr" "$command" "$var=$value" -- "$@"
  run "$tag.opt" "$folder" "$azr" "$command" -- "$option" "$@"
  if ! filled "$tag.var" || ! filled "$tag.opt" || ! filled "$tag.base"; then
    bad "$var: a run failed (status $(cat "$WORK/$tag.base/status") / $(cat "$WORK/$tag.var/status") / $(cat "$WORK/$tag.opt/status"))"
    return
  fi
  local moves="the same files as without"
  same "$tag.var" "$tag.base" || moves="differs from the run without ($(diff -rq "$WORK/$tag.var/output" "$WORK/$tag.base/output" | wc -l) files)"
  if same "$tag.var" "$tag.opt"; then ok "$var=$value = $option on $folder: identical files; $moves"
  else bad "$var=$value and $option on $folder: the files differ ($(diff -rq "$WORK/$tag.var/output" "$WORK/$tag.opt/output" | sed -n '1,2p'))"; fi
  local d
  d="$(cd "$WORK/$tag.opt" && "$AZURE2_BIN" --describe "$azr" "$option" 2>/dev/null | grep -F "($var) = ")"
  case "$d" in *"(option)") ;; *) bad "$option: --describe says '$d'";; esac
  d="$(cd "$WORK/$tag.var" && env "$var=$value" "$AZURE2_BIN" --describe "$azr" 2>/dev/null | grep -F "($var) = ")"
  case "$d" in *"(environment)") ;; *) bad "$var: --describe says '$d'";; esac
}

echo "1-2. variable = option, and the effect"
[ -n "$TABLE" ] && [ -f "$TABLE" ] || { echo "  FAIL  set AZURE2_COULOMB_TABLE to the Coulomb table file"; exit 1; }
switch AZURE2_COULOMB_TABLE "$TABLE" "--coulomb-table=$TABLE" upgrade/T1 T1.azr "$CALC"
printf -- '--no-brune\n' > "$WORK/options.txt"
switch AZURE_OPTIONS_FILE "$WORK/options.txt" "--options-file=$WORK/options.txt" upgrade/T1 T1.azr "$CALC"
SRIM="$(cd "$TESTS/../stopping/data/srim" && pwd)"
switch AZURE2_SRIM_DATA_DIR "$SRIM" "--srim-data-dir=$SRIM" stopping_yield/models stopping_yield.azr "$CALC"
switch AZURE_REFINE_QUADRATURE 2 --refine-quadrature=2 upgrade/T2 T2.azr "$CALC"
switch AZURE2_EC_TABLE_NODES 80 --ec-table-nodes=80 13N 13N.azr "$CALC"
switch AZURE2_SWAP_COUPLING_ORDER 1 --swap-coupling-order 7Li_p_ay 7Li_p_ay.azr "$CALC"
switch AZURE_BOUND_SHIFT_NUMERICAL 1 --bound-shift-numerical upgrade/T1 T1.azr "$CALC"
switch AZURE_NO_NEUTRAL_CLOSED_FORM 1 --no-neutral-closed-form reaction_rate/models n12C_truth.azr "$EXTRAP"
switch AZURE_SHIFT_DE_1EV 1 --shift-derivative-1ev upgrade/T1 T1.azr "$CALC"
switch AZURE_LEVEL_SPACE 1 --level-space upgrade/T2 T2.azr "$CALC"
switch AZURE_GRAD_ADJOINT 1 --gradient-adjoint upgrade/T2 T2.azr "$FIT" --use-gradient --max-calls=40
switch AZURE_GRAD_DENSE 1 --gradient-dense upgrade/T2 T2.azr "$FIT" --use-gradient --max-calls=40 --gradient-adjoint
switch AZURE_TE_WEIGHTS_FD 1 --fold-weights-by-differences upgrade/T2 T2.azr "$FIT" --use-gradient --max-calls=40
switch AZURE2_THM_FD 1 --thm-by-differences 18O_p_a_thm 18O_p_a_thm.azr "$FIT" --max-calls=20

echo "3. precedence"
printf -- '--refine-quadrature=3\n' > "$WORK/refine3.txt"
desc() { (cd "$WORK/azure_refine_quadrature.base" && env "$@" --describe T2.azr 2>/dev/null | grep -F "(AZURE_REFINE_QUADRATURE) = "); }
d1="$(desc AZURE_REFINE_QUADRATURE=2 AZURE_OPTIONS_FILE="$WORK/refine3.txt" "$AZURE2_BIN" --refine-quadrature=4)"
d2="$(desc AZURE_REFINE_QUADRATURE=2 AZURE_OPTIONS_FILE="$WORK/refine3.txt" "$AZURE2_BIN")"
d3="$(desc AZURE_REFINE_QUADRATURE=2 "$AZURE2_BIN")"
d4="$(desc AZURE_OPTIONS_FILE=/no/such/file "$AZURE2_BIN" "--options-file=$WORK/refine3.txt")"
if [ "$d1" = "  --refine-quadrature=K (AZURE_REFINE_QUADRATURE) = 4 (option)" ] &&
   [ "$d2" = "  --refine-quadrature=K (AZURE_REFINE_QUADRATURE) = 3 (options file)" ] &&
   [ "$d3" = "  --refine-quadrature=K (AZURE_REFINE_QUADRATURE) = 2 (environment)" ] &&
   [ "$d4" = "  --refine-quadrature=K (AZURE_REFINE_QUADRATURE) = 3 (options file)" ]; then
  ok "option 4 over options file 3 over variable 2; --options-file= over AZURE_OPTIONS_FILE"
else
  bad "precedence: '$d1' / '$d2' / '$d3' / '$d4'"
fi
run prec.all upgrade/T2 T2.azr "$CALC" AZURE_REFINE_QUADRATURE=2 AZURE_OPTIONS_FILE="$WORK/refine3.txt" -- --refine-quadrature=4
run prec.four upgrade/T2 T2.azr "$CALC" -- --refine-quadrature=4
filled prec.all && same prec.all prec.four && ok "the run with 4 / 3 / 2 given computes what --refine-quadrature=4 alone does" ||
  bad "the run with option, file and variable differs from the option alone"

echo "4. refusals alike"
refused() {  # refused WHAT TEXT VAR=VALUE... -- FLAGS...
  local what="$1" text="$2"; shift 2
  local vars=()
  while [ $# -gt 0 ] && [ "$1" != -- ]; do vars+=("$1"); shift; done
  [ $# -gt 0 ] && shift
  local out; out="$(cd "$WORK" && env "${vars[@]}" "$AZURE2_BIN" --no-gui "$@" no_such_project.azr < /dev/null 2>&1)"
  local rc=$?
  if [ "$rc" = 2 ] && grep -qF -- "$text" <<< "$out"; then ok "$what: exit 2, $text"
  else bad "$what: exit $rc, $(printf '%s' "$out" | sed -n '1,2p')"; fi
}
refused "--refine-quadrature=65" "ERROR: --refine-quadrature needs a whole number from 1 to 64, not '65'" -- --refine-quadrature=65
refused "AZURE_REFINE_QUADRATURE=65" "ERROR: AZURE_REFINE_QUADRATURE needs a whole number from 1 to 64, not '65'" \
  AZURE_REFINE_QUADRATURE=65 --
printf -- '--refine-quadrature=65\n' > "$WORK/refine65.txt"
refused "--refine-quadrature=65 in the options file" "ERROR: --refine-quadrature needs a whole number from 1 to 64" \
  AZURE_OPTIONS_FILE="$WORK/refine65.txt" --
refused "--srim-data-dir= without a value" "ERROR: --srim-data-dir= needs a value: --srim-data-dir=DIR" -- --srim-data-dir=
refused "--srim-data-dir alone" "ERROR: Unknown option --srim-data-dir." -- --srim-data-dir
refused "--ec-table-nodes= without a value" "ERROR: --ec-table-nodes= needs a value" -- --ec-table-nodes=
mkdir -p "$WORK/empty"
run srim.var stopping_yield/models stopping_yield.azr "$CALC" AZURE2_SRIM_DATA_DIR="$WORK/empty" --
run srim.opt stopping_yield/models stopping_yield.azr "$CALC" -- "--srim-data-dir=$WORK/empty"
if [ "$(cat "$WORK/srim.var/status")" != 0 ] && [ "$(cat "$WORK/srim.var/status")" = "$(cat "$WORK/srim.opt/status")" ] &&
   grep -q ERROR "$WORK/srim.var/log" &&
   [ "$(grep ERROR "$WORK/srim.var/log" | sed -n '1,3p')" = "$(grep ERROR "$WORK/srim.opt/log" | sed -n '1,3p')" ]; then
  ok "an empty SRIM directory refuses a project with a target alike by option and variable (exit $(cat "$WORK/srim.var/status"))"
else
  bad "an empty SRIM directory: exit $(cat "$WORK/srim.var/status") / $(cat "$WORK/srim.opt/status")"
fi

echo "5. --threads"
refused "--threads=0" "ERROR: --threads needs a whole number from 1 to 4096, not '0'" -- --threads=0
refused "--threads=x" "ERROR: --threads needs a whole number" -- --threads=x
refused "--threads alone" "ERROR: --threads needs a value" -- --threads
run thr.opt upgrade/T1 T1.azr "$CALC" OMP_NUM_THREADS=3 -- --threads=1
run thr.var upgrade/T1 T1.azr "$CALC" OMP_NUM_THREADS=1 --
if grep -qxF "Threads: 1 (--threads, over OMP_NUM_THREADS=3)" "$WORK/thr.opt/log" &&
   grep -qxF "Threads: 1 (OMP_NUM_THREADS)" "$WORK/thr.var/log" && filled thr.opt && same thr.opt thr.var; then
  ok "the header: the option over OMP_NUM_THREADS, said; --threads=1 = OMP_NUM_THREADS=1 (identical files)"
else
  bad "--threads: $(grep Threads "$WORK/thr.opt/log") / $(grep Threads "$WORK/thr.var/log")"
fi
d="$(cd "$WORK/thr.opt" && OMP_NUM_THREADS=3 "$AZURE2_BIN" --describe T1.azr --threads=2 2>/dev/null | grep -F "  threads = ")"
[ "$d" = "  threads = 2 (--threads, over OMP_NUM_THREADS=3)" ] && ok "--describe: $d" || bad "--describe: '$d'"
dj="$(cd "$WORK/thr.opt" && env -u OMP_NUM_THREADS "$AZURE2_BIN" --describe T1.azr --json --threads=2 2>/dev/null)"
if printf '%s' "$dj" | python3 -c 'import json,sys; d=json.load(sys.stdin); sys.exit(0 if d["threads"]=={"value": 2, "source": "--threads"} and len(d["switches"])==14 else 1)'; then
  ok "--describe --json: threads and the fourteen switches"
else
  bad "--describe --json: no threads or switches member"
fi

echo
if [ "$failures" -gt 0 ]; then echo "FAIL: $failures check(s)"; exit 1; fi
echo "all checks passed"
