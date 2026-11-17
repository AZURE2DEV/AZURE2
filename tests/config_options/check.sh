#!/usr/bin/env bash
#
# The options of the calculation as <config> keys (docs/upgrade/project_format3.md section 3.12), in the program.
#
#   1. key = flag: for every option and both of its values, tests/upgrade/T1 with the key in <config> writes the same
#      output files, byte for byte, as T1 with the program's flag (or, for the default value of an option that has no
#      flag for it, as T1 alone).  A value that moves T1's numbers is shown to (its run differs from T1 alone).
#   2. flag over key over default: --use-brune on a project with use_brune=0 gives T1 alone's files.
#   3. refusals: use_brune=yes (E_CONFIG_OPTION_VALUE), use_park=1 with use_brune=0 and use_rmc=1 alone
#      (E_CONFIG_OPTIONS_CONTRADICT, both keys named), each by --validate (exit 2) and by the run (exit 1, no output);
#      --use-rmc on a project with use_brune=1 (a flag against a key: refused by the run, the sources named).
#      Contradicting flags without a key (--use-rmc alone, --use-park with --no-brune in either order,
#      --no-ec-convection-current alone): refused by the run with what, allowed and the fix in flags.
#      The removed option use_gsl_coul: the key (either value) by --validate and the run, the flag by the run.
#   4. --describe: each option's value and source (default, project, command line), for every project, text and --json.
#
#   ./tests/config_options/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
. "$HERE/../lib/edit_guard.sh"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
T1="$HERE/../upgrade/T1"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/config_options.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
cp -R "$T1/data" "$WORK/data"  # --validate reads the data files beside the project

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# with_keys FILE KEY=VALUE...: T1.azr with these lines before </config>, written to FILE (the edit must bite).
add_keys() { awk -v keys="$KEYS" '/^<\/config>/ { n = split(keys, k, " "); for (i = 1; i <= n; i++) print k[i] } { print }' "$1"; }
with_keys() {
  local dest="$1"; shift
  KEYS="$*" edited "$T1/T1.azr" "$dest" add_keys
}
# run NAME PROJECT [FLAGS...]: mode 1 of PROJECT (a T1.azr variant) in $WORK/NAME; sets $rc.
run() {
  local dir="$WORK/$1" proj="$2"; shift 2
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$T1/data" "$dir/data"
  cp "$proj" "$dir/T1.azr"
  (cd "$dir" && "$AZURE2_BIN" calculate --no-gui "$@" T1.azr > run.log 2>&1)
  rc=$?
}
same() { diff -r "$WORK/$1/output" "$WORK/$2/output" > /dev/null 2>&1; }

run base "$T1/T1.azr"
[ "$rc" = 0 ] && [ -s "$WORK/base/output/chiSquared.out" ] && ok "T1 alone runs" || bad "T1 alone: exit $rc"

# 1. key = flag, both values of every option.  A row: key value flags-for-that-value moves-T1 (1: the run must differ
# from T1 alone).  Options that need a partner (use_rmc needs Brune off, ec_convection_current=0 the exact EL form)
# carry it in both runs.
n=0
while read -r key value flags moves partner; do
  [ -z "$key" ] && continue
  n=$((n + 1))
  partner_key=""; partner_flag=""
  if [ "$partner" != "-" ]; then partner_key="${partner%%:*}"; partner_flag="${partner#*:}"; fi
  with_keys "$WORK/k$n.azr" "$key=$value" $partner_key || { bad "$key=$value: the key was not added"; continue; }
  run "k$n" "$WORK/k$n.azr"
  krc=$rc
  if [ "$flags" = "-" ] && [ -z "$partner_flag" ]; then ref=base; else
    run "f$n" "$T1/T1.azr" $( [ "$flags" = "-" ] || echo "$flags" ) $partner_flag
    ref="f$n"
  fi
  if [ "$krc" = 0 ] && [ -s "$WORK/k$n/output/chiSquared.out" ] && same "k$n" "$ref"; then
    label="$flags"; [ "$flags" = "-" ] && label="no flag"
    ok "$key=$value ${partner_key:+($partner_key) }= $label${partner_flag:+ $partner_flag}: output files identical"
  else
    bad "$key=$value: exit $krc, files differ from the run with $flags ($(diff -rq "$WORK/k$n/output" "$WORK/$ref/output" 2>&1 | sed -n '1,3p'))"
  fi
  if [ "$moves" = 1 ]; then
    same "k$n" base && bad "$key=$value: T1's files did not move (the key did nothing)" || ok "$key=$value: T1's files move"
  fi
done <<'ROWS'
use_brune 1 --use-brune 0 -
use_brune 0 --no-brune 1 -
use_park 1 --use-park 1 -
use_park 0 - 0 -
transform 1 - 0 -
transform 0 --no-transform 1 -
ignore_externals 1 --ignore-externals 0 -
ignore_externals 0 --no-ignore-externals 0 -
use_long_wavelength 1 - 0 -
use_long_wavelength 0 --no-long-wavelength 1 -
ec_convection_current 1 --ec-convection-current 0 -
ec_convection_current 0 --no-ec-convection-current 1 use_long_wavelength=0:--no-long-wavelength
use_rmc 0 - 0 -
use_rmc 1 --use-rmc 1 use_brune=0:--no-brune
use_coulomb_table 1 --coulomb-table 0 -
use_coulomb_table 0 --no-coulomb-table 0 -
use_ec_table 1 --ec-table 0 -
use_ec_table 0 --no-ec-table 1 -
replace_grids 1 --grid-replacement 0 -
replace_grids 0 --no-grid-replacement 0 -
ROWS
[ "$n" = 20 ] && ok "20 cases: every option, both values" || bad "$n cases, expected 20"

# 2. a flag wins over the key
with_keys "$WORK/brune0.azr" "use_brune=0"
run over "$WORK/brune0.azr" --use-brune
[ "$rc" = 0 ] && same over base && ok "--use-brune over use_brune=0: T1 alone's files" || bad "--use-brune over use_brune=0: exit $rc or files differ"

# 3. refusals
refused() {  # refused NAME CODE TEXT... : --validate exit 2 and the run exit 1 without output, the code and texts printed
  local name="$1" code="$2"; shift 2
  local v; v="$(cd "$WORK" && "$AZURE2_BIN" --validate "$name.azr" 2>&1)"; local vrc=$?
  run "$name" "$WORK/$name.azr"
  local all=1 t
  for t in "$code" "$@"; do
    grep -qF -- "$t" <<< "$v" || all=0
    grep -qF -- "$t" "$WORK/$name/run.log" || all=0
  done
  if [ "$vrc" = 2 ] && [ "$rc" = 1 ] && [ ! -s "$WORK/$name/output/chiSquared.out" ] && [ "$all" = 1 ]; then
    ok "$name: --validate exit 2, the run exit 1 without output, $code $*"
  else
    bad "$name: --validate exit $vrc, run exit $rc, texts found $all: $(printf '%s' "$v" | grep -F "$code" | sed -n 1p)"
  fi
}
with_keys "$WORK/badvalue.azr" "use_brune=yes" && refused badvalue E_CONFIG_OPTION_VALUE "use_brune=yes"
with_keys "$WORK/parkbrune.azr" "use_park=1" "use_brune=0" && \
  refused parkbrune E_CONFIG_OPTIONS_CONTRADICT "use_park=1 (this file)" "use_brune=0 (this file)"
with_keys "$WORK/rmc.azr" "use_rmc=1" && refused rmc E_CONFIG_OPTIONS_CONTRADICT "use_rmc=1 (this file)" "use_brune=1 (the default)"
# The removed option use_gsl_coul (2026-11-16; until then two rows of the table above, key = flag for both values): the
# key with either value, and the flag, are refused with the one message that says why and what to do.
GSL_WHY="wrong by up to 5 % in the penetrability inside the outer turning point"
for value in 0 1; do
  with_keys "$WORK/gslkey$value.azr" "use_gsl_coul=$value" && \
    refused "gslkey$value" E_CONFIG_KEY_REMOVED "use_gsl_coul is no longer a key of <config>: delete the line." "$GSL_WHY"
done
run gslflag "$T1/T1.azr" --gsl-coul
if [ "$rc" = 2 ] && [ ! -s "$WORK/gslflag/output/chiSquared.out" ] \
   && grep -qF -- "--gsl-coul is no longer an option of AZURE2: leave it out." "$WORK/gslflag/run.log" \
   && grep -qF -- "$GSL_WHY" "$WORK/gslflag/run.log"; then
  ok "--gsl-coul: refused (exit 2, no output), saying why and to leave it out"
else
  bad "--gsl-coul: exit $rc, $(sed -n 1,3p "$WORK/gslflag/run.log")"
fi
described="$("$AZURE2_BIN" --describe "$T1/T1.azr" 2>&1)"
grep -q "use_gsl_coul" <<< "$described" && bad "--describe still lists use_gsl_coul" || ok "--describe does not list use_gsl_coul"
with_keys "$WORK/brune1.azr" "use_brune=1"
v="$(cd "$WORK" && "$AZURE2_BIN" --validate brune1.azr 2>&1)"; vrc=$?
run flagkey "$WORK/brune1.azr" --use-rmc
if [ "$vrc" = 0 ] && [ "$rc" = 1 ] && grep -qF "E_CONFIG_OPTIONS_CONTRADICT" "$WORK/flagkey/run.log" \
   && grep -qF "use_rmc=1 (command line)" "$WORK/flagkey/run.log" && grep -qF "use_brune=1 (this file)" "$WORK/flagkey/run.log"; then
  ok "--use-rmc on a project with use_brune=1: valid alone, the run refused naming both and their sources"
else
  bad "--use-rmc with use_brune=1: --validate exit $vrc, run exit $rc"
fi

# Contradicting flags without a key: refused by the run (exit 1, no output) as a session refuses the keywords: the
# same code and the same three parts (what, allowed, fix), the fix in flags.
flags_refused() {  # flags_refused NAME TEXT -- FLAGS...
  local name="$1" text="$2"; shift 3
  run "$name" "$T1/T1.azr" "$@"
  if [ "$rc" = 1 ] && [ ! -s "$WORK/$name/output/chiSquared.out" ] && grep -qF "E_CONFIG_OPTIONS_CONTRADICT" "$WORK/$name/run.log" \
     && grep -qF -- "$text" "$WORK/$name/run.log"; then
    ok "$*: the run refused (exit 1, no output): $text"
  else
    bad "$*: exit $rc, $(grep -F ERROR "$WORK/$name/run.log" | sed -n 1p)"
  fi
}
flags_refused rmcflag "use_rmc=1 (command line) with use_brune=1 (the default) cannot hold together: the Reich-Moore \
capture is computed with constant boundary conditions only (Brune is on by default); allowed: use_rmc=1 with use_brune=0, \
or use_rmc=0; fix: add --no-brune (or leave out --use-rmc)" -- --use-rmc
flags_refused parkfirst "use_park=1 (command line) with use_brune=0 (command line) cannot hold together: Park's observed \
amplitudes are defined with Brune's shift functions; allowed: use_park=1 with use_brune=1, or use_park=0; fix: --use-brune \
in place of --no-brune (or leave out --use-park)" -- --use-park --no-brune
flags_refused brunefirst "use_park=1 (command line) with use_brune=0 (command line) cannot hold together" -- --no-brune --use-park
flags_refused convection "ec_convection_current=0 (command line) with use_long_wavelength=1 (the default) cannot hold \
together" -- --no-ec-convection-current

# 4. --describe: value and source of each option
d="$(cd "$WORK/k2" && "$AZURE2_BIN" --describe T1.azr --no-transform 2>&1)"
if grep -qF "  use_brune = 0 (project)" <<< "$d" && grep -qF "  transform = 0 (command line)" <<< "$d" \
   && grep -qF "  use_coulomb_table = 1 (default)" <<< "$d"; then
  ok "--describe: use_brune 0 from the project, transform 0 from the command line, use_coulomb_table 1 by default"
else
  bad "--describe: $(printf '%s' "$d" | grep -A12 'options of the calculation')"
fi
d0="$(cd "$WORK/base" && "$AZURE2_BIN" --describe T1.azr 2>&1)"
if [ "$(printf '%s\n' "$d0" | grep -A10 -F "options of the calculation (value, source):" | grep -c ' (default)$')" = 10 ]; then
  ok "--describe of T1 alone: the options block, all ten at their default"
else
  bad "--describe of T1 alone: no options block with the ten defaults"
fi
dj="$(cd "$WORK/k2" && "$AZURE2_BIN" --describe T1.azr --json 2>/dev/null)"
printf '%s' "$dj" | python3 -c 'import json,sys; d=json.load(sys.stdin); o=d["options"]; sys.exit(0 if o["use_brune"]=={"value": False, "source": "project"} and len(o)==10 else 1)' \
  && ok "--describe --json: options with value and source" || bad "--describe --json: no options member"

echo
if [ "$failures" -gt 0 ]; then echo "FAIL: $failures check(s)"; exit 1; fi
echo "all checks passed"
