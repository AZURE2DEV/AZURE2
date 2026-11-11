#!/usr/bin/env bash
# The audit's stand-in for sed, awk, perl and python3 (tools/edit_audit.sh puts copies named after them first on
# PATH).  It runs the real tool unchanged and appends one line per call to $EDIT_AUDIT_LOG saying whether the call
# changed what it was given:
#   INPLACE  <same|changed> <file>       an in-place edit (sed -i, perl -i, awk -i inplace)
#   FILTER   <same|differs> <input>      a filter whose input is a file argument, a redirected file or a pipe:
#                                        "same" means its output is its input byte for byte (an edit that hit nothing)
#   PYTHON   <n changed> <files>         a python3 call: the files under $EDIT_AUDIT_SCOPE it changed
# followed by a tab and the arguments.  A call the shim cannot judge (no input it can see) is logged as OTHER.
# Nothing here changes what the tool does; the comparison runs after it.

self="$(basename "$0")"
real_var="EDIT_AUDIT_REAL_$self"
real="${!real_var}"
[ -n "$real" ] || { echo "edit_audit shim: no real $self" >&2; exit 127; }
log="${EDIT_AUDIT_LOG:-/dev/null}"
tmp="$(mktemp -d "${EDIT_AUDIT_TMP:-/tmp}/shim.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT
args="$(printf '%q ' "$@")"
args="${args//$'\n'/\\n}"
emit() { printf '%s\t%s\t%s\t%s\n' "${EDIT_AUDIT_CHECK:-?}" "$self" "$1" "$args" >> "$log"; }
sha() { sha1sum "$1" 2>/dev/null | cut -c1-40; }

# python3: compare the project-like files of the check's scratch before and after.
if [ "$self" = python3 ]; then
  scope="${EDIT_AUDIT_SCOPE:-}"
  list() {
    [ -n "$scope" ] && [ -d "$scope" ] || return 0
    find "$scope" -path '*/output' -prune -o -path '*/checks' -prune -o -path '*/shim.*' -prune -o -type f \
      \( -name '*.azr' -o -name '*.dat' -o -name '*.par' -o -name '*.sav' -o -name '*.txt' -o -name '*.csv' \
         -o -name '*.json' -o -name '*.extrap' \) -size -20M -print 2>/dev/null | sort | while read -r f; do
      printf '%s %s\n' "$(sha "$f")" "$f"
    done
  }
  list > "$tmp/before"
  "$real" "$@"
  st=$?
  list > "$tmp/after"
  changed="$(diff "$tmp/before" "$tmp/after" | "$EDIT_AUDIT_REAL_awk" '/^>/ { print $3 }' | sort -u)"
  n=0; [ -n "$changed" ] && n="$(printf '%s\n' "$changed" | wc -l)"
  emit "PYTHON $n $(printf '%s' "$changed" | tr '\n' ' ')"
  exit $st
fi

# Which arguments are input files, and is the call in place?
inplace=0
inputs=()
skipnext=0
prog_seen=0
program=""
case "$self" in
  sed)
    for a in "$@"; do
      if [ $skipnext = 1 ]; then skipnext=0; continue; fi
      case "$a" in
        -i*|--in-place*) inplace=1 ;;
        -e|-f) skipnext=1; prog_seen=1 ;;
        -*) ;;
        *) if [ $prog_seen = 0 ]; then prog_seen=1; elif [ -f "$a" ]; then inputs+=("$a"); fi ;;
      esac
    done ;;
  awk)
    for a in "$@"; do
      if [ $skipnext = 1 ]; then
        skipnext=0
        [ "$prev" = "-i" ] && [ "$a" = inplace ] && inplace=1
        [ "$prev" = "-f" ] && prog_seen=1
        continue
      fi
      case "$a" in
        -v|-f|-F|-i) skipnext=1; prev="$a" ;;
        -v*|-F*) ;;
        -*) ;;
        *) if [ $prog_seen = 0 ]; then prog_seen=1; program="$a"
           # an operand var=value is an assignment only when var is a name (AZUREOut_aa=1_R=2.out is a file)
           elif [[ "$a" =~ ^[A-Za-z_][A-Za-z0-9_]*= ]]; then :
           elif [ -f "$a" ]; then inputs+=("$a"); fi ;;
      esac
    done ;;
  perl)
    reads=0
    for a in "$@"; do
      if [ $skipnext = 1 ]; then skipnext=0; continue; fi
      case "$a" in
        -e|-E|-M*) [ "$a" = -e ] || [ "$a" = -E ] && skipnext=1; prog_seen=1 ;;
        -*)
          [[ "$a" =~ ^-[a-zA-Z0-9]*i ]] && inplace=1
          [[ "$a" =~ ^-[a-zA-Z0-9]*[np] ]] && reads=1 ;;
        *) if [ $prog_seen = 0 ]; then prog_seen=1; elif [ -f "$a" ]; then inputs+=("$a"); fi ;;
      esac
    done ;;
esac

if [ $inplace = 1 ]; then
  for f in "${inputs[@]}"; do cp "$f" "$tmp/before.$(printf '%s' "$f" | sha1sum | cut -c1-8)"; done
  "$real" "$@"
  st=$?
  for f in "${inputs[@]}"; do
    b="$tmp/before.$(printf '%s' "$f" | sha1sum | cut -c1-8)"
    if cmp -s "$b" "$f"; then emit "INPLACE same $f"; else emit "INPLACE changed $f"; fi
  done
  exit $st
fi

# A filter: its input is the file arguments, else standard input.  A pipe is read here and handed on (not when the
# program reads no input: an awk BEGIN program, perl without -n/-p).
src=""
if [ ${#inputs[@]} -gt 0 ]; then
  cat "${inputs[@]}" > "$tmp/in" 2>/dev/null
  src="${inputs[*]}"
else
  fd0="$(readlink /proc/$$/fd/0 2>/dev/null)"
  noread=0
  [ "$self" = awk ] && [[ "$program" =~ ^[[:space:]]*BEGIN ]] && noread=1
  [ "$self" = perl ] && [ "${reads:-0}" = 0 ] && noread=1
  if [ $noread = 0 ] && [ -f "$fd0" ]; then
    cat > "$tmp/in"; src="<$fd0"
  elif [ $noread = 0 ] && [[ "$fd0" == pipe:* ]]; then
    cat > "$tmp/in"; src="<pipe"
  fi
fi
if [ -z "$src" ]; then
  "$real" "$@"
  st=$?
  emit "OTHER"
  exit $st
fi
if [ ${#inputs[@]} -gt 0 ]; then "$real" "$@" > "$tmp/out"; else "$real" "$@" < "$tmp/in" > "$tmp/out"; fi
st=$?
cat "$tmp/out"
kind=text
grep -qE '^[[:space:]]*(<(config|levels|data|targets|thm|segmentsData)>|dataset[[:space:]])' "$tmp/in" 2>/dev/null && kind=project
[[ "$src" =~ \.(azr|dat|par|sav|extrap)( |$) ]] && kind=project
if [ -s "$tmp/in" ] && cmp -s "$tmp/in" "$tmp/out"; then emit "FILTER same $kind $src"; else emit "FILTER differs $kind $src"; fi
exit $st
