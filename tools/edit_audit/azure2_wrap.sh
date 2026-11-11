#!/usr/bin/env bash
# The audit's stand-in for the AZURE2 binary (tools/edit_audit.sh passes it to a check as its binary).  It runs the
# real binary ($EDIT_AUDIT_REAL_AZURE2) with the same arguments, input and environment, and appends one line per run
# to $EDIT_AUDIT_LOG:
#   RUN <input key> <output key> <exit status> <project>   <tab> arguments
# The input key is a hash of what decides the run: the arguments, the standard input, the AZURE*/OMP* environment,
# and every file of the project's folder except output/ and checks/ (the project, its data, a parameter file).
# The output key is a hash of the exit status, the standard output and the files the run wrote or changed in that
# folder.  Two runs of one check with
# different input keys and one output key are variants that made no difference.

real="${EDIT_AUDIT_REAL_AZURE2:?}"
log="${EDIT_AUDIT_LOG:-/dev/null}"
tmp="$(mktemp -d "${EDIT_AUDIT_TMP:-/tmp}/azwrap.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT
args="$(printf '%q ' "$@")"
args="${args//$'\n'/\\n}"
awk_real="${EDIT_AUDIT_REAL_awk:-awk}"

# The project: the first argument that is an existing file ending in .azr; its folder is what the run reads.
proj=""
for a in "$@"; do
  case "$a" in *.azr) [ -f "$a" ] && { proj="$a"; break; } ;; esac
done
dir=""
[ -n "$proj" ] && dir="$(cd "$(dirname "$proj")" && pwd)"

fdhash() {  # what the run reads: the project, the other files beside it (not other projects), its data/ folder and
            # every argument that is a file: "sha path" lines
  {
    [ -n "$dir" ] && (cd "$dir" && find . -maxdepth 1 -type f ! -name '*.azr' ! -name '*.log' -print; \
                       [ -d data ] && find data -type f -print) 2>/dev/null | "${EDIT_AUDIT_REAL_sed:-sed}" "s|^|$dir/|"
    for a in "$@"; do [ -f "$a" ] && printf '%s\n' "$a"; done
  } | LC_ALL=C sort -u | while read -r f; do
    [ -f "$f" ] && printf '%s %s\n' "$(sha1sum < "$f" | cut -c1-40)" "$(basename "$f")"
  done
}

# Standard input is read here and handed on (the menu answers a check pipes in); a terminal or /dev/null is left.
fd0="$(readlink /proc/$$/fd/0 2>/dev/null)"
stdin_used=0
if [[ "$fd0" == pipe:* ]] || [ -f "$fd0" ]; then cat > "$tmp/stdin"; stdin_used=1; else : > "$tmp/stdin"; fi

{
  # The arguments with every file named by its content (below), not its path: copies of one project in two folders
  # are the same input.
  for a in "$@"; do if [ -f "$a" ]; then echo "arg FILE"; else printf 'arg %s\n' "$a"; fi; done
  env | LC_ALL=C sort | grep -E '^(AZURE|OMP_)' | grep -vE '^(AZURE2_COULOMB_TABLE)='
  sha1sum < "$tmp/stdin"
  fdhash "$@"
} > "$tmp/inkey"
inkey="$(sha1sum < "$tmp/inkey" | cut -c1-12)"

touch "$tmp/stamp"
sleep 0.01
# Standard output is copied (tee) for the key; standard error is left alone, so the order of the two is kept where a
# check merges them.  A pipe, not a process substitution: tee has written everything when this script ends.
if [ $stdin_used = 1 ]; then "$real" "$@" < "$tmp/stdin" | tee "$tmp/stdout"; else "$real" "$@" | tee "$tmp/stdout"; fi
st=${PIPESTATUS[0]}

# The output key: the exit status, the standard output (lines with a time, a path or a thread count left out) and
# the files the run wrote.
outkey="-"
{
  echo "status $st"
  "$awk_real" '!/[0-9][0-9]:[0-9][0-9]:[0-9][0-9]|[Ee]lapsed|seconds|[0-9] s$|[0-9]s\)|wall|time:|thread|\//' "$tmp/stdout"
} > "$tmp/outkey"
if [ -n "$dir" ]; then
  (cd "$dir" && find . -type f -newer "$tmp/stamp" ! -name 'run.log' ! -name '*.log' -size -50M -print 2>/dev/null \
     | LC_ALL=C sort | while read -r f; do
       # Lines that carry a time or a date are left out of the key.
       printf '%s %s\n' "$("$awk_real" '!/[0-9][0-9]:[0-9][0-9]:[0-9][0-9]|[Ee]lapsed|seconds|[0-9] s$|wall|time:/' "$f" | sha1sum | cut -c1-40)" "$f"
     done) >> "$tmp/outkey"
fi
outkey="$(sha1sum < "$tmp/outkey" | cut -c1-12)"
printf '%s\tAZURE2\tRUN %s %s %s %s\t%s\n' "${EDIT_AUDIT_CHECK:-?}" "$inkey" "$outkey" "$st" "${proj:--}" "$args" >> "$log"
if [ -n "${EDIT_AUDIT_KEEP:-}" ]; then
  mkdir -p "$EDIT_AUDIT_KEEP"
  cp "$tmp/inkey" "$EDIT_AUDIT_KEEP/$inkey.in"
  [ -n "$proj" ] && cp "$proj" "$EDIT_AUDIT_KEEP/$inkey.azr"
  [ -f "$tmp/outkey" ] && cp "$tmp/outkey" "$EDIT_AUDIT_KEEP/$inkey.$outkey.out"
fi
exit $st
