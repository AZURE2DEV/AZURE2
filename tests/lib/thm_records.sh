# Sourced by the THM check.sh scripts: the THM part of a format-3 project (project_format3.md section 3.9) edited on
# its records.  A data set joins an experiment by its own field `thm <experiment>` (`thm global`: none); the <thm>
# section holds the global options and each experiment's keys, without a member list.
#
#   thm_section FILE BLOCK          FILE's <thm> section replaced by BLOCK (lines between <thm> and </thm>), or
#                                   appended when FILE has none; an empty BLOCK removes the section
#   thm_join FILE EXP DATASET...    the named data sets of FILE name experiment EXP (`thm EXP`); EXP = global
#                                   leaves them in none
#   thm_set FILE BLOCK DATASET...   thm_section FILE BLOCK, then the named data sets join the one experiment BLOCK
#                                   defines (thm_join); a BLOCK without an experiment leaves them in none, one with
#                                   several is an error (join them one by one)
#   thm_weight FILE RECORD TABLE    the data set or extrapolation RECORD of FILE weighted by the table TABLE
#                                   (`thm_weight TABLE`, a file relative to the project)
#   thm_datasets FILE               the names of FILE's data sets, in the order of <data>, one per line
#   THM_NOFOLD                      an awk program (awk "$THM_NOFOLD" FILE): FILE without its folds, i.e. without the
#                                   records of <targets> and the beam and quadrature fields of its data sets and
#                                   extrapolations (the THM projects' only target records are their beams)
#
# FILE is edited in place.  A data set that is not in FILE, or has no `thm` field, is an error (exit status 1 and a
# message on stderr): a misspelt name must not leave a data set out of its experiment unnoticed.

thm_section() {
  local f="$1" block="$2"
  awk -v block="$block" '
    /^<thm>/ { skip = 1; if (block != "") printf "<thm>\n%s\n</thm>\n", block; done = 1; next }
    skip && /^<\/thm>/ { skip = 0; next }
    !skip { print }
    END { if (!done && block != "") printf "<thm>\n%s\n</thm>\n", block }' "$f" > "$f.thm_tmp" && mv "$f.thm_tmp" "$f"
}

thm_join() {
  local f="$1" exp="$2"
  shift 2
  awk -v x="$exp" -v names="$*" '
    BEGIN { n = split(names, a, " "); for (i = 1; i <= n; i++) want[a[i]] = 1 }
    /^<data>/ { data = 1 } /^<\/data>/ { data = 0 }
    data && $1 == "dataset" { cur = $2 }
    data && cur in want && $1 == "thm" { sub(/thm[ \t]+[^ \t#]+/, sprintf("%-14s %s", "thm", x)); joined[cur] = 1 }
    { print }
    END { for (k in want) if (!(k in joined)) { print "thm_join: no THM data set " k > "/dev/stderr"; bad = 1 } exit bad }' \
    "$f" > "$f.thm_tmp" && mv "$f.thm_tmp" "$f" || { rm -f "$f.thm_tmp"; return 1; }
}

thm_set() {
  local f="$1" block="$2" names
  shift 2
  thm_section "$f" "$block" || return 1
  names="$(printf '%s\n' "$block" | sed -n 's/^[[:space:]]*experiment\[\([^]]*\)\].*/\1/p' | sort -u)"
  [ -n "$names" ] && [ $# -gt 0 ] || return 0
  [ "$(printf '%s\n' "$names" | wc -l)" -eq 1 ] || { echo "thm_set: the block defines several experiments" >&2; return 1; }
  thm_join "$f" "$names" "$@"
}

THM_NOFOLD='/^<targets>/ { T = 1; print; next } /^<\/targets>/ { T = 0 } T { next } $1 == "beam" || $1 == "quadrature" { next } { print }'

thm_weight() {
  awk -v r="$2" -v t="$3" '
    { print }
    ($1 == "dataset" || $1 == "extrapolation") && $2 == r { printf "    %-14s %s\n", "thm_weight", t; found = 1 }
    END { if (!found) { print "thm_weight: no data set or extrapolation " r > "/dev/stderr"; exit 1 } }' \
    "$1" > "$1.thm_tmp" && mv "$1.thm_tmp" "$1" || { rm -f "$1.thm_tmp"; return 1; }
}

thm_datasets() { awk '/^<data>/ { d = 1 } /^<\/data>/ { d = 0 } d && $1 == "dataset" { print $2 }' "$1"; }
