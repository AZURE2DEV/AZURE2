# Sourced by tests/*/check.sh that make a variant of a project, data or parameter file by editing its text: the edit
# fails when it changes nothing, so a pattern that no longer matches (a line written differently after a conversion)
# stops the check instead of letting it compare a run with itself (docs/upgrade/format3_tree.md section 7).
#
#   edited SRC DEST CMD [ARGS...]   DEST = the output of `CMD ARGS... SRC` (awk, sed, perl, a shell function that
#                                   reads its last argument); exit status 1 with a message on stderr when CMD fails or
#                                   when DEST equals SRC byte for byte.  DEST may be SRC (the file is replaced).
#   must_differ A B [WHAT]          exit status 1 with a message when the files A and B are equal (a variant made some
#                                   other way, e.g. by a heredoc or pyazr, that must not be the base).
#
# Use it as `edited ... || bad "..."` or under `set -e`; the message names the command, so the failing edit is found.

edited() {
  local src="$1" dest="$2" tmp
  shift 2
  tmp="$(mktemp "${TMPDIR:-/tmp}/edited.XXXXXX")" || return 1
  if ! "$@" "$src" > "$tmp"; then
    echo "edited: '$*' failed on $src" >&2
    rm -f "$tmp"
    return 1
  fi
  if cmp -s "$src" "$tmp"; then
    echo "edited: '$*' changed nothing in $src (the edit hit nothing)" >&2
    rm -f "$tmp"
    return 1
  fi
  mv "$tmp" "$dest"
}

must_differ() {
  if cmp -s "$1" "$2"; then
    echo "must_differ: ${3:-$2} is the same as $1 (the variant is the base)" >&2
    return 1
  fi
}
