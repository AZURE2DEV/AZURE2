#!/bin/sh
# compare_erya.sh BUILD_DIR OUT.tsv: the table of this module's stopping against ERYA's (docs/upgrade/stopping.md
# section 9).  Needs a build with USE_ERYA (target stopping_erya_compare); runs from external/, where ERYA finds its
# coefficient file, and prints the largest differences per projectile and material.
set -e
[ $# -eq 2 ] || { echo "usage: $0 BUILD_DIR OUT.tsv" >&2; exit 2; }
here=$(cd "$(dirname "$0")" && pwd)
tool="$1/stopping/stopping_erya_compare"
[ -x "$tool" ] || { echo "$tool not built (configure with USE_ERYA)" >&2; exit 2; }
out=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
(cd "$here/../../external" && "$tool") > "$out"
awk -F'\t' 'NR > 1 { k = $1 " " $2; d = $7 < 0 ? -$7 : $7; if (d > m[k]) { m[k] = d; at[k] = $3; s[k] = $7 } }
            END { for (k in m) printf "%-12s max |rel diff| %.3f at %s MeV (%+.3f)\n", k, m[k], at[k], s[k] }' "$out" | sort
